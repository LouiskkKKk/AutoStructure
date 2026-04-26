"""
validation_refine.py
--------------------
Post-pipeline structure auto-fix loop using Phenix MolProbity validation
and Coot real-space refinement.

Workflow
--------
For each iteration (up to *max_iterations*, default 10):

  1. Run ``phenix.molprobity`` on the current PDB.
  2. Parse stdout + result files to extract problematic residues in five
     categories: Cbeta deviations, cis/twisted peptides, rotamer outliers,
     Ramachandran outliers, and geometry-restraints violations.
  3. If no *new* issues since the previous round (or zero issues) → stop.
  4. Build ±5-residue windows around every problem residue.
  5. Merge overlapping windows to minimise the number of Coot calls.
  6. Inject the original AF2 PDB's HELIX/SHEET records into the current PDB
     so that secondary-structure restraints are correct.
  7. Run Coot real-space refine (with secondary-structure restraints) on
     every merged window in a single headless session.
  8. Log per-round statistics and save the iteration PDB + validation report
     under *output_dir*/iteration_XX/.

Public API
----------
run_auto_fix(pdb_path, map_path, af2_pdb, output_dir, ...) -> str
    Entry point.  Returns path to the final fixed PDB.

check_dependencies(phenix_molprobity_exe, coot_exe) -> None
    Raises RuntimeError with a clear message if an executable is missing.
"""

from __future__ import annotations

import logging
import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple

log = logging.getLogger("autostructure.validation_refine")

# ---------------------------------------------------------------------------
# Issue category constants
# ---------------------------------------------------------------------------

ISSUE_CBETA        = "cbeta"
ISSUE_CIS_TWISTED  = "cis_twisted"
ISSUE_ROTAMER      = "rotamer"
ISSUE_RAMACHANDRAN = "ramachandran"
ISSUE_GEOMETRY     = "geometry"

_ALL_CATEGORIES = (
    ISSUE_CBETA,
    ISSUE_CIS_TWISTED,
    ISSUE_ROTAMER,
    ISSUE_RAMACHANDRAN,
    ISSUE_GEOMETRY,
)

# ---------------------------------------------------------------------------
# Coot RSR script template
# ---------------------------------------------------------------------------

_COOT_RSR_SCRIPT_TEMPLATE = """\
; Auto-generated Coot real-space refine script – do not edit manually
(let* ((mol  (read-pdb "{pdb_path}"))
       (map  (read-ccp4-map "{map_path}" 0)))
  (set-imol-refinement-map map)
  (make-and-draw-secondary-structure mol)
  (set-secondary-structure-restraints-type 2)
{refine_calls}  (write-pdb-file mol "{output_path}")
  (coot-real-exit 0))
"""

_COOT_REFINE_ZONE_TMPL = (
    '  (refine-zone mol "{chain_id}" {start_res} {end_res} "")\n'
    '  (accept-regularization)\n'
)

# ---------------------------------------------------------------------------
# Dependency checking
# ---------------------------------------------------------------------------

def check_dependencies(phenix_molprobity_exe: str, coot_exe: str) -> None:
    """
    Verify that the required external executables are reachable.

    Raises
    ------
    RuntimeError
        With a human-readable message naming the missing executable and a
        hint about how to configure it.
    """
    for exe, label, cfg_key in [
        (phenix_molprobity_exe, "Phenix MolProbity", "tools.phenix_molprobity"),
        (coot_exe,              "Coot",               "tools.coot"),
    ]:
        if shutil.which(exe) is None and not os.path.isfile(exe):
            raise RuntimeError(
                f"{label} executable not found: '{exe}'.\n"
                f"Install {label} and either add it to PATH or set "
                f"'{cfg_key}' in config.yaml."
            )


# ---------------------------------------------------------------------------
# Run phenix.molprobity
# ---------------------------------------------------------------------------

def run_molprobity(
    pdb_path: str,
    output_dir: str,
    phenix_molprobity_exe: str = "phenix.molprobity",
    log_path: Optional[str] = None,
) -> Tuple[str, str, int]:
    """
    Execute ``phenix.molprobity`` on *pdb_path*.

    Parameters
    ----------
    pdb_path              : PDB file to validate
    output_dir            : directory for Phenix output files
    phenix_molprobity_exe : name/path of the ``phenix.molprobity`` executable
    log_path              : optional file to capture combined stdout/stderr

    Returns
    -------
    (stdout, stderr, returncode)
    """
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    stem = Path(pdb_path).stem
    cmd = [
        phenix_molprobity_exe,
        os.path.abspath(pdb_path),
        f"output.prefix={stem}",
        f"output_dir={os.path.abspath(output_dir)}",
    ]

    log.debug("Running: %s", " ".join(cmd))
    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        check=False,
        cwd=output_dir,
    )

    # Also try to append any generated text result files to stdout for parsing
    combined_stdout = result.stdout
    for candidate in Path(output_dir).glob(f"{stem}*.txt"):
        try:
            combined_stdout += "\n" + candidate.read_text(errors="replace")
        except OSError:
            pass

    if log_path:
        Path(log_path).parent.mkdir(parents=True, exist_ok=True)
        with open(log_path, "w") as fh:
            fh.write("=== CMD ===\n")
            fh.write(" ".join(cmd) + "\n\n")
            fh.write("=== STDOUT ===\n")
            fh.write(result.stdout)
            fh.write("\n=== STDERR ===\n")
            fh.write(result.stderr)

    return combined_stdout, result.stderr, result.returncode


# ---------------------------------------------------------------------------
# Parse molprobity output
# ---------------------------------------------------------------------------

# Flexible regex matching a residue-data row.
# Captures: chain(1), resno(int), icode(optional), resname(3-letter), rest
_RESIDUE_ROW_RE = re.compile(
    r"""
    ^\s*
    ([A-Za-z0-9 ])          # chain id (1 char, may be space)
    \s+
    (-?\d+)                 # residue sequence number
    \s+
    ([A-Za-z ]{0,4})        # optional insertion code / empty field
    \s*
    ([A-Z]{3})              # residue name (3 uppercase letters)
    (?:\s+(.*))?            # rest of line (optional)
    \s*$
    """,
    re.VERBOSE,
)

# Section-header keywords (case-insensitive) that introduce each category
_SECTION_HEADERS: List[Tuple[str, str]] = [
    ("ramachandran",                ISSUE_RAMACHANDRAN),
    ("rotamer",                     ISSUE_ROTAMER),
    ("c-beta",                      ISSUE_CBETA),
    ("cbeta",                       ISSUE_CBETA),
    ("cis",                         ISSUE_CIS_TWISTED),
    ("twisted",                     ISSUE_CIS_TWISTED),
    ("geometry restraints",         ISSUE_GEOMETRY),
    ("geometry_restraints",         ISSUE_GEOMETRY),
    ("bond length",                 ISSUE_GEOMETRY),
    ("bond angle",                  ISSUE_GEOMETRY),
    ("bond deviation",              ISSUE_GEOMETRY),
    ("non-bonded",                  ISSUE_GEOMETRY),
    ("dihedral",                    ISSUE_GEOMETRY),
    ("chirality",                   ISSUE_GEOMETRY),
    ("planarity",                   ISSUE_GEOMETRY),
]

# Keywords that indicate a line is a section *header*, not data
_HEADER_KEYWORDS = re.compile(
    r"(summary|statistic|analysis|score|result|total|residues analyzed|"
    r"outlier count|percent|chain\s+resid|chain\s+resseq|^-+$|^=+$|"
    r"^#|^>|favored|allowed|not favored)",
    re.IGNORECASE,
)


def _detect_section(line: str) -> Optional[str]:
    """Return an issue-category constant if *line* looks like a section header."""
    lower = line.lower()
    for keyword, category in _SECTION_HEADERS:
        if keyword in lower:
            return category
    return None


def _is_data_row(line: str) -> bool:
    """Return True if the line looks like a tabular residue-data row."""
    stripped = line.strip()
    if not stripped or len(stripped) < 8:
        return False
    if _HEADER_KEYWORDS.search(stripped):
        return False
    return bool(_RESIDUE_ROW_RE.match(line))


def parse_molprobity_output(text: str) -> List[dict]:
    """
    Parse ``phenix.molprobity`` text output (stdout or result file).

    Returns
    -------
    list of dicts, each containing:
        chain        : chain ID (str, single character)
        resseq       : residue sequence number (int)
        icode        : insertion code (str, may be empty)
        resname      : residue name (str, 3 letters)
        atom_name    : atom name if available (str or None)
        issue_type   : one of the ISSUE_* constants
        description  : raw tail text from the parsed line
    """
    issues: List[dict] = []
    current_section: Optional[str] = None

    for line in text.splitlines():
        stripped = line.strip()
        if not stripped:
            continue

        # If the line looks like a residue-data row, parse it immediately
        # (before checking for section headers, to avoid misclassifying data
        # lines that happen to contain keywords like "CIS" or "TWISTED").
        m = _RESIDUE_ROW_RE.match(line)
        if m and current_section is not None:
            chain   = m.group(1).strip() or " "
            resseq  = int(m.group(2))
            icode   = m.group(3).strip()
            resname = m.group(4).strip()
            rest    = (m.group(5) or "").strip()

            # Skip obvious column-header rows
            if resname.upper() in ("RES", "NAME", "RESD", "CHN"):
                continue

            # Geometry section: try to extract an atom name from the tail
            atom_name = None
            if current_section == ISSUE_GEOMETRY:
                atom_m = re.search(r'\b([A-Z][A-Z0-9]{0,3})\b', rest)
                if atom_m:
                    atom_name = atom_m.group(1)

            issues.append({
                "chain":       chain,
                "resseq":      resseq,
                "icode":       icode,
                "resname":     resname,
                "atom_name":   atom_name,
                "issue_type":  current_section,
                "description": rest,
            })
            continue

        # Check whether this line introduces a new section
        detected = _detect_section(stripped)
        if detected is not None:
            current_section = detected
            continue

        if current_section is None:
            continue

        # A non-matching, non-empty line after a section header resets
        # section only if it looks like a top-level header itself
        if re.match(r'^[A-Z][A-Z ]{5,}$', stripped):
            current_section = None

    return issues


# ---------------------------------------------------------------------------
# Collect residues from PDB
# ---------------------------------------------------------------------------

def collect_chain_residues(pdb_path: str) -> Dict[str, List[int]]:
    """
    Return a dict mapping each chain ID to a sorted list of residue sequence
    numbers present in the PDB's ATOM/HETATM records.
    """
    chains: Dict[str, set] = {}
    with open(pdb_path) as fh:
        for line in fh:
            if not (line.startswith("ATOM  ") or line.startswith("HETATM")):
                continue
            chain = line[21]
            try:
                resseq = int(line[22:26].strip())
            except ValueError:
                continue
            chains.setdefault(chain, set()).add(resseq)
    return {ch: sorted(res) for ch, res in chains.items()}


# ---------------------------------------------------------------------------
# Build and merge refinement windows
# ---------------------------------------------------------------------------

def build_refine_windows(
    issues: List[dict],
    chain_residues: Dict[str, List[int]],
    window_half: int = 5,
) -> List[Tuple[str, int, int]]:
    """
    For every issue, create a window of up to *window_half* residues on each
    side of the problem residue (clamped to chain boundaries).

    Parameters
    ----------
    issues          : list of issue dicts (from ``parse_molprobity_output``)
    chain_residues  : {chain_id: sorted residue list} (from ``collect_chain_residues``)
    window_half     : half-window size (default 5 → at most 11 residues)

    Returns
    -------
    list of (chain_id, start_resseq, end_resseq) tuples (unsorted, may overlap)
    """
    windows: List[Tuple[str, int, int]] = []
    for issue in issues:
        chain  = issue["chain"]
        resseq = issue["resseq"]
        res_list = chain_residues.get(chain, [])
        if not res_list:
            continue

        try:
            idx = res_list.index(resseq)
        except ValueError:
            # Residue not in list; find the closest one
            idx = min(range(len(res_list)), key=lambda i: abs(res_list[i] - resseq))

        lo = max(0, idx - window_half)
        hi = min(len(res_list) - 1, idx + window_half)
        windows.append((chain, res_list[lo], res_list[hi]))

    return windows


def merge_windows(
    windows: List[Tuple[str, int, int]],
) -> List[Tuple[str, int, int]]:
    """
    Merge overlapping or adjacent windows on the same chain.

    Overlapping / touching windows are merged into a single larger window to
    reduce the number of Coot calls.

    Parameters
    ----------
    windows : list of (chain_id, start_resseq, end_resseq)

    Returns
    -------
    Sorted, deduplicated, merged list of (chain_id, start_resseq, end_resseq).
    """
    if not windows:
        return []

    # Sort by chain then start residue
    sorted_wins = sorted(windows, key=lambda w: (w[0], w[1]))
    merged: List[Tuple[str, int, int]] = [sorted_wins[0]]

    for chain, start, end in sorted_wins[1:]:
        last_chain, last_start, last_end = merged[-1]
        if chain == last_chain and start <= last_end + 1:
            # Overlapping or adjacent: extend
            merged[-1] = (last_chain, last_start, max(last_end, end))
        else:
            merged.append((chain, start, end))

    return merged


# ---------------------------------------------------------------------------
# Copy secondary-structure records from source PDB to target PDB
# ---------------------------------------------------------------------------

_SS_RECORD_PREFIXES = ("HELIX ", "SHEET ")


def copy_ss_records(
    source_pdb: str,
    target_pdb: str,
    output_pdb: str,
) -> str:
    """
    Replace the HELIX/SHEET records in *target_pdb* with those from *source_pdb*,
    writing the result to *output_pdb*.

    Logic
    -----
    * All non-HELIX/SHEET lines from *target_pdb* are preserved in order.
    * HELIX/SHEET lines from *target_pdb* are discarded.
    * HELIX/SHEET lines from *source_pdb* are inserted before the first ATOM
      record in the output (standard PDB placement).

    Parameters
    ----------
    source_pdb : PDB whose HELIX/SHEET records provide the reference SS
    target_pdb : PDB to update (ATOM records come from here)
    output_pdb : destination path for the merged PDB

    Returns
    -------
    str  The *output_pdb* path.
    """
    # Read source SS records
    ss_lines: List[str] = []
    with open(source_pdb) as fh:
        for line in fh:
            if any(line.startswith(p) for p in _SS_RECORD_PREFIXES):
                ss_lines.append(line if line.endswith("\n") else line + "\n")

    # Read target lines, skipping its own SS records
    target_lines: List[str] = []
    with open(target_pdb) as fh:
        for line in fh:
            if not any(line.startswith(p) for p in _SS_RECORD_PREFIXES):
                target_lines.append(line if line.endswith("\n") else line + "\n")

    # Find insertion point: before the first ATOM/HETATM record
    insert_pos = next(
        (i for i, l in enumerate(target_lines)
         if l.startswith("ATOM  ") or l.startswith("HETATM")),
        0,
    )

    output_lines = (
        target_lines[:insert_pos]
        + ss_lines
        + target_lines[insert_pos:]
    )

    Path(output_pdb).parent.mkdir(parents=True, exist_ok=True)
    with open(output_pdb, "w") as fh:
        fh.writelines(output_lines)

    return output_pdb


# ---------------------------------------------------------------------------
# Coot real-space refine
# ---------------------------------------------------------------------------

def _write_coot_rsr_script(
    pdb_path: str,
    map_path: str,
    windows: List[Tuple[str, int, int]],
    output_path: str,
) -> str:
    """Write a Coot Scheme script that refines every window in *windows*."""
    refine_calls = ""
    for chain_id, start_res, end_res in windows:
        # Escape chain id in case it contains a double-quote (shouldn't happen
        # with standard PDB chain IDs, but guard anyway)
        safe_chain = chain_id.replace('"', '\\"')
        refine_calls += _COOT_REFINE_ZONE_TMPL.format(
            chain_id=safe_chain,
            start_res=start_res,
            end_res=end_res,
        )

    script = _COOT_RSR_SCRIPT_TEMPLATE.format(
        pdb_path=os.path.abspath(pdb_path),
        map_path=os.path.abspath(map_path),
        refine_calls=refine_calls,
        output_path=os.path.abspath(output_path),
    )

    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".scm", delete=False, prefix="autostructure_rsr_"
    ) as tmp:
        tmp.write(script)
        return tmp.name


def run_coot_real_space_refine_windows(
    pdb_path: str,
    map_path: str,
    windows: List[Tuple[str, int, int]],
    output_path: str,
    coot_exe: str = "coot",
    log_path: Optional[str] = None,
) -> str:
    """
    Run Coot headless real-space refinement on the given *windows*.

    A single Coot session is launched that:
    1. Loads *pdb_path* and *map_path*.
    2. Enables secondary-structure restraints (type 2: helices + strands).
    3. Calls ``(refine-zone …)`` and ``(accept-regularization)`` for every
       window.
    4. Writes the refined model to *output_path*.

    Parameters
    ----------
    pdb_path    : current PDB file (SS records should be pre-injected)
    map_path    : EM density map
    windows     : list of (chain_id, start_resseq, end_resseq)
    output_path : destination for the refined PDB
    coot_exe    : path/name of the Coot executable
    log_path    : optional file to capture Coot stdout/stderr

    Returns
    -------
    str  Absolute path to *output_path*.

    Raises
    ------
    RuntimeError  if Coot exits with a non-zero return code.
    """
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)

    script_path = _write_coot_rsr_script(
        pdb_path=pdb_path,
        map_path=map_path,
        windows=windows,
        output_path=output_path,
    )

    cmd = [coot_exe, "--no-graphics", "--script", script_path]

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=False,
        )

        if log_path:
            Path(log_path).parent.mkdir(parents=True, exist_ok=True)
            with open(log_path, "w") as fh:
                fh.write("=== CMD ===\n")
                fh.write(" ".join(cmd) + "\n\n")
                fh.write("=== Coot Script ===\n")
                try:
                    fh.write(Path(script_path).read_text())
                except OSError:
                    pass
                fh.write("\n=== STDOUT ===\n")
                fh.write(result.stdout)
                fh.write("\n=== STDERR ===\n")
                fh.write(result.stderr)

        if result.returncode != 0:
            raise RuntimeError(
                f"Coot real-space refine exited with code {result.returncode}.\n"
                f"STDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
            )
    finally:
        try:
            os.unlink(script_path)
        except OSError:
            pass

    return os.path.abspath(output_path)


# ---------------------------------------------------------------------------
# Issue-set helpers
# ---------------------------------------------------------------------------

def _issue_key(issue: dict) -> tuple:
    """Return a hashable key that identifies a unique residue/type combination."""
    return (
        issue["chain"],
        issue["resseq"],
        issue["icode"],
        issue["resname"],
        issue["issue_type"],
    )


def _count_by_category(issues: List[dict]) -> Dict[str, int]:
    """Return per-category issue counts."""
    counts: Dict[str, int] = {c: 0 for c in _ALL_CATEGORIES}
    for iss in issues:
        counts[iss["issue_type"]] = counts.get(iss["issue_type"], 0) + 1
    return counts


# ---------------------------------------------------------------------------
# Main auto-fix loop
# ---------------------------------------------------------------------------

def run_auto_fix(
    pdb_path: str,
    map_path: str,
    af2_pdb: str,
    output_dir: str,
    phenix_molprobity_exe: str = "phenix.molprobity",
    coot_exe: str = "coot",
    max_iterations: int = 10,
    window_half: int = 5,
) -> str:
    """
    Iteratively validate and refine a PDB until all MolProbity issues are
    resolved or the maximum number of iterations is reached.

    Parameters
    ----------
    pdb_path              : starting PDB (output of the main pipeline)
    map_path              : EM density map used throughout
    af2_pdb               : original AlphaFold2 PDB (source of SS records)
    output_dir            : top-level directory for iteration subdirectories
    phenix_molprobity_exe : ``phenix.molprobity`` executable name/path
    coot_exe              : Coot executable name/path
    max_iterations        : maximum refinement rounds (default 10)
    window_half           : residues on each side of a problem site (default 5)

    Returns
    -------
    str  Path to the final (best available) fixed PDB.
    """
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    log.info("Auto-fix: checking dependencies …")
    check_dependencies(phenix_molprobity_exe, coot_exe)

    current_pdb = os.path.abspath(pdb_path)
    previous_issue_keys: Optional[set] = None

    for iteration in range(1, max_iterations + 1):
        iter_dir = os.path.join(output_dir, f"iteration_{iteration:02d}")
        Path(iter_dir).mkdir(parents=True, exist_ok=True)

        log.info("Auto-fix iteration %d/%d", iteration, max_iterations)
        log.info("  Input PDB: %s", current_pdb)

        # ---- Step 1: Inject AF2 secondary-structure records ----------------
        ss_injected_pdb = os.path.join(iter_dir, "input_with_ss.pdb")
        copy_ss_records(
            source_pdb=af2_pdb,
            target_pdb=current_pdb,
            output_pdb=ss_injected_pdb,
        )
        log.debug("  SS records injected → %s", ss_injected_pdb)

        # ---- Step 2: Run MolProbity validation -----------------------------
        validation_dir = os.path.join(iter_dir, "validation")
        validation_log = os.path.join(iter_dir, "molprobity.log")
        log.info("  Running phenix.molprobity …")
        stdout, stderr, returncode = run_molprobity(
            pdb_path=ss_injected_pdb,
            output_dir=validation_dir,
            phenix_molprobity_exe=phenix_molprobity_exe,
            log_path=validation_log,
        )

        if returncode != 0:
            log.warning(
                "  phenix.molprobity exited with code %d (non-fatal; continuing).",
                returncode,
            )

        # Save validation report
        report_path = os.path.join(iter_dir, "molprobity_report.txt")
        with open(report_path, "w") as fh:
            fh.write(stdout)
        log.info("  Validation report: %s", report_path)

        # ---- Step 3: Parse issues ------------------------------------------
        issues = parse_molprobity_output(stdout)
        current_keys = {_issue_key(i) for i in issues}
        counts = _count_by_category(issues)

        log.info(
            "  Issues found: %d  (rama=%d  rota=%d  cbeta=%d  cis=%d  geo=%d)",
            len(issues),
            counts[ISSUE_RAMACHANDRAN],
            counts[ISSUE_ROTAMER],
            counts[ISSUE_CBETA],
            counts[ISSUE_CIS_TWISTED],
            counts[ISSUE_GEOMETRY],
        )

        # Write per-iteration issue summary
        summary_path = os.path.join(iter_dir, "issues_summary.txt")
        with open(summary_path, "w") as fh:
            fh.write(f"Iteration {iteration}\n")
            fh.write(f"Total issues: {len(issues)}\n")
            for cat in _ALL_CATEGORIES:
                fh.write(f"  {cat}: {counts[cat]}\n")
            fh.write("\nDetails:\n")
            for iss in issues:
                fh.write(
                    f"  {iss['issue_type']:15s}  chain={iss['chain']}  "
                    f"resseq={iss['resseq']:5d}  resname={iss['resname']}  "
                    f"{iss['description']}\n"
                )

        # ---- Step 4: Convergence check -------------------------------------
        if len(issues) == 0:
            log.info("  No issues found – converged after %d iteration(s).", iteration)
            return os.path.abspath(current_pdb)

        if previous_issue_keys is not None:
            new_issues = current_keys - previous_issue_keys
            if not new_issues:
                log.info(
                    "  No *new* issues since previous round – "
                    "converged after %d iteration(s).",
                    iteration,
                )
                return os.path.abspath(current_pdb)

        previous_issue_keys = current_keys

        if iteration == max_iterations:
            log.warning(
                "  Reached maximum %d iterations; %d issue(s) remain.",
                max_iterations,
                len(issues),
            )
            return os.path.abspath(current_pdb)

        # ---- Step 5: Build and merge refinement windows --------------------
        chain_residues = collect_chain_residues(ss_injected_pdb)
        raw_windows    = build_refine_windows(issues, chain_residues, window_half)
        merged_windows = merge_windows(raw_windows)
        log.info(
            "  %d raw window(s) → %d merged window(s) to refine.",
            len(raw_windows),
            len(merged_windows),
        )

        if not merged_windows:
            log.warning("  No valid windows – skipping refinement.")
            return os.path.abspath(current_pdb)

        # ---- Step 6: Run Coot real-space refine ----------------------------
        refined_pdb = os.path.join(iter_dir, "refined.pdb")
        coot_log    = os.path.join(iter_dir, "coot_rsr.log")
        log.info("  Running Coot RSR on %d window(s) …", len(merged_windows))
        run_coot_real_space_refine_windows(
            pdb_path=ss_injected_pdb,
            map_path=map_path,
            windows=merged_windows,
            output_path=refined_pdb,
            coot_exe=coot_exe,
            log_path=coot_log,
        )

        log.info("  Refined PDB: %s", refined_pdb)
        current_pdb = refined_pdb

    # Should not reach here (the loop always returns inside)
    return os.path.abspath(current_pdb)
