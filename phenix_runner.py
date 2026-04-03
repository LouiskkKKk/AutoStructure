"""
phenix_runner.py
----------------
Wrappers for running Phenix real-space refinement in two modes:

  • rigid_body_refine  – for structured pieces (helices, sheets)
  • morphing_refine    – for disordered/loop pieces

Both functions call ``phenix.real_space_refine`` via subprocess and return
the path to the refined PDB file.
"""

from __future__ import annotations

import glob as _glob
import os
import subprocess
from pathlib import Path
from typing import Optional


# ---------------------------------------------------------------------------
# Internal helper
# ---------------------------------------------------------------------------

def _find_output_pdb(output_dir: str, stem: str) -> Optional[str]:
    """
    Locate the refined PDB written by Phenix.

    Phenix appends ``_real_space_refined.pdb`` or ``_real_space_refined_XXX.pdb``
    to the input file stem.
    """
    patterns = [
        os.path.join(output_dir, f"{stem}_real_space_refined.pdb"),
        os.path.join(output_dir, f"{stem}_real_space_refined_*.pdb"),
        os.path.join(output_dir, "*.pdb"),
    ]
    for pat in patterns:
        found = sorted(_glob.glob(pat))
        if found:
            # Return the most recently modified file
            return max(found, key=os.path.getmtime)
    return None


def _run_phenix_real_space_refine(
    pdb_path: str,
    map_path: str,
    output_dir: str,
    resolution: float,
    extra_args: list,
    phenix_exe: str,
    log_path: Optional[str],
) -> str:
    """
    Core helper that builds and executes the phenix.real_space_refine command.

    Returns
    -------
    str  Path to the refined PDB file.

    Raises
    ------
    RuntimeError  if Phenix exits non-zero or the output PDB cannot be found.
    """
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    stem = Path(pdb_path).stem

    cmd = [
        phenix_exe,
        os.path.abspath(pdb_path),
        os.path.abspath(map_path),
        f"resolution={resolution}",
        f"output.file_name_prefix={stem}",
        f"output.directory={os.path.abspath(output_dir)}",
    ] + extra_args

    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        check=False,
        cwd=output_dir,
    )

    if log_path:
        Path(log_path).parent.mkdir(parents=True, exist_ok=True)
        with open(log_path, "w") as lf:
            lf.write("=== CMD ===\n")
            lf.write(" ".join(cmd) + "\n\n")
            lf.write("=== STDOUT ===\n")
            lf.write(result.stdout)
            lf.write("\n=== STDERR ===\n")
            lf.write(result.stderr)

    if result.returncode != 0:
        raise RuntimeError(
            f"phenix.real_space_refine exited with code {result.returncode}.\n"
            f"STDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
        )

    refined_pdb = _find_output_pdb(output_dir, stem)
    if refined_pdb is None:
        raise RuntimeError(
            f"Phenix finished but no output PDB was found in {output_dir}"
        )

    return os.path.abspath(refined_pdb)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def rigid_body_refine(
    pdb_path: str,
    map_path: str,
    output_dir: str,
    resolution: float = 3.0,
    phenix_exe: str = "phenix.real_space_refine",
    max_iterations: int = 50,
    log_path: Optional[str] = None,
) -> str:
    """
    Run Phenix real-space rigid-body refinement on a structured piece.

    Parameters
    ----------
    pdb_path       : input piece PDB (after ChimeraX fitmap)
    map_path       : EM density map
    output_dir     : directory for Phenix output files
    resolution     : map resolution in Å
    phenix_exe     : path/name of phenix.real_space_refine executable
    max_iterations : maximum macro-cycles
    log_path       : optional path to save Phenix stdout/stderr

    Returns
    -------
    str  Path to the rigid-body refined PDB.
    """
    extra_args = [
        "rigid_body=True",
        "minimization_global=False",
        "adp=False",
        f"macro_cycles={max_iterations}",
    ]
    return _run_phenix_real_space_refine(
        pdb_path=pdb_path,
        map_path=map_path,
        output_dir=output_dir,
        resolution=resolution,
        extra_args=extra_args,
        phenix_exe=phenix_exe,
        log_path=log_path,
    )


def morphing_refine(
    pdb_path: str,
    map_path: str,
    output_dir: str,
    resolution: float = 3.0,
    phenix_exe: str = "phenix.real_space_refine",
    morphing_weight: float = 1.0,
    log_path: Optional[str] = None,
) -> str:
    """
    Run Phenix real-space morphing refinement on a disordered/loop piece.

    Parameters
    ----------
    pdb_path         : input piece PDB (after ChimeraX fitmap)
    map_path         : EM density map
    output_dir       : directory for Phenix output files
    resolution       : map resolution in Å
    phenix_exe       : path/name of phenix.real_space_refine executable
    morphing_weight  : weight for morphing restraints (higher → more conformational
                       change allowed)
    log_path         : optional path to save Phenix stdout/stderr

    Returns
    -------
    str  Path to the morphing-refined PDB.
    """
    extra_args = [
        "morphing=True",
        "rigid_body=False",
        "minimization_global=True",
        "adp=True",
        f"morphing_weight={morphing_weight}",
    ]
    return _run_phenix_real_space_refine(
        pdb_path=pdb_path,
        map_path=map_path,
        output_dir=output_dir,
        resolution=resolution,
        extra_args=extra_args,
        phenix_exe=phenix_exe,
        log_path=log_path,
    )
