"""
coot_runner.py
--------------
Wrapper for running Coot in headless (--script) mode to perform
fit-protein refinement on disordered/loop pieces.

Public API
----------
fit_protein(pdb_path, map_path, output_path, coot_exe, log_path)
    -> str  (path to the fitted PDB)
"""

from __future__ import annotations

import os
import subprocess
import tempfile
from pathlib import Path
from typing import Optional


# ---------------------------------------------------------------------------
# Coot script template
# ---------------------------------------------------------------------------

_COOT_SCRIPT_TEMPLATE = """\
; Auto-generated Coot script – do not edit manually
(let* ((mol  (read-pdb "{pdb_path}"))
       (map  (read-ccp4-map "{map_path}" 0)))
  (set-imol-refinement-map map)
  (fit-protein mol)
  (write-pdb-file mol "{output_path}")
  (coot-real-exit 0))
"""


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _write_coot_script(pdb_path: str, map_path: str, output_path: str) -> str:
    """Write a Coot Scheme script to a temp file and return its path."""
    script = _COOT_SCRIPT_TEMPLATE.format(
        pdb_path=os.path.abspath(pdb_path),
        map_path=os.path.abspath(map_path),
        output_path=os.path.abspath(output_path),
    )
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".scm", delete=False, prefix="autostructure_coot_"
    ) as tmp:
        tmp.write(script)
        return tmp.name


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def fit_protein(
    pdb_path: str,
    map_path: str,
    output_path: str,
    coot_exe: str = "coot",
    log_path: Optional[str] = None,
) -> str:
    """
    Run Coot's fit-protein function on a disordered piece in headless mode.

    Coot is invoked as::

        coot --no-graphics --script <script.scm>

    The script reads the PDB and map, calls ``(fit-protein mol)``, and saves
    the result to *output_path*.

    Parameters
    ----------
    pdb_path    : input piece PDB (after ChimeraX fitmap)
    map_path    : EM density map (CCP4/MRC format)
    output_path : destination for the fitted PDB
    coot_exe    : path/name of the Coot executable
    log_path    : optional path to save Coot stdout/stderr

    Returns
    -------
    str  Absolute path to the fitted piece PDB.

    Raises
    ------
    RuntimeError  if Coot exits with a non-zero return code.
    """
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)

    script_path = _write_coot_script(
        pdb_path=pdb_path,
        map_path=map_path,
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
            with open(log_path, "w") as lf:
                lf.write("=== CMD ===\n")
                lf.write(" ".join(cmd) + "\n\n")
                lf.write("=== STDOUT ===\n")
                lf.write(result.stdout)
                lf.write("\n=== STDERR ===\n")
                lf.write(result.stderr)

        if result.returncode != 0:
            raise RuntimeError(
                f"Coot exited with code {result.returncode}.\n"
                f"STDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
            )
    finally:
        try:
            os.unlink(script_path)
        except OSError:
            pass

    return os.path.abspath(output_path)
