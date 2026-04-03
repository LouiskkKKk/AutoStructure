"""
chimerax_runner.py
------------------
Wrappers for calling ChimeraX in headless (--nogui) mode to perform
fit-in-map (fitmap) alignment.

Public API
----------
fit_whole_model(pdb_path, map_path, output_path, resolution, chimerax_exe, **kwargs)
    -> str  (path to the fitted PDB)

fit_piece(piece_pdb, map_path, output_path, resolution, chimerax_exe, **kwargs)
    -> str  (path to the fitted PDB)

Both functions write a temporary ChimeraX command script, execute ChimeraX,
and return the path to the saved output PDB.
"""

from __future__ import annotations

import os
import subprocess
import tempfile
from pathlib import Path
from typing import Optional


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

_CHIMERAX_SCRIPT_TEMPLATE = """\
# Auto-generated ChimeraX script – do not edit manually
open {pdb_path}
open {map_path}
fitmap #{model_id} inMap #{map_id} metric {metric} resolution {resolution} maxSteps {max_steps}
save {output_path} #{model_id}
exit
"""


def _build_script(
    pdb_path: str,
    map_path: str,
    output_path: str,
    resolution: float = 3.0,
    metric: str = "correlation",
    max_steps: int = 2000,
    model_id: int = 1,
    map_id: int = 2,
) -> str:
    """Return a ChimeraX command script as a string."""
    return _CHIMERAX_SCRIPT_TEMPLATE.format(
        pdb_path=os.path.abspath(pdb_path),
        map_path=os.path.abspath(map_path),
        output_path=os.path.abspath(output_path),
        resolution=resolution,
        metric=metric,
        max_steps=max_steps,
        model_id=model_id,
        map_id=map_id,
    )


def _run_chimerax(script: str, chimerax_exe: str, log_path: Optional[str] = None) -> None:
    """Write *script* to a temp file and execute ChimeraX in headless mode."""
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".cxc", delete=False, prefix="autostructure_"
    ) as tmp:
        tmp.write(script)
        script_path = tmp.name

    cmd = [chimerax_exe, "--nogui", "--exit", "--script", script_path]

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
                lf.write("=== STDOUT ===\n")
                lf.write(result.stdout)
                lf.write("\n=== STDERR ===\n")
                lf.write(result.stderr)

        if result.returncode != 0:
            raise RuntimeError(
                f"ChimeraX exited with code {result.returncode}.\n"
                f"STDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
            )
    finally:
        try:
            os.unlink(script_path)
        except OSError:
            pass


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def fit_whole_model(
    pdb_path: str,
    map_path: str,
    output_path: str,
    resolution: float = 3.0,
    chimerax_exe: str = "chimerax",
    metric: str = "correlation",
    max_steps: int = 2000,
    log_path: Optional[str] = None,
) -> str:
    """
    Fit the entire predicted PDB structure into the EM map using ChimeraX fitmap.

    Parameters
    ----------
    pdb_path     : path to the input PDB (e.g. AlphaFold2 prediction)
    map_path     : path to the EM density map (MRC/CCP4/etc.)
    output_path  : where to write the fitted PDB
    resolution   : nominal map resolution (Å) — informational for ChimeraX
    chimerax_exe : path/name of ChimeraX executable
    metric       : fitmap optimisation metric ('correlation' | 'cam')
    max_steps    : maximum optimisation iterations
    log_path     : optional path to save ChimeraX stdout/stderr

    Returns
    -------
    str  Absolute path to the fitted PDB file.
    """
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    script = _build_script(
        pdb_path=pdb_path,
        map_path=map_path,
        output_path=output_path,
        resolution=resolution,
        metric=metric,
        max_steps=max_steps,
    )
    _run_chimerax(script, chimerax_exe, log_path=log_path)
    return os.path.abspath(output_path)


def fit_piece(
    piece_pdb: str,
    map_path: str,
    output_path: str,
    resolution: float = 3.0,
    chimerax_exe: str = "chimerax",
    metric: str = "correlation",
    max_steps: int = 2000,
    log_path: Optional[str] = None,
) -> str:
    """
    Fit a single piece PDB into the EM map using ChimeraX fitmap.

    Parameters
    ----------
    piece_pdb    : path to the piece PDB file
    map_path     : path to the EM density map
    output_path  : where to write the fitted piece PDB
    resolution   : nominal map resolution (Å)
    chimerax_exe : path/name of ChimeraX executable
    metric       : fitmap optimisation metric
    max_steps    : maximum optimisation iterations
    log_path     : optional path to save ChimeraX stdout/stderr

    Returns
    -------
    str  Absolute path to the fitted piece PDB file.
    """
    return fit_whole_model(
        pdb_path=piece_pdb,
        map_path=map_path,
        output_path=output_path,
        resolution=resolution,
        chimerax_exe=chimerax_exe,
        metric=metric,
        max_steps=max_steps,
        log_path=log_path,
    )
