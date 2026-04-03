"""
pipeline.py
-----------
AutoStructure main pipeline.

Workflow
--------
1. Global fitmap  – fit the full AlphaFold2 PDB into the EM map (ChimeraX)
2. SS parse       – extract secondary-structure annotation from the fitted PDB
3. Split          – divide the structure into 'structured' and 'loop' pieces
4. Per-piece fit  – fit each piece independently into the EM map (ChimeraX)
5. Refine         – rigid-body refine structured pieces (Phenix);
                    fit-protein (Coot) or morphing (Phenix) for loop pieces
6. Merge          – combine all refined pieces into a single output PDB

Usage (command line)
--------------------
    python pipeline.py --pdb model.pdb --map density.mrc [--config config.yaml]

Usage (Python API)
------------------
    from pipeline import run_pipeline
    run_pipeline(pdb_path="model.pdb", map_path="density.mrc")
"""

from __future__ import annotations

import argparse
import concurrent.futures
import logging
import os
from pathlib import Path
from typing import Optional

import yaml

from chimerax_runner import fit_piece, fit_whole_model
from coot_runner import fit_protein as coot_fit_protein
from phenix_runner import morphing_refine, rigid_body_refine
from piece_splitter import merge_pdbs, split_pdb
from secondary_structure import get_pieces

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger("autostructure")


# ---------------------------------------------------------------------------
# Config helpers
# ---------------------------------------------------------------------------

_DEFAULT_CONFIG = {
    "tools": {
        "chimerax": "chimerax",
        "phenix": "phenix.real_space_refine",
        "coot": "coot",
    },
    "fitmap": {
        "resolution": 3.0,
        "metric": "correlation",
        "max_steps": 2000,
    },
    "secondary_structure": {
        "min_helix_length": 4,
        "min_sheet_length": 3,
        "loop_min_length": 3,
    },
    "splitting": {
        "overlap_residues": 2,
    },
    "refinement": {
        "nproc": 4,
        "rigid_body_max_iterations": 50,
        "morphing_weight": 1.0,
        "coot_fit_protein": True,
    },
    "output": {
        "work_dir": "autostructure_output",
    },
}


def _load_config(config_path: Optional[str]) -> dict:
    cfg = dict(_DEFAULT_CONFIG)
    if config_path and os.path.isfile(config_path):
        with open(config_path) as fh:
            user = yaml.safe_load(fh) or {}
        for section, values in user.items():
            if isinstance(values, dict):
                cfg.setdefault(section, {})
                cfg[section].update(values)
            else:
                cfg[section] = values
    return cfg


# ---------------------------------------------------------------------------
# Per-piece refinement task (runs inside a process pool)
# ---------------------------------------------------------------------------

def _refine_piece(args: tuple) -> dict:
    """
    Fit and refine one piece.  Designed to run in a worker process.

    Parameters
    ----------
    args : (piece_dict, map_path, work_dir, cfg_dict)

    Returns
    -------
    Annotated piece dict with 'fitted_pdb_path' and 'refined_pdb_path' set.
    """
    piece, map_path, work_dir, cfg = args

    tools     = cfg["tools"]
    fitmap    = cfg["fitmap"]
    refine    = cfg["refinement"]

    tag = f"piece_{piece['chain']}_{piece['start']}_{piece['end']}"
    fit_dir    = os.path.join(work_dir, "fitted_pieces")
    refine_dir = os.path.join(work_dir, "refined_pieces", tag)
    log_dir    = os.path.join(work_dir, "logs")

    os.makedirs(fit_dir,    exist_ok=True)
    os.makedirs(refine_dir, exist_ok=True)
    os.makedirs(log_dir,    exist_ok=True)

    # ---- Step: fit piece into map ----
    fitted_path = os.path.join(fit_dir, f"{tag}_fitted.pdb")
    fit_piece(
        piece_pdb=piece["pdb_path"],
        map_path=map_path,
        output_path=fitted_path,
        resolution=fitmap["resolution"],
        chimerax_exe=tools["chimerax"],
        metric=fitmap["metric"],
        max_steps=fitmap["max_steps"],
        log_path=os.path.join(log_dir, f"{tag}_fitmap.log"),
    )
    piece["fitted_pdb_path"] = fitted_path

    # ---- Step: refine ----
    if piece["type"] == "structured":
        refined_path = rigid_body_refine(
            pdb_path=fitted_path,
            map_path=map_path,
            output_dir=refine_dir,
            resolution=fitmap["resolution"],
            phenix_exe=tools["phenix"],
            max_iterations=refine["rigid_body_max_iterations"],
            log_path=os.path.join(log_dir, f"{tag}_rigid_body.log"),
        )
    else:
        # Loop / disordered piece
        if refine["coot_fit_protein"]:
            refined_path = os.path.join(refine_dir, f"{tag}_coot_fitted.pdb")
            coot_fit_protein(
                pdb_path=fitted_path,
                map_path=map_path,
                output_path=refined_path,
                coot_exe=tools["coot"],
                log_path=os.path.join(log_dir, f"{tag}_coot.log"),
            )
        else:
            refined_path = morphing_refine(
                pdb_path=fitted_path,
                map_path=map_path,
                output_dir=refine_dir,
                resolution=fitmap["resolution"],
                phenix_exe=tools["phenix"],
                morphing_weight=refine["morphing_weight"],
                log_path=os.path.join(log_dir, f"{tag}_morphing.log"),
            )

    piece["refined_pdb_path"] = refined_path
    return piece


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def run_pipeline(
    pdb_path: str,
    map_path: str,
    config_path: Optional[str] = None,
    output_path: Optional[str] = None,
) -> str:
    """
    Run the full AutoStructure pipeline.

    Parameters
    ----------
    pdb_path    : AlphaFold2 (or similar) predicted PDB
    map_path    : experimental EM density map (MRC/CCP4)
    config_path : path to a YAML config file (optional; falls back to defaults)
    output_path : path for the final merged PDB  (optional; auto-generated if omitted)

    Returns
    -------
    str  Path to the final refined and merged PDB.
    """
    cfg      = _load_config(config_path)
    work_dir = cfg["output"]["work_dir"]
    tools    = cfg["tools"]
    fitmap   = cfg["fitmap"]
    ss_cfg   = cfg["secondary_structure"]
    split    = cfg["splitting"]
    nproc    = cfg["refinement"]["nproc"]

    Path(work_dir).mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # Step 1 – Global fitmap
    # ------------------------------------------------------------------
    log.info("Step 1: Global fitmap – fitting %s into %s", pdb_path, map_path)
    globally_fitted_pdb = os.path.join(work_dir, "globally_fitted.pdb")
    fit_whole_model(
        pdb_path=pdb_path,
        map_path=map_path,
        output_path=globally_fitted_pdb,
        resolution=fitmap["resolution"],
        chimerax_exe=tools["chimerax"],
        metric=fitmap["metric"],
        max_steps=fitmap["max_steps"],
        log_path=os.path.join(work_dir, "logs", "global_fitmap.log"),
    )
    log.info("  -> %s", globally_fitted_pdb)

    # ------------------------------------------------------------------
    # Step 2 – Parse secondary structure & split into pieces
    # ------------------------------------------------------------------
    log.info("Step 2: Parsing secondary structure from %s", globally_fitted_pdb)
    pieces = get_pieces(
        globally_fitted_pdb,
        min_helix=ss_cfg["min_helix_length"],
        min_sheet=ss_cfg["min_sheet_length"],
        loop_min_length=ss_cfg["loop_min_length"],
    )
    n_structured = sum(1 for p in pieces if p["type"] == "structured")
    n_loop       = len(pieces) - n_structured
    log.info("  -> %d pieces total (%d structured, %d loop)", len(pieces), n_structured, n_loop)

    pieces_dir = os.path.join(work_dir, "pieces")
    pieces = split_pdb(
        globally_fitted_pdb,
        pieces,
        pieces_dir,
        overlap=split["overlap_residues"],
    )
    log.info("  -> pieces written to %s", pieces_dir)

    # ------------------------------------------------------------------
    # Step 3-5 – Per-piece fitmap + refinement (parallelised)
    # ------------------------------------------------------------------
    log.info("Step 3-5: Per-piece fitmap and refinement (%d workers)", nproc)

    worker_args = [
        (piece, map_path, work_dir, cfg)
        for piece in pieces
    ]

    refined_pieces = []
    if nproc > 1:
        with concurrent.futures.ProcessPoolExecutor(max_workers=nproc) as pool:
            futures = {pool.submit(_refine_piece, a): a[0] for a in worker_args}
            for future in concurrent.futures.as_completed(futures):
                original = futures[future]
                tag = f"chain {original['chain']} {original['start']}-{original['end']}"
                try:
                    result = future.result()
                    log.info("  [OK]  %s  -> %s", tag, result.get("refined_pdb_path"))
                    refined_pieces.append(result)
                except Exception as exc:
                    log.error("  [FAIL] %s: %s", tag, exc)
                    raise
    else:
        for args in worker_args:
            result = _refine_piece(args)
            tag = f"chain {result['chain']} {result['start']}-{result['end']}"
            log.info("  [OK]  %s  -> %s", tag, result.get("refined_pdb_path"))
            refined_pieces.append(result)

    # Restore original ordering (ProcessPoolExecutor may reorder)
    refined_pieces.sort(key=lambda p: (p["chain"], p["start"]))

    # ------------------------------------------------------------------
    # Step 6 – Merge all pieces into the final PDB
    # ------------------------------------------------------------------
    log.info("Step 6: Merging %d pieces", len(refined_pieces))
    if output_path is None:
        stem        = Path(pdb_path).stem
        output_path = os.path.join(work_dir, f"{stem}_autostructure_refined.pdb")

    final_pdb = merge_pdbs(refined_pieces, output_path)
    log.info("Pipeline complete.  Final PDB: %s", final_pdb)
    return final_pdb


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="AutoStructure: automatically refine protein structures "
                    "using ChimeraX, Phenix, and Coot.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--pdb",    required=True,  help="Input PDB file (AlphaFold2 prediction)")
    p.add_argument("--map",    required=True,  help="EM density map (MRC/CCP4)")
    p.add_argument("--config", default=None,   help="Path to config.yaml")
    p.add_argument("--output", default=None,   help="Output PDB path")
    return p


if __name__ == "__main__":
    args = _build_parser().parse_args()
    final = run_pipeline(
        pdb_path=args.pdb,
        map_path=args.map,
        config_path=args.config,
        output_path=args.output,
    )
    print(f"Done: {final}")
