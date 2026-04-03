"""
secondary_structure.py
----------------------
Parse secondary-structure annotations from a PDB file and partition the chain
into contiguous pieces labelled as "structured" (helix / sheet) or "loop".

Public API
----------
parse_secondary_structure(pdb_path) -> dict
    Returns a dict: {chain_id: {resseq: 'H'|'E'|'L'}}

get_pieces(pdb_path, min_helix=4, min_sheet=3, loop_min_length=3)
    -> list[dict]
    Each element: {
        'type'  : 'structured' | 'loop',
        'chain' : str,
        'start' : int,   # first residue sequence number
        'end'   : int,   # last residue sequence number (inclusive)
        'residues': list[int]
    }
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Dict, List, Tuple

try:
    from Bio import PDB as biopdb
    from Bio.PDB import DSSP
    _BIOPYTHON = True
except ImportError:
    _BIOPYTHON = False


# ---------------------------------------------------------------------------
# Low-level PDB record parsers
# ---------------------------------------------------------------------------

def _parse_helix_records(pdb_path: str) -> List[Tuple[str, int, int]]:
    """Return list of (chain_id, start_resseq, end_resseq) from HELIX records."""
    results = []
    with open(pdb_path) as fh:
        for line in fh:
            if not line.startswith("HELIX "):
                continue
            chain = line[19]
            try:
                start = int(line[21:25].strip())
                end   = int(line[33:37].strip())
            except ValueError:
                continue
            results.append((chain, start, end))
    return results


def _parse_sheet_records(pdb_path: str) -> List[Tuple[str, int, int]]:
    """Return list of (chain_id, start_resseq, end_resseq) from SHEET records."""
    results = []
    with open(pdb_path) as fh:
        for line in fh:
            if not line.startswith("SHEET "):
                continue
            chain = line[21]
            try:
                start = int(line[22:26].strip())
                end   = int(line[33:37].strip())
            except ValueError:
                continue
            results.append((chain, start, end))
    return results


def _collect_residues_from_pdb(pdb_path: str) -> Dict[str, List[int]]:
    """Return {chain_id: sorted list of residue sequence numbers} from ATOM records."""
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
# Public: parse secondary structure
# ---------------------------------------------------------------------------

def parse_secondary_structure(pdb_path: str) -> Dict[str, Dict[int, str]]:
    """
    Parse secondary structure annotations from a PDB file.

    Returns
    -------
    dict  {chain_id: {resseq: 'H' (helix) | 'E' (sheet/strand) | 'L' (loop)}}

    Strategy
    --------
    1. Try HELIX/SHEET records embedded in the PDB (written by AlphaFold2, Rosetta, etc.)
    2. Fall back to BioPython DSSP if the records are absent and DSSP is installed.
    3. If neither is available, all residues are labelled 'L'.
    """
    pdb_path = str(pdb_path)
    all_residues = _collect_residues_from_pdb(pdb_path)

    # Initialise everything as loop
    ss: Dict[str, Dict[int, str]] = {}
    for chain, residues in all_residues.items():
        ss[chain] = {r: "L" for r in residues}

    helix_records = _parse_helix_records(pdb_path)
    sheet_records = _parse_sheet_records(pdb_path)

    has_records = bool(helix_records or sheet_records)

    if has_records:
        for chain, start, end in helix_records:
            if chain not in ss:
                continue
            for r in range(start, end + 1):
                if r in ss[chain]:
                    ss[chain][r] = "H"

        for chain, start, end in sheet_records:
            if chain not in ss:
                continue
            for r in range(start, end + 1):
                if r in ss[chain]:
                    ss[chain][r] = "E"

        return ss

    # ---- Fall back to DSSP ----
    if _BIOPYTHON:
        try:
            parser = biopdb.PDBParser(QUIET=True)
            structure = parser.get_structure("tmp", pdb_path)
            model = next(structure.get_models())
            dssp = DSSP(model, pdb_path)
            _dssp_map = {"H": "H", "G": "H", "I": "H",
                         "B": "E", "E": "E",
                         "T": "L", "S": "L", "-": "L", "C": "L"}
            for key, value in dssp.property_dict.items():
                chain_id = key[0]
                resseq   = key[1][1]
                raw_ss   = value[2]
                mapped   = _dssp_map.get(raw_ss, "L")
                if chain_id in ss and resseq in ss[chain_id]:
                    ss[chain_id][resseq] = mapped
            return ss
        except Exception:
            pass  # DSSP failed; return all-loop labels

    return ss


# ---------------------------------------------------------------------------
# Public: get_pieces
# ---------------------------------------------------------------------------

def get_pieces(
    pdb_path: str,
    min_helix: int = 4,
    min_sheet: int = 3,
    loop_min_length: int = 3,
) -> List[dict]:
    """
    Partition every chain in *pdb_path* into contiguous pieces.

    Parameters
    ----------
    pdb_path        : path to the PDB file
    min_helix       : minimum helix residues to classify a run as 'structured'
    min_sheet       : minimum sheet residues to classify a run as 'structured'
    loop_min_length : minimum residues for a loop to become a standalone piece
                      (shorter loops are merged with adjacent structured pieces)

    Returns
    -------
    list of dicts, each containing:
        type      : 'structured' | 'loop'
        chain     : chain ID (str)
        start     : first residue seqnum (int)
        end       : last residue seqnum (int, inclusive)
        residues  : sorted list of residue seqnums
    """
    ss_map = parse_secondary_structure(pdb_path)
    all_residues = _collect_residues_from_pdb(pdb_path)
    pieces: List[dict] = []

    for chain in sorted(ss_map.keys()):
        residues = all_residues.get(chain, [])
        if not residues:
            continue

        # Build a run-length encoded list of (ss_label, residue_list)
        runs: List[Tuple[str, List[int]]] = []
        current_label = ss_map[chain].get(residues[0], "L")
        current_run: List[int] = [residues[0]]

        for r in residues[1:]:
            label = ss_map[chain].get(r, "L")
            if label == current_label:
                current_run.append(r)
            else:
                runs.append((current_label, current_run))
                current_label = label
                current_run = [r]
        runs.append((current_label, current_run))

        # Apply minimum-length thresholds
        for label, run in runs:
            is_structured = False
            if label == "H" and len(run) >= min_helix:
                is_structured = True
            elif label == "E" and len(run) >= min_sheet:
                is_structured = True

            if not is_structured and len(run) < loop_min_length:
                # Merge tiny loops into the previous piece if one exists
                if pieces and pieces[-1]["chain"] == chain:
                    pieces[-1]["residues"] = sorted(pieces[-1]["residues"] + run)
                    pieces[-1]["end"] = pieces[-1]["residues"][-1]
                    continue

            pieces.append({
                "type"    : "structured" if is_structured else "loop",
                "chain"   : chain,
                "start"   : run[0],
                "end"     : run[-1],
                "residues": list(run),
            })

    return pieces


# ---------------------------------------------------------------------------
# CLI helper
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import sys
    import json

    if len(sys.argv) < 2:
        print("Usage: python secondary_structure.py <pdb_file>")
        sys.exit(1)

    pieces = get_pieces(sys.argv[1])
    for i, p in enumerate(pieces):
        print(f"Piece {i:3d} | chain={p['chain']} | {p['type']:10s} | "
              f"res {p['start']:5d}-{p['end']:5d} | length={len(p['residues'])}")
