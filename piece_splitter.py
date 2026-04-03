"""
piece_splitter.py
-----------------
Split a PDB file into sub-PDB files (pieces) by residue range, and later
merge refined pieces back into a single combined PDB.

Public API
----------
split_pdb(pdb_path, pieces, output_dir, overlap=2)
    -> list[dict]  (pieces annotated with 'pdb_path' key)

merge_pdbs(pieces, output_path)
    -> str  (path to the merged PDB)
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Dict, List, Optional

try:
    from Bio import PDB as biopdb
    from Bio.PDB import PDBIO, Select
    _BIOPYTHON = True
except ImportError:
    _BIOPYTHON = False

    class Select:  # noqa: N801 – dummy base so class body is always valid
        """Fallback stub used when BioPython is not installed."""


# ---------------------------------------------------------------------------
# BioPython-based implementation
# ---------------------------------------------------------------------------

class _ResidueRangeSelect(Select):
    """BioPython Select subclass that keeps only a specified residue range."""

    def __init__(self, chain_id: str, residue_set: set):
        self.chain_id    = chain_id
        self.residue_set = residue_set

    def accept_chain(self, chain):
        return chain.get_id() == self.chain_id

    def accept_residue(self, residue):
        return residue.get_id()[1] in self.residue_set


def _expand_with_overlap(residues: List[int], all_residues: List[int], overlap: int) -> List[int]:
    """Expand a residue list by *overlap* neighbours on each side."""
    if overlap <= 0 or not residues:
        return list(residues)
    full_set = set(all_residues)
    expanded = set(residues)
    sorted_all = sorted(all_residues)
    idx_map = {r: i for i, r in enumerate(sorted_all)}

    for r in list(residues):
        if r not in idx_map:
            continue
        pos = idx_map[r]
        for delta in range(1, overlap + 1):
            lo = pos - delta
            hi = pos + delta
            if lo >= 0:
                expanded.add(sorted_all[lo])
            if hi < len(sorted_all):
                expanded.add(sorted_all[hi])

    return sorted(expanded)


def _split_biopython(
    pdb_path: str,
    pieces: List[dict],
    output_dir: str,
    overlap: int,
) -> List[dict]:
    parser = biopdb.PDBParser(QUIET=True)
    structure = parser.get_structure("model", pdb_path)
    io = PDBIO()
    io.set_structure(structure)

    # Build per-chain all-residue list for overlap expansion
    chain_residues: Dict[str, List[int]] = {}
    for model in structure:
        for chain in model:
            cid = chain.get_id()
            chain_residues[cid] = sorted({r.get_id()[1] for r in chain.get_residues()})
        break  # use first model only

    annotated: List[dict] = []
    for idx, piece in enumerate(pieces):
        chain   = piece["chain"]
        residues = piece["residues"]
        all_res  = chain_residues.get(chain, residues)

        expanded = _expand_with_overlap(residues, all_res, overlap)
        piece_path = os.path.join(output_dir, f"piece_{idx:04d}_{piece['type']}_chain{chain}.pdb")
        sel = _ResidueRangeSelect(chain, set(expanded))
        io.save(piece_path, sel)

        annotated.append({**piece, "pdb_path": piece_path, "expanded_residues": expanded})

    return annotated


# ---------------------------------------------------------------------------
# Plain-text fallback (no BioPython)
# ---------------------------------------------------------------------------

def _split_plaintext(
    pdb_path: str,
    pieces: List[dict],
    output_dir: str,
    overlap: int,
) -> List[dict]:
    """Extract ATOM/HETATM lines by residue range using plain string parsing."""

    # Read all ATOM/HETATM lines grouped by (chain, resseq)
    atom_lines: Dict[tuple, List[str]] = {}
    with open(pdb_path) as fh:
        for line in fh:
            if not (line.startswith("ATOM  ") or line.startswith("HETATM")):
                continue
            chain  = line[21]
            try:
                resseq = int(line[22:26].strip())
            except ValueError:
                continue
            atom_lines.setdefault((chain, resseq), []).append(line)

    chain_residues: Dict[str, List[int]] = {}
    for chain, resseq in atom_lines:
        chain_residues.setdefault(chain, []).append(resseq)
    chain_residues = {c: sorted(set(r)) for c, r in chain_residues.items()}

    annotated: List[dict] = []
    for idx, piece in enumerate(pieces):
        chain    = piece["chain"]
        residues = piece["residues"]
        all_res  = chain_residues.get(chain, residues)

        expanded = _expand_with_overlap(residues, all_res, overlap)
        piece_path = os.path.join(output_dir, f"piece_{idx:04d}_{piece['type']}_chain{chain}.pdb")

        with open(piece_path, "w") as out:
            for r in expanded:
                for line in atom_lines.get((chain, r), []):
                    out.write(line)
            out.write("END\n")

        annotated.append({**piece, "pdb_path": piece_path, "expanded_residues": expanded})

    return annotated


# ---------------------------------------------------------------------------
# Public: split_pdb
# ---------------------------------------------------------------------------

def split_pdb(
    pdb_path: str,
    pieces: List[dict],
    output_dir: str,
    overlap: int = 2,
) -> List[dict]:
    """
    Write each piece in *pieces* to a separate PDB file under *output_dir*.

    Parameters
    ----------
    pdb_path   : source PDB file (full model)
    pieces     : list of piece dicts as returned by secondary_structure.get_pieces()
    output_dir : directory where piece PDB files are written (created if absent)
    overlap    : number of residues added to each side of a piece for context

    Returns
    -------
    The same list of dicts, each extended with:
        pdb_path          : path to the written piece PDB file
        expanded_residues : residue seqnums actually written (incl. overlap)
    """
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    if _BIOPYTHON:
        return _split_biopython(pdb_path, pieces, output_dir, overlap)
    return _split_plaintext(pdb_path, pieces, output_dir, overlap)


# ---------------------------------------------------------------------------
# Public: merge_pdbs
# ---------------------------------------------------------------------------

def merge_pdbs(pieces: List[dict], output_path: str) -> str:
    """
    Merge refined piece PDB files back into a single PDB.

    Only the *core* residues (piece["residues"]) are taken from each piece's
    refined PDB; overlap residues are discarded.  The result is written to
    *output_path*.

    Parameters
    ----------
    pieces      : annotated piece dicts (must have 'pdb_path' and 'residues')
    output_path : destination file path for the merged PDB

    Returns
    -------
    str  The *output_path* argument.
    """
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w") as out:
        for piece in pieces:
            core_set = set(piece["residues"])
            chain    = piece["chain"]
            piece_pdb = piece.get("refined_pdb_path") or piece.get("pdb_path")

            if not piece_pdb or not os.path.isfile(piece_pdb):
                print(f"[merge_pdbs] WARNING: piece file not found: {piece_pdb}")
                continue

            with open(piece_pdb) as fh:
                for line in fh:
                    if not (line.startswith("ATOM  ") or line.startswith("HETATM")):
                        continue
                    if line[21] != chain:
                        continue
                    try:
                        resseq = int(line[22:26].strip())
                    except ValueError:
                        continue
                    if resseq in core_set:
                        out.write(line)

        out.write("END\n")

    return output_path


# ---------------------------------------------------------------------------
# CLI helper
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import sys
    from secondary_structure import get_pieces

    if len(sys.argv) < 3:
        print("Usage: python piece_splitter.py <pdb_file> <output_dir> [overlap]")
        sys.exit(1)

    pdb  = sys.argv[1]
    odir = sys.argv[2]
    ovlp = int(sys.argv[3]) if len(sys.argv) > 3 else 2

    ps = get_pieces(pdb)
    annotated = split_pdb(pdb, ps, odir, overlap=ovlp)
    for p in annotated:
        print(f"  -> {p['pdb_path']}  ({p['type']}, chain {p['chain']}, "
              f"res {p['start']}-{p['end']})")
