"""
tests/test_validation_refine.py
--------------------------------
Unit tests for validation_refine.py.

All tests use mock objects / in-memory data so that neither Phenix nor Coot
needs to be installed in the test environment.
"""

from __future__ import annotations

import os
import sys
import textwrap
import tempfile
from pathlib import Path
from unittest import mock

import pytest

# Make sure the project root is on sys.path
sys.path.insert(0, str(Path(__file__).parent.parent))

from validation_refine import (
    ISSUE_CBETA,
    ISSUE_CIS_TWISTED,
    ISSUE_GEOMETRY,
    ISSUE_RAMACHANDRAN,
    ISSUE_ROTAMER,
    build_refine_windows,
    check_dependencies,
    collect_chain_residues,
    copy_ss_records,
    merge_windows,
    parse_molprobity_output,
    run_coot_real_space_refine_windows,
    run_molprobity,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _write_tmp_pdb(content: str) -> str:
    """Write *content* to a temp PDB file and return its path."""
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".pdb", delete=False, prefix="test_vr_"
    ) as f:
        f.write(content)
        return f.name


_SIMPLE_PDB = textwrap.dedent("""\
    HELIX    1   1 ALA A    1  ALA A    5  1                                   5
    SHEET    1   A 2 VAL A   10  PHE A  12  0
    ATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00 10.00           N
    ATOM      2  CA  ALA A   1       2.000   3.000   4.000  1.00 10.00           C
    ATOM      3  N   GLY A   2       3.000   4.000   5.000  1.00 10.00           N
    ATOM      4  CA  GLY A   2       4.000   5.000   6.000  1.00 10.00           C
    ATOM      5  N   VAL A   3       5.000   6.000   7.000  1.00 10.00           N
    ATOM      6  CA  VAL A   4       6.000   7.000   8.000  1.00 10.00           C
    ATOM      7  N   LEU A   5       7.000   8.000   9.000  1.00 10.00           N
    ATOM      8  CA  LEU B  10       8.000   9.000  10.000  1.00 10.00           C
    ATOM      9  CA  PHE B  11       9.000  10.000  11.000  1.00 10.00           C
    ATOM     10  CA  TRP B  12      10.000  11.000  12.000  1.00 10.00           C
    END
""")


# ---------------------------------------------------------------------------
# parse_molprobity_output
# ---------------------------------------------------------------------------

class TestParseMolprobityOutput:

    SAMPLE_OUTPUT = textwrap.dedent("""\
        MolProbity Analysis
        ===================

        Ramachandran not favored:
          chain  resseq  icode  resname  phi    psi
          A      10             MET      -60.0  130.0
          A      20             GLY      -50.0  120.0

        Rotamer outliers:
          chain  resseq  icode  resname  rotamer
          A      15             LEU      ...

        C-beta deviations (>0.25 A):
          chain  resseq  icode  resname  deviation
          A      25             PHE      0.35

        Cis and twisted peptides:
          chain  resseq  icode  resname  omega  type
          A      30             GLY      5.0    CIS

        Geometry restraints:
          chain  resseq  icode  resname  type     deviation
          B      50             LYS      bond     0.08
    """)

    def test_ramachandran_parsed(self):
        issues = parse_molprobity_output(self.SAMPLE_OUTPUT)
        rama = [i for i in issues if i["issue_type"] == ISSUE_RAMACHANDRAN]
        assert len(rama) == 2
        assert rama[0]["chain"] == "A"
        assert rama[0]["resseq"] == 10
        assert rama[0]["resname"] == "MET"

    def test_rotamer_parsed(self):
        issues = parse_molprobity_output(self.SAMPLE_OUTPUT)
        rota = [i for i in issues if i["issue_type"] == ISSUE_ROTAMER]
        assert len(rota) == 1
        assert rota[0]["resseq"] == 15

    def test_cbeta_parsed(self):
        issues = parse_molprobity_output(self.SAMPLE_OUTPUT)
        cb = [i for i in issues if i["issue_type"] == ISSUE_CBETA]
        assert len(cb) == 1
        assert cb[0]["resseq"] == 25
        assert cb[0]["resname"] == "PHE"

    def test_cis_twisted_parsed(self):
        issues = parse_molprobity_output(self.SAMPLE_OUTPUT)
        ct = [i for i in issues if i["issue_type"] == ISSUE_CIS_TWISTED]
        assert len(ct) == 1
        assert ct[0]["resseq"] == 30

    def test_geometry_parsed(self):
        issues = parse_molprobity_output(self.SAMPLE_OUTPUT)
        geo = [i for i in issues if i["issue_type"] == ISSUE_GEOMETRY]
        assert len(geo) == 1
        assert geo[0]["chain"] == "B"
        assert geo[0]["resseq"] == 50
        assert geo[0]["resname"] == "LYS"

    def test_empty_output_returns_empty_list(self):
        assert parse_molprobity_output("") == []

    def test_no_issues_section(self):
        text = "MolProbity score: 0.5\nAll good!\n"
        assert parse_molprobity_output(text) == []

    def test_insertion_code_parsed(self):
        text = textwrap.dedent("""\
            Ramachandran not favored:
              A      10    A    MET      -60.0  130.0
        """)
        issues = parse_molprobity_output(text)
        assert len(issues) == 1
        assert issues[0]["icode"] == "A"


# ---------------------------------------------------------------------------
# collect_chain_residues
# ---------------------------------------------------------------------------

class TestCollectChainResidues:

    def test_basic(self):
        pdb = _write_tmp_pdb(_SIMPLE_PDB)
        try:
            result = collect_chain_residues(pdb)
            assert "A" in result
            assert "B" in result
            assert result["A"] == [1, 2, 3, 4, 5]
            assert result["B"] == [10, 11, 12]
        finally:
            os.unlink(pdb)

    def test_empty_pdb(self):
        pdb = _write_tmp_pdb("END\n")
        try:
            assert collect_chain_residues(pdb) == {}
        finally:
            os.unlink(pdb)


# ---------------------------------------------------------------------------
# build_refine_windows
# ---------------------------------------------------------------------------

class TestBuildRefineWindows:

    def _chain_residues(self):
        return {"A": list(range(1, 21))}  # residues 1–20

    def test_centered_window(self):
        issues = [{"chain": "A", "resseq": 10, "issue_type": ISSUE_RAMACHANDRAN}]
        wins = build_refine_windows(issues, self._chain_residues(), window_half=5)
        assert len(wins) == 1
        chain, lo, hi = wins[0]
        assert chain == "A"
        assert lo == 5    # 10 - 5
        assert hi == 15   # 10 + 5

    def test_clamp_at_chain_start(self):
        issues = [{"chain": "A", "resseq": 2, "issue_type": ISSUE_RAMACHANDRAN}]
        wins = build_refine_windows(issues, self._chain_residues(), window_half=5)
        assert wins[0][1] == 1  # clamped to start

    def test_clamp_at_chain_end(self):
        issues = [{"chain": "A", "resseq": 19, "issue_type": ISSUE_RAMACHANDRAN}]
        wins = build_refine_windows(issues, self._chain_residues(), window_half=5)
        assert wins[0][2] == 20  # clamped to end

    def test_unknown_chain_skipped(self):
        issues = [{"chain": "Z", "resseq": 5, "issue_type": ISSUE_CBETA}]
        wins = build_refine_windows(issues, self._chain_residues(), window_half=5)
        assert wins == []

    def test_multiple_issues(self):
        issues = [
            {"chain": "A", "resseq": 3,  "issue_type": ISSUE_CBETA},
            {"chain": "A", "resseq": 18, "issue_type": ISSUE_ROTAMER},
        ]
        wins = build_refine_windows(issues, self._chain_residues(), window_half=5)
        assert len(wins) == 2


# ---------------------------------------------------------------------------
# merge_windows
# ---------------------------------------------------------------------------

class TestMergeWindows:

    def test_no_overlap(self):
        wins = [("A", 1, 5), ("A", 10, 15)]
        merged = merge_windows(wins)
        assert merged == [("A", 1, 5), ("A", 10, 15)]

    def test_overlapping(self):
        wins = [("A", 1, 10), ("A", 8, 20)]
        merged = merge_windows(wins)
        assert merged == [("A", 1, 20)]

    def test_adjacent(self):
        wins = [("A", 1, 5), ("A", 6, 10)]
        merged = merge_windows(wins)
        assert merged == [("A", 1, 10)]

    def test_different_chains_not_merged(self):
        wins = [("A", 1, 10), ("B", 5, 15)]
        merged = merge_windows(wins)
        assert len(merged) == 2

    def test_empty(self):
        assert merge_windows([]) == []

    def test_single(self):
        assert merge_windows([("A", 1, 5)]) == [("A", 1, 5)]

    def test_unsorted_input(self):
        wins = [("A", 10, 20), ("A", 1, 5), ("A", 3, 12)]
        merged = merge_windows(wins)
        # Should be sorted and merged: (1,5) + (3,12) + (10,20) → (1,20)
        assert merged == [("A", 1, 20)]

    def test_three_chains(self):
        wins = [("A", 1, 5), ("B", 1, 5), ("C", 1, 5)]
        merged = merge_windows(wins)
        assert len(merged) == 3


# ---------------------------------------------------------------------------
# copy_ss_records
# ---------------------------------------------------------------------------

class TestCopySsRecords:

    _SOURCE = textwrap.dedent("""\
        HELIX    1   1 ALA A    1  ALA A    5  1                                   5
        SHEET    1   A 2 VAL A  10  PHE A  12  0
        ATOM      1  N   ALA A   1       0.0     0.0     0.0   1.00 10.00           N
        END
    """)

    _TARGET = textwrap.dedent("""\
        HELIX    1   1 OLD A    1  OLD A    3  1                                   3
        REMARK  old remark
        ATOM      1  N   ALA A   1       0.0     0.0     0.0   1.00 10.00           N
        ATOM      2  CA  ALA A   1       1.0     1.0     1.0   1.00 10.00           C
        END
    """)

    def test_ss_records_replaced(self):
        src = _write_tmp_pdb(self._SOURCE)
        tgt = _write_tmp_pdb(self._TARGET)
        with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as out:
            out_path = out.name
        try:
            copy_ss_records(src, tgt, out_path)
            content = Path(out_path).read_text()
            # Old HELIX record should be gone
            assert "OLD" not in content
            # New HELIX/SHEET from source should be present
            assert "HELIX    1   1 ALA" in content
            assert "SHEET    1   A 2 VAL" in content
            # Remark and ATOM records from target should be preserved
            assert "old remark" in content
            assert "ATOM      1" in content
        finally:
            os.unlink(src)
            os.unlink(tgt)
            os.unlink(out_path)

    def test_target_without_ss(self):
        src = _write_tmp_pdb(self._SOURCE)
        tgt_content = "ATOM      1  N   ALA A   1       0.0     0.0     0.0   1.00 10.00           N\nEND\n"
        tgt = _write_tmp_pdb(tgt_content)
        with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as out:
            out_path = out.name
        try:
            copy_ss_records(src, tgt, out_path)
            content = Path(out_path).read_text()
            assert "HELIX" in content
            assert "ATOM" in content
        finally:
            os.unlink(src)
            os.unlink(tgt)
            os.unlink(out_path)

    def test_source_without_ss(self):
        """When source has no SS records, target HELIX/SHEET are stripped."""
        src_content = "ATOM      1  N   ALA A   1       0.0     0.0     0.0   1.00 10.00           N\nEND\n"
        src = _write_tmp_pdb(src_content)
        tgt = _write_tmp_pdb(self._TARGET)
        with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as out:
            out_path = out.name
        try:
            copy_ss_records(src, tgt, out_path)
            content = Path(out_path).read_text()
            assert "HELIX" not in content
        finally:
            os.unlink(src)
            os.unlink(tgt)
            os.unlink(out_path)


# ---------------------------------------------------------------------------
# check_dependencies
# ---------------------------------------------------------------------------

class TestCheckDependencies:

    def test_missing_phenix_raises(self):
        with mock.patch("shutil.which", return_value=None), \
             mock.patch("os.path.isfile", return_value=False):
            with pytest.raises(RuntimeError, match="Phenix MolProbity"):
                check_dependencies("phenix.molprobity", "coot")

    def test_missing_coot_raises(self):
        def fake_which(exe):
            return "/usr/bin/phenix.molprobity" if "phenix" in exe else None

        with mock.patch("shutil.which", side_effect=fake_which), \
             mock.patch("os.path.isfile", return_value=False):
            with pytest.raises(RuntimeError, match="Coot"):
                check_dependencies("phenix.molprobity", "coot")

    def test_both_present_no_error(self):
        with mock.patch("shutil.which", return_value="/usr/bin/tool"):
            check_dependencies("phenix.molprobity", "coot")


# ---------------------------------------------------------------------------
# run_molprobity (mocked subprocess)
# ---------------------------------------------------------------------------

class TestRunMolprobity:

    def test_returns_stdout_stderr_returncode(self, tmp_path):
        fake_result = mock.MagicMock()
        fake_result.stdout = "MolProbity output"
        fake_result.stderr = ""
        fake_result.returncode = 0

        pdb = tmp_path / "test.pdb"
        pdb.write_text("ATOM  ...\n")

        with mock.patch("subprocess.run", return_value=fake_result):
            stdout, stderr, rc = run_molprobity(
                pdb_path=str(pdb),
                output_dir=str(tmp_path),
            )

        assert "MolProbity output" in stdout
        assert rc == 0

    def test_log_path_written(self, tmp_path):
        fake_result = mock.MagicMock()
        fake_result.stdout = "stdout text"
        fake_result.stderr = "stderr text"
        fake_result.returncode = 0

        pdb = tmp_path / "model.pdb"
        pdb.write_text("ATOM  ...\n")
        log_file = tmp_path / "mol.log"

        with mock.patch("subprocess.run", return_value=fake_result):
            run_molprobity(
                pdb_path=str(pdb),
                output_dir=str(tmp_path),
                log_path=str(log_file),
            )

        assert log_file.exists()
        content = log_file.read_text()
        assert "stdout text" in content
        assert "stderr text" in content


# ---------------------------------------------------------------------------
# run_coot_real_space_refine_windows (mocked subprocess)
# ---------------------------------------------------------------------------

class TestRunCootRsrWindows:

    def test_coot_called_with_no_graphics(self, tmp_path):
        pdb = tmp_path / "input.pdb"
        pdb.write_text("ATOM  ...\n")
        map_ = tmp_path / "map.mrc"
        map_.write_text("")
        out = tmp_path / "refined.pdb"

        fake_result = mock.MagicMock()
        fake_result.returncode = 0
        fake_result.stdout = ""
        fake_result.stderr = ""

        captured_cmd = []

        def fake_run(cmd, **kwargs):
            captured_cmd.extend(cmd)
            return fake_result

        with mock.patch("subprocess.run", side_effect=fake_run):
            run_coot_real_space_refine_windows(
                pdb_path=str(pdb),
                map_path=str(map_),
                windows=[("A", 1, 11)],
                output_path=str(out),
            )

        assert "--no-graphics" in captured_cmd
        assert "--script" in captured_cmd

    def test_nonzero_returncode_raises(self, tmp_path):
        pdb = tmp_path / "input.pdb"
        pdb.write_text("ATOM  ...\n")
        map_ = tmp_path / "map.mrc"
        map_.write_text("")
        out = tmp_path / "refined.pdb"

        fake_result = mock.MagicMock()
        fake_result.returncode = 1
        fake_result.stdout = "error"
        fake_result.stderr = "bad"

        with mock.patch("subprocess.run", return_value=fake_result):
            with pytest.raises(RuntimeError, match="Coot real-space refine"):
                run_coot_real_space_refine_windows(
                    pdb_path=str(pdb),
                    map_path=str(map_),
                    windows=[("A", 1, 11)],
                    output_path=str(out),
                )

    def test_multiple_windows_all_in_script(self, tmp_path):
        """Verify that all windows appear in the generated Coot script."""
        pdb = tmp_path / "input.pdb"
        pdb.write_text("ATOM  ...\n")
        map_ = tmp_path / "map.mrc"
        map_.write_text("")
        out = tmp_path / "refined.pdb"

        script_content = []

        fake_result = mock.MagicMock()
        fake_result.returncode = 0
        fake_result.stdout = ""
        fake_result.stderr = ""

        def fake_run(cmd, **kwargs):
            # Read the script that was written
            script_idx = cmd.index("--script") + 1
            script_path = cmd[script_idx]
            try:
                with open(script_path) as f:
                    script_content.append(f.read())
            except OSError:
                pass
            return fake_result

        windows = [("A", 1, 11), ("A", 50, 60), ("B", 5, 15)]
        with mock.patch("subprocess.run", side_effect=fake_run):
            run_coot_real_space_refine_windows(
                pdb_path=str(pdb),
                map_path=str(map_),
                windows=windows,
                output_path=str(out),
            )

        assert script_content, "Script was not read"
        script = script_content[0]
        assert '"A" 1 11' in script
        assert '"A" 50 60' in script
        assert '"B" 5 15' in script
        assert "set-secondary-structure-restraints-type 2" in script
