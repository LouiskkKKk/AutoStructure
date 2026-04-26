# AutoStructure

Automated protein structure refinement pipeline using **ChimeraX**, **Phenix**, and **Coot**.

## Overview

AutoStructure takes an AlphaFold2 (or similar) predicted PDB and an experimental
cryo-EM density map, then performs the following steps automatically:

1. **Global fit-in-map** – fit the full predicted model into the EM map (ChimeraX `fitmap`)
2. **Secondary-structure parsing** – extract helix/sheet annotations from the fitted PDB
3. **Piece splitting** – divide the structure into *structured* (helix/sheet) and *loop* pieces
4. **Per-piece fit-in-map** – independently fit each piece into the EM map (ChimeraX `fitmap`)
5. **Refinement**
   - Structured pieces → **Phenix rigid-body real-space refinement**
   - Loop/disordered pieces → **Coot fit-protein** (or **Phenix morphing**, configurable)
6. **Merge** – reassemble all refined pieces into a single output PDB
7. **Auto-fix** *(optional)* – iterative MolProbity validation + Coot real-space refine loop

## Requirements

| Software | Version tested |
|----------|---------------|
| Python   | ≥ 3.8         |
| ChimeraX | ≥ 1.0         |
| Phenix   | ≥ 1.20        |
| Coot     | ≥ 0.9         |

Python dependencies:

```bash
pip install -r requirements.txt
```

## Usage

### Windows (Recommended: One-Command Launcher)

```powershell
PS> .\run_windows.ps1 -PdbFile "E:\data\model.pdb" -MapFile "E:\data\map.mrc"
```

The launcher script will auto-detect ChimeraX, Phenix, and Coot, then run the full pipeline.
For detailed setup, tool paths, and troubleshooting, see [WINDOWS_SETUP.md](WINDOWS_SETUP.md).

### Command line (All platforms)

```bash
python pipeline.py --pdb alphafold_model.pdb --map density.mrc
```

Optional arguments:

```
--config config.yaml   Custom configuration file (tool paths, resolution, etc.)
--output refined.pdb   Output PDB path (default: <work_dir>/<stem>_autostructure_refined.pdb)
--auto-fix             Enable the post-pipeline auto-fix loop (MolProbity + Coot RSR)
--no-auto-fix          Disable the auto-fix loop even if config has auto_fix.enabled: true
```

### Python API

```python
from pipeline import run_pipeline

final_pdb = run_pipeline(
    pdb_path="alphafold_model.pdb",
    map_path="density.mrc",
    config_path="config.yaml",   # optional
    auto_fix=True,               # optional; overrides config auto_fix.enabled
)
print(f"Refined model: {final_pdb}")
```

## Configuration

Copy and edit `config.yaml` to customise tool paths and parameters:

```yaml
tools:
  chimerax: "/path/to/chimerax"
  phenix:   "phenix.real_space_refine"
  coot:     "/path/to/coot"

fitmap:
  resolution: 3.0       # EM map resolution in Å
  metric: correlation   # correlation | cam
  max_steps: 2000

refinement:
  nproc: 4              # parallel workers for per-piece refinement
  coot_fit_protein: false  # false → use Phenix morphing for loops (recommended for Windows)

output:
  work_dir: autostructure_output
```

### Windows-Specific Configuration

On Windows, use forward slashes (`/`) in paths and prefer **Phenix morphing** for loop refinement (more stable than Coot headless):

```yaml
tools:
  chimerax: "D:/ChimeraX 1.8/bin/ChimeraX-console.exe"
  phenix: "E:/Phenix/phenix_bin/phenix.real_space_refine.bat"
  phenix_molprobity: "E:/Phenix/phenix_bin/phenix.molprobity.bat"
  coot: "E:/wincoot/bin/coot-bin.exe"

refinement:
  nproc: 1              # Use 1 for initial run; increase after verification
  coot_fit_protein: false  # Use Phenix morphing on Windows
```

For automated tool detection and config generation, use the provided launcher:
```powershell
.\ run_windows.ps1 -PdbFile "model.pdb" -MapFile "map.mrc"
```

See [WINDOWS_SETUP.md](WINDOWS_SETUP.md) for Windows-specific paths, tool locations, and troubleshooting.

## Auto-fix: MolProbity Validation + Coot Real-Space Refine Loop

The optional **auto-fix** step runs *after* the main pipeline finishes.  It
repeatedly validates the merged PDB with `phenix.molprobity`, identifies
problematic residues, and fixes them with Coot real-space refinement until
the structure is clean or the iteration limit is reached.

### What it fixes

| MolProbity category       | Description                                    |
|---------------------------|------------------------------------------------|
| Ramachandran outliers      | φ/ψ angles in disallowed regions               |
| Rotamer outliers           | side-chain rotamer in poor conformation        |
| Cβ deviations              | Cβ atom position > 0.25 Å from ideal          |
| Cis / twisted peptides     | ω angle deviating from 0° / 180°              |
| Geometry restraints        | bond-length / bond-angle / chirality outliers  |

### Enabling auto-fix

**Option A – CLI flag:**
```bash
python pipeline.py --pdb model.pdb --map density.mrc --auto-fix
```

**Option B – config file (`config.yaml`):**
```yaml
auto_fix:
  enabled: true
  max_iterations: 10   # stop after 10 rounds even if issues remain
  window_half: 5       # refine ±5 residues around each problem site
```

**Option C – Python API:**
```python
from pipeline import run_pipeline
run_pipeline("model.pdb", "density.mrc", auto_fix=True)
```

### Required dependencies

| Tool | Used for |
|------|----------|
| `phenix.molprobity` | Structure validation (Ramachandran, rotamers, Cbeta, geometry) |
| `coot` | Headless real-space refinement with secondary-structure restraints |

Configure the executables in `config.yaml`:
```yaml
tools:
  phenix_molprobity: "phenix.molprobity"   # or absolute path
  coot: "coot"                              # or absolute path
```

### Output directory layout

```
autostructure_output/
├── globally_fitted.pdb
├── ...  (main pipeline outputs)
├── <stem>_autostructure_refined.pdb   ← main pipeline final PDB
└── auto_fix/
    ├── iteration_01/
    │   ├── input_with_ss.pdb          ← current PDB with AF2 SS records injected
    │   ├── validation/                ← phenix.molprobity output files
    │   ├── molprobity.log             ← combined stdout/stderr log
    │   ├── molprobity_report.txt      ← full validation text report
    │   ├── issues_summary.txt         ← per-category issue counts
    │   ├── refined.pdb                ← Coot RSR output
    │   └── coot_rsr.log               ← Coot stdout/stderr log
    ├── iteration_02/
    │   └── ...
    └── ...
```

### Iteration logic

1. Inject original AF2 HELIX/SHEET records into the current PDB (to ensure
   correct secondary-structure restraints in Coot).
2. Run `phenix.molprobity` and parse the five categories listed above.
3. **Stop** if there are no issues, or if no *new* issues appeared since the
   previous round (converged).
4. Build a ±`window_half`-residue window around each problem residue, then
   merge overlapping windows.
5. Run a single Coot headless session that refines all merged windows with
   secondary-structure restraints enabled.
6. Repeat from step 1 with the new PDB, up to `max_iterations` rounds.

### Running tests

```bash
python -m pytest tests/test_validation_refine.py -v
```

Tests use `unittest.mock` to simulate Phenix/Coot calls and exercise all
core logic (parser, window builder, window merger, SS-record injection,
subprocess wrappers) without requiring the external tools to be installed.

## Module Overview

| File | Purpose |
|------|---------|
| `pipeline.py` | Main entry point – orchestrates the full workflow |
| `chimerax_runner.py` | ChimeraX CLI wrapper (`fit_whole_model`, `fit_piece`) |
| `secondary_structure.py` | Parse HELIX/SHEET records; partition chains into pieces |
| `piece_splitter.py` | Split a PDB into sub-files; merge refined pieces back |
| `phenix_runner.py` | Phenix rigid-body and morphing refinement wrappers |
| `coot_runner.py` | Coot headless fit-protein wrapper |
| `validation_refine.py` | MolProbity validation parser + Coot RSR auto-fix loop |
| `config.yaml` | Default configuration |
| `config.local.yaml` | Windows-compatible config (auto-generated by `run_windows.ps1`) |
| `run_windows.ps1` | Windows launcher script (auto-detects tools, creates config, runs pipeline) |
| `WINDOWS_SETUP.md` | Windows-specific setup, tool paths, troubleshooting guide |
| `tests/test_validation_refine.py` | Unit tests for the auto-fix module (mock-based) |

## Platform-Specific Notes

### Windows

- **Recommended:** Use the `run_windows.ps1` launcher for one-command execution
- **Tool paths:** Use `/` in YAML config files (not `\`)
- **Loop refinement:** Use Phenix morphing (`coot_fit_protein: false`) for stability
- **Parallelization:** Start with `nproc: 1` for easier debugging; increase later
- See [WINDOWS_SETUP.md](WINDOWS_SETUP.md) for detailed tool installation and troubleshooting

### Linux / macOS

- Standard usage applies; see **Usage** and **Configuration** sections above
- Tool paths typically in `/usr/local/bin` or similar standard locations
- Coot headless mode generally works well
