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

### Command line

```bash
python pipeline.py --pdb alphafold_model.pdb --map density.mrc
```

Optional arguments:

```
--config config.yaml   Custom configuration file (tool paths, resolution, etc.)
--output refined.pdb   Output PDB path (default: <work_dir>/<stem>_autostructure_refined.pdb)
```

### Python API

```python
from pipeline import run_pipeline

final_pdb = run_pipeline(
    pdb_path="alphafold_model.pdb",
    map_path="density.mrc",
    config_path="config.yaml",   # optional
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
  coot_fit_protein: true  # false → use Phenix morphing for loops

output:
  work_dir: autostructure_output
```

## Module Overview

| File | Purpose |
|------|---------|
| `pipeline.py` | Main entry point – orchestrates the full workflow |
| `chimerax_runner.py` | ChimeraX CLI wrapper (`fit_whole_model`, `fit_piece`) |
| `secondary_structure.py` | Parse HELIX/SHEET records; partition chains into pieces |
| `piece_splitter.py` | Split a PDB into sub-files; merge refined pieces back |
| `phenix_runner.py` | Phenix rigid-body and morphing refinement wrappers |
| `coot_runner.py` | Coot headless fit-protein wrapper |
| `config.yaml` | Default configuration |

