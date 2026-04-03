# AutoStructure on Windows

This guide covers setup and troubleshooting for Windows systems.

## Quick Start

### Option 1: One-Command Launch (Recommended)

```powershell
PS> .\run_windows.ps1 -PdbFile "E:\data\model.pdb" -MapFile "E:\data\map.mrc"
```

The script will:
- Auto-detect ChimeraX, Phenix, and Coot installations
- Create a Windows-compatible config file
- Run the full pipeline
- Report success or failure with color-coded output

### Option 2: Manual Setup

1. **Create Python virtual environment:**
   ```powershell
   python -m venv .venv
   .\.venv\Scripts\Activate.ps1
   pip install -r requirements.txt
   ```

2. **Copy and customize config:**
   ```powershell
   Copy-Item config.yaml config.local.yaml
   ```
   Edit `config.local.yaml` with your tool paths (see below).

3. **Run pipeline:**
   ```powershell
   python pipeline.py --pdb "path\to\model.pdb" --map "path\to\map.mrc" --config config.local.yaml
   ```

## Windows Configuration Template

Create `config.local.yaml` with paths matching your installation:

```yaml
tools:
  chimerax: "D:/ChimeraX 1.8/bin/ChimeraX-console.exe"
  phenix: "E:/Phenix/phenix_bin/phenix.real_space_refine.bat"
  coot: "E:/wincoot/bin/coot-bin.exe"

fitmap:
  resolution: 3.0
  metric: "correlation"
  max_steps: 2000

refinement:
  nproc: 1
  coot_fit_protein: false    # ← Use Phenix morphing (more stable on Windows)

output:
  work_dir: "autostructure_output"
```

### Finding Tool Paths

#### ChimeraX
- Navigate to your ChimeraX installation, typically:
  - `C:\Program Files\ChimeraX` or
  - `D:\ChimeraX 1.8`
- Look for `bin\ChimeraX-console.exe`

#### Phenix
- Check your Phenix root, typically:
  - `E:\Phenix` or
  - `C:\Phenix`
- Use `phenix_bin\phenix.real_space_refine.bat` (not `.exe`)

#### Coot
- Check typically:
  - `E:\wincoot\bin`
- Use `coot-bin.exe` (not just `coot`, which may be a shell script)

## Windows-Specific Compatibility Notes

### 1. Coot Headless Mode Issue

**Problem:** Coot may hang in headless script mode on some Windows builds, leaving the Python subprocess waiting indefinitely.

**Solution:** Use Phenix morphing refinement for loop pieces instead:
```yaml
refinement:
  coot_fit_protein: false  # Use Phenix morphing
```

This is the default and recommended setting for Windows.

### 2. Phenix Parameter Compatibility

**Problem:** Older Phenix parameter syntax like `adp=True` or `rigid_body=True` causes ambiguity errors in recent versions.

**Solution:** The code has been updated to use explicit parameter paths:
- ✓ Use: `refinement.run=rigid_body` or `refinement.run=morphing`
- ✗ Don't use: `rigid_body=True`, `morphing=True`, `adp=True`

### 3. Path Separators

PowerShell paths can use `/` or `\`, but for consistency with YAML:
- Use `/` in `config.yaml` files
- PowerShell scripts handle both automatically

Example:
```yaml
tools:
  chimerax: "D:/ChimeraX/bin/ChimeraX-console.exe"  # Use forward slash in YAML
  phenix: "E:/Phenix/phenix_bin/phenix.real_space_refine.bat"
```

### 4. Performance: Single Worker for First Run

The default config uses `nproc: 1` (single worker). This simplifies initial debugging on Windows.

For faster processing on multi-core systems, increase after verification:
```yaml
refinement:
  nproc: 4  # Adjust to your CPU count
```

## Troubleshooting

### Error: "WinError 193 %1 is not a valid Win32 application"

**Cause:** Executable path points to a script wrapper instead of a binary.

**Solution:**
- For Coot: Use `coot-bin.exe`, not `coot` (no extension)
- For Phenix: Use `phenix.real_space_refine.bat`, not the Python module

### Error: "Unrecognized PHIL parameters"

**Cause:** Phenix version incompatibility with parameter syntax.

**Check:** The `phenix_runner.py` has been updated. If you modified it, verify:
```python
# ✓ Correct
extra_args = [
    "refinement.run=morphing",
    "refinement.max_iterations=100",
    "refinement.macro_cycles=5",
]

# ✗ Wrong
extra_args = [
    "morphing=True",
    "adp=True",
    "output.file_name_prefix=...",  # Use output_dir instead
]
```

### Error: "Cannot find ChimeraX/Phenix/Coot"

**Solution:** Use explicit paths in the launcher:
```powershell
.\run_windows.ps1 -PdbFile "path\model.pdb" `
                  -MapFile "path\map.mrc" `
                  -ChimeraXExe "D:\ChimeraX 1.8\bin\ChimeraX-console.exe" `
                  -PhenixExe "E:\Phenix\phenix_bin\phenix.real_space_refine.bat" `
                  -CootExe "E:\wincoot\bin\coot-bin.exe"
```

### Pipeline Takes Too Long

**Check:** Most time is spent in per-piece refinement (normal for 90+ pieces).

**Optimize:**
1. Reduce `fitmap.max_steps` (default 2000, try 500-1000)
2. Reduce `refinement.macro_cycles` in `phenix_runner.py`
3. Increase `refinement.nproc` if using `pheniex_runner.py`

## Advanced Usage

### Custom Configuration File

```powershell
.\run_windows.ps1 -PdbFile "model.pdb" `
                  -MapFile "map.mrc" `
                  -ConfigFile "my_config.yaml" `
                  -OutputDir "my_output"
```

### Using Coot for Loop Refinement (if stable)

If your Coot build is headless-friendly:
```yaml
refinement:
  coot_fit_protein: true
```

### Running from Batch File

Create `run_autostructure.bat`:
```batch
@echo off
powershell -NoProfile -ExecutionPolicy Bypass -File "%~dp0run_windows.ps1" ^
    -PdbFile %1 -MapFile %2
pause
```

Then drag-and-drop files or use:
```cmd
run_autostructure.bat "model.pdb" "map.mrc"
```

## File Structure After Run

```
autostructure_output/
├── globally_fitted.pdb          # Step 1: Global fit result
├── pieces/                       # Step 3: Split pieces
│   ├── piece_0000_loop_chainA.pdb
│   ├── piece_0001_structured_chainA.pdb
│   └── ...
├── fitted_pieces/               # Step 4: Per-piece fits
│   ├── piece_A_1_3_fitted.pdb
│   └── ...
├── refined_pieces/              # Step 5: Refined pieces
│   ├── piece_A_1_3/
│   │   └── piece_A_1_3_fitted_real_space_refined_000.pdb
│   └── ...
├── logs/                         # Step logs
│   ├── global_fitmap.log
│   ├── piece_A_1_3_fitmap.log
│   ├── piece_A_1_3_morphing.log
│   └── ...
└── inmapfold_autostructure_refined.pdb  # FINAL OUTPUT
```

## Next Steps

- Validate the output PDB with MolProbity or similar tools
- Compare quality metrics (correlation, geometry) before/after refinement
- Adjust parameters and re-run if needed (intermediate files are retained)
