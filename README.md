# STB-Suite — Siesta ToolBox Suite

> **⚠️ Development Version — Not Yet Released**
> This is an unreleased, work-in-progress version of STB-Suite. Features may be incomplete, APIs may change without notice, and not all tools have been exhaustively tested. Use at your own risk in production workflows.

A comprehensive Python toolkit for pre- and post-processing [SIESTA](https://siesta-project.org) DFT simulations, developed at the University of Brasília.

---

## Features

### Pre-processing (Calculation Setup)
| Tool | Command | Description |
|---|---|---|
| Input File Generator | `stb-inputfile` | Generate `calc.fdf` from a structure file (supports total energy, relax, AIMD, bands, +D3) |
| K-Grid Generator | `stb-kgrid` | Suggest a Monkhorst-Pack k-point mesh from a target density |
| K-Path Generator | `stb-kpath` | Generate high-symmetry k-paths for band structure calculations |
| Strain Generator | `stb-strain` | Create strained structures for stress-strain curves |
| Elastic Inputs | `stb-elasticInputs` | Generate deformed structures for elastic tensor calculations |
| 2D Monolayer Stacker | `stb-2Dstacking` | Stack two monolayers into a heterostructure (ZSL algorithm) |
| Cohesive Energy Setup | `stb-cohesive` | Prepare input folders for cohesive energy calculations |
| Phonon Displacement Generator | `stb-phononsCreate` | Generate Phonopy displacement folders for SIESTA phonon runs |

### Post-processing (Analysis)
| Tool | Command | Description |
|---|---|---|
| Bands Analyzer | `stb-bands` | Parse `.bands` files, calculate band gaps, and plot band structures |
| PDOS Parser | `stb-dos` | Extract projected DOS from `PDOS.xml` (total, per atom, per species) |
| DOS Convolution | `stb-convdos` | Apply Gaussian smoothing to DOS data files |
| Structure Analyzer | `stb-structural` | Calculate ECN and structural properties |
| Symmetry Analyzer | `stb-symmetry` | Determine space group and crystal symmetry |
| Strain Analysis | `stb-strainAnalysis` | Extract stress-strain curves from `strain_*/` folders |
| Elastic Analysis | `stb-elasticAnalysis` | Compute stiffness tensor, Young's modulus, and mechanical stability |
| Bader Charge Analysis | `stb-bader` | Compute atomic charges via the Bader AIM method (`.RHO` + `.XV`) |
| Work Function Calculator | `stb-workfunction` | Calculate the work function from the electrostatic potential (`.VT`) |
| Charge Density Plotter | `stb-density` | Export 2D planar cuts or 3D charge density clouds from `.RHO` |
| Cohesive Energy Analysis | `stb-cohesiveAnalysis` | Process SIESTA outputs and compute cohesive energy per atom |
| Phonon Post-processing | `stb-phononsPos` | Extract forces, build `FORCE_SETS`, and compute thermal properties |

### Utilities
| Tool | Command | Description |
|---|---|---|
| File Translator | `stb-translate` | Convert between structure formats (FDF, POSCAR, CIF, XYZ, FHI, XSF, DFTB+) |
| Grid to Cube | `stb-cube` | Convert SIESTA grid files (`.RHO`, `.VT`, `.VH`) to Gaussian `.cube` |
| Clean Tool | `stb-clean` | Remove calculation artifacts from a directory |
| Wantibexos Interface | `stb-siesta2wtb` | Convert SIESTA Hamiltonian to WantiBEXOS format |

All tools are also accessible through the interactive menu:
```
stb-suite
```

---

## Installation

### Requirements
- Python **3.12**
- [Conda](https://docs.conda.io/en/latest/) (recommended)

### Recommended: Conda environment

```bash
# Create a dedicated environment with Python 3.12
conda create -n stb python=3.12
conda activate stb

# Install dependencies available on conda-forge
conda install -c conda-forge numpy matplotlib scipy ase pymatgen spglib phonopy

# Install the remaining dependencies and the package via pip
pip install sisl pybader
pip install .
```

### Alternative: pip only

```bash
pip install .
```

> **Note:** `pybader` may require additional system dependencies depending on your OS. If the pip install fails, refer to the [pybader documentation](https://github.com/adam-kerrigan/pybader).

---

## Quick Start

Launch the interactive suite menu:
```bash
stb-suite
```

Or run individual tools directly, for example:
```bash
# Generate a SIESTA input file from a structure
stb-inputfile structure.fdf --type relax

# Suggest a k-point mesh
stb-kgrid --file structure.fdf --type fdf --density 0.2

# Plot band structure
stb-bands --file siesta.bands --shift fermi
```

Each tool supports `--help` for full usage details:
```bash
stb-bands --help
```

---

## Typical Workflow

```
[Structure: CIF / POSCAR / FDF]
        │
        ├── stb-translate      → convert to FDF if needed
        ├── stb-kgrid          → choose k-point mesh
        ├── stb-kpath          → generate k-path for bands
        └── stb-inputfile      → build calc.fdf
                │
                ▼
         *** SIESTA runs ***
                │
        ├── stb-bands          → band structure & gap
        ├── stb-dos            → projected DOS
        ├── stb-bader          → atomic charges
        ├── stb-workfunction   → work function
        ├── stb-density        → charge density maps
        └── stb-elasticAnalysis / stb-strainAnalysis / stb-cohesiveAnalysis
```

---

## Status

> This package is under active development. The tools listed above are functional but **have not been exhaustively tested** across all use cases, crystal systems, or SIESTA versions. Bug reports and feedback are welcome via the [issue tracker](https://github.com/stb-suite/stb-suite/issues).

---

## License

MIT License — see [LICENSE](LICENSE) for details.

## Author

Developed by **Dr. Carlos M. O. Bastos** — University of Brasília, 2025.
