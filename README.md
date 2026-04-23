```markdown
# CarboKitten Extension

[![License](https://img.shields.io/badge/license-GPLv3-blue.svg)](LICENSE)
[![Julia](https://img.shields.io/badge/julia-1.10+-9558B2.svg)](https://julialang.org)
[![Docs](https://img.shields.io/badge/docs-stable-blue.svg)](docs/)
[![DOI](https://shields.io)](https://doi.org/10.xxxx/joss.xxxxx)

## Overview
**CarboKitten** is an open-source, Julia-based Stratigraphic Forward Model (SFM) derived from **CarboCAT (Burgess, 2012)**. It simulates carbonate platform evolution through sediment diffusion and ecological competition.

This repository provides an **extended implementation** specifically enhanced for reservoir-oriented modeling workflows, introducing physics-based forcing and high-resolution stratigraphic tracking.

![Model overview](docs/images/model_ext_overview.png)

---

## Key Features

The extension introduces several critical modifications to the core CarboKitten engine:

*   **Dynamic Carbonate Production:** Production rates vary through user-defined depth-dependent and time-dependent curves.
*   **Spatially Variable Subsidence:** Prescription of spatially distributed subsidence fields with optional temporal scaling.
*   **Layer-Based Stratigraphy & Compaction:** Stores stratigraphy as a discrete stack of layers. A factory-dependent engine calculates porosity reduction at every timestep, ensuring preserved thickness reflects cumulative burial history.
*   **Physics-Based Wave Model:** Multi-component hydrodynamic forcing that derives orbital velocity from wave amplitude, period, and direction.
*   **Automated Facies Classification:** Deposits are classified based on sediment proportions, wave energy, and water depth.
*   **Enhanced Visualization:** Integrated routines for fence diagrams, map views, and chronostratigraphic (Wheeler) plots.
*   **Quantitative Metrics Export:** Automated export of platform statistics, patch size analysis, and factory proportions to CSV for sensitivity analysis.

---

## Installation

This version is tested with **Julia 1.10+**.

1. **Clone the repository:**

   ```bash
   git clone https://github.com/username/CarboKitten.jl
   cd CarboKitten.jl
   ```

2. **Instantiate environment:**

   ```bash
   julia --project -e 'using Pkg; Pkg.instantiate()'
   ```
---

# Quick Start : Running the Example

We provide a fully reproducible workflow that generates all figures and metrics presented in the associated JOSS paper.
Run the example from the root directory:


```bash
julia --project run_extension_example.jl
```

Outputs will be saved to `examples/extension/output/` and include:
* Fence diagrams and map-view visualizations.
* CSV files containing thickness statistics and factory distribution metrics.

---

# Outputs

The model produces several outputs including:

* stratigraphic grids
* facies classifications
* platform thickness statistics
* chronostratigraphic visualizations
* CSV files for quantitative analysis

---

# Repository structure

```
CarboKitten.jl
├── src/                # Core model implementation
├── examples/           # Workflow examples
│   └── extension/      # Extension-specific simulation
├── data/               # Input datasets
├── docs/               # Original model documentation
├── docs_extension/     # Extension-specific documentation
├── test/               # Unit tests
├── ext/                # Visualization tools & optional extensions
└── Project.toml        # Julia dependencies
```

---
# Documentation
Comprehensive guides on model inputs and workflows are located in the `/docs_extension` directory. For the base engine logic, refer to the original `/docs`.

---

# Citation

If you use this software, please cite:

Salmon et al. (2025)
*CarboKitten: An open-source stratigraphic forward model for carbonate platform evolution.*
Journal of Open Source Software.

---

# License

This project is distributed under the terms of the **GNU General Public License v3.0**.
See `LICENSE` for details.

---

# Acknowledgements

This work was supported by the Swiss State Secretariat for Education, Research and Innovation (SERI) and the European Union Horizon Europe programme (Grant No. 101147618).

---
```



