[![Python 3.7+](https://img.shields.io/badge/python-3.7+-blue.svg)](https://www.python.org/) 
[![License: LGPL v3](https://img.shields.io/badge/License-LGPL%20v3-blue.svg)](https://www.gnu.org/licenses/lgpl-3.0)
[![DOI](https://zenodo.org/badge/276040050.svg)](https://doi.org/10.5281/zenodo.17945808)

<a href="https://github.com/MABeeskow/GebPy">
<img src="https://raw.githubusercontent.com/MABeeskow/GebPy/master/documents/readme_images/readme_header.png" 
width="67%">
</a>

---

# GebPy

**GebPy** is an open-source, Python-based **engine for the synthetic generation of geological data**, with a focus 
on **minerals, rocks, and stratigraphic sequences**.

GebPy is designed to generate **physically, chemically, and crystallographically consistent** datasets under controlled 
assumptions. The core idea is that bulk rock properties are primarily determined by the underlying **mineral 
assemblage**, while structural and tectonic effects are intentionally excluded.

The project targets **research, teaching, and data-driven exploration workflows**, including the generation of 
high-quality training data for machine learning applications.

---

## 🚀 Project Status

* **Current version:** 1.1.5
* **License:** LGPL-3.0
* **Development:** Active
* **API stability:** evolving (pre-2.0)

The core engine of GebPy is considered stable. The public API may evolve until a future 2.0.0 release, where long-term 
API stability is intended.

---

## 🧱 Design Philosophy

GebPy is designed as a **scientific engine**, not as a monolithic application.

Key principles:

* Engine-first architecture
* Explicit physical, chemical, and crystallographic assumptions
* Reproducibility over visual polish
* Notebook-first workflows

A graphical user interface (GUI) exists for exploratory use but is considered **legacy** and will be deprecated in 
future releases. The recommended way to use GebPy is via **Python scripts and Jupyter notebooks**.

---

## 🚀 Installation

Install the latest stable release from PyPI:

```bash
pip install gebpy
```

For development or research use:

```bash
git clone https://github.com/MABeeskow/GebPy.git
cd GebPy
pip install .
```

---

## 🧪 Quick Example

```python
from gebpy import Carbonates

calcite = Carbonates(name="Calcite", random_seed=42)
data = calcite.generate_dataset(number=1)

print(data["rho"])
```

---

## 💎 Synthetic Mineral Data

Synthetic mineral data are generated based on:

* elemental properties derived from mineral formulas
* crystallographic cell parameters defining the crystal lattice
* elastic parameters (bulk and shear modulus) taken from literature or databases

From these inputs, GebPy derives, among others:

* molar mass $M$
* molar and cell volume $V_m$ and $V$
* density $\varrho$
* seismic velocities $v_P$ and $v_S$
* elastic properties $E$ and $\nu$
* electromagnetic properties $GR$, $PE$, $U$

In addition, the full **chemical composition** of each mineral is calculated.

---

## 🪨 Synthetic Rock Data

Based on synthetic mineral datasets, GebPy computes bulk rock properties, including chemical composition.

Calculated properties include:

* bulk density $\varrho_B$
* porosity $\phi$
* seismic velocities $v_P$ and $v_S$
* elastic properties $E$ and $\nu$
* electromagnetic properties $GR$, $PE$, $U$

The current database contains **more than 100 minerals**, allowing the modeling of common sedimentary, igneous, 
metamorphic, and ore-related rocks.

---

## 🏞️ Stratigraphic Data

Since GebPy can generate synthetic rock data, it also supports the construction of **stratigraphic sequences** by 
stacking lithological units.

Please note that GebPy intentionally does **not** consider effects related to:

* structural geology
* tectonics
* fluid flow
* geothermal gradients
* lithostatic pressure

These simplifications are deliberate and serve to maintain conceptual clarity and reproducibility.

---

## 🎯 Intended Use Cases

* Teaching mineralogy, petrology, and geophysics
* Synthetic data generation for machine learning
* Benchmarking and validation of exploration algorithms
* Conceptual testing of geological hypotheses

GebPy is **not** intended to replace field observations or laboratory measurements.

---

## 💭 Citation

If you use GebPy in academic or industrial work, please cite:

> Beeskow, M. A. (2025).  
> *GebPy – A Python engine for synthetic geological data generation*.

A DOI will be provided via **Zenodo** for stable releases.

---

## ⚖️ License

GebPy is released under the **LGPL-3.0 License**.

You are free to use, modify, and redistribute the software under the terms of this license. Modifications to the GebPy 
source code must remain open.

---

## 👤 Author

**Maximilian A. Beeskow**
Geoscientist (Geochemistry, Mineral Systems, Hydrothermal Processes)
RWTH Aachen University

* ORCID: https://orcid.org/0000-0003-2456-5228
* GitHub: [https://github.com/MABeeskow](https://github.com/MABeeskow)
* Zenodo: https://doi.org/10.5281/zenodo.17945809

---

## 🗂️ Legacy Code Notice

Folders such as `gebpy_legacy/` and `modules_archive/` contain archived or transitional code from earlier development 
stages. These components are **not part of the public API** and are preserved solely for reproducibility and historical 
reference.

---

*Last updated: 16.12.2025*
