# 🧭 Valley Scanner Research Repository

## Datasets, Validation Results, and Reproducibility Resources

**Author:** Jacob Orellana Real

**Contact:** [jacoboreore@gmail.com](mailto:jacoboreore@gmail.com)

---

# 📘 Overview

The **Valley Scanner** project investigates numerical and geometric approaches for exploring the landscape of the Riemann zeta function along the critical line.

This repository contains the research materials associated with the project, including:

* Published datasets
* Validation results
* Visualization utilities
* Reproducibility resources
* Manuscript and preprint materials

The repository is intended to support independent verification, reproducibility, and further exploration of the numerical observations described in the accompanying publication.

---

# 🎯 Project Motivation

The Valley Scanner project explores the geometric structure of the Hardy (Z(t)) function through continuous sampling along the critical line.

Rather than relying on a predefined sequence of special points, the method investigates the local geometry of the function by examining the evolution of ( |Z(t)| ) across successive samples.

This approach provides a visual and computational framework for locating candidate minima, which may subsequently be refined and numerically validated as nontrivial zeros of the Riemann zeta function.

The objective of the project is not to replace established methods, but to provide an alternative geometric perspective for numerical exploration and validation.

---

# 🧩 Repository Contents

## 📄 Academic Materials

This repository includes manuscript and preprint materials associated with the Valley Scanner project.

The publication discusses:

* Numerical exploration of the Hardy (Z(t)) landscape
* Geometric interpretation of local minima
* Validation against known zeta zeros
* Statistical comparisons with expected zero spacing
* Reproducibility procedures

---

## 📂 Datasets

The repository contains datasets generated during the experimental validation process.

Examples include:

* Refined zero samples
* Spacing measurements
* Validation datasets
* High-(t) numerical experiments

Typical dataset structure:

```text
t, absZ, spacing
```

Where:

* **t** = imaginary component on the critical line
* **absZ** = ( |Z(t)| )
* **spacing** = distance to the subsequent validated zero

These datasets are provided to support independent verification and reproduction of reported observations.

---

## 🧪 Visualization Playground

The file:

```bash
playground.py
```

provides a simple visualization environment for exploring the local geometry of ( |Z(t)| ).

Example:

```bash
python3 playground.py
```

The script produces graphical representations of local maxima and minima within selected ranges of the critical line.

This utility is intended for educational and exploratory purposes.

---

## ⚙️ Technical Notes

Selected experiments were performed using:

* MPFR arbitrary-precision arithmetic
* Multi-core CPU execution
* Distributed computation workflows for large-scale evaluations

The datasets included in this repository represent the published outputs of those experiments.

---

# 📈 Reproducibility

The Valley Scanner project emphasizes reproducibility.

This repository therefore provides:

1. Public datasets used in validation.
2. Supporting visualization tools.
3. Documentation of experimental procedures.
4. References to published numerical results.

Researchers are encouraged to independently verify, reproduce, and critique the presented results.

---

# 💡 Quick Start

Clone the repository:

```bash
git clone https://github.com/zjore/z-research.git
cd z-research
```

Install dependencies:

```bash
pip install -r requirements.txt
```

Run the visualization playground:

```bash
python3 playground.py
```

---

### 🌐 Web Interface

The same functionality is available through the web interface:

👉 [**https://p56yzukrvv.us-east-1.awsapprunner.com/**](https://p56yzukrvv.us-east-1.awsapprunner.com/)

This interface allows users to launch the same computational modes from a browser-based environment with live progress tracking and dataset export capabilities.

---

# 🛰️ Licensing and Citation

## Repository Components

| Component                           | License                                                        |
| ----------------------------------- | -------------------------------------------------------------- |
| Datasets                            | CC BY 4.0                                                      |
| Python visualization utilities      | MIT License                                                    |
| Manuscripts and preprints           | Copyright © Jacob Orellana                                     |
| Valley Scanner computational engine | Separate repository under Valley Scanner Research License v1.0 |

---

## Citation

If you use the datasets, visualizations, numerical results, or supporting materials contained in this repository, please cite:

> Orellana, J.
>
> *Valley Scanner: A Numerical Framework for Exploring the Riemann Zeta Landscape.*

When available, please include the corresponding DOI and repository URL.

---

## Attribution

Academic and research use of the materials contained in this repository is encouraged.

Appropriate attribution to the Valley Scanner project and its author is required when reproducing, extending, or redistributing the published materials.

---

## Computational Engine

The Valley Scanner computational engine is maintained in a separate repository under the **Valley Scanner Research License v1.0**.

This repository contains research artifacts and reproducibility resources only.

---

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17566257.svg)](https://doi.org/10.5281/zenodo.17566257)

---

*"Walk the mountains, rest at the valleys. All is revealed with symmetry."*

— *Valley Scanner Project*
