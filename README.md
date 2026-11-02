# iCG875, eciCG875 and MPECMM

This repository provides the genome-scale metabolic model **iCG875**, its
enzyme-constrained counterpart **eciCG875**, and **MPECMM**, the
**multi-enzyme-constrained metabolic model construction method**
(多重酶约束模型构建方法) used to build, calibrate and analyze the two models.

Both models are distributed in **SBML Level 3 Version 1** with the **fbc**
(Flux Balance Constraints) package and can be opened directly by COBRA-based
tools (`cobrapy`, the COBRA Toolbox, memote, etc.).

---

## 1. Overview

| Name | Type | Description |
|---|---|---|
| `iCG875` | Genome-scale metabolic model (GEM) | Updated base GEM reconstructed for the studied strain; defines the metabolic reaction network, stoichiometry and flux bounds. |
| `eciCG875` | Enzyme-constrained GEM (ecModel) | Enzyme-constrained version of `iCG875`, integrating enzyme kinetic (kcat) and enzyme-usage (proteome allocation) constraints on top of the stoichiometric network. |
| `MPECMM` | Construction method / workflow | Multi-enzyme-constrained metabolic model construction method (多重酶约束模型构建方法) used to build, calibrate and analyze the models above. |

---

## 2. Repository contents

```
.
├── README.md            # this file
├── models/              # model files (SBML L3V1 + fbc)
│   ├── iCG875.xml       #   updated genome-scale metabolic model
│   └── eciCG875.xml     #   enzyme-constrained model
├── doc/                 # raw source data collected for model construction
└── test/                # testing and simulation results
```

| Path | Contents | Opened by |
|---|---|---|
| `models/iCG875.xml` | Genome-scale metabolic model (SBML L3V1 + fbc) | `cobra.io.read_sbml_model` / COBRA Toolbox |
| `models/eciCG875.xml` | Enzyme-constrained model (SBML L3V1 + fbc) | `cobra.io.read_sbml_model` (the enzyme parameters are preserved by the MPECMM workflow) |
| `doc/` | Raw source data used during model reconstruction and enzyme-constraint calibration |
| `test/` | Model validation, simulation and analysis results |

---

## 3. Model summary

| Model | Metabolites | Reactions | Genes |
|---|---|---|---|
| `iCG875` | ~1042 | ~1480 | 875 |
| `eciCG875` | ~1042 | ~1480 (+ enzyme-usage constraints) | 875 |

> Please confirm the exact metabolite/reaction/gene counts from the returned
> values after loading the released models, and update the table if needed.

The enzyme-constrained model `eciCG875` extends `iCG875` by coupling each
gene-associated reaction to an enzymatic capacity term derived from kcat values
and the available enzyme (proteome) budget.

---

## 4. Getting started

### 4.1 Environment

```shell
$ conda create -n icg875 python=3.9
$ conda activate icg875
$ pip install cobra libsbml pandas numpy
```

### 4.2 Load a model

```python
from cobra.io import read_sbml_model

model = read_sbml_model("models/iCG875.xml")
print(len(model.metabolites), "metabolites,",
      len(model.reactions), "reactions,",
      len(model.genes), "genes")

# inspect the objective value
print("Objective:", model.optimize().objective_value)
```

The same call loads `models/eciCG875.xml`; if the enzyme constraints are stored
as additional parameters, load them through the MPECMM interface so the kcat /
enzyme-usage constraints are preserved.

### 4.3 Run the MPECMM workflow

MPECMM is the multi-enzyme-constrained metabolic model construction method used
to derive `eciCG875` from `iCG875`. The typical pipeline is:

1. Collect and curate raw data in `doc/` (genome annotation, enzyme identifiers
   and kinetic parameters, proteome constraints).
2. Build the base GEM `iCG875` and add the enzyme (kcat and enzyme-usage)
   constraints to obtain `eciCG875`.
3. Calibrate the enzyme parameters of `eciCG875`.
4. Run validation and simulation; results are written to `test/`.

---

## 5. Data sources

The models are reconstructed / curated from public and in-house resources
(collected under `doc/`), including (non-exhaustive): BiGG (metabolite
namespace), UniProt / BRENDA / SABIO-RK (enzyme identifiers and kinetic
parameters), and the reference genome annotation of the strain. Please see the
**Materials & Methods** of the accompanying manuscript for the full provenance
of every data source.

---

## 6. Testing & simulation results

The `test/` folder contains the model validation and simulation/analysis
results (e.g. growth-rate predictions, flux distributions and the comparison
between `iCG875` and `eciCG875`), supporting the assessment reported in the
manuscript.

---

## 7. File format

Both models follow **SBML Level 3 Version 1** with the **fbc** package, ensuring
portability across COBRA toolchains. They are also prepared for deposition in a
public model repository (e.g. BioModels); an accession number, once assigned,
will be added here and in the manuscript.

---

## 8. License & contact

The models are released for academic use. For questions, please open an issue
in this repository or contact the corresponding author.

---

内容由 AI 生成
