# In Silico Drug Discovery Pipeline — LRRK2 Inhibitors for Parkinson's Disease

This repository contains an end-to-end computational drug discovery pipeline written in R. The pipeline identifies novel, potent, and safe lead compounds targeting **Leucine-Rich Repeat Kinase 2 (LRRK2)**, a heavily validated genetic target for **Parkinson's Disease**.

By combining API-driven data curation, cheminformatics, machine learning (QSAR), and 3D molecular docking, this project circumvents traditional wet-lab limitations and algorithmically explicitly hunts for **structural novelty**.

---

## 🚀 Final Identified Lead Compound

Out of thousands of compounds screened, the pipeline identified **CHEMBL2152707 (Cmpd_001)** as the clear winner through our multi-parameter priority scoring formula. 

* **ChEMBL ID:** CHEMBL2152707
* **Max Phase:** Preclinical
* **Molecular Formula:** C18H21N5O
* **Molecular Weight:** 323.40
* **Molecule Type:** Small molecule
* **Canonical SMILES:** `Cn1cc(/C=C/c2n[nH]c3ccnc(OC4CCCCC4)c23)cn1`

This compound successfully achieved an ideal balance of predicted high potency, excellent structural novelty, optimal docking affinity in the LRRK2 kinase pocket, and perfect compliance with CNS (Blood-Brain Barrier) penetration parameters.

---

## 🧪 Pipeline Overview

The computational pipeline executes across 9 sequential phases:

1. **Target Selection:** Validates LRRK2 using the Open Targets GraphQL API and PubMed literature.
2. **Data Collection:** Fetches bioactivity data (IC50) dynamically via the ChEMBL REST API.
3. **Featurization:** Calculates 1D/2D physicochemical properties and ECFP4 Circular Fingerprints using the `rcdk` package.
4. **QSAR Modeling:** Trains a Random Forest machine learning model on a hierarchical scaffold-split to predict pIC50 based on molecular structure.
5. **Chemical Space:** Projects multi-dimensional chemical data into 2D space using PCA (global) and t-SNE (local neighborhoods) to validate SAR.
6. **Virtual Screening:** Evaluates the library using a custom algorithm that rewards both *predicted potency* and *structural novelty*.
7. **3D Preparation & Docking:** Utilizes **AlphaFold2** to generate the LRRK2 kinase domain receptor, and AutoDock Vina / PyRx to predict binding affinity (kcal/mol).
8. **ADMET Filtering:** Applies Lipinski's Rule of 5, Veber filters, and the Pfizer CNS Multi-Parameter Optimization (MPO) boundaries to ensure safety and blood-brain barrier penetration.
9. **Priority Scoring:** Ranks the final candidates using a composite mathematical formula, identifying the final lead compound.

---

## 💻 Reproducibility Instructions

To reproduce the analysis locally, follow these steps:

### 1. Prerequisites
You will need **R** and **RStudio** installed.
Because the pipeline relies on the Chemistry Development Kit (`rcdk`), you **must** have a Java Development Kit (JDK) installed (e.g., JDK 21). 

Ensure your `JAVA_HOME` path is correctly set in your R environment before running. Inside the R script, adjust the following line to match your system path:
```R
Sys.setenv(JAVA_HOME = "C:/Program Files/Java/jdk-21.0.11") # Adjust to your JDK path
```

### 2. Install Required Packages
The script is designed to automatically install missing packages. It relies on:
`httr2`, `jsonlite`, `dplyr`, `tibble`, `purrr`, `tidyr`, `ggplot2`, `DT`, `stringr`, `scales`, `rentrez`, `knitr`, `rcdk`, `fingerprint`, `caret`, `randomForest`, `Rtsne`.

### 3. Run the Pipeline
Open the provided `.Rmd` file in RStudio and click **Knit**, or run the chunks sequentially. This will recreate the data curation, QSAR modeling, chemical space visualization, and 2D virtual screening.

### 4. AlphaFold2 & 3D Molecular Docking
Because LRRK2 is massive (>2500 amino acids), computing its 3D structure locally is resource-prohibitive.
1. We utilized the Google Colab environment to fold the Kinase Domain (residues 1875–2132). 
2. **Colab Link Used:** [AlphaFold2.ipynb (sokrypton/ColabFold)](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb)
3. Download the resulting `.pdb` file.
4. Export the top docking candidates' SMILES from the R script.
5. Use **PyRx / AutoDock Vina** locally to prepare the ligands (PDBQT format) and dock them into the ATP-binding pocket of the AlphaFold receptor.

---

## 📚 Glossary

* **pIC50:** The negative base-10 logarithm of the IC50 value (in Molar). It linearizes dose-response data. A higher pIC50 means higher potency (e.g., pIC50 9 = 1 nM).
* **QSAR (Quantitative Structure-Activity Relationship):** A machine learning approach that models the relationship between a molecule's mathematical features (descriptors) and its biological activity.
* **ECFP4 (Extended-Connectivity Fingerprints):** A method to digitize molecules into a binary string by hashing the topological neighborhoods (up to radius 2) of every atom.
* **t-SNE:** A statistical method for visualizing high-dimensional data by giving each datapoint a location in a two or three-dimensional map, preserving local data clusters.
* **AlphaFold2:** An AI program developed by DeepMind that performs highly accurate predictions of protein 3D structures directly from amino acid sequences.
* **Molecular Docking:** A computational simulation of a candidate ligand binding to a protein receptor to predict the preferred orientation and calculate binding affinity (kcal/mol).
* **ADMET:** Absorption, Distribution, Metabolism, Excretion, and Toxicity. The pharmacokinetic profile of a drug.
* **CNS MPO (Central Nervous System Multi-Parameter Optimization):** An algorithm developed by Pfizer that uses properties like LogP, Molecular Weight, and TPSA to predict if a drug can successfully cross the Blood-Brain Barrier.
