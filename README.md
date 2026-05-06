<div align="center">

# 🧬 LRRK2 In Silico Drug Discovery
**End-to-end R pipeline identifying novel, CNS-penetrant LRRK2 inhibitors for Parkinson's Disease.**

[![R](https://img.shields.io/badge/R-276DC3?style=flat-square&logo=r&logoColor=white)](#)
[![ChEMBL](https://img.shields.io/badge/ChEMBL-Data-00A68A?style=flat-square)](#)
[![AutoDock Vina](https://img.shields.io/badge/AutoDock_Vina-Docking-FF6F00?style=flat-square)](#)

[Overview](#overview) · [Setup](#setup) · [Steps](#steps) · [Results](#results) · [Glossary](#glossary) ·[Reference](#reference)

</div>

---

## Overview
Computational pipeline predicting LRRK2 kinase domain inhibitors using QSAR modeling, structural novelty scoring, and molecular docking. Data is LRRK2 bioactivity measurements sourced from the[ChEMBL Database](https://www.ebi.ac.uk/chembl/).

## Setup
Install dependencies and configure the environment.
```R
# Install cheminformatics and modeling packages
required_packages <- c("httr2", "dplyr", "rcdk", "fingerprint", "randomForest", "ggplot2", "Rtsne")
to_install <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(to_install) > 0) install.packages(to_install)
invisible(lapply(required_packages, library, character.only = TRUE))

# Point rcdk to Java installation
Sys.setenv(JAVA_HOME = "C:/Program Files/Java/jdk-21.0.11")
```

## Steps

Query Open Targets and ChEMBL APIs for raw target validation and bioactivity data.
```R
# Fetch LRRK2 (CHEMBL1075104) IC50 data
url <- "https://www.ebi.ac.uk/chembl/api/data/activity.json?target_chembl_id=CHEMBL1075104&standard_type__in=IC50,Ki&limit=500"
req <- httr2::request(url) %>% httr2::req_headers(Accept = "application/json")
raw_data <- httr2::resp_body_json(httr2::req_perform(req), simplifyVector = TRUE)

# Filter exact measurements and convert to pIC50 scale
curated <- raw_data %>%
  filter(standard_relation == "=", standard_units %in% c("nM", "uM", "M")) %>%
  mutate(pIC50 = -log10(standard_value * 1e-9)) # assuming nM base for example
```

Convert SMILES to 3D molecules, calculate properties, and generate topological fingerprints.
```R
# Parse SMILES strings into molecular objects
mols <- parse.smiles(curated$canonical_smiles)

# Calculate physicochemical descriptors
desc_vals <- eval.desc(mols, c("org.openscience.cdk.qsar.descriptors.molecular.WeightDescriptor", "org.openscience.cdk.qsar.descriptors.molecular.XLogPDescriptor"))

# Generate Extended-Connectivity Fingerprints (ECFP4) for machine learning
fps <- lapply(mols, get.fingerprint, type = "circular")
fp_matrix <- do.call(rbind, lapply(fps, function(x) as.numeric(x@bits)))
```

Train a Random Forest model on scaffold-split chemical data to predict pIC50.
```R
# Train Random Forest to predict potency based on descriptors
rf_model <- randomForest(x = train_x, y = train_y, ntree = 500, mtry = 3, importance = TRUE)

# Predict activity for test compounds
all_predicted <- predict(rf_model, newdata = all_desc)
```

Filter candidates by structural novelty and export for physical docking validation.
```R
# Combine Novelty (1 - Tanimoto) and QSAR predictions into a screening score
curated <- curated %>%
  mutate(Combined_Score = 0.5 * Novelty_Norm + 0.5 * QSAR_Norm) %>%
  arrange(desc(Combined_Score))

# Export top 8 candidates
write.csv(head(curated, 8), "docking_candidates.csv", row.names = FALSE)
```

Fold the target protein and run docking simulation outside of R.
```text
# 1. Provide LRRK2 Kinase domain sequence (Residues 1875–2132) to AlphaFold2.
# 2. Prepare AlphaFold .pdb receptor and ligand .pdbqt files in PyRx.
# 3. Execute AutoDock Vina against ATP-binding pocket (Met1947 gatekeeper).
```

Calculate composite score factoring potency, novelty, safety, and CNS penetration.
```R
# Rank leads using a normalized multi-parameter matrix
prioritisation <- curated %>%
  mutate(
    Priority_Score = 0.20 * pIC50_Norm + 0.15 * Novelty_Norm + 
                     0.25 * Docking_Norm + 0.25 * ADMET_Score + 0.15 * CNS_Bonus
  ) %>%
  arrange(desc(Priority_Score))
```

## Results

| Visual / Output | Description |
| :--- | :--- |
| **Feature Importance Plot** | Horizontal bar chart of `%IncMSE` showing `LogP` and `MW` drive inhibition. |
| **Chemical Space Map** | t-SNE scatter plot clustering structurally similar active compounds. |
| **CNS MPO Radar Chart** | Overlay of candidate footprint vs optimal Blood-Brain Barrier limits. |
| **Priority Score Bar Chart** | Final ranking of docked compounds based on multi-parameter matrix. |
| **Lead Candidate** | **`CHEMBL2152707`** emerges as the top novel, safe, CNS-penetrant inhibitor. |

## Glossary
* **pIC50**: Negative log of the IC50 value. Linearizes dose-response for machine learning.
* **ECFP4**: Extended-Connectivity Fingerprints. Converts molecular topology into binary arrays.
* **QSAR**: Quantitative Structure-Activity Relationship. Statistical model linking chemical structure to biological activity.
* **t-SNE**: Algorithm used to visualize high-dimensional chemical fingerprints in 2D space.
* **ADMET**: Absorption, Distribution, Metabolism, Excretion, and Toxicity.
* **CNS MPO**: Multi-Parameter Optimization score specifically for predicting Blood-Brain Barrier permeability.

## Reference
Protein structure folding powered by AlphaFold2 Colab:[https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb)
