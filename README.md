# Gene Expression Modeling with multiome data (scATAC-Seq + scRNA-Seq)

This repository contains code and analysis related to my BMI 206 (BIOSTAT) individual project. Built lasso and glms to predict gene expression
from nearby enhancer accessibility. This was done in comparison to SCENT https://github.com/immunogenomics/SCENT/tree/main 

**Key Analyses:**

* **SCENT Analysis:** Exploration of enhancer-gene relationships using the SCENT tool, this runs individual gene x peak poisson regressions
* **GLM Modeling:** Fitting Generalized Linear Models (GLM) to predict gene expression: gene_expr ~ peak_covariate_matrix
* **Lasso Regression:** Implementing Lasso regression to identify sparse sets of influential enhancers.

**Files:**

* `SCENT_1_qsub_fibroblast_not_bin....`: Script for running SCENT analysis (wynton (SGE) submission).
* `SCENT_parallelization.R`: Script for parallelizing SCENT analyses.
* `glms_single_core.R`: R script for running GLMs/Lasso on a single core (for loops)
* `BMI_206_individual_analysis.ipynb`: Jupyter Notebook for individual-level analysis (SEE FINAL FIGURES HERE)
* `Group_sensitivity_analyses.ipynb`: Jupyter Notebook for group project, looking at sensitivity of SCENT derived coefficients from with
  inclusion of nUMI as a covariate
* `getting_bw.ipynb`: Jupyter Notebook for greating "pseudo" bigwigs from the anndata objects so I can plot coverage
