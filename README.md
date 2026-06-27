
# Dynamic Local Average Treatment Effects in Time Series

**Reproducibility Package for "Dynamic Local Average Treatment Effects in Time Series"**

* **Date Assembled**: June 2026
* **Author(s)**: Alessandro Casini, Luca Mattia Rolla, Adam McCloskey, Raimondo Pala
* **Contact**: 
  * Alessandro Casini (alessandro.casini@uniroma2.it, University of Rome Tor Vergata, Italy)
  * Luca Mattia Rolla (lmrolla92@gmail.com, University of Rome Tor Vergata, Italy)
  * Adam McCloskey (adam.mccloskey@colorado.edu, University of Colorado at Boulder, USA)
  * Raimondo Pala (raimondo.pala@uniroma2.it, University of Rome Tor Vergata, Italy)
* **Repository URL**: [Insert GitHub URL]

## Overview

This repository contains the code and data to reproduce the simulations and empirical results for the paper "Dynamic Local Average Treatment Effects in Time Series" (Casini, Rolla, McCloskey, & Pala, 2025). 

The paper addresses the pervasive issue of weak identification in instrumental variable (IV) models utilizing time-series data. It establishes that strong identification often holds only in specific, contiguous subsamples where the instrument actively induces changes in treatment status. The methodology relies on dynamic programming algorithms to select the most strongly identified subsample, facilitating robust structural inference.

The core methodological contributions are divided into three sections:
1. **Section 4**: Introduces a supremum-type pre-test ($F_T^*$) for identification failure across all possible sample partitions.
2. **Section 5**: Develops the estimation theory for the strongly identified subsample ($\hat{S}_T$), providing algorithms for Two-Stage Least Squares (2SLS) and Feasible Generalized Least Squares (FGLS) under heteroskedasticity and autocorrelation.
3. **Section 6**: Establishes identification-robust inference procedures by extending the Anderson-Rubin (AR), Lagrange Multiplier (LM), and Conditional Likelihood Ratio (CLR) tests to the dynamic LATE framework.

## Repository Structure

The project is organized into modular directories mapping directly to the paper's econometric sections, alongside core algorithmic dependencies and a unified user interface.

### `sez4`: Testing for Identification Failure
Contains scripts to reproduce the Monte Carlo simulations and theoretical critical values for the $F_T^*$ test (Section 4).
* `runSection4_Simul.m`: Master script to execute the simulation suite for Section 4 exhibits (e.g., Table 4, Figure 5, Figure 6). Incorporates options for exact legacy replication (fixed random seeds, specific variance scaling).
* `piLATE_Ftest.m`: Core function computing the supremum $\pi$-LATE $F^*$ statistic across possible structural breaks.
* `Compute_CVs.m` / `RunMc_CVs_Fstat.m`: Subroutines located in `cvals/` to simulate the non-standard asymptotic distributions (Wiener processes) required for the critical values.

### `sez5`: Subsample Estimation
Contains scripts to replicate the estimation properties (Bias and MSE) detailed in Section 5.
* `runSection5_Simul.m`: Master script to run the simulations evaluating the point estimators.
* `est_2sls_pilate.m`: Implementation of the 2SLS estimator on the selected subsample, featuring options for White (HC1) or Newey-West (HAC) standard errors via exact matrix projection.
* `est_fgls_pilate.m`: Implementation of the FGLS estimator for seemingly unrelated regression (SUR) representations, relying on multivariate break detection algorithms.

### `sez6`: Robust Inference
Contains the codebase for computing size and size-adjusted power of identification-robust tests (Section 6).
* `runSection6_Simul.m`: Master script executing the computational suite for Tables 9-10 and Figures 7-8.
* `robust_LATE_inference.m`: Unified engine generating the Anderson-Rubin, Lagrange Multiplier, and Conditional Likelihood Ratio statistics (AMS and ST variants) conditioned on the estimated optimal partition.

### Core Dependencies
* **`m_break`**: Contains the fundamental univariate dynamic programming routines for identifying structural changes and maximizing the first-stage F-statistic.
* **`Mbreak_Multivariate`**: Contains the multivariate extensions (Qu and Perron, 2007) utilized primarily for the FGLS/SUR estimation protocols in `sez5`.

## Unified Routine: `Dynamic_LATE_Analysis.m`

For researchers intending to apply the dynamic LATE methodology to proprietary datasets without engaging with the low-level simulation architecture, the repository provides a high-level, unified interface in the root directory: `Dynamic_LATE_Analysis.m`.

### Functionality
This function acts as a comprehensive analytical pipeline:
1. **Input Parsing & Sandboxing**: Accepts structural variables (`Y`, `D`, `Z`, `X`) and an `options` structure (defining $\pi_L$, max breaks, HAC parameters, etc.). It safely modifies the MATLAB path to include dependencies during execution and restores the global environment upon completion.
2. **Structural Extraction**: Iterates the identification process over a grid of minimal sample fractions ($\pi_L$), extracting the optimal break dates and maximizing the identification strength.
3. **Automated Inference**: Computes the robust AR, LM, and LR statistics on the selected subsamples and maps them directly against the built-in asymptotic critical values (for $\alpha = 0.10, 0.05, 0.01$) derived in the paper.
4. **Output**: Returns a structured `Results_Table` ready for econometric post-processing.

**Basic Usage:**
```matlab
% options.piLow = 0.6; options.m_breaks = 5; options.Hac_option = true;
Results = Dynamic_LATE_Analysis(Y, D, Z, X, options);
