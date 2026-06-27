# Dynamic Local Average Treatment Effects in Time Series

**Reproducibility Package for "Dynamic Local Average Treatment Effects in Time Series"**

* **Date Assembled**: June 2026
* **Author(s)**: Alessandro Casini, Adam McCloskey, Luca Rolla, Raimondo Pala
* **Contact**: 
  * Alessandro Casini (alessandro.casini@uniroma2.it, University of Rome Tor Vergata, Italy)
  * Adam McCloskey (adam.mccloskey@colorado.edu, University of Colorado at Boulder, USA)
  * Luca Rolla (lmrolla92@gmail.com, University of Rome Tor Vergata, Italy)
  * Raimondo Pala (raimondo.pala@uniroma2.it, University of Rome Tor Vergata, Italy)
* **Repository URL**: [Insert GitHub URL]

## Overview

This repository contains the code and data to reproduce the simulations and empirical results for the paper "Dynamic Local Average Treatment Effects in Time Series" (Casini, McCloskey, Rolla, & Pala, 2025). 

This paper discusses identification, estimation, and inference on dynamic local average treatment effects (LATEs) in instrumental variables (IVs) settings. We document pervasive weak identification in applied work using IVs with time series data. However, we find that strong identification often holds in large subsamples for which the instrument induces changes in the treatment. Motivated by this, we introduce a method based on dynamic programming to detect the most strongly-identified subsample and show how to use this subsample to improve estimation and inference. We also develop new identification-robust inference procedures.

The core methodological contributions are divided into three sections:
1. **Section 4**: Proposes a pre-test for identification failure based on the supremum of the first-stage F-statistic over all possible data partitions ($F_T^*$). The test seeks to determine whether there is any subsample such that identification is strong.
2. **Section 5**: Discusses the estimation of $\beta$ for the most strongly-identified subsample $\hat{S}_T$. We study the properties of the two-stage least-squares (2SLS) estimator and the feasible generalized least-squares (FGLS) estimator under autocorrelated errors.
3. **Section 6**: Develops identification-robust inference for $\beta$ when the subsample $S_{0,T}$ is unknown. We construct identification-robust tests—namely, Anderson-Rubin (AR), Lagrange multiplier (LM), and conditional likelihood ratio (CLR) tests—evaluated at the estimated partition $\hat{S}_T$.

## Repository Structure

The project is organized into modular directories mapping directly to the paper's sections, alongside core algorithmic dependencies and a unified user interface.

### `sez4`: Testing for Identification Failure
Contains scripts to reproduce the Monte Carlo simulations and critical values of $F_T^*$ (Section 4).
* `runSection4_Simul.m`: Master script to execute the simulation suite for Section 4 exhibits (e.g., Table 4, Table 5, Figure 5, Figure 6). Incorporates options for exact legacy replication (fixed random seeds, specific variance scaling).
* `piLATE_Ftest.m`: Core function computing the supremum $\pi$-LATE $F_T^*$ statistic across possible structural breaks.
* `Compute_CVs.m` / `RunMc_CVs_Fstat.m`: Subroutines located in `cvals/` to simulate the non-standard asymptotic distributions required for the critical values.

### `sez5`: Subsample Estimation
Contains scripts to replicate the estimation properties detailed in Section 5.
* `runSection5_Simul.m`: Master script to run the simulations evaluating the point estimators.
* `est_2sls_pilate.m`: Implementation of the 2SLS estimator for the most strongly-identified subsample $\hat{S}_T$, featuring options for White (HC1) or Newey-West (HAC) standard errors.
* `est_fgls_pilate.m`: Implementation of the FGLS estimator for seemingly unrelated regression (SUR) representations, relying on multivariate break detection algorithms.

### `sez6`: Robust Inference
Contains the codebase for computing size and size-adjusted power of identification-robust tests (Section 6).
* `runSection6_Simul.m`: Master script executing the computational suite for Tables 9-10 and Figures 7-8.
* `robust_LATE_inference.m`: Unified engine generating the Anderson-Rubin, Lagrange Multiplier, and Conditional Likelihood Ratio statistics evaluated at the estimated partition $\hat{S}_T$.

### Core Dependencies
* **`m_break`**: Contains the fundamental univariate dynamic programming routines for identifying structural changes and maximizing the first-stage F-statistic.
* **`Mbreak_Multivariate`**: Contains the multivariate extensions (Qu and Perron, 2007) utilized primarily for the FGLS/SUR estimation protocols in `sez5`.

## Unified Routine: `Dynamic_LATE_Analysis.m`

For researchers intending to apply the dynamic LATE methodology to proprietary datasets, the repository provides a unified interface in the root directory: `Dynamic_LATE_Analysis.m`.

### Functionality
This function implements the sequential procedure described in the paper:
1. **Data Input**: Accepts structural variables (`Y`, `D`, `Z`, `X`) and an `options` structure.
2. **Subsample Detection**: Detects the most strongly-identified subsample by maximizing the first-stage F-statistic over a grid of minimal sample fractions ($\pi_L$).
3. **Inference**: Computes the robust AR, LM, and LR statistics on the selected subsamples and reports the built-in asymptotic critical values derived in the paper.
4. **Output**: Returns a structured `Results_Table`.

**Basic Usage:**
```matlab
% options.piLow = 0.6; options.m_breaks = 5; options.Hac_option = true;
Results = Dynamic_LATE_Analysis(Y, D, Z, X, options);

## System and Software Environment

The results in the paper were obtained using the following computational environment.

* **Software**:
  * **MATLAB Version**: R2025b
* **MATLAB Toolboxes Required**:
  * Curve Fitting Toolbox
  * Deep Learning Toolbox
  * Econometrics Toolbox
  * Financial Toolbox
  * Global Optimization Toolbox
  * Optimization Toolbox
  * Parallel Computing Toolbox
  * Statistics and Machine Learning Toolbox
  * Symbolic Math Toolbox

## References

* Andrews, D. W. K., Moreira, M., & Stock, J. H. (2006). Optimal Two-Sided Invariant Similar Tests for Instrumental Variables Regression. *Econometrica*, 74(3), 715-752.
* Bai, J., & Perron, P. (1998). Estimating and Testing Linear Models with Multiple Structural Changes. *Econometrica*, 66(1), 47-78.
* Casini, A., McCloskey, A., Rolla, L., & Pala, R. (2025). Dynamic Local Average Treatment Effects in Time Series. *Unpublished manuscript*, University of Rome Tor Vergata.
* Imbens, G. W., & Angrist, J. D. (1994). Identification and Estimation of Local Average Treatment Effects. *Econometrica*, 62(2), 467-475.
* Lewis, D. J. (2022). Robust Inference in Models Identified via Heteroskedasticity. *Review of Economics and Statistics*, 104(3), 512-525.
* Nakamura, E., & Steinsson, J. (2018). High-Frequency Identification of Monetary Non-Neutrality: The Information Effect. *The Quarterly Journal of Economics*, 133(3), 1283-1330.
* Newey, W. K., & West, K. D. (1987). A Simple, Positive Semi-Definite, Heteroskedasticity and Autocorrelation Consistent Covariance Matrix. *Econometrica*, 55(3), 703-708.
* Qu, Z., & Perron, P. (2007). Estimating and Testing Structural Changes in Multivariate Regressions. *Econometrica*, 75(2), 459-502.
* Rigobon, R. (2003). Identification through Heteroskedasticity. *Review of Economics and Statistics*, 85(4), 777-792.
