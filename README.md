# Dynamic Local Average Treatment Effects in Time Series

**Package for "Dynamic Local Average Treatment Effects in Time Series"**

* **Date Assembled**: June 2026
* **Author(s)**: Alessandro Casini, Adam McCloskey, Luca Rolla, Raimondo Pala
* **Contact**: 
  * Alessandro Casini (alessandro.casini@uniroma2.it, University of Rome Tor Vergata, Italy)
  * Adam McCloskey (adam.mccloskey@colorado.edu, University of Colorado at Boulder, USA)
  * Luca Rolla (lmrolla92@gmail.com, University of Rome Tor Vergata, Italy)
  * Raimondo Pala (raimondo.pala@uniroma2.it, University of Rome Tor Vergata, Italy)
* **Repository URL**: [Insert GitHub URL]

## Overview

This repository implements testing, estimation and inference based on IV in time series. It also contains the code and data to reproduce the simulations and empirical results based on the methods of Sections 4-6 in the paper. 

The core methodological contributions are divided into three sections:
1. **Section 4**: Proposes a pre-test for identification failure based on the supremum of the first-stage F-statistic over all possible data partitions ($F_T^*$). The test seeks to determine whether there is any subsample such that identification is strong.
2. **Section 5**: Estimates the most strongly-identified subsample and discusses the estimation of $\beta$ for the most strongly-identified subsample $\hat{S}_T$. 
3. **Section 6**: Develops identification-robust inference for $\beta$. We construct identification-robust tests—namely, Anderson-Rubin (AR), Lagrange multiplier (LM), and conditional likelihood ratio (CLR) tests—evaluated at the estimated partition $\hat{S}_T$.

## Repository Structure

The project is organized into separate directories mapping directly to the paper's sections, alongside core algorithmic dependencies and a unified user interface.

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
```

### Example: `test_Dynamic_LATE_all_DGPs.m`
The `test_Dynamic_LATE_all_DGPs.m` script executes a sequential test of the `Dynamic_LATE_Analysis.m` pipeline. It demonstrates the utility of the methodology across two Data Generating Processes (DGPs) defined in the paper's supplement to evaluate identification, estimation, and inference performance.

#### Scenario 1 (Subsample Estimation)
This scenario implements the linear IV regression model with exogenous covariates ($X_t$) defined by Equations (S.C.2) and (S.C.3).  
* **Codebase Mapping**: Algorithmically, this is an adaptation of the `generate_data_tab67` subroutine from `runSection5_Simul.m` in the `Section5` directory. It replicates the exact setup employed to evaluate the Bias and MSE of the point estimators in **Tables 6 and 7**.
* **Equations**:
  $Y_t = \beta D_t + \gamma_1 X_t + u_t$
  $D_t = \theta Z_t + \gamma_2 X_t + e_t$
* **Mechanism**: It evaluates the properties of the 2SLS and FGLS estimators on the strongly-identified subsample $\hat{S}_T$. In regimes 1 and 3, identification is local-to-zero ($\theta = d T^{-1/2}$). In Regime 2 (the anomaly), the instrument is irrelevant ($\theta = 0$), allowing the test to verify the consistency of $\hat{\beta}(\hat{S}_T)$ when identification fails over parts of the sample.

#### Scenario 2 (Robust Inference / Perverse Identification)
This scenario employs the baseline Model M1 (without exogenous covariates) but modifies the first-stage coefficient $\theta_t$ to create a "perverse identification" setting.  
* **Codebase Mapping**: Computationally, this DGP is drawn directly from `runSection6_Simul.m` in the `Section6` directory. It implements the specific subroutine designed for the stress-test used to generate the size-adjusted power curves under perverse identification shown in **Figure 8**.
* **Equations**:
  $Y_t = \beta D_t + u_t$
  $D_t = \theta_t Z_t + e_t$
* **Mechanism**: We set $\theta_1 = \theta_3 = d T^{-1/2}$ (strong identification) and $\theta_2 = -0.5 \theta_1$ (perverse identification in the anomaly regime).
* **Objective**: This tests the robustness of the identification-robust AR, LM, and CLR tests. It demonstrates that the proposed procedures remain valid even when the first-stage coefficient changes sign across regimes, a setting where full-sample IV estimators typically suffer from severe bias.

**How to run the test:**
To execute the full validation suite, run the following command from the MATLAB root directory:
```matlab
run test_Dynamic_LATE_all_DGPs.m
```

## Replication of the Results

### `Section4`: Testing for Identification Failure
Contains scripts to reproduce the Monte Carlo simulations and critical values of $F_T^*$ (Section 4).

#### `Folder 1: ftest`
* `runSection4_Simul.m`: Master script to execute the simulation suite for Section 4 exhibits (e.g., Table 4, Table 5, Figure 5, Figure 6). Incorporates options for exact legacy replication (fixed random seeds, specific variance scaling).
* `piLATE_Ftest.m`: Core function computing the supremum $\pi$-LATE $F_T^*$ statistic across possible structural breaks.

#### `Folder 2: cvals`
* `Compute_CVs.m` / `RunMc_CVs_Fstat.m`: Subroutines located in `cvals/` to simulate the non-standard asymptotic distributions required for the critical values.

### `Section5`: Subsample Estimation
Contains scripts to replicate the estimation properties detailed in Section 5.
* `runSection5_Simul.m`: Master script to run the simulations evaluating the point estimators.
* `est_2sls_pilate.m`: Implementation of the 2SLS estimator for the most strongly-identified subsample $\hat{S}_T$, featuring options for White (HC1) or Newey-West (HAC) standard errors.
* `est_fgls_pilate.m`: Implementation of the FGLS estimator for seemingly unrelated regression (SUR) representations, relying on multivariate break detection algorithms.

### `Section6`: Robust Inference
Contains the codebase for computing size and size-adjusted power of identification-robust tests (Section 6).
* `runSection6_Simul.m`: Master script executing the computational suite for Tables 9-10 and Figures 7-8.
* `robust_LATE_inference.m`: Unified engine generating the Anderson-Rubin, Lagrange Multiplier, and Conditional Likelihood Ratio statistics evaluated at the estimated partition $\hat{S}_T$.

