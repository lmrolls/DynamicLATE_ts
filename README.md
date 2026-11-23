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

#### Subfolder 1: `ftest`
* `runSection4_Simul.m`: Master script to execute the simulation suite for Section 4 exhibits (e.g., Table 4, Table 5, Figure 5, Figure 6). Incorporates options for exact legacy replication (fixed random seeds, specific variance scaling).
* `piLATE_Ftest.m`: Core function computing the supremum $\pi$-LATE $F_T^*$ statistic across possible structural breaks.

**Example of usage:**
```matlab 
run runSection4_Simul.m
```

#### Subfolder 2: `cvals`
* `Compute_CVs.m` / `RunMc_CVs_Fstat.m`: Subroutines located in `cvals/` to simulate the non-standard asymptotic distributions required for the critical values.

### `Section5`: Subsample Estimation
Contains scripts to replicate the estimation properties detailed in Section 5.
* `runSection5_Simul.m`: Master script to run the simulations evaluating the point estimators.
* `est_2sls_pilate.m`: Implementation of the 2SLS estimator for the most strongly-identified subsample $\hat{S}_T$, featuring options for White (HC1) or Newey-West (HAC) standard errors.
* `est_fgls_pilate.m`: Implementation of the FGLS estimator for seemingly unrelated regression (SUR) representations, relying on multivariate break detection algorithms.

**Example of usage:**
```matlab 
run runSection5_Simul.m
```

### `Section6`: Robust Inference
Contains the codebase for computing size and size-adjusted power of identification-robust tests (Section 6).
* `runSection6_Simul.m`: Master script executing the computational suite for Tables 9-10 and Figures 7-8, and for the supplement (Section N.D): `'tabAR'` / `'figAR'` (Table 2 and Figure 1 of the supplement, AR(1) errors), `'tabMIV'` / `'figMIV'` (Table 3 and Figure 2 of the supplement, three instruments), `'suppAR'`, `'suppMIV'`, `'supp'`. The supplement targets use the seeds of the legacy scripts; the published runs used 10 000 replications for the AR(1) exhibits and 5 000 for the three-instrument ones.
* `robust_LATE_inference.m`: Unified engine generating the Anderson-Rubin, Lagrange Multiplier, and Conditional Likelihood Ratio statistics evaluated at the estimated partition $\hat{S}_T$.

**Example of usage:**
```matlab 
run runSection6_Simul.m                      % Tables 9-10, Figures 7-8
runSection6_Simul(10000, 'suppAR', true);    % supplement, AR(1) errors (Table 2, Figure 1)
runSection6_Simul(5000, 'suppMIV', true);    % supplement, three instruments (Table 3, Figure 2)
```

### `Section7`: Empirical Application
Reproduces the exhibits of Section 7 (monetary policy, Nakamura-Steinsson data in `m_break/dataForMatlabAltEstTick.csv`). Each driver is the legacy script of the replication package turned into a function; results are written to a `Paper_Exhibits_Results_SEC7_*` folder.
* `runSection7_Table1.m`: Table 1, $F_T^*$ tests for identification failure (30-minute policy news and 1-day 2-Year yields, $\pi_L \in \{0.6,\dots,1\}$).
* `runSection7_Table2.m`: Table 2 (TSLS and GMM estimates of $\beta$, full sample and $\hat{\pi}$-sample) and Figure 3 (strongly-identified subsample). The GMM rows use the replication code of Lewis (2022), which is not distributed here: copy its `estimation` folder into `Section7/external/Lewis2022/` (or pass its path as second argument); without it the table is produced with the TSLS rows only.
* `runSection7_Table3.m`: Table 3, weak identification-robust inference on $\beta$ (LM and CLR, full sample and $\hat{S}_T$, three significance levels).
* `runSection7_Compliers.m`: Section 7.4, identification of compliers (Figure 4, share of compliers) and test of the exclusion restriction.
* `runSection7_Figure1.m`: Figure 1, survey of first-stage $F$ and $F^*$ statistics (data in `Section7/data/`).
* `runSection7_Figure2.m`: Figure 2, break dates in the 2-Year yields from the Casini-Perron (2024) All-Inside procedure, an external dependency: copy the folder `All-Inside_Matlab` into `Section7/external/AllInside/` (or pass its path); without it the published dates are plotted.
* `legacy_IRI/`, `legacy_mbreak/`: helper functions of the replication package used verbatim (`S_Dates_m_M2_Heter`, `LR_critical_value`, `F_ST_Stat`, `TSLS_Robust`).

**Example of usage:**
```matlab 
cd Section7
runSection7_Table1(); runSection7_Table2(); runSection7_Table3(); runSection7_Compliers(); runSection7_Figure1(); runSection7_Figure2();
```

### Exhibit numbering: arXiv v1 vs. v2
The exhibit numbers used throughout this repository (folder READMEs, targets such as `'tab10'`, export headers, result file names) are those of arXiv v1. The submitted version (arXiv v2) moved some exhibits and renumbered them; the correspondence is:

| Repository (v1) | arXiv v2 | Where in v2 |
|---|---|---|
| Tables 1-3, Figures 1-4 (Sections 1 and 7) | unchanged | main text |
| Table 4 (critical values of $F_T^*$) | Table 4 | online supplement, Section S.A |
| Table 5, Figures 5-6 (Section 4) | unchanged | online supplement, Section S.C.1 |
| Tables 6-8 (Section 5) | Tables 2-4 | supplement not for online publication, Section N.F.1 |
| Tables 9-10, Figures 7-8 (Section 6) | Tables 6-7, Figures 7-8 | online supplement, Section S.C.2 |
| "Supplement Table 2", "Supplement Figure 1" (AR(1) errors; code comments say Section N.D) | Table 5, Figure 1 | supplement not for online publication, Section N.F.2 |
| "Supplement Table 3", "Supplement Figure 2" (three instruments, model (N.1)-(N.2)) | Table 6, Figure 2 | supplement not for online publication, Section N.F.2 |
