<br />
<div align="center">
  <a href="https://github.com/cfjimenezv07/Network_Spectral_Analysis">
    <img src="UoY.png" alt="York Logo" height="150">
    <img src="Nest.png" alt="Nest Logo" height="150">
  </a>

<h3 align="center">Frequency-Domain Analysis of Time Series with Network-Structured Dependence: Application to Global Bank Connectedness</h3>
</div>

## Abstract
<p align="justify">
Financial spillovers in interconnected systems, such as global banking networks, require tools that capture temporal and frequency dynamics, while incorporating the underlying network topology. While current network time series models are developed in the time-domain, frequency-domain approaches, which reveal how cross-nodal dependencies vary across different cycles, remain under-explored. This paper develops a spectral analysis framework that accommodates flexible forms of network dependence, including interactions mediated through intermediate nodes. This ensures that inter-nodal relationships are not restricted to direct connections, a feature crucial for capturing indirect financial spillovers. We define the network time series spectral density, alongside coherence and partial coherence, and propose both parametric and network-constrained nonparametric methods for their estimation. Simulations and theoretical results demonstrate the strong performance of the parametric approach when the data-generating process aligns with the model structure, whereas the nonparametric alternative provides robustness against model misspecification. An application to global bank connectedness shows that the proposed spectral measures capture inter-bank frequency-specific spillover effects, yielding results consistent with existing measures while additionally uncovering richer patterns of volatility transmission that are intimately connected to the network topology.
</p>

### Main Results
The R script files in the `R Code` folder should be used in the following order:

#### 1. Core Functions & Methodologies
* **Aux_GNAR_spec.R**: Auxiliary functions for the computation of spectral quantities such as parametric and nonparametric estimators, as well as tailor-made functions for the global bank network connectedness application.

#### 2. Network Analysis & Empirical Application
* **Bank_N_conn.R**: Handles the bank network connectedness setup and structural modeling framework.
* **Threshold_estimating.R**: Threshold estimation for the global bank network connectedness using the network of Demirer et al. (2018) published in the *Journal of Applied Econometrics* (JAE).

#### 3. Simulation Frameworks
* **Simulation_5N_N.R**: Simulation framework implemented for a 5-node network configuration under a correctly specified data-generating process.
* **Simulation_5N_N_MM.R**: Simulation framework for a 5-node network incorporating model misspecification, with model selection optimized via BIC.
* **Simulation_10N_N.R**: Extended simulation setup implemented for a larger 10-node network structure.
* **Simulation_10N_N_MM.R**: Extended 10-node simulation framework accounting for model misspecification.

#### 4. Plotting & Visualization
* **Plots_Estimation_methods.R**: Generates empirical figures comparing the performance of the parametric and nonparametric spectral estimation methodologies.
* **Plots_r_spec.R**: Visualizes the $r$-stage neighborhood spectral density responses and multi-stage frequency-domain connectedness patterns.

## Contact
**arXiv link:** https://arxiv.org/abs/2510.06157

**Cristian F. Jiménez-Varón** - cristian.jimenezvaron@york.ac.uk

**Marina I. Knight** - marina.knight@york.ac.uk

## Acknowledgements
The authors gratefully acknowledge support from the EPSRC NeST Programme Grant EP/X002195/1.

<br />
<div align="center">
  <a href="https://github.com/cfjimenezv07/GNAR_Spectral_Analysis"><strong>Explore R code »</strong></a>
</div>
