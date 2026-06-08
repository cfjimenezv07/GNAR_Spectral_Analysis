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

#### 1. Setup, Data Preprocessing, and Network Topology Construction
* [cite_start]**Aux_GNAR_spec.R**: Installs and loads required packages and libraries for managing time series with network-structured dependence (TS-NSD)[cite: 8, 31].
* [cite_start]**read_and_format_data.R**: Loads the raw global banking log-volatility datasets [cite: 220, 222] [cite_start]and processes them into structured network time series objects[cite: 63].
* [cite_start]**construct_network.R**: Constructs the underlying network adjacency matrices using empirical features such as Generalized Forecast Error Variance Decomposition (GFEVD) via a Lasso-VAR framework[cite: 222, 224, 225].

#### 2. Spectral Density Estimation
* [cite_start]**net_spec_density.R**: Implements the main function to evaluate network time series spectral density matrices across Fourier frequencies[cite: 104, 124].
* [cite_start]**parametric_estimation.R**: Computes the parametric spectrum by embedding model parameters directly into a $d$-dimensional network VAR coefficient representation (e.g., global GNAR or NAR settings)[cite: 78, 83, 102].
* **nonparametric_estimation.R**: Constructs raw and smoothed Fourier periodogram matrices to offer data-driven, model-agnostic spectral alternatives[cite: 121, 124, 125, 156].
* **constrained_optimization.R**: Enforces known physical topological or multi-hop adjacency constraints onto real-valued augmented spectral representations using matrix optimizations[cite: 159, 172, 174].

#### 3. Network Coherence and Interdependence Metrics
* **network_coherence.R**: Calculates squared coherence to determine frequency-specific cross-nodal dependencies across pairs of institutions[cite: 106].
* **partial_coherence.R**: Evaluates partial coherence by inverting the spectral density matrix to separate direct network interactions from multi-stage/multi-hop pathways[cite: 3, 106].

#### 4. Empirical Application & Volatility Transmission Visualizations
* **global_bank_analysis.R**: Executes the complete financial empirical study monitoring systemic risk and frequency-dependent volatility co-movements among 57 major international banks[cite: 14, 230].
* [cite_start]**plot_networks.R**: Generates visualizations of regional clusters, multi-stage neighborhood structures, phase relationships, and localized volatility transmission plots across variable time horizons[cite: 229, 237].

## Contact
[cite_start]**arXiv link:** https://arxiv.org/abs/2510.06157 

[cite_start]**Cristian F. Jiménez-Varón** - cristian.jimenezvaron@york.ac.uk [cite: 1, 9]

[cite_start]**Marina I. Knight** - marina.knight@york.ac.uk [cite: 1, 9]

## Acknowledgements
[cite_start]The authors gratefully acknowledge support from the EPSRC NeST Programme Grant EP/X002195/1.

<br />
<div align="center">
  <a href="https://github.com/cfjimenezv07/Network_Spectral_Analysis"><strong>Explore R code »</strong></a>
</div>
