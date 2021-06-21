# MetapopulationTrade-Model
Simulation and visualization code for manuscript by Fair, Bauch, and Anand (2021). Model describes a system where agricultural land area, food supply, and human population growth are coupled. Associated preprint is posted at ???.

To run a set of simulations make all required option selections (detailed within each script) and launch the appropriate script.

All data needed to run scripts are included in the "InputFiles" folder.

## Scripts

* **Metapop_Globalfit.m** runs global paramter fitting under low & high yield scenarios.

* **Metapop_Globalfit_Viz.m** visualizes output from Metapop_Globalfit.m

* **Metapop_ParameterPlanes.m** runs simulations sampling different points in the parameter plane for either (beta, gamma) or (beta_A, beta_I)

* **Metapop_ParameterPlanes_Viz.m** visualizes output from Metapop_ParameterPlanes.m

* **Metapop_NetworkExp.m** runs simulations on networks with different densities and rewiring probabilities.

* **Metapop_NetworkExp_Viz.m** visualizes output from Metapop_NetworkExp.m

* **Metapop_HeterogeneityResponsiveness.m** runs simulations with heterogeneous patch responsiveness (i.e. heterogeneous gamma-values)

* **Metapop_HeterogeneityImportThreshold.m** runs simulations with heterogeneous patch import demand thresholds (i.e. heterogeneous beta_I-values)

* **Metapop_HeterogeneityExp_Viz.m** visualizes output from Metapop_HeterogeneityResponsiveness.m and Metapop_HeterogeneityImportThreshold.m

* **Metapop_HeterogeneityExp_betaprop_Viz.m** visualizes output from Metapop_HeterogeneityImportThreshold.m across a broader range of values for the percent of patches with high (vs. low) import demand thresholds

## Input files

* **Meta_FittingTS_Global_Feb1.csv** contains data from the United Nations on global population size (https://population.un.org/wpp/Download/Standard/Population/), food supply, agricultural land area, total land area, etc. (http://www.fao.org/faostat/en/#home)

* **Ky_35_params_NOurban_Feb5.csv** contains parameters from global fitting under a low yield scenario

* **Ky_70_params_NOurban_Feb5.csv** contains parameters from global fitting under a high yield scenario

* **Networks_forParameterPlanes** folder contains networks required to run Metapop_ParameterPlanes.m

* **Networks_forNetworkExperiment** folder contains networks required to run Metapop_NetworkExperiment.m

* **Networks_forHeterogeneityExperiments** folder contains networks required to run Metapop_HeterogeneityEfficiency.m and Metapop_HeterogeneityImportDemand.m

* **NetworkStatistics_forNetworkExperiment** folder contains statistics for characterising networks, required to run Metapop_NetworkExp_Viz.m

* **NetworkStatistics_forHeterogeneityExperiments** folder contains statistics for characterising networks, required to run Metapop_HeterogeneityExp_Viz.m and Metapop_HeterogeneityExp_betaprop_Viz.m

## System & hardware requirements

 * Windows 10 Pro Version 2004

 * Sufficient RAM to support storage of data structures during simulations

## Installation guides

All software should install within a few minutes on a standard computer, the versions listed here are those the scripts have been tested on.

 * MATLAB R2016b https://www.mathworks.com/products/matlab.html
