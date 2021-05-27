# MetapopulationTrade-Model
Simulation and visualization code for manuscript by Fair, Bauch, and Anand (2021). Model of SARS-CoV-2 transmission used to explore how physical factors, human behavior, and non-pharmaceutical interventions (NPIs) impact the persistence of novel pathogens. Associated preprint is posted at https://doi.org/10.1101/2021.05.03.21256551.

To run a set of simulations or generate a figure, make all required option selections (detailed within each script) and launch the appropriate script.

All data needed to run visualization scripts are included in the "InputFiles" folder, except for inputs to COVIDviz_singlepatch_birthdeathexp.R, due to file size restrictions. These can be generated using covidHier0.35Lite_singlepatch_birthdeathexp.R.

All code has a GNU GPLv3 license (see CODE_LICENSE for details), and all data has a CC-BY-4.0 license (see DATA_LICENSE for details).

## Scripts

* **covidHier0.35Lite_singlepatch_scenarios.R** runs simulations for different NPI scenarios, visualizes results and outputs figure

* **covidHier0.35Lite_singlepatch_birthdeathexp.R** runs simulations for experiment on the impact of accounting for births/deaths and outputs results of these simulations

* **COVIDviz_singlepatch_birthdeathexp.R** visualizes results from covidHier0.35Lite_singlepatch_birthdeathexp.R and outputs figure

* **covidHier0.35Lite_singlepatch_parametervary.R** runs simulations for experiment varying a single parameter at a time and outputs results of these simulations

* **COVIDviz_singlepatch_parametervary.R** visualizes results from covidHier0.35Lite_singlepatch_parametervary.R and outputs figure


## Input files

* **COVIDsim_CCS_fadeoutdat_X.R**, where X is a placeholder for the parameter being varied, contains data on fade-outs and infections from simulations generated using covidHier0.35Lite_singlepatch_parametervary.R

* **COVIDsim_CCS_ccsdat_X.R**, where X is a placeholder for the parameter being varied, contains data on critical community size values from simulations generated using covidHier0.35Lite_singlepatch_parametervary.R


## System & hardware requirements

 * Windows 10 Pro Version 2004

 * Sufficient RAM to support storage of data structures during simulations

## Installation guides

All software should install within a few minutes on a standard computer, the versions listed here are those the scripts have been tested on.

 * R Version 4.0.3 https://www.r-project.org/

 * R Studio Version 1.2.5019 (IDE for R) https://rstudio.com/ 

