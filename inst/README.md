This folder contains supplemental code for the paper `Robust marginal functional principal component analysis for repeated measurements data" (Brune, Radojicic, Filzmoser, 2024).

The provided code requires the `robLFDA ` R-package as available on GitHub. It can be installed by running `devtools::install_github("b-brune/robLFDA")` in the R console.

# Reproduction code for data examples in the paper:

* reproduction_file_mortality.R; reproduces the analysis and figures for the mortality application
* reproduction_file_uvf.R; reproduces the analysis and figures for the analysis of photovoltaic spectra

# Simulation results

The simulations were run using the `batchtools` framework in R. If you aim to rerun the simulation study, this can be achieved using the following two files. However, the simulation will run for multiple days on > 80 Cores, thus we also provide the reduced RData files with the results.

* simulation_registry_setup.R; contains the code to set up the experiment registry for the simulation study
* preprocess_batchtools_output.R; processes and `unwraps` the batchtools output (nested dataframes)

Additionally, those two files load 

* DATA_GENERATION.R; contains the code to generate the datasets used for evaluation of the method
* evaluation_functions.R; utility functions to evaluate the output
 
The **preprocessed simulation results** can be found in 

* simulation_amplitude_outliers.Rdata
* simulation_structural_outliers.Rdata

By loading those two files into the file `plots_paper.R` the graphics reported in the paper can be reproduced.