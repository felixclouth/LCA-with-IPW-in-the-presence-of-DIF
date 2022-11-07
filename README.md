# LCA-with-IPW-in-the-presence-of-DIF
Accompanying code for the paper: Three-step latent class analysis with inverse propensity weighting in the presence of differential item functioning

Results for the simulated data example can be reproduced using the .R file. Within this script, syntax files (.lgs) are generated. Latent GOLD (https://www.statisticalinnovations.com/latent-gold-6-0/) is called from the R script to run these syntax files using the example.sav file as input to simulate data and analyze it.

The data for the real-life example can be downloaded here (https://www.nlsinfo.org/content/cohorts/nlsy79). The Data Selection.NLSY79 file can be uploaded on this website to obtain the selection of variables needed for the analysis. All necessary step to prepare the data for analysis are done in the .R file. For the analysis, first, the Step one.lgs files needs to be run in Latent GOLD. This file estimates the measurement model and generates a dataset (.sav file) containing the posterior class membership probabilities and classifications. The Step three.lgs file then calls this new data file to estimate the structural part of the model.
