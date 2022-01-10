# tensile-strength-model
This repository contains the data code to reproduce the results in the manuscript ``Knots and their effect on the tensile strength of lumber'' by Shuxian Fan, Samuel WK Wong, and James V Zidek.

- To process the experimental data, run `data_process.R`
- To generate simulated pieces and strengths and fit the Bayesian model, run `fit_simulated.R`
- To fit the Bayesian model to the Douglas Fir dataset, run `fit_realdata.R`
- To apply K-fold cross validation using the Douglas Fir dataset, run `fit_realdata_cv.R`. Then `predict-cv-fits.R` will compute prediction metrics for comparison with baseline regression models
- To summarize the model fits and make tables, run `make_tables.R`
