# Prediction and Inference With Missing Data in Patient Alert Systems

# Author Contributions Checklist Form

## Data

### Abstract

34923 observations on 171 predictor variables (vital signs, lab results, nursing assessments, demographics,…).  Response is time until adverse event, which is right censored for most observations.  Many variables are missing for many observations.

### Availability

Data cannot be made available.  It is electronic health data and thus must be protected.
Code to generate data from the simulation cases in the paper has been made available.


## Code

### Abstract 

R code with the main function “weibull.MCMC” to conduct the estimation of the model described in the paper, along with many other supporting functions.

### Description 

The file sim_cases.R is an R script that can produce data from the simulation cases 1-4 in the paper, and shows the use of the estimation function “weibull.MCMC”.

R session information for the version of R and packages used by the authors
can be found below

> sessionInfo()
R version 3.4.3 (2017-11-30)
Platform: x86_64-apple-darwin15.6.0 (64-bit)
Running under: macOS Sierra 10.12.6

Matrix products: default
BLAS: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRblas.0.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] splines   parallel  stats     graphics  grDevices utils     datasets 
[8] methods   base     

other attached packages:
 [1] truncnorm_1.0-8     MASS_7.3-49         missForest_1.4     
 [4] itertools_0.1-3     randomForest_4.6-12 glmnet_2.0-13      
 [7] Matrix_1.2-12       gbm_2.1.3           lattice_0.20-35    
[10] doMC_1.3.5          iterators_1.0.9     foreach_1.4.4      
[13] mclust_5.4          survival_2.41-3    

loaded via a namespace (and not attached):
[1] codetools_0.2-15 grid_3.4.3       compiler_3.4.3  


## Instructions for use

Below we provide instructions on how to (i) analyze a synthetic data set that
is similar to the BPR ata set used in the manuscript, and (ii) run the
simulation cases used in the manuscript and summarize the results.


### Example Analysis Results

In a directory that houses "analyze_example.R", "weibull_response_sparse_DPM.R",
and "jasa_data.RData", source the file "analyze_example.R" to load in the
synthetic example data and fit the sDPM_MNAR model.

NOTE: The results will be similar but not identical to those in the manuscript
because the actual BPR data can not be made available publicly.  Instead a
synthetic data set intended to mimic that of the BPR data has been generated.
It has mostly the same columns as the original data, variables are scaled to
have mean 0, varinace 1, but otherwise represent the distirbution of the actual
BPR data.

### Simulation Results 

To run the simulation cases used in the manuscript, do the following...
At the command line in the directory that houses "sim_cases.R" and
"weibull_response_sparse_DPM.R", execute the command

mkdir Results

and then the command

./sim_cases.R A B

for all combinations of A=1,2,3,4, and B=0,2,...,95.  'A' designates the
simulation case, and 'B' designates the seed to set at the beginning of the
sim_cases.R file before it generates the data.

NOTE: The runs for cases 5,6,7,8, will be similar but not identical to those in
the manuscript because the actual BPR data can not be made available publicly.
Instead a synthetic data set intended to mimic that of the BPR data has been
generated

NOTE: Each run takes several hours.  You will need to parallelize these runs
to make this feasible.  We ran all 8 x 96 = 768 combinations with the
following commands

seq 0 96 | xargs -P32 -I% -n1 ./sim_cases.R 1 %
seq 0 96 | xargs -P32 -I% -n1 ./sim_cases.R 2 %
seq 0 96 | xargs -P32 -I% -n1 ./sim_cases.R 3 %
seq 0 96 | xargs -P32 -I% -n1 ./sim_cases.R 4 %
seq 0 96 | xargs -P32 -I% -n1 ./sim_cases.R 5 %
seq 0 96 | xargs -P32 -I% -n1 ./sim_cases.R 6 %
seq 0 96 | xargs -P32 -I% -n1 ./sim_cases.R 7 %
seq 0 96 | xargs -P32 -I% -n1 ./sim_cases.R 8 %

Each of these 768 executions of the "/sim_cases.R" code will generate a file
named "CaseSeed_A_B.RData" in the local folder ./Results. You can then source
the file "summarize_sims.R" to create the resulting tables and plots from the
simulation study section of the manuscript.
