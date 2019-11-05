
########################
### Read in the data ###
########################

## The file "jasa_data.RData" contains the object "jasa.data" which is a
## list of the data objects
load("jasa_data.RData")

## X is the design matrix.  The first 89 variables are covariates, the remaining
## 34 are indicators of missingness for the 34 of those 89 variables that have
## missing values.
X <- jasa.data$X

## y is the time to event (or censoring)
y <- jasa.data$y

## y is the outcome (event or not)
outcome <- jasa.data$outcome

## bin is an indicator vector to indicate whether each column of X
## is binary variable or not
bin <- jasa.data$bin

## discrete is an indicator vector to indicate whether each column of X
## is a discrete (but not binary) variable or not
discrete <- jasa.data$discrete

## limits is a matrix with two columns (min, max) that contains hard limits
## for some of the variables. The default is to just set these all to (-Inf, Inf).
limits <- jasa.data$limits


## Hold out some data for validation ##
set.seed(22)
n <- nrow(X)
np <- 2000
ind.p <- sample(1:n, np, replace=FALSE)

## Held out observations
Xp <- X[ind.p,]
yp <- y[ind.p]
outcomep <- outcome[ind.p]

## Training observations
Xo <- X[-ind.p,]
yo <- y[-ind.p]
outcomeo <- outcome[-ind.p]


## Make a training data set that includes Xp but with y's censored at 0.
## That is, these observations will provide no information for training.
## We only want them included in order to make predictions at their X values.
Xa <- rbind(Xo, Xp)
ya <- c(yo,rep(0,np))
outcomea <- c(outcomeo,rep(0,np))




############################################
### Defining settings for MCMC algorithm ###
############################################

## Number of MCMC iterations to run.  Set here to only 2000 to save time.
## Better to use more (e.g., ~10,000+) for an actual analysis.
N.mcmc <- 2000

## Begin recording (stop burning in) after this many iterations
begin <- 1000

## Record in a posterior sample every 'every' iterations
every <- 2

## XZ updates take much longer than the other updates.  Update XZ only
## every 'everyXZ' iterations.
everyXZ <- 2

## Plot MCMC progress, trace plots, etc., every 'nplot' iterations.
nplot <- 20


#############################################
### Fit MVN-MAR model and sDPM-MNAR model ###
#############################################

## Load packages
library("survival")
library("mclust")
library("doMC")
library("glmnet")
library("truncnorm")

source("weibull_response_sparse_DPM.R")

## Recommend using multiple cores to speed up the XZ updates
registerDoMC()
options(cores=8)

## Set 'rho', a vector of the prior probabilities of variable inclusion in the
## regression model.  Defaults to intercept~1 and all variables equal 0.5.
rho <- .5

## There are pr=89 variables that are not missing indicators.  The remaining
## 34 variables (missing indicators) are not used in this MAR model.
pr <- 89

## Run the MVN-MAR model.  This may take a few hours to complete.
set.seed(12)
ans.MVN.MAR <- weibull.MCMC(X=Xa[,1:pr], y=ya, outcome=outcomea, rho=rho, N.mcmc=N.mcmc,
every=every, nplot=nplot, bin.X=bin[1:pr], discrete.X=discrete[1:pr],
limits.X=limits[1:pr,], everyXZ=everyXZ, M=1, yp=yp, outcomep=outcomep)



## Run the sDPM-MNAR model.  This may take a few hours to complete.

## Set the maximum number of allowed mixture components for the sparse DPM
Mmod <- 20

## Set the possible number of components for initialization (defaults to 1:Mmod)
G.vec <- 3:8

## Here we will use all 89+34 columns of X, but we set the rho elements for
## indicators for variable missingness to all have ~0 probability of being part
## of the regression.
rho <- c(1-1E-15,rep(.5,pr),rep(1E-200,ncol(Xa)-pr))

## Use beta estimate from MVN-MAR to initialize beta
beta.init <- c(apply(ans.MVN.MAR$beta,2,mean)*(apply(abs(ans.MVN.MAR$beta)>0,2,mean)>.1),
rep(0,ncol(Xa)-pr))

## Fit model
set.seed(12)
ans.sDPM.MNAR <- weibull.MCMC(X=Xa, y=ya, outcome=outcomea, rho=rho, N.mcmc=N.mcmc,
every=every, nplot=nplot, bin.X=bin, discrete.X=discrete, limits.X=limits,
everyXZ=everyXZ, M=Mmod, G=G.vec, begin=begin, yp=yp, outcomep=outcomep,
beta.init=beta.init)



## Evaluate hold-out (posterior mean) predictions 
pred.sDPM.MNAR <- apply(ans.sDPM.MNAR$risk[,(length(yo)+1):length(y)],2,mean)

## calculate survival concordance
survConcordance(Surv(yp, event=outcomep)~pred.sDPM.MNAR)$concord

## print out posterior summary of regression coefficients
beta <- ans.sDPM.MNAR$beta[,2:90]

beta.post <- as.data.frame(cbind(apply(beta,2,mean),t(apply(beta,2,quantile,c(.025,.975)))))
rownames(beta.post) <- colnames(Xa)[1:89]
colnames(beta.post) <- c("post mean", "LL.95", "UL.95")
print(beta.post)


