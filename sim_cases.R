#! /usr/local/bin/Rscript


## Description of Cases ##
#
# All cases 50/50 continuous/binary, (n,p,q)=(2500,30)
#
# Case 1 - x ~ MVN, 60% MAR x1,x16,  20% MAR rest
#
# Case 2 - x ~ DPM, 60% MAR x1,x16,  20% MAR rest
#
# Case 3 - x ~ MVN, 60% MNAR x1,x16, 20% MNAR x2,x3,x17,x18, 20% MAR rest
#
# Case 4 - x ~ DPM, 60% MNAR x1,x16, 20% MNAR x2,x3,x17,x18, 20% MAR rest
#
# Case 5 - x ~ BPR, 60% MAR x1,x16,  20% MAR rest
#
# Case 6 - x ~ BPR, 2 predictors removed, 60% MAR x1,x16, 20% MAR rest
#
# Case 7 - x ~ BPR, 60% MNAR x1,x16, 20% MNAR x2,x3,x17,x18, 20% MAR rest
#
# Case 8 - x ~ BPR, 2 predictors removed, 60% MNAR x1,x16, ...
#
##########


library("survival")
library("mclust")
library("doMC")
library("gbm")
library("glmnet")
library("missForest")
library("MASS")
library("truncnorm")
library("methods")

source("weibull_response_sparse_DPM.R")
registerDoMC()
options(cores=8)
nb <- options()$cores


args <- commandArgs(TRUE)
case <- as.numeric(args[1])
seed <- as.numeric(args[2])


set.seed(seed)


beta <- c(0,2.0,1.0,-1.0,rep(0,12),3,1.0,-1.0,rep(0,12))*.5
rho <- .2
p <- length(beta)-1
bin.X <- rep(0,p)
bin.X[1:(p/2)] <- 1
discrete.X <- rep(0,p)
limits.X <- cbind(rep(-Inf,p),rep(Inf,p))
p.miss <- rep(.2,p)
p.outcome <- .2
kappa <- 2
df.wish <- p+5
lam.mu <- 1
gamma <- c(rep(1,3),rep(0,p/2-3),rep(1,3),rep(0,p/2-3))


if(case==1){
  omega <- NULL
  muX <- NULL
  inform.col <- c()
  Mmax <- 1
  p.miss[c(1,(p/2)+1)] <- .6
}

if(case==2){
  omega <- NULL
  muX <- NULL
  inform.col <- c()
  Mmax <- 10
  p.miss[c(1,(p/2)+1)] <- .6
}

if(case==3){
  inform.col <- c(1:3,(p/2)+1:3)
  omega <- NULL
  muX <- NULL
  Mmax <- 1
  p.miss[inform.col] <- .25
  p.miss[1] <- .4
  p.miss[(p/2)+1] <- .6
  p.miss[(p/2)+2] <- .25
  df.wish <- p+15
}

if(case==4){
  inform.col <- c(1:3,(p/2)+1:3)
  omega <- NULL
  muX <- NULL
  Mmax <- 10
  p.miss[inform.col] <- .25
  p.miss[1] <- .4
  p.miss[(p/2)+1] <- .6
  p.miss[(p/2)+2] <- .25
  df.wish <- p+15
}

if(case==5 || case==6){
  inform.col <- c()
  omega <- NULL
  muX <- NULL
  Mmax <- 10
  p.miss <- rep(.1,p)
  p.miss[1] <- .4
  p.miss[(p/2)+1] <- .6
  p.miss[c(2,3,(p/2)+2:3)] <- .2
}

if(case==7 || case==8){
  inform.col <- c(1:3,(p/2)+1:3)
  omega <- NULL
  muX <- NULL
  Mmax <- 1
  p.miss <- rep(.1,p)
  p.miss[1] <- .4
  p.miss[(p/2)+1] <- .6
  p.miss[c(2,3,(p/2)+2:3)] <- .2
}


#############################
####### Generate Data #######
#############################


nt <- 2500
np <- 500


## Fixed Params
delta <- 2
b.miss <- 1.5*rep(1,length(inform.col))
s.miss <- rep(1,length(inform.col))
Mmod <- 10
n <- nt + np



## Generate Data
omega <- 1
M <- 1
while(M<Mmax && M<3){
  v <- 1-c(rbeta(Mmax-1,1,delta),0)
  v <- c(rbeta(Mmax-1,1.5,.2),0)
  omega <- makeprobs(v)
  omega <- sort(omega[omega>.025], decreasing=TRUE)
  omega <- omega/sum(omega)
  M <- length(omega)
}
M
omega





#### Generate mu and Sigma from varible selection prior ####

c.vars <- (gamma==1)
nc.vars <- (gamma==0)
p1 <- sum(c.vars)
p2 <- sum(nc.vars)

if(is.null(muX)){
  muX <- nuX <- matrix(0,M,p)
}

Psi <- (df.wish-p-1)*diag(1,p)
Psi[c(p/2+1,p/2+2),c(p/2+1,p/2+2)] <- (df.wish-p-1)*matrix(c(1,-.6,-.6,1),2,2)
SigmaX <- QX <- RX <- array(0,c(M,p,p))
Q <- rwish(df.wish,ginv.gp(Psi)$inv)
Q22 <- Q[nc.vars,nc.vars]
Q21 <- Q[nc.vars,c.vars]
nu <- rmvnorm(1, rep(0,p), 1/lam.mu*ginv.gp(Q)$inv)
b2 <- Q21%*%nu[c.vars]+Q22%*%nu[nc.vars]
ans.Q22 <- ginv.gp(Q22)

M.x <- rep(0,p)
for(m in 1:M){
  Psi.m <- Psi
  if(m==3)
    Psi.m[c(p/2+1,p/2+2),c(p/2+1,p/2+2)] <- (df.wish-p-1)*matrix(c(1,.6,.6,1),2,2)
  S.m11.inv <- rwish(df.wish-p2,ginv.gp(Psi.m[c.vars,c.vars])$inv)
  Q.m11 <- S.m11.inv + t(Q21)%*%ans.Q22$inv%*%Q21
  QX[m,c.vars,c.vars] <- Q.m11
  QX[m,c.vars,nc.vars] <- t(Q21)
  QX[m,nc.vars,c.vars] <- Q21
  QX[m,nc.vars,nc.vars] <- Q22
  SX.m <- ginv.gp(QX[m,,])$inv
  M.x[c(p/2+1,p/2+2)] <- (m==1)*c(3*sqrt(SX.m[p/2+1,p/2+1]),-3*sqrt(SX.m[p/2+2,p/2+2])) +
                         (m==2)*c(-3*sqrt(SX.m[p/2+1,p/2+1]),0*sqrt(SX.m[p/2+2,p/2+2])) +
			 (m==3)*c(6*sqrt(SX.m[p/2+1,p/2+1]),4*sqrt(SX.m[p/2+2,p/2+2]))
  nuX[m,c.vars] <- rmvnorm(1,M.x[c.vars], 1/lam.mu*ginv.gp(S.m11.inv)$inv)
  nuX[m,nc.vars] <- M.x[nc.vars]+ans.Q22$inv%*%(b2-Q21%*%(nuX[m,c.vars]-1/lam.mu*M.x[c.vars]))
  scale.m <- diag(1/sqrt(diag(QX[m,,])))
  diag(scale.m)[bin.X==0] <- 1
  RX[m,,] <- scale.m%*%QX[m,,]%*%t(scale.m)
  muX[m,] <- nuX[m,]*(1/diag(scale.m))
  muX[m,bin.X==1] <- nuX[m,bin.X==1]*.8
  SigmaX[m,,] <- ginv.gp(RX[m,,])$inv
}
for(m in 1:M)
print(summary(abs(cov2cor(SigmaX[m,,])[upper.tri(SigmaX[1,,])])))
print(max(abs(cov2cor(SigmaX[m,,]))[-1,1]))
print(max(abs(cov2cor(SigmaX[m,,]))[upper.tri(SigmaX[1,,])]))


## generate X
phi <- sample(1:M, n, replace=TRUE, prob=omega)
Z <- matrix(0,n,p)
for(m in 1:M){
  ind.m <- which(phi==m)
  Z[ind.m,] <- rmvnorm(length(ind.m),muX[m,],SigmaX[m,,])
}

if(0){
par(mfrow=c(floor(sqrt(p)),ceiling(sqrt(p))))
cor(Z[,c(p/2+2,p/2+1)])
for(j in 1:p){
  if(j==(p/2+2)){
    mu.foo <- apply(Z[,c(p/2+2,p/2+1)],2,mean)
    cov.foo <- cov(Z[,c(p/2+2,p/2+1)])
    Z.foo <- rmvnorm(10000,mu.foo,cov.foo)
    plot(rbind(Z.foo,Z[,c(p/2+2,p/2+1)]),col=0)
    points(Z.foo,col=rgb(.5,.5,.5,alpha=.4))
    points(Z[,c(p/2+2,p/2+1)],col=phi)
  }
  else
    plot(Z[,c(j,p/2+1)],col=phi, main=paste("X",j,"by X",p/2+1,sep=""))
}
}





## create dichotomous and scale the X's to mean 0, sd .5
Xc <- Z
for(j in which(bin.X==1))
  Xc[,j] <- as.numeric(Z[,j]>0)
apply(Xc[,bin.X==1],2,mean)

Xbar <- sdX <- rep(0,p)
mult <- rep(1,p)
for(j in which(bin.X==0)){
  Xbar[j] <- mean(Xc[,j])
  sdX[j] <- sd(Xc[,j])
  Xc[,j] <- (Xc[,j] - Xbar[j])/sdX[j]*.5
  muX[,j] <- (muX[,j] - Xbar[j])/sdX[j]*.5
  mult[j] <- 1/sdX[j]*.5
}
for(m in 1:M)
  SigmaX[m,,] <- diag(mult)%*%SigmaX[m,,]%*%diag(mult)
Psi <- diag(mult)%*%Psi%*%diag(mult)



if(case>=5){
  load("jasa_data.RData")
  X.jasa <- jasa.data$X
  bin.jasa <- jasa.data$bin
  discrete.jasa <- jasa.data$discrete
  limits.jasa <- jasa.data$limits
  avail.cols <- 1:89
  ind.90 <- apply(is.na(X.jasa[,avail.cols]), 2, mean) < .9
  avail.cols <- avail.cols[ind.90]
  repeat{
    cols <- sample(avail.cols, p, replace=FALSE)
    X.temp <- X.jasa[,cols[c(1:3,p/2+1:3)]]
    avail.rows <- which(apply(!is.na(X.temp),1,prod)==1)
    if(length(avail.rows) < 1000)
      next
    rows <- sample(avail.rows, n, replace=TRUE)
    Xc <- X.jasa[rows, cols]
    bin.X <- bin.jasa[cols]
    discrete.X <- discrete.jasa[cols]
    limits.X <- limits.jasa[cols,]
    break
  }
  rownames(Xc) <- NULL
  Xc[is.na(Xc)] <- NA
length(avail.rows)
}





## create some missing entries
X <- Xc
for(j in 1:p){
  if(any(inform.col==j)){
    k <- which(inform.col==j)
    z.j <- qnorm(p.miss[j]) + (b.miss[k]*scale(Z[,j]) +
           s.miss[k]*rnorm(n,0,1))/sqrt(b.miss[k]^2 + s.miss[k]^2)
    X[z.j>0,j] <- NA 
  }    
  else{
    ind.na.j <- rbinom(n,1,p.miss[j])==1
    X[ind.na.j,j] <- NA
  }
}

apply(is.na(X),2,mean)



## generate y
Xc[is.na(Xc)] <- 0
risk <- exp(as.numeric(cbind(1,Xc[,1:p])%*%beta))
y <- yc <- rweibull(n, kappa, 1/risk)
T <- quantile(y,p.outcome)
outcome <- as.numeric(y<T)
y[!outcome] <- T
mean(outcome)
survConcordance(Surv(y, event=outcome)~risk)$concord


if(case==6 || case==8){
  X <- X[,-c(2,17)]
  p <- p-2
  bin.X <- bin.X[-c(2,17)]
}




if(0){
par(mfrow=c(2,2))
for(jj in 1:4){
  boxplot(Z[,inform.col[jj]]~is.na(X[,inform.col[jj]]))
  print(get.ROC(is.na(X[,inform.col[jj]]), Z[,inform.col[jj]])$auc)
  print(table(is.na(X[,inform.col[jj]]), outcome))
}
}



## add missingness indicator columns
for(j in 1:p){
  foo.j <- as.numeric(is.na(X[,j]))
  if(any(foo.j==1))
    X <- cbind(X,foo.j)
}
bin.X <- c(bin.X, rep(1,ncol(X)-p))
discrete.X <- c(discrete.X, rep(0,ncol(X)-p))
limits.X <- rbind(limits.X, cbind(rep(-Inf,ncol(X)-p),rep(Inf,ncol(X)-p)))

## create missing ind and levels for Factors for glm
Xd <- as.data.frame(X[,1:p])
for(j in 1:(p/2))
  Xd[,j] <- as.factor(Xd[,j])
Xf <- add.missing.ind.cols(Xd)

## hold out np obs for prediction
Xp <- X[(n-np+1):n,]
Xfp <- Xf[(n-np+1):n,]
Xcp <- Xc[(n-np+1):n,]
yp <- y[(n-np+1):n]
ycp <- yc[(n-np+1):n]
outcomep <- outcome[(n-np+1):n]
outcomecp <- rep(1,np)
riskp <- risk[(n-np+1):n]
X <- X[1:(n-np),]
Xf <- Xf[1:(n-np),]
y <- y[1:(n-np)]
yc <- yc[1:(n-np)]
outcome <- outcome[1:(n-np)]



## Make training data set including Xp with censored at 0 y's 
Xa <- rbind(X,Xp)
ya <- c(y,rep(0,np))
outcomea <- c(outcome,rep(0,np))

print(survConcordance(Surv(ycp, event=outcomecp)~riskp)$concord)


###############################
####### Fit Regressions #######
###############################


N.mcmc <- 10000
begin <- 6000
every <- 4
everyXZ <- 1
nplot <- 100
Ly <- 50
beta.init <- NULL

N.mcmc <- 10
begin <- 5
every <- 1
everyXZ <- 1
nplot <- 100
Ly <- 2
beta.init <- NULL


## Fit MVN, MAR (~1 hour) ##
G.vec <- 1
Mmod <- 1

set.seed(12)
system.time(ans.MVN.MAR <- weibull.MCMC(X=Xa[,1:p], y=ya, outcome=outcomea, rho=rho,
rho.prop.01=.3, rho.prop.10=.3, A.tau2=2, B.tau2=.1, rho.kappa=.999, A.kappa=1, B.kappa=1,
prop.sd.kappa=.05, kappa.prop.01=.1, kappa.prop.10=.1, rhoX=0.1, p.ad=.3, p.swap=.2, M.nu=0,
S.nu=100, A.lambda=1, B.lambda=1, prop.sd.lambda=.5, A.Psi=2, B.Psi=2, prop.sd.Psi=.1,
A.eta=2, B.eta=2, prop.sd.eta=.25, M.Beta=0, A.delta=1, B.delta=.1, N.mcmc=N.mcmc,
every=every, nplot=nplot, nback=1000, dfp=20, groups=NULL, mult1=2, mult2=2,
bin.X=bin.X[1:p], discrete.X=discrete.X[1:p], limits.X=limits.X[1:p,], everyXZ=everyXZ,
maxplot=30, beta.init=beta.init, phi.init=NULL, omega.init=NULL, M=Mmod, G=G.vec,
begin=begin, nb=nb, etagp=5, Ly=Ly, yp=ycp, outcomep=outcomecp))

pred.MVN.MAR <- apply(ans.MVN.MAR$risk[,(nt+1):n],2,mean)
risk.CI.MVN.MAR <- t(apply(ans.MVN.MAR$risk[,(nt+1):n],2,quantile,prob=c(.025,.975)))

survConcordance(Surv(ycp, event=outcomecp)~pred.MVN.MAR)$concord
mean(log(riskp) >= risk.CI.MVN.MAR[,1] & log(riskp) <= risk.CI.MVN.MAR[,2]) 





## Fit sparse DPM, MAR (~1 hour) ##

Mmod <- 10
G.vec <- 1:5

set.seed(12)
system.time(ans.sDPM.MAR <- weibull.MCMC(X=Xa[,1:p], y=ya, outcome=outcomea, rho=rho,
rho.prop.01=.3, rho.prop.10=.3, A.tau2=2, B.tau2=.1, rho.kappa=.999, A.kappa=1, B.kappa=1,
prop.sd.kappa=.05, kappa.prop.01=.1, kappa.prop.10=.1, rhoX=0.1, p.ad=.3, p.swap=.2, M.nu=0,
S.nu=100, A.lambda=1, B.lambda=1, prop.sd.lambda=.5, A.Psi=2, B.Psi=2, prop.sd.Psi=.1,
A.eta=2, B.eta=2, prop.sd.eta=.25, M.Beta=0, A.delta=1, B.delta=.1, N.mcmc=N.mcmc,
every=every, nplot=nplot, nback=1000, dfp=20, groups=NULL, mult1=2, mult2=2,
bin.X=bin.X[1:p], discrete.X=discrete.X[1:p], limits.X=limits.X[1:p,], everyXZ=everyXZ,
maxplot=30, beta.init=beta.init, phi.init=NULL, omega.init=NULL, M=Mmod, G=G.vec,
begin=begin, nb=nb, etagp=5, Ly=Ly, yp=ycp, outcomep=outcomecp))

pred.sDPM.MAR <- apply(ans.sDPM.MAR$risk[,(nt+1):n],2,mean)
risk.CI.sDPM.MAR <- t(apply(ans.sDPM.MAR$risk[,(nt+1):n],2,quantile,prob=c(.025,.975)))

survConcordance(Surv(ycp, event=outcomecp)~pred.sDPM.MAR)$concord
mean(log(riskp) >= risk.CI.sDPM.MAR[,1] & log(riskp) <= risk.CI.sDPM.MAR[,2]) 





## Fit MVN, MNAR (~1 hour) ##
G.vec <- 1
Mmod <- 1

set.seed(12)
system.time(ans.MVN.MNAR <- weibull.MCMC(X=Xa, y=ya, outcome=outcomea,
rho=c(1-1E-15,rep(rho,p),rep(1E-200,ncol(Xa)-p)), rho.prop.01=.3, rho.prop.10=.3,
A.tau2=2, B.tau2=.1, rho.kappa=.999, A.kappa=1, B.kappa=1, prop.sd.kappa=.05,
kappa.prop.01=.1, kappa.prop.10=.1, rhoX=0.1, p.ad=.3, p.swap=.2, M.nu=0, S.nu=100,
A.lambda=1, B.lambda=1, prop.sd.lambda=.5, A.Psi=2, B.Psi=2, prop.sd.Psi=.1,
A.eta=2, B.eta=2, prop.sd.eta=.25, M.Beta=0, A.delta=1, B.delta=.1, N.mcmc=N.mcmc,
every=every, nplot=nplot, nback=1000, dfp=20, groups=NULL, mult1=2, mult2=2,
bin.X=bin.X, discrete.X=discrete.X, limits.X=limits.X, everyXZ=everyXZ, maxplot=30,
beta.init=beta.init, phi.init=NULL, omega.init=NULL, M=Mmod, G=G.vec, begin=begin,
nb=nb, etagp=5, Ly=Ly, yp=ycp, outcomep=outcomecp))

pred.MVN.MNAR <- apply(ans.MVN.MNAR$risk[,(nt+1):n],2,mean)
risk.CI.MVN.MNAR <- t(apply(ans.MVN.MNAR$risk[,(nt+1):n],2,quantile,prob=c(.025,.975)))

survConcordance(Surv(ycp, event=outcomecp)~pred.MVN.MNAR)$concord
mean(log(riskp) >= risk.CI.MVN.MNAR[,1] & log(riskp) <= risk.CI.MVN.MNAR[,2]) 





## Fit sparse DPM, MNAR (~1 hour) ##
Mmod <- 10
G.vec <- 1:5

set.seed(12)
system.time(ans.sDPM.MNAR <- weibull.MCMC(X=Xa, y=ya, outcome=outcomea,
rho=c(1-1E-15,rep(rho,p),rep(1E-200,ncol(Xa)-p)), rho.prop.01=.3, rho.prop.10=.3,
A.tau2=2, B.tau2=.1, rho.kappa=.999, A.kappa=1, B.kappa=1, prop.sd.kappa=.05,
kappa.prop.01=.1, kappa.prop.10=.1, rhoX=0.1, p.ad=.3, p.swap=.2, M.nu=0, S.nu=100,
A.lambda=1, B.lambda=1, prop.sd.lambda=.5, A.Psi=2, B.Psi=2, prop.sd.Psi=.1, A.eta=2,
B.eta=2, prop.sd.eta=.25, M.Beta=0, A.delta=1, B.delta=.1, N.mcmc=N.mcmc, every=every,
nplot=nplot, nback=1000, dfp=20, groups=NULL, mult1=2, mult2=2, bin.X=bin.X,
discrete.X=discrete.X, limits.X=limits.X, everyXZ=everyXZ, maxplot=30,
beta.init=beta.init, phi.init=NULL, omega.init=NULL, M=Mmod, G=G.vec, begin=begin,
nb=nb, etagp=5, Ly=Ly, yp=ycp, outcomep=outcomecp))

pred.sDPM.MNAR <- apply(ans.sDPM.MNAR$risk[,(nt+1):n],2,mean)
risk.CI.sDPM.MNAR <- t(apply(ans.sDPM.MNAR$risk[,(nt+1):n],2,quantile,prob=c(.025,.975)))

survConcordance(Surv(ycp, event=outcomecp)~pred.sDPM.MNAR)$concord
mean(log(riskp) >= risk.CI.sDPM.MNAR[,1] & log(riskp) <= risk.CI.sDPM.MNAR[,2]) 





## Fit GBM (~1 min) ##
set.seed(11)
n.trees <- seq(200,800,by=10)
shrinkage <- 10^(seq(-3,-1,length=8))
ans.cv <- cv.gbm(X[,1:p], Surv(y,outcome), n.trees=n.trees, interaction.depth=2,
n.minobsinnode = 5, shrinkage=shrinkage, bag.fraction = .5, distribution="coxph",
foldid=NULL, nfolds=5, seed=220)
ans.cv

ans.gbm <- gbm.fit(X[,1:p], Surv(y,outcome), n.trees=ans.cv$n.trees, interaction.depth=2,
n.minobsinnode=5, shrinkage=ans.cv$shrinkage, bag.fraction=.5, distribution="coxph")
lamhat.gbm <- predict(ans.gbm, Xp[,1:p], n.trees=ans.cv$n.trees)
survConcordance(Surv(ycp, event=outcomecp)~lamhat.gbm)$concord





## Fit depth 1 GBM (~1 min) ##
set.seed(11)
n.trees <- seq(200,800,by=10)
shrinkage <- 10^(seq(-3,-1,length=8))
ans.cv1 <- cv.gbm(X[,1:p], Surv(y,outcome), n.trees=n.trees, interaction.depth=1,
n.minobsinnode = 5, shrinkage=shrinkage, bag.fraction = .5, distribution="coxph",
foldid=NULL, nfolds=5, seed=220)
ans.cv1

ans.gbm1 <- gbm.fit(X[,1:p], Surv(y,outcome), n.trees=ans.cv1$n.trees, interaction.depth=1,
n.minobsinnode=5, shrinkage=ans.cv1$shrinkage, bag.fraction=.5, distribution="coxph")
lamhat.gbm1 <- predict(ans.gbm1, Xp[,1:p], n.trees=ans.cv1$n.trees)
survConcordance(Surv(ycp, event=outcomecp)~lamhat.gbm1)$concord



## Fit glmnet after RF interpolation (~10 min) ##

ya2 <- ya
ya2[ya2==0] <- 1E-8
ans.RFEN1 <- RFEN(Xa[,1:p], ya2, outcomea, bin.X[1:p], maxiter=3, family="cox")

lamhat.RFEN1 <- as.numeric(predict(ans.RFEN1$fit, ans.RFEN1$Xh[(nt+1):n,], s=ans.RFEN1$fit$lambda.min))
survConcordance(Surv(ycp, event=outcomecp)~lamhat.RFEN1)$concord




## Fit glmnet after RF interpolation using missing indicators (~10 min) ##

ans.RFEN2 <- RFEN(Xa, ya2, outcomea, bin.X, maxiter=3, family="cox", use=1:p)

lamhat.RFEN2 <- as.numeric(predict(ans.RFEN2$fit, ans.RFEN2$Xh[(nt+1):n,], s=ans.RFEN2$fit$lambda.min))
survConcordance(Surv(ycp, event=outcomecp)~lamhat.RFEN2)$concord




## Fit stepwise cox using missing indicator vars and filling in 0 for missing vals (~1 min) ##

yS <- Surv(y,outcome)
X0 <- X
X0[is.na(X0)] <- 0
X0p <- as.data.frame(Xp)
X0p[is.na(X0p)] <- 0
Sdata <- as.data.frame(X0)
names(Sdata) <- paste0("X.",1:ncol(Sdata))
names(X0p) <- paste0("X.",1:ncol(Sdata))
Sdata$yS <- yS
formu <- as.formula(paste("~ ",paste(names(Sdata)[1:ncol(X0)], collapse=" + ")))
fit1 <- coxph(yS ~ 1, data=Sdata)
foo0 <- stepAIC(fit1, scope=list(lower=~1, upper=formu),data=Sdata,direction="both")
lamhat.cont.sw <- as.numeric(predict(foo0, newdata=X0p))
survConcordance(Surv(ycp, event=outcomecp)~lamhat.cont.sw)$concord
lamhat.fac.sw <- lamhat.cont.sw
survConcordance(Surv(ycp, event=outcomecp)~lamhat.fac.sw)$concord










#################################
### Save Workspace to a file ####
#################################


if(seed<10)
  fname <- paste("Results/CaseSeed_",case,"_0",seed,".RData",sep="")
if(seed>=10)
  fname <- paste("Results/CaseSeed_",case,"_",seed,".RData",sep="")
save(list=ls(all.names=TRUE), file=fname)









