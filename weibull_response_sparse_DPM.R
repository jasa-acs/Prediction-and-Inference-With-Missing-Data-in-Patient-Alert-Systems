


#############################################################################################
############## Main Function to fit Weibull Regression models via MCMC  #####################
#############################################################################################

weibull.MCMC <- function(X, y, outcome=NULL, rho=.5, rho.prop.01=.3, rho.prop.10=.3, A.tau2=2, B.tau2=.1, rho.kappa=.999, A.kappa=1, B.kappa=1, prop.sd.kappa=.05, kappa.prop.01=.1, kappa.prop.10=.1, rhoX=0.1, p.ad=.3, p.swap=.2, M.nu=0, S.nu=100, A.lambda=1, B.lambda=.1, prop.sd.lambda=.2, A.Psi=NULL, B.Psi=NULL, prop.sd.Psi=.02, A.eta=2, B.eta=2, prop.sd.eta=1, M.Beta=0, A.delta=1, B.delta=.1, N.mcmc=5000, every=1, nplot=10, nback=1000, dfp=20, groups=NULL, mult1=2, mult2=2, bin.X=NULL, discrete.X=NULL, limits.X=NULL, everyXZ=1, maxplot=30, beta.init=NULL, gamma.init=NULL, phi.init=NULL, omega.init=NULL, M=20, G=NULL, begin=0, nb=getDoParWorkers(), etagp=10, Ly=50, yp=NULL, outcomep=NULL){

### INPUTS ###
# X - matrix of predictor variables
# y - vector of event (or censor) times
# outcome - indicator for event or not (i.e., 1=event, 0= right censored)
# rho - a vector of the prior probabilities of variable inclusion in the regression model.
#       Defaults to intercept~1 and all variables rho=0.5.
# rho.prop.01 - tuning parameter for rho proposals; probability rho going from 0 -> 1,
#               defaults to .3
# rho.prop.10 - tuning parameter for rho proposals; probability rho going from 1 -> 0,
#               defaults to .3
# A.tau2, B.tau2 - prior for tau2 ~ IG(A.tau2, B.tau2)
# rho.kappa - Not used.
# A.kappa, B.kappa - prior for kappa ~ Gamma(A.kappa, B.kappa)
# prop.sd.kappa - scale parameter for MH proposal for kappa
# kappa.prop.01 - Not used
# kappa.prop.10 - Not used
# rhoX - prior probabilities of variable inclusion in the sDPM for the model
#        for the X distribution.
# p.ad - tuning parameter for variable selection in sDPM proposal.
#        Probability of proposing an additional variable
# p.swap - tuning parameter for variable selection in sDPM proposal.
#          Probablity of swapping variables between informative/noninformative
# M.nu, S.nu - prior for nu ~ N(M.nu, S.nu)
# A.lambda, B.lambda - prior for lambda ~ 
# prop.sd.lambda - prior on lambda ~ Gamma(A.lambda, B.lambda)
# A.Psi, B.Psi - prior on psi ~ Gamma(A.Psi, B.Psi) defaults to gamma(2,2).
# prop.sd.Psi - scale parameter for MH proposal for Psi
# A.eta, B.eta - prior on eta ~ p + etagp + Gamma(A.eta, B.eta)
#                where p is number of Z variables and etagp is a predefined constant
# etagp - See above; defaults to 10.
# prop.sd.eta - scale parameter for MH proposal for eta
# M.Beta - Not used
# A.delta, B.delta - prior on delta ~ Gamma(A.delta, B.delta)
# N.mcmc - Number of MCMC iterations to run.
# begin - Start recording values of parameters after 'begin' burn-in iterations
# every - Record current value of parameters every 'every' iterations
# nplot - Plot a summary of MCMC progress every 'nplot' iterations
# nback - Plot traces in the summary from current iteration back through nback iterations.
# dfp - degrees of freedom for the t-proposal for random walk MH proposals.
# groups - which variables (columns of X) should be grouped together
#          (and thus updated together) for variable selection purposes.
# mult1 - Tuning parameter (variance inflation factor) for beta proposals.
# mult2 - Tuning parameter (variance inflation factor) for XZ proposals.
# bin.X -  bin.X is an indicator vector to indicate whether each column of X is binary
#          variable or not
# discrete.X - indicator vector to indicate whether each column of X is a discrete
#              (but not binary) variable or not.
# limits.X - a matrix with two columns (min, max) that contains hard limits for some
#            of the variables. The default is to just set these all to (-Inf, Inf).
# everyXZ - XZ updates take much longer than the other updates.  Update XZ only
#           every 'everyXZ' iterations.
# maxplot - Maximum number of plots to attempt to display during MCMC progress updates
# beta.init - Starting value for beta. defaults to all 0's
# gamma.init - starting value for gamma, defaults to an Elastic Net fit after Ly
#              iterations of missing value updates assuming 1 component normal distn.
# phi.init - Starting values for phi.  Defaults to a kmeans clustering on the continuous
#            variables after Ly iterations of missing value updates.
# omega.init - starting values for omega.  Defaults to proportions determined from the
#              kmeans fit described above for phi.init.
# M - Maximum number of components to allow in the DPM approximation/truncation.
# G - vector of number of components to try for the initialization of kmeans clustering
#     Defaults to 1:M
# nb - Number of blocks to use for parallel updates of XZ.  Defaults to number of
#      available cores
# Ly - See above for gamma.init, phi.init
# yp - Optional vector of event/censor times for the last length(yp) observations that
#      are being held out.  Must set y=0 and outcome=0 for these elements to effectively
#      hold them out of the estimation process (i.e., give them zero influence on the
#      likelihood)
# outcomep - Optional vecotr of outcomes (event or censored) for the last length(yp)
#            observations that are being held out.
####################

### OUTPUTS ###
# beta - An N x p matrix of N MCMC samples of the p elements of beta, where
#        N = (N.mcmc-begin)%/%every
# tau2 - An N x nG matrix of N MCMC samples of the nG=length(groups) elements of tau2
# kappa - A N vector of MCMC samples of kappa
# muX - An N x M x p array of N MCMC samples of the M mean vectors for the sDPM of X/Z
# RX - An N x M x p x p array of N MCMC samples of the M covariance matrices
#      for the sDPM of X/Z.
# omega - An N x M matrix of MCMC samples of the omega vector (mixture probabilities).
# X - The input X matrix
# y - The input y vector
# outcome - The input outcome vector
# bin.X - The input bin.X vector
# risk - a N x n matrix of N MCMC samples of the predicted log-risk for each observation.
##################


 ## Create needed variables ##
  n <- length(y)
  p <- ncol(X)
  if(is.null(groups))
    groups <- as.list(1:(p+1))
  else
    groups <- c(list(1),sapply(groups,"+",1))
  nG <- length(groups)
  if(is.null(outcome))
    outcome <- rep(1,n)
  y.Surv <- Surv(y, event=outcome)

  if(length(rho.prop.01)==1)
    rho.prop.01 <- rep(rho.prop.01, nG)
  if(length(rho.prop.10)==1)
    rho.prop.10 <- rep(rho.prop.10, nG)
  if(length(rho)==1)
    rho <- rep(rho, nG)
  rho.prop.01[1] <- .9
  rho.prop.10[1] <- .1
  rho[1] <- 1-1E-15
  if(length(rhoX)==1)
    rhoX <- rep(rhoX, p)

  if(is.null(G))
    G <- 1:M
  if(length(mult1)==1)
    mult1 <- rep(mult1, nG)

  if(length(M.nu)==1)
    M.nu <- rep(M.nu,p)
  if(length(S.nu)==1)
    S.nu <- diag(S.nu,p)
  if(is.null(A.Psi))
    A.Psi <- 2
  if(is.null(B.Psi))
    B.Psi <- 2


  if(is.null(bin.X)){
    bin.X <- rep(0,p)
    for(j in 1:p){
      bin.X[j] <- as.numeric(length(unique(X[!is.na(X[,j]),j]))<3)
    }
  }
  if(is.null(discrete.X)){
    discrete.X <- rep(0,p)
    for(j in 1:p){
      foo <- length(unique(X[!is.na(X[,j]),j]))
      discrete.X[j] <- as.numeric(foo>=3 && foo<(M+1))
    }
  }
  if(is.null(limits.X))
    limits.X <- cbind(rep(-Inf,p),rep(Inf,p))
  ind.scale <- rep(FALSE,p)
  for(j in 2:nG){
    bin.j <- bin.X[groups[[j]][1]-1]
    if(bin.j)
      ind.scale[groups[[j]][1]-1] <- TRUE
  }

  gammap.1 <<- gamma2(1)
  gammapp.1 <<- gamma3(1)

 ## From here on out X is a p+1 column matrix (1st col is intercept) ##
  X <- cbind(1,X)

## initialize kappa
  kappa.now <- 1

  Psi.now <- diag(A.Psi/B.Psi,p)
  eta.now <- A.eta/B.eta + p + etagp
  lambda.now <- A.lambda/B.lambda

 ## initialize Xc, Z, cluster probs, means, Precisions, and beta
cat("\nInitializing Clusters \n")
  ans.XZ <- initialize.XZ(X, y, outcome, bin.X, discrete.X, limits.X, M, G, phi.init, beta.init, gamma.init, kappa.now, groups, rho, rhoX, ind.scale, lambda.now, eta.now, Psi.now, Ly=Ly, nb=nb)

  Xc.now <- ans.XZ$Xc
  Z.now <- ans.XZ$Z
  muX.now <- ans.XZ$muX
  QX.now <- ans.XZ$QX
  ldet.now <- ans.XZ$ldet.now
  ldet11.now <- ans.XZ$ldet11.now
  gamma.now <- ans.XZ$gamma.now
  phi.now <- ans.XZ$phi.now
  beta.now <- ans.XZ$beta.now
  Beta.now <- M.Beta
  if(length(Beta.now)==1)
    Beta.now <- rep(Beta.now,p+1)
  nu.now <- apply(muX.now,2,mean)


  delta.now <- A.delta/B.delta
  ans.vee <- update.vee(phi.now, delta.now, M)
  vee.now <- ans.vee$vee
  omega.now <- ans.vee$omega
  
  tau2.now <- rep(B.tau2/A.tau2, nG)
  risk.now <- rep(0,n)

cat("\nAllocating Memory for Posterior Objects \n")

 ## Allocate posterior objects for which to store MCMC samples
  beta <- array(0, c((N.mcmc-begin)%/%every, p+1))
  kappa <- rep(0, (N.mcmc-begin)%/%every)
  lambda <- rep(0, (N.mcmc-begin)%/%every)
  eta <- rep(0, (N.mcmc-begin)%/%every)
  psi <- rep(0, (N.mcmc-begin)%/%every)
  tau2 <- matrix(0, (N.mcmc-begin)%/%every, nG)
  muX <- array(0, c((N.mcmc-begin)%/%every, M, p))
  RX <- array(0, c((N.mcmc-begin)%/%every,M,p,p))
  vee <- matrix(0, (N.mcmc-begin)%/%every, M)
  omega <- matrix(0, (N.mcmc-begin)%/%every, M)
  risk <- matrix(0, (N.mcmc-begin)%/%every, n)
  
  accept.beta <- matrix(0, N.mcmc, nG)
  accept.kappa <- rep(0, N.mcmc)
  accept.eta <- rep(0, N.mcmc)
  accept.lambda <- rep(0, N.mcmc)
  accept.Psi <- rep(0, N.mcmc)
  accept.mQg <- rep(0, N.mcmc)
  accept.XZ <- rep(0, sum(is.na(X[,-1])))

  concord1.avg <- .8
  concord2.avg <- .8
  h.avg <- .05

 
 ################
 ## Begin MCMC ##
 ################
 cat("\n")
  for(it in 1:N.mcmc){

if(it%%10==0)
  cat("\nIteration", it, "out of", N.mcmc)

  ## y update
    yc.now <- update.y(y, outcome, Xc.now, beta.now, phi.now, kappa.now)


  ## beta update
    ans.beta <- update.beta(y, yc.now, outcome, Xc.now, beta.now, kappa.now, phi.now, Beta.now, tau2.now, groups, mult1, rho, rho.prop.01, rho.prop.10)
    beta.now <- ans.beta$beta
    accept.beta[it,] <- ans.beta$accept


  ## tau^2 update
    tau2.now <- update.tau2(beta.now, groups, Beta.now, A.tau2, B.tau2)


  ## kappa update
    ans.kappa <- update.kappa(y, outcome, Xc.now, beta.now, kappa.now, phi.now, A.kappa, B.kappa, rho.kappa, prop.sd.kappa, kappa.prop.01, kappa.prop.10, dfp)
    kappa.now <- ans.kappa$kappa
    accept.kappa[it] <- ans.kappa$accept

  ## X and Z updates only every everyXZ iterations
    if(it%%everyXZ==0){
      ans.XZ <- update.XZ(X, Xc.now, Z.now, bin.X, discrete.X, limits.X, muX.now, QX.now, ldet11.now, y, outcome, beta.now, kappa.now, mult2, omega.now, gamma.now, nb)
      Xc.now <- ans.XZ$Xc
      Z.now <- ans.XZ$Z
      phi.now <- ans.XZ$phi.now
      accept.XZ <- accept.XZ + ans.XZ$accept
    }


  ## Update vee/omega
    ans.vee <- update.vee(phi.now, delta.now, M)
    vee.now <- ans.vee$vee
    omega.now <- ans.vee$omega


  ## Update delta
    delta.now <- update.delta(vee.now, A.delta, B.delta)



  ## Psi update
    ans.Psi <- update.Psi(Psi.now, muX.now, QX.now, A.Psi, B.Psi, eta.now, nu.now, lambda.now, gamma.now, prop.sd.Psi, dfp, ind.scale)
    Psi.now <- ans.Psi$Psi
    accept.Psi[it] <- ans.Psi$accept
    

  ## lambda update
    ans.lambda <- update.lambda(muX.now, QX.now, ldet11.now, gamma.now, nu.now, lambda.now, eta.now, Psi.now, A.lambda, B.lambda, prop.sd.lambda, dfp, ind.scale)
    lambda.now <- ans.lambda$lambda
    accept.lambda[it] <- ans.lambda$accept


  ## eta update 
    ans.eta <- update.eta(muX.now, QX.now, gamma.now, nu.now, lambda.now, eta.now, Psi.now, A.eta, B.eta, prop.sd.eta, dfp, etagp, ind.scale)
    eta.now <- ans.eta$eta
    accept.eta[it] <- ans.eta$accept

  ## muX, QX, gamma update
    ans.mQg <- update.muX.QX.gamma(gamma.now, muX.now, QX.now, Z.now, p.ad, p.swap, rhoX, phi.now, nu.now, lambda.now, eta.now, Psi.now, ind.scale, nb, it)

    muX.now <- ans.mQg$muX.now
    QX.now <- ans.mQg$QX.now
    gamma.now <- ans.mQg$gamma.now
    ldet.now <- ans.mQg$ldet.now
    ldet11.now <- ans.mQg$ldet.now
    accept.mQg[it] <- ans.mQg$accept


  ## scale mu and Q

   ## scale muX and QX
    RX.now <- QX.now
    nuX.now <- muX.now
    for(m in 1:M){
      scale.m <- diag(1/sqrt(diag(QX.now[m,,])))
      diag(scale.m)[!ind.scale] <- 1
      RX.now[m,,] <- scale.m%*%QX.now[m,,]%*%t(scale.m)
      nuX.now[m,] <- muX.now[m,]*(1/diag(scale.m))
    }

  ## End of Updates

    risk.now <- as.numeric(Xc.now%*%beta.now)

   ## record params.now in posterior sample
    if(it>begin && (it-begin)%%every==0){

      beta[(it-begin)/every,] <- beta.now
      tau2[(it-begin)/every,] <- tau2.now
      kappa[(it-begin)/every] <- kappa.now
      lambda[(it-begin)/every] <- lambda.now
      eta[(it-begin)/every] <- eta.now
      psi[(it-begin)/every] <- Psi.now[1,1]
      omega[(it-begin)/every,] <- omega.now
      vee[(it-begin)/every,] <- vee.now
      muX[(it-begin)/every,,] <- nuX.now
      RX[(it-begin)/every,,,] <- RX.now
      risk[(it-begin)/every,] <- risk.now
    }
    

   ## Summarize and Plot posterior
    if((it-begin)%%nplot==0){
    
      concord1 <- survConcordance(y.Surv~risk.now)$concord
      concord1.avg <- (1-h.avg)*concord1.avg + h.avg*concord1
      if(is.null(yp)){
        concord2 <- survConcordance(Surv(yc.now, event=rep(1,n))~risk.now)$concord
        concord2.avg <- (1-h.avg)*concord2.avg + h.avg*concord2
      }
      else{
        np <- length(yp)
        concord2 <- survConcordance(Surv(yp, event=outcomep)~risk.now[(n-np+1):n])$concord
        concord2.avg <- (1-h.avg)*concord2.avg + h.avg*concord2
      }
      pp <- min(maxplot,p)
      pcv <- sum(gamma.now==1)
      max.pw.Z <- 10
      N.plots <- (pp+1)+1+2+3 + min(max.pw.Z,choose(pcv,2))+1
      cols <- min(13, ceiling(sqrt(N.plots)))
      rows <- min(13, ceiling(N.plots/cols))
      max.pw.Z <- 10 + rows*cols-N.plots
      it.e <- floor((it-begin)/every)
      ind.now <- max(1,floor(it.e/2),it.e-nback+1):max(1,it.e)
      par(mfrow=c(rows,cols), mar=c(2,2,2,1))

     ## Print Rsq and recent average
      cat("\n\nConcordance  =",concord1)
      cat("\nConcordance (recent avg) =",concord1.avg)
      if(is.null(yp)){
        cat("\nConcordance2  =",concord2)
        cat("\nConcordance2 (recent avg) =",concord2.avg,"\n")
      }
      else{
        cat("\nConcordance-CV  =",concord2)
        cat("\nConcordance-CV (recent avg) =",concord2.avg,"\n")
      }
      like.y <- get.like(y, outcome, Xc.now, beta.now, kappa.now, phi.now)
      cat("\nlog(like) =",like.y,"\n")

      cat("\ngamma.now =",gamma.now,"\n")
      cat("\npsi.now =",Psi.now[1,1],"\n")

     ## Print acceptance %
      cat("\nmu Q gamma acceptance = ", mean(accept.mQg[1:it]))
      cat("\nlambda acceptance = ", mean(accept.lambda[1:it]))
      cat("\neta acceptance = ", mean(accept.eta[1:it]))
      cat("\nPsi acceptance = ", mean(accept.Psi[1:it]))
      cat("\nkappa acceptance = ", mean(accept.kappa[1:it]))
      cat("\nbeta acceptance = \n")
      print(summary(colMeans(accept.beta[1:it,,drop=FALSE])))
      cat("\nnX/Z acceptance = \n")
      print(summary(accept.XZ/floor(it/everyXZ)))
      cat("\n\n\n")
    }
    if((it-begin)%%nplot==0 && it>begin){
     ## Plot Betas
      Beta.ord <- c(1,order(-abs(colMeans(beta[ind.now,-1,drop=FALSE])))+1)
      for(k in 1:(pp+1))
        plot(beta[ind.now,Beta.ord[k]],ylab="",main=paste("Beta_",Beta.ord[k]-1,sep=""),cex=.5)

#     ## Plot tau2's
#     plot(tau2[ind.now,1],ylab="",main="tau2",cex=.5)
     
     ## Plot kappa
      plot(kappa[ind.now],ylab="",main="kappa",cex=.5)

     ## Plot lambda
      plot(lambda[ind.now],ylab="",main="lambda",cex=.5)

     ## Plot eta
      plot(eta[ind.now],ylab="",main="eta",cex=.5)

     ## Plot Psi
      plot(psi[ind.now],ylab="",main="Psi",cex=.5)

    ## Bar plot of omega.now
      barplot(omega.now, main="omega")
      abline(h=exp(-6), col=2)
      abline(h=exp(-4), col=4)

      lomega <- -sort(-log(omega.now+1E-300))
      barplot(lomega, main="omega", yaxt='n')
      at <- ceiling(1/10*round(seq(min(lomega),0,length=4),0))*10
      labels <- format(exp(at),nsmall=1,digits=1)
      axis(2,at=at, labels=labels)
      abline(h=-6, col=2)
      abline(h=-4, col=4)

     ## Plot Z.now w/phi.now for c.vars
      ind.cv <- which(gamma.now==1)
      names.cv <- paste("x",ind.cv,sep="")
      n.cv <- length(ind.cv)
      indp <- 1
      if(n.cv==1){
        ind2 <- ifelse(ind.cv==1,2,1)
        plot(Z.now[,ind.cv], Z.now[,ind2], col=phi.now, main=paste(names.cv[1],"by",names.cv[ind2]),pch=16)
      }
      else{
        for(j in 1:(n.cv-1)){
          for(k in (j+1):n.cv){
            if(indp <= max.pw.Z){
	      plot(Z.now[,ind.cv[j]], Z.now[,ind.cv[k]], col=phi.now, pch=16, cex=.7, cex.axis=.8, mgp=c(2,.35,0))
	      title(main=paste(names.cv[j],"by",names.cv[k]),line=.5, cex.main=.8)
              indp <- indp + 1	  	      
	    }
	    else
	      break
	  }
        }
      }           
    }
  }
  return(list(beta=beta, tau2=tau2, kappa=kappa, muX=muX, RX=RX, omega=omega, vee=vee, X=X, y=y, outcome=outcome, bin.X=bin.X, risk=risk))
}







############################################################################################
## All other functions below are lower level functions mean to serve the main function    ##
## above and/or the comparison methods in the simulation                                  ##
############################################################################################


######################################################################################
################### MCMC Updating Functions used in weibull.MCMC######################
######################################################################################


#### Positive proposal generation and densitiy functions for rho and kappa ####
dlogt <- function(x, df=1, sigma=1, mu=0){
  num <- gamma((df+1)/2)
  denom <- x*sigma*sqrt(pi*df)*gamma(df/2)*(1+1/df*((log(x)-mu)/sigma)^2)^((df+1)/2)
  return(num/denom)
}

ldlogt <- function(x, df=1, sigma=1, mu=0){
  num <- lgamma((df+1)/2)
  denom <- log(x)+log(sigma)+.5*log(pi*df)+lgamma(df/2)+((df+1)/2)*log(1+1/df*((log(x)-mu)/sigma)^2)
  return(num-denom)
}


rlogt <- function(n, df=1, sigma=1, mu=0){
  return(exp(rt(n,df)*sigma+mu))
}


r.pos.proposal <- function(prev, df, sigma){
  return(rlogt(1, df=df, sigma=sigma, mu=log(prev)))
}

d.pos.proposal <- function(prop, prev, df, sigma){
  return(ldlogt(prop, df=df, sigma=sigma, mu=log(prev)))
}



## Evaluate complete data model likelihood
get.like <- function(y, outcome, X, beta.now, kappa.now, phi.now, lambda.now=NULL){

  n <- length(y)
  if(is.null(lambda.now)){
    lambda.now <- exp(as.numeric(X%*%beta.now))
  }
  ind.cens <- (outcome==0)
  ans <- 0
  if(any(ind.cens))
    ans <- ans + sum(pweibull(y[ind.cens], kappa.now, 1/lambda.now[ind.cens], lower.tail=FALSE, log.p=TRUE))
  if(any(!ind.cens))
    ans <- ans + sum(dweibull(y[!ind.cens], kappa.now, 1/lambda.now[!ind.cens], log=TRUE))
  return(ans)
}



## Evaluate likelihood for a single observation; used in XZ update
get.1obs.like <- function(y, outcome, x, beta.now, kappa.now){

  lambda.now <- exp(sum(beta.now*x))
  if(outcome==0)
    ans <- pweibull(y, kappa.now, 1/lambda.now, lower.tail=FALSE, log.p=TRUE)
  else
    ans <- dweibull(y, kappa.now, 1/lambda.now, log=TRUE)
  return(ans)
}




## generate complete y's for censored observations
update.y <- function(y, outcome, Xc.now, beta.now, phi.now, kappa.now){

  ind.cens <- which(outcome==0)
  if(length(ind.cens)>0){
    y.cens <- y[ind.cens]
    lambda.cens <- exp(as.numeric(Xc.now[ind.cens,,drop=FALSE]%*%beta.now))
    lprob <- log(runif(length(ind.cens)))+pweibull(y.cens, kappa.now, 1/lambda.cens, lower.tail=FALSE,log.p=TRUE)
    y.cens <- qweibull(lprob, kappa.now, 1/lambda.cens, lower.tail=FALSE,log.p=TRUE)
    y[ind.cens] <- y.cens
  }
  return(y)
}





## Kappa MH Update
update.kappa <- function(y, outcome, Xc.now, beta.now, kappa.now, phi.now, A.kappa, B.kappa, rho.kappa, prop.sd.kappa, kappa.prop.01, kappa.prop.10, dfp){

  accept <- 0
  B.k.now <- as.numeric(kappa.now!=1)
  prob.p.g.now <- ifelse(B.k.now==0, kappa.prop.01, 1-kappa.prop.10)
  B.k.p <- rbinom(1,1,prob.p.g.now)
  prob.now.g.p <- ifelse(B.k.p==0, kappa.prop.01, 1-kappa.prop.10)

  if(B.k.p==1)
    kappa.p <- r.pos.proposal(kappa.now, df=dfp, sigma=prop.sd.kappa)
  else
    kappa.p <- 1
  like.p <- get.like(y, outcome, Xc.now, beta.now, kappa.p, phi.now)
  like.now <- get.like(y, outcome, Xc.now, beta.now, kappa.now, phi.now)
  pi.now <- B.k.now*(log(rho.kappa)+ dgamma(kappa.now, A.kappa, B.kappa, log=TRUE)) + (1-B.k.now)*log(1-rho.kappa)
  pi.p <- B.k.p*(log(rho.kappa)+ dgamma(kappa.p, A.kappa, B.kappa, log=TRUE)) + (1-B.k.p)*log(1-rho.kappa)
  dprop.now <- B.k.p*(log(prob.p.g.now)+d.pos.proposal(kappa.now, kappa.p, df=dfp, sigma=prop.sd.kappa)) + (1-B.k.p)*log(1-prob.p.g.now)
  dprop.p <- B.k.now*(log(prob.now.g.p)+d.pos.proposal(kappa.p, kappa.now, df=dfp, sigma=prop.sd.kappa)) + (1-B.k.now)*log(1-prob.now.g.p)
  MH.ratio <- exp((like.p + pi.p - dprop.p) - (like.now + pi.now - dprop.now))
  if(runif(1) < MH.ratio){
    kappa.now <- kappa.p
    accept <- 1
  }
  return(list(kappa=kappa.now, accept=accept))
}




## Gibbs update for tau
update.tau2 <- function(beta.now, groups, Beta.now, A.tau2, B.tau2){

  nG <- length(groups)
  p <- length(Beta.now)-1
  indneq0 <- which(beta.now!=0)
  nneq0 <- length(indneq0)
  A.s <- A.tau2 + nneq0/2
  B.s <- B.tau2 + sum(beta.now[indneq0]^2)/2      
  tau2.now <- rep(1/rgamma(1,A.s,B.s), nG)
  return(tau2.now)
}



## MH update for beta
update.beta <- function(y, yc.now, outcome, Xc.now, beta.now, kappa.now, phi.now, Beta.now, tau2.now, groups, mult1, rho, rho.prop.01, rho.prop.10){

  n <- length(y)
  ys <- -log(yc.now)
  nG <- length(groups)
  accept <- rep(0,nG)
  lambda.now <- exp(as.numeric(Xc.now%*%beta.now))

  for(k in 1:nG){
    sig2 <- 1/kappa.now^2*(gammapp.1 - gammap.1^2)*mult1[k]
    J.k <- length(groups[[k]])
    B.k.now <- as.numeric(beta.now[groups[[k]]][1]!=0)
    prob.p.g.now <- ifelse(B.k.now==0, rho.prop.01[k], 1-rho.prop.10[k])
    B.k.p <- rbinom(1,1,prob.p.g.now)
    prob.now.g.p <- ifelse(B.k.p==0, rho.prop.01[k], 1-rho.prop.10[k])
    Q0 <- diag(1/tau2.now[k], J.k)
    M0 <- Beta.now[groups[[k]]]
    Q0M0 <- Q0%*%M0
    dp.g.now <- dnow.g.p <- pi.p <- pi.now <- 0
    beta.p <- beta.now
    lambda.p <- lambda.now

    beta.k.now <- beta.now[groups[[k]]]
    X.k <- Xc.now[,groups[[k]],drop=FALSE]
    log.lambda.mk <- as.numeric(log(lambda.now) - X.k%*%beta.k.now)
    yt <- as.numeric(ys + 1/kappa.now*gammap.1 - log.lambda.mk)
    Q.s <- as.matrix(1/sig2*(crossprod(X.k) + sig2*Q0))
    ans.inv <- ginv.gp(Q.s)
    S.s <- ans.inv$inv
    mu.s <- as.numeric(1/sig2*S.s%*%( sig2*Q0M0 + crossprod(X.k,yt) ))
      
    if(B.k.p==1)
      beta.k.p <- rmvnorm(1,mu=mu.s, S.sqrt=ans.inv$sqrt.inv)
    else
      beta.k.p <- rep(0,J.k)

    dp.g.now <- dp.g.now + B.k.p*(log(prob.p.g.now)+dmvnorm(beta.k.p, mu.s, S.inv=Q.s, log.det=-ans.inv$log.det)) + (1-B.k.p)*log(1-prob.p.g.now)
    dnow.g.p <- dnow.g.p + B.k.now*(log(prob.now.g.p)+dmvnorm(beta.k.now, mu.s, S.inv=Q.s, log.det=-ans.inv$log.det)) + (1-B.k.now)*log(1-prob.now.g.p)

    pi.p <- pi.p + B.k.p*(log(rho[k])+sum(dnorm(beta.k.p,M0,sqrt(tau2.now[k]),log=TRUE))) + (1-B.k.p)*log(1-rho[k])
    pi.now <- pi.now + B.k.now*(log(rho[k])+sum(dnorm(beta.k.now,M0,sqrt(tau2.now[k]),log=TRUE))) + (1-B.k.now)*log(1-rho[k])
    beta.p[groups[[k]]] <- beta.k.p
    lambda.p <- as.numeric(exp(log.lambda.mk + X.k%*%beta.k.p))

    like.p <- get.like(y, outcome, kappa.now=kappa.now, lambda.now=lambda.p)
    like.now <- get.like(y, outcome, kappa.now=kappa.now, lambda.now=lambda.now)
    MH.ratio <- exp(like.p+pi.p+dnow.g.p-like.now-pi.now-dp.g.now)
    if(runif(1) < MH.ratio){
      beta.now <- beta.p
      lambda.now <- lambda.p
      accept[k] <- 1
    }
  }
  return(list(beta=beta.now, accept=accept))
}




## identify reasonable starting values for missing X's and latent Z's
initialize.XZ <- function(X, y, outcome, bin.X, discrete.X, limits.X, M, G=1:M, phi.init, beta.init, gamma.init, kappa.now, groups, rho, rhoX, ind.scale, lambda.now, eta.now, Psi.now, Ly=10, nb, use.mclust=FALSE){

  Xc.now <- X
  Z.now <- X[,-1]
  p <- ncol(X)-1
  n <- nrow(X)

  for(j in 1:p){
    if(bin.X[j]==1)
      Z.now[,j] <- 2*X[,j+1]-1
  }

 ## initialize Xc.now, Z.now
  muX <- apply(Z.now,2,mean,na.rm=TRUE)
  SigmaX <- cov(Z.now,use="pairwise.complete.obs")
  QX <- try(ginv.gp(SigmaX)$inv,silent=TRUE)
  while(is.character(QX[1])){
    diag(SigmaX) <- 1.05*diag(SigmaX)
    QX <- try(ginv.gp(SigmaX)$inv,silent=TRUE)
  }

  for(j in 1:p){
    ind.miss.j <- which(is.na(X[,j+1]))
    Z.now[ind.miss.j,j] <- muX[j]
    if(bin.X[j]==1)
      Xc.now[,j+1] <- as.numeric(Z.now[,j]>0)
    else
      Xc.now[,j+1] <- Z.now[,j]
  }    
  gamma.now <- rep(1,p)
  if(is.null(beta.init)){
    rhofoo <- rep(rho, sapply(groups, length))
    indkeep <- which(rhofoo > .0001)[-1] - 1
    yS <- cbind(y,outcome)
    yS[yS[,1]==0,] <- 1E-6
    colnames(yS) <- c("time", "status")
    beta.1 <- rep(0,p+1)
    beta.1[1] <- -mean(log(y[y>0])) + 1/kappa.now*gammap.1
    ans.fit <- cv.glmnet(Xc.now[,indkeep], yS, family="cox", alpha=.5, nfolds=5)
    ind <- which(ans.fit$glmnet.fit$lambda==ans.fit$lambda.1se)
    beta.1[indkeep+1] <- ans.fit$glmnet.fit$beta[,ind]
    beta.now <- beta.1
  }
  else
    beta.now <- beta.init
  QX.now <- array(QX,c(1,p,p))
  muX.now <- array(muX,c(1,p))
  for(i in 1:Ly){
    ans.i <- update.XZ(X, Xc.now, Z.now, bin.X, discrete.X, limits.X, muX.now, QX.now, ldet11.now=0, y, outcome, beta.now, kappa.now=1, mult2=1, omega.now=1, gamma.now, nb=nb)
    Xc.now <- ans.i$Xc.now
    Z.now <- ans.i$Z.now
  }
    

 ## initial beta
  if(is.null(beta.init)){
    rhofoo <- rep(rho, sapply(groups, length))
    indkeep <- which(rhofoo > .0001)[-1] - 1
    yS <- cbind(y,outcome)
    yS[yS[,1]==0,] <- 1E-6
    colnames(yS) <- c("time", "status")
    beta.1 <- rep(0,p+1)
    beta.1[1] <- -mean(log(y[y>0])) + 1/kappa.now*gammap.1
    ans.fit <- cv.glmnet(Xc.now[,indkeep], yS, family="cox", alpha=.5, nfolds=5)
    ind <- which(ans.fit$glmnet.fit$lambda==ans.fit$lambda.1se)
    beta.1[indkeep+1] <- ans.fit$glmnet.fit$beta[,ind]
    beta.now <- beta.1
  }
  else{
    beta.now <- beta.init
  }
  
 ## initial gamma
  if(is.null(gamma.init))
    gamma.now <- ((beta.now[-1]!=0 | rhoX >=.99)*1)
  else
    gamma.now <- gamma.init
  if(all(gamma.now==1))
    gamma.now[p] <- 0
  if(all(gamma.now==0)){
    if(any(bin.X==1))
      gamma.now[min(which(bin.X==0))] <- 1
    else
      gamma.now[1] <- 1
  }
  
 ## initial clustering
  if(sum(gamma.now)>1)
    c.vars <- (gamma.now==1 & bin.X==0)
  else
    c.vars <- (gamma.now==1)
  nc.vars <- !c.vars


  if(is.null(phi.init)){
    if(use.mclust){
      foo2 <- Mclust(Z.now[,c.vars,drop=FALSE],G,modelNames="VVV")
      phi.now <- foo2$classification
    }
    else{
      ans.kmeans <- mykmeans(Z.now[,c.vars,drop=FALSE],G)
      phi.now <- ans.kmeans$cluster
    }
    ord.m <- order(-table(phi.now))
    phi.now <- match(phi.now, ord.m)
  }
  else{
    phi.now <- phi.init
  }
  K.foo <- max(phi.now)
  muX.now <- matrix(0,M,p)
  QX.now <- SigmaX.now <- array(0,c(M,p,p))
  bs <- matrix(0,M,p)
  nvec <- ldet.now <- ldet11.now <- rep(0,M)
  prec.avg <- 0


  ans.rmuQ <- rmuQ(Z.now, phi.now, gamma.now, nu.now, lambda.now, eta.now, Psi.now, M, ind.scale, muX.now=NULL, QX.now=NULL)
  muX.now <- ans.rmuQ$mu
  QX.now <- ans.rmuQ$Q
  ldet.now <- ans.rmuQ$ldet.vec
  ldet11.now <- ans.rmuQ$ldet11.vec


  return(list(Xc.now=Xc.now, Z.now=Z.now, muX.now=muX.now, QX.now=QX.now, phi.now=phi.now, beta.now=beta.now, ldet.now=ldet.now, ldet11.now=ldet11.now, gamma.now=gamma.now))
}










## MH update for missing X's and latent Z's

update.XZ <- function(X, Xc.now, Z.now, bin.X, discrete.X, limits.X, muX.now, QX.now, ldet11.now, y, outcome, beta.now, kappa.now, mult2, omega.now, gamma.now, nb=getDoParWorkers()){

  n <- nrow(X)
  p <- ncol(X)-1
  udis <- list()
  for(j in 1:p){
    if(discrete.X[j]==1){
      udis[[j]] <- sort(unique(Xc.now[,j+1]))
      udis[[j]] <- c(udis[[j]], Inf)
    }
    else
      udis[[j]] <- c()
  }
  
  M <- length(omega.now)
  c.vars <- which(gamma.now==1)
  nc.vars <- which(gamma.now==0)
  p1 <- length(c.vars)
  S11.inv <- array(0,c(M,p1,p1))
  ans.Q22 <- ginv.gp(QX.now[1,nc.vars,nc.vars])
  Q12Q22.invQ21 <- QX.now[1,c.vars,nc.vars]%*%ans.Q22$inv%*%QX.now[1,nc.vars,c.vars]
  for(m in 1:M)
    S11.inv[m,,] <- QX.now[m,c.vars,c.vars]-Q12Q22.invQ21
  
## for each missing j, draw Zij (and effectively Xij) | rest, and accept/reject
  blocks <- list()
  ind.now <- 1
  inc <- ceiling(n/nb)
  for(b in 1:nb){
    blocks[[b]] <- ind.now:min(ind.now+inc-1, n)
    ind.now <- ind.now+inc
  }

  XZ <- foreach(b=1:nb, .combine=rbind)%dopar%{
    Xc.now.b <- matrix(0, length(blocks[[b]]), p+1)
    Z.now.b <- matrix(0, length(blocks[[b]]), p)
    accept.b <- matrix(NA, length(blocks[[b]]), p)
    phi.now.b <- rep(0, length(blocks[[b]]))
    for(k in 1:length(blocks[[b]])){
      i <- blocks[[b]][k]
      
     ## Gibbs sample a new phi.now[i]
      xc.now.i <- Xc.now[i,]
      z.now.i <- Z.now[i,]
      probs.i <- dmvmixnorm(z.now.i[c.vars], omega.now, muX.now[,c.vars,drop=FALSE], S11.inv, ldet11.now)
      phi.i <- sample(1:M, 1, prob=probs.i)

     ## for each missing j, draw Xij, Zij | rest 
      ind.miss.i <- which(is.na(X[i,-1]))
      accept.i <- rep(0,p)
      if(length(ind.miss.i)>0){
        for(j in ind.miss.i){
          mu.miss.ij <- muX.now[phi.i,j] - 1/QX.now[phi.i,j,j]*sum(QX.now[phi.i,j,-j]*(z.now.i[-j] - muX.now[phi.i,-j]))
          z.miss.p <- rnorm(1, mu.miss.ij, sqrt(mult2/QX.now[phi.i,j,j]))
	  if(bin.X[j]==1)
	    xc.miss.p <- as.numeric(z.miss.p > 0)
	  else if(discrete.X[j]==1)
	    xc.miss.p <- udis[[j]][max(1,sum(udis[[j]]<=z.miss.p))]
	  else
	    xc.miss.p <- z.miss.p
	  if(xc.miss.p < limits.X[j,1])
	    xc.miss.p <- limits.X[j,1]	    
	  if(xc.miss.p > limits.X[j,2])
	    xc.miss.p <- limits.X[j,2]
          xc.p.i <- xc.now.i
          z.p.i <- z.now.i
          xc.p.i[j+1] <- xc.miss.p
          z.p.i[j] <- z.miss.p
	  if(mult2==1){
	    dZ.p <- dZ.now <- pi.p <- pi.now <- 0
	  }
	  else{
            dZ.p <- dnorm(z.miss.p, mu.miss.ij, sqrt(mult2/QX.now[phi.i,j,j]), log=TRUE)
            dZ.now <- dnorm(z.now.i[j], mu.miss.ij, sqrt(mult2/QX.now[phi.i,j,j]), log=TRUE)
            pi.p <- dnorm(z.miss.p, mu.miss.ij, sqrt(1/QX.now[phi.i,j,j]), log=TRUE)
            pi.now <- dnorm(z.now.i[j], mu.miss.ij, sqrt(1/QX.now[phi.i,j,j]), log=TRUE)
	  }
          likey.p <- get.1obs.like(y[i], outcome[i], xc.p.i, beta.now, kappa.now)
          likey.now <- get.1obs.like(y[i], outcome[i], xc.now.i, beta.now, kappa.now)
          MH.ratio <- exp(likey.p + pi.p - dZ.p - likey.now - pi.now + dZ.now)
          if(runif(1) < MH.ratio){
            xc.now.i <- xc.p.i
	    z.now.i <- z.p.i
            accept.i[j] <- 1
          }
	  else
	    accept.i[j] <- 0
        }
      }
    ## for each binary, non-missing j, draw Zij | rest 
      ind.bin.nomiss.i <- which(!is.na(X[i,-1]) & bin.X==1)
      if(length(ind.bin.nomiss.i)>0){
        for(j in ind.bin.nomiss.i){
          mu.ij <- muX.now[phi.i,j] - 1/QX.now[phi.i,j,j]*sum(QX.now[phi.i,j,-j]*(z.now.i[-j]-muX.now[phi.i,-j]))
          s2.ij <- 1/QX.now[phi.i,j,j]
          z.now.i[j] <- myrtruncnorm(1, mu.ij, sqrt(s2.ij), T=0, greater=as.logical(xc.now.i[j+1]))
        }
      }
    ## for each discrete, ordinal, non-missing j, draw Zij | rest 
      ind.dis.nomiss.i <- which(!is.na(X[i,-1]) & discrete.X==1)
      if(length(ind.dis.nomiss.i)>0){
        for(j in ind.dis.nomiss.i){
          mu.ij <- muX.now[phi.i,j] - 1/QX.now[phi.i,j,j]*sum(QX.now[phi.i,j,-j]*(z.now.i[-j]-muX.now[phi.i,-j]))
          s2.ij <- 1/QX.now[phi.i,j,j]
          indu <- which(udis[[j]]==xc.now.i[j+1])
          aa <- ifelse(indu==1, -Inf, xc.now.i[j+1])
          bb <- udis[[j]][indu+1]
          z.now.i[j] <- rtruncnorm(1,a=aa, b=bb, mu.ij, sqrt(s2.ij))
        }
      }

     ## sample latent Zij for observed, but left or right censored X's
      ind.lcens <- X[i,-1]<=limits.X[,1]
      ind.rcens <- X[i,-1]>=limits.X[,2]
      ind.lrobs.i <- which(!discrete.X & (ind.lcens | ind.rcens))
      nlrobs.i <- length(ind.lrobs.i)
      for(l in index(1,nlrobs.i)){
        j <- ind.lrobs.i[l]
        mu.ij <- muX.now[phi.i,j] - 1/QX.now[phi.i,j,j]*sum(QX.now[phi.i,j,-j]*(z.now.i[-j]-muX.now[phi.i,-j]))
        s2.ij <- 1/QX.now[phi.i,j,j]
	aa <- ifelse(ind.rcens[j], limits.X[j,2], -Inf)
	bb <- ifelse(ind.lcens[j], limits.X[j,1], Inf)
        z.now.i[j] <- rtruncnorm(1,a=aa, b=bb, mu.ij, sqrt(s2.ij))
      }

      Xc.now.b[k,] <- xc.now.i
      Z.now.b[k,] <- z.now.i
      accept.b[k,] <- accept.i
      phi.now.b[k] <- phi.i
    }
    cbind(Xc.now.b, Z.now.b, accept.b, phi.now.b)
  }
  Xc.now <- XZ[,1:(p+1)]
  Z.now <- XZ[,(p+2):(2*p+1)]
  accept <- XZ[,(2*p+2):(3*p+1)]
  phi.now <- XZ[,3*p+2]
  accept <- accept[is.na(X[,-1])]

  return(list(Xc.now=Xc.now, Z.now=Z.now, phi.now=phi.now, accept=accept))
}






## Generate a Proposal for MH update of gamma
rgamma.prop <- function(gamma.now, p.ad, p.swap){

  p <- length(gamma.now)
  groups <- rep(1,p)
  gamma.p <- gamma.now
  probs <- c(p.ad, p.swap, 1-p.ad-p.swap)
  probs[probs<0] <- 0
  move <- sample(1:3, 1, prob=probs)
  ng.vec <- as.numeric(table(groups))
  G <- length(ng.vec)
  g <- sample(1:G, 1, prob=ng.vec)
  ind.g <- which(groups==g)
  ncv <- sum(gamma.now==1)
  ind0.g <- which(gamma.now==0 & groups==g)
  ind1.g <- which(gamma.now==1 & groups==g)
  if(move==1){
    if(length(ind1.g)==1)
      ind <- sample(ind0.g, 1)
    else if(ncv==(p-1))
      ind <- sample(ind1.g, 1)
    else
      ind <- sample(ind.g, 1)
    gamma.p[ind] <- 1-gamma.now[ind]
  }
  if(move==2){
    ind120 <- ifelse(length(ind1.g)==1, ind1.g, sample(ind1.g, 1))
    ind021 <- ifelse(length(ind0.g)==1, ind0.g, sample(ind0.g, 1))
    gamma.p[ind120] <- 0
    gamma.p[ind021] <- 1
  }
  return(gamma.p)
}



## evaluate density of the proposal for gamma
dgamma.prop <- function(gamma.p, gamma.now, p.ad, p.swap){

  p <- length(gamma.now)
  groups <- rep(1,p)
  ndiff <- sum(abs(gamma.p-gamma.now))
  if(ndiff==0)
    ans <- log(1-p.ad-p.swap)
  else{
    ng.vec <- as.numeric(table(groups))
    G <- length(ng.vec)
    g <- groups[which(gamma.now!=gamma.p)][1]
    p.g <- ng.vec[g]/sum(ng.vec)
    ind.g <- which(groups==g)
    ncv <- sum(gamma.now==1)
    ind0.g <- which(gamma.now==0 & groups==g)
    ind1.g <- which(gamma.now==1 & groups==g)
    if(ndiff==2)
      ans <- log(p.swap*p.g/length(ind1.g)/length(ind0.g))
    if(ndiff==1){
      if(length(ind1.g)==1 || length(ind0.g)==1)
        ans <- p.ad*p.g/(ng.vec[g]-1)
      else
        ans <- p.ad*p.g/ng.vec[g]
    }
  }
  return(ans)
}



## Generate a Proposal for MH update of mu and Q
rmuQ <- function(Z.now, phi.now, gamma.now, nu.now, lambda.now, eta.now, Psi.now, M, ind.scale, muX.now=NULL, QX.now=NULL, use.now=NULL){

  n <- nrow(Z.now)
  p <- ncol(Z.now)
  c.vars <- (gamma.now==1)
  c1.vars <- c.vars & !ind.scale
  c2.vars <- c.vars & ind.scale
  ind1.1 <- c1.vars[c.vars]
  ind1.2 <- c2.vars[c.vars]
  nc.vars <- !c.vars
  y1 <- Z.now[,c.vars,drop=FALSE]
  y1.1 <- Z.now[,c1.vars,drop=FALSE]
  y1.2 <- Z.now[,c2.vars,drop=FALSE]
  y2 <- Z.now[,nc.vars,drop=FALSE]
  p1 <- ncol(y1)
  p1.1 <- ncol(y1.1)
  p1.2 <- ncol(y1.2)
  p2 <- ncol(y2)
  P11 <- Psi.now[c.vars,c.vars,drop=FALSE]
  P11.11 <- Psi.now[c1.vars,c1.vars,drop=FALSE]
  ans.P11.11 <- ginv.gp(P11.11)
  P11.22 <- Psi.now[c2.vars,c2.vars,drop=FALSE]
  P11.12 <- Psi.now[c1.vars,c2.vars,drop=FALSE]
  P11.22.1 <- P11.22 - t(P11.12)%*%ans.P11.11$inv%*%P11.12  
  ans.P11.22.1 <- ginv.gp(P11.22.1)
  P12 <- Psi.now[c.vars,nc.vars,drop=FALSE]
  P22 <- Psi.now[nc.vars,nc.vars,drop=FALSE]
  ans.P11 <- ginv.gp(P11)
  P22.1 <- P22 - t(P12)%*%ans.P11$inv%*%P12
  ans.P22.1 <- ginv.gp(P22.1)
  y1bar <- apply(y1,2,mean)
  y1barM <- matrix(y1bar,n,p1,byrow=TRUE)
  y2bar <- apply(y2,2,mean)
  y2barM <- matrix(y2bar,n,p2,byrow=TRUE)
  V11 <- crossprod(y1-y1barM) + n*lambda.now/(n+lambda.now)*y1bar%*%t(y1bar) + P11
  V22 <- crossprod(y2-y2barM) + n*lambda.now/(n+lambda.now)*y2bar%*%t(y2bar) + P22
  V12 <- crossprod(y1-y1barM,y2-y2barM) + n*lambda.now/(n+lambda.now)*y1bar%*%t(y2bar) + P12
  ans.V11 <- ginv.gp(V11)
  V22.1 <- V22 - t(V12)%*%ans.V11$inv%*%V12
  ans.V22.1 <- ginv.gp(V22.1)
  ldet.P11 <- ans.P11$log.det
  ldet.V11 <- ans.V11$log.det
  ldet.P22.1 <- ans.P22.1$log.det
  ldet.V22.1 <- ans.V22.1$log.det
  ldet.vec <- ldet11.vec <- rep(0,M)
  if(is.null(use.now))
    use.now <- rep(!is.null(muX.now),5)
  
 ## Draw Q22
  if(use.now[1])
    Q22 <- QX.now[1,nc.vars,nc.vars]
  else
    Q22 <- rwish(n+eta.now, ans.V22.1$inv)
  ans.Q22 <- ginv.gp(Q22)
  dens <- dwish(Q22, n+eta.now, ans.V22.1$inv)
  dprior <- dwish(Q22, eta.now, ans.P22.1$inv)

  if(use.now[2])
    Q21 <- matrix(QX.now[1,nc.vars,c.vars],sum(nc.vars),sum(c.vars))
  else
    Q21 <- rmatnorm(-Q22%*%t(V12)%*%ans.V11$inv, U.sqrt=ans.Q22$sqrt, V.sqrt=ans.V11$sqrt.inv)
  dens <- dens + dmatnorm(Q21, -Q22%*%t(V12)%*%ans.V11$inv, ans.U=ans.Q22, ans.V=list(inv=V11, log.det=-ans.V11$log.det))
  dprior <- dprior + dmatnorm(Q21, -Q22%*%t(P12)%*%ans.P11$inv, ans.U=ans.Q22, ans.V=list(inv=P11, log.det=-ans.P11$log.det))

  mu.b2.s <- as.numeric(n/(n+lambda.now)*(Q22%*%y2bar+Q21%*%y1bar))
  S.b2.sqrt <- 1/sqrt(n+lambda.now)*ans.Q22$sqrt
  S.b2.inv <- (n+lambda.now)*ans.Q22$inv
  log.det.b2 <- -p2*log(n+lambda.now) + ans.Q22$log.det
  Sp.b2.inv <- (lambda.now)*ans.Q22$inv
  log.detp.b2 <- -p2*log(lambda.now) + ans.Q22$log.det
  if(use.now[3])
    b2 <- as.numeric(Q21%*%muX.now[1,c.vars] + Q22%*%muX.now[1,nc.vars])
  else
    b2 <- rmvnorm(1,mu.b2.s,S.sqrt=S.b2.sqrt)
  dens <- dens + dmvnorm(b2, mu.b2.s, S.inv=S.b2.inv, log.det=log.det.b2)
  dprior <- dprior + dmvnorm(b2, rep(0,p2), S.inv=Sp.b2.inv, log.det=log.detp.b2)

  Q <- array(0,c(M,p,p))
  mu <- matrix(0,M,p)
  V11.11 <- array(0,c(M,p1.1,p1.1))
  V11.12 <- array(0,c(M,p1.1,p1.2))
  V11.22.1 <- array(0,c(M,p1.2,p1.2))
  ym1.1bar <- matrix(0,M,p1.1)
  ym1.2bar <- matrix(0,M,p1.2)
  ans.V11.11 <- list()
  nvec <- rep(0,M)
  for(m in 1:M){
    ym1.1 <- y1.1[phi.now==m,,drop=FALSE]
    ym1.2 <- y1.2[phi.now==m,,drop=FALSE]
    nm <- nvec[m] <- nrow(ym1.1)
    if(nm > 0){
      ym1.1bar[m,] <- apply(ym1.1,2,mean)
      ym1.1barM <- matrix(ym1.1bar[m,],nm,p1.1,byrow=TRUE)
      ym1.2bar[m,] <- apply(ym1.2,2,mean)
      ym1.2barM <- matrix(ym1.2bar[m,],nm,p1.2,byrow=TRUE)
      V11.11[m,,] <- crossprod(ym1.1-ym1.1barM) + nm*lambda.now/(nm+lambda.now)*ym1.1bar[m,]%*%t(ym1.1bar[m,]) + P11.11
      Vm11.22 <- crossprod(ym1.2-ym1.2barM) + nm*lambda.now/(nm+lambda.now)*ym1.2bar[m,]%*%t(ym1.2bar[m,]) + P11.22
      V11.12[m,,] <- crossprod(ym1.1-ym1.1barM,ym1.2-ym1.2barM)+nm*lambda.now/(nm+lambda.now)*ym1.1bar[m,]%*%t(ym1.2bar[m,])+P11.12
    }
    else{
      ym1.1bar[m,] <- rep(0,p1.1)
      ym1.2bar[m,] <- rep(0,p1.2)
      V11.11[m,,] <- P11.11
      Vm11.22 <- P11.22
      V11.12[m,,] <- P11.12
    }
    ans.V11.11[[m]] <- ginv.gp(V11.11[m,,])
    V11.22.1[m,,] <- Vm11.22 - t(matrix(V11.12[m,,],p1.1,p1.2))%*%ans.V11.11[[m]]$inv%*%matrix(V11.12[m,,],p1.1,p1.2)
  }

  P11.22.s <- apply(V11.22.1,c(2,3),sum) - (M-1)*P11.22.1
  ans.P11.22.s <- ginv.gp(P11.22.s)
  if(use.now[4]){
    S.m11.inv <- QX.now[1,c.vars,c.vars] - t(Q21)%*%ans.Q22$inv%*%Q21
    Q11.22 <- S.m11.inv[ind1.2,ind1.2,drop=FALSE]
  }
  else{
    Q11.22 <- rwish(n+eta.now-p2, ans.P11.22.s$inv)
  }
  ans.Q11.22 <- ginv.gp(Q11.22)
  dens <- dens + dwish(Q11.22, n+eta.now-p2, ans.P11.22.s$inv)
  dprior <- dprior + dwish(Q11.22, eta.now-p2, ans.P11.22.1$inv)
  
  for(m in 1:M){
    nm <- nvec[m]
    lambda.s <- nm + lambda.now
    if(use.now[5]){
      S.m11.inv <- QX.now[m,c.vars,c.vars] - t(Q21)%*%ans.Q22$inv%*%Q21
      Qm11.21 <- S.m11.inv[ind1.2,ind1.1,drop=FALSE]
      Qm11.11 <- S.m11.inv[ind1.1,ind1.1,drop=FALSE]
      mu.bm1.2.s <- as.numeric(nm/(nm+lambda.now)*(Q11.22%*%ym1.2bar[m,] + Qm11.21%*%ym1.1bar[m,]))
      bm1.2 <- as.numeric(Qm11.21%*%muX.now[m,c1.vars] + Q11.22%*%muX.now[m,c2.vars])
      mu.m1.1.s <- nm/(nm+lambda.now)*ym1.1bar[m,]
      mu.m1.1 <- muX.now[m,c1.vars]
      S.m11.11.inv <- Qm11.11 - t(Qm11.21)%*%ans.Q11.22$inv%*%Qm11.21
    }
    else{
      Qm11.21 <- rmatnorm(-Q11.22%*%t(matrix(V11.12[m,,],p1.1,p1.2))%*%ans.V11.11[[m]]$inv, U.sqrt=ans.Q11.22$sqrt, V.sqrt=ans.V11.11[[m]]$sqrt.inv)
      mu.bm1.2.s <- as.numeric(nm/(nm+lambda.now)*(Q11.22%*%ym1.2bar[m,] + Qm11.21%*%ym1.1bar[m,]))
      S.bm1.2.sqrt <- 1/sqrt(nm+lambda.now)*ans.Q11.22$sqrt
      bm1.2 <- rmvnorm(1,mu.bm1.2.s,S.sqrt=S.bm1.2.sqrt)
      mu.m1.1.s <- nm/(nm+lambda.now)*ym1.1bar[m,]
      ans.m <- rniw(mu.m1.1.s,lambda.s,matrix(V11.11[m,,],p1.1,p1.1), nm + eta.now-p2-p1.2)
      mu.m1.1 <- ans.m$mu
      Qm11.11 <- ans.m$Q + t(Qm11.21)%*%ans.Q11.22$inv%*%Qm11.21
      S.m11.11.inv <- ans.m$Q
      S.m11.inv <- matrix(0,p1,p1)
      S.m11.inv[ind1.1,ind1.1] <- Qm11.11
      S.m11.inv[ind1.1,ind1.2] <- t(Qm11.21)
      S.m11.inv[ind1.2,ind1.1] <- Qm11.21
      S.m11.inv[ind1.2,ind1.2] <- Q11.22
    }
    S.bm1.2.inv <- (nm+lambda.now)*ans.Q11.22$inv
    log.det.bm1.2 <- -p1.2*log(nm+lambda.now) + ans.Q11.22$log.det
    Sp.bm1.2.inv <- (lambda.now)*ans.Q11.22$inv
    log.detp.bm1.2 <- -p1.2*log(lambda.now) + ans.Q11.22$log.det
    
    dens <- dens + dmatnorm(Qm11.21, -Q11.22%*%t(matrix(V11.12[m,,],p1.1,p1.2))%*%ans.V11.11[[m]]$inv, ans.U=ans.Q11.22, ans.V=list(inv=matrix(V11.11[m,,],p1.1,p1.1), log.det=-ans.V11.11[[m]]$log.det))
    dprior <- dprior + dmatnorm(Qm11.21, -Q11.22%*%t(P11.12)%*%ans.P11.11$inv, ans.U=ans.Q11.22, ans.V=list(inv=P11.11, log.det=-ans.P11.11$log.det))
    
    dens <- dens + dmvnorm(bm1.2, mu.bm1.2.s, S.inv=S.bm1.2.inv, log.det=log.det.bm1.2)
    dprior <- dprior + dmvnorm(bm1.2, rep(0,p1.2), S.inv=Sp.bm1.2.inv, log.det=log.detp.bm1.2)
    
    dens <- dens + dniw(mu.m1.1, S.m11.11.inv, mu.m1.1.s, lambda.s, matrix(V11.11[m,,],p1.1,p1.1), nm + eta.now-p2-p1.2)
    dprior <- dprior + dniw(mu.m1.1, S.m11.11.inv, rep(0,p1.1), lambda.now, P11.11, eta.now-p2-p1.2)
    
    mu.m1.2 <- ans.Q11.22$inv%*%(bm1.2-Qm11.21%*%mu.m1.1)
    mu[m,c1.vars] <- mu.m1.1
    mu[m,c2.vars] <- mu.m1.2
    mu.m2 <- ans.Q22$inv%*%(b2-Q21%*%mu[m,c.vars])
    mu[m,nc.vars] <- mu.m2

    Q[m,c.vars,c.vars] <- S.m11.inv + t(Q21)%*%ans.Q22$inv%*%Q21
    Q[m,c.vars,nc.vars] <- t(Q21)
    Q[m,nc.vars,c.vars] <- Q21
    Q[m,nc.vars,nc.vars] <- Q22

    ldet11.vec[m] <- -ginv.gp(S.m11.11.inv)$log.det - ans.Q11.22$log.det
    ldet.vec[m] <- ldet11.vec[m] - ans.Q22$log.det
  }
  return(list(mu=mu, Q=Q, dens=dens, dprior=dprior, ldet.vec=ldet.vec, ldet11.vec=ldet11.vec))
}



## evaluate prior density for mu and Q
get.pi.muQ <- function(gamma.now, nu.now, lambda.now, eta.now, Psi.now, M, ind.scale, muX.now, QX.now){

  p <- length(nu.now)
  c.vars <- (gamma.now==1)
  c1.vars <- c.vars & !ind.scale
  c2.vars <- c.vars & ind.scale
  ind1.1 <- c1.vars[c.vars]
  ind1.2 <- c2.vars[c.vars]
  nc.vars <- !c.vars
  p1 <- sum(c.vars)
  p1.1 <- sum(c1.vars)
  p1.2 <- sum(c2.vars)
  p2 <- sum(nc.vars)
  P11 <- Psi.now[c.vars,c.vars,drop=FALSE]
  P11.11 <- Psi.now[c1.vars,c1.vars,drop=FALSE]
  ans.P11.11 <- ginv.gp(P11.11)
  P11.22 <- Psi.now[c2.vars,c2.vars,drop=FALSE]
  P11.12 <- Psi.now[c1.vars,c2.vars,drop=FALSE]
  P11.22.1 <- P11.22 - t(P11.12)%*%ans.P11.11$inv%*%P11.12  
  ans.P11.22.1 <- ginv.gp(P11.22.1)
  P12 <- Psi.now[c.vars,nc.vars,drop=FALSE]
  P22 <- Psi.now[nc.vars,nc.vars,drop=FALSE]
  ans.P11 <- ginv.gp(P11)
  P22.1 <- P22 - t(P12)%*%ans.P11$inv%*%P12
  ans.P22.1 <- ginv.gp(P22.1)
  ldet.P11 <- ans.P11$log.det
  ldet.P22.1 <- ans.P22.1$log.det
  
 ## Draw Q22
  Q22 <- QX.now[1,nc.vars,nc.vars]
  ans.Q22 <- ginv.gp(Q22)
  dprior <- dwish(Q22, eta.now, ans.P22.1$inv)

  Q21 <- matrix(QX.now[1,nc.vars,c.vars],sum(nc.vars),sum(c.vars))
  dprior <- dprior + dmatnorm(Q21, -Q22%*%t(P12)%*%ans.P11$inv, ans.U=ans.Q22, ans.V=list(inv=P11, log.det=-ans.P11$log.det))

  Sp.b2.inv <- (lambda.now)*ans.Q22$inv
  log.detp.b2 <- -p2*log(lambda.now) + ans.Q22$log.det
  b2 <- as.numeric(Q21%*%muX.now[1,c.vars] + Q22%*%muX.now[1,nc.vars])
  dprior <- dprior + dmvnorm(b2, rep(0,p2), S.inv=Sp.b2.inv, log.det=log.detp.b2)

  S.m11.inv <- QX.now[1,c.vars,c.vars] - t(Q21)%*%ans.Q22$inv%*%Q21
  Q11.22 <- S.m11.inv[ind1.2,ind1.2,drop=FALSE]
  ans.Q11.22 <- ginv.gp(Q11.22)
  dprior <- dprior + dwish(Q11.22, eta.now-p2, ans.P11.22.1$inv)
  
  for(m in 1:M){
    S.m11.inv <- QX.now[m,c.vars,c.vars] - t(Q21)%*%ans.Q22$inv%*%Q21
    Qm11.21 <- S.m11.inv[ind1.2,ind1.1,drop=FALSE]
    Qm11.11 <- S.m11.inv[ind1.1,ind1.1,drop=FALSE]
    bm1.2 <- as.numeric(Qm11.21%*%muX.now[m,c1.vars] + Q11.22%*%muX.now[m,c2.vars])
    mu.m1.1 <- muX.now[m,c1.vars]
    S.m11.11.inv <- Qm11.11 - t(Qm11.21)%*%ans.Q11.22$inv%*%Qm11.21

    Sp.bm1.2.inv <- (lambda.now)*ans.Q11.22$inv
    log.detp.bm1.2 <- -p1.2*log(lambda.now) + ans.Q11.22$log.det
    
    dprior <- dprior + dmatnorm(Qm11.21, -Q11.22%*%t(P11.12)%*%ans.P11.11$inv, ans.U=ans.Q11.22, ans.V=list(inv=P11.11, log.det=-ans.P11.11$log.det))
    
    dprior <- dprior + dmvnorm(bm1.2, rep(0,p1.2), S.inv=Sp.bm1.2.inv, log.det=log.detp.bm1.2)

    dprior <- dprior + dniw(mu.m1.1, S.m11.11.inv, rep(0,p1.1), lambda.now, P11.11, eta.now-p2-p1.2)    
  }
  return(dprior)
}




## Gibbs update for muX and QX, and gammaX
update.muX.QX.gamma <- function(gamma.now, muX.now, QX.now, Z.now, p.ad, p.swap, rhoX, phi.now, nu.now, lambda.now, eta.now, Psi.now, ind.scale, nb, it){

  M <- nrow(muX.now)
  if(it==1)
    gamma.p <- gamma.now
  else
    gamma.p <- rgamma.prop(gamma.now, p.ad, p.swap)
  d.gamma.p <- dgamma.prop(gamma.p, gamma.now, p.ad, p.swap)
  d.gamma.now <- dgamma.prop(gamma.now, gamma.p, p.ad, p.swap)
  pi.gamma.p <- sum(gamma.p*log(rhoX))+sum((1-gamma.p)*log(1-rhoX))
  pi.gamma.now <- sum(gamma.now*log(rhoX))+sum((1-gamma.now)*log(1-rhoX))

  ans.muQ.p <- rmuQ(Z.now, phi.now, gamma.p, nu.now, lambda.now, eta.now, Psi.now, M, ind.scale, muX.now=NULL, QX.now=NULL)
  muX.p <- ans.muQ.p$mu
  QX.p <- ans.muQ.p$Q
  d.muQ.p <- ans.muQ.p$dens
  pi.muQ.p <- ans.muQ.p$dprior
  ldet.p <- ans.muQ.p$ldet.vec
  ldet11.p <- ans.muQ.p$ldet11.vec

  ans.muQ.now <- rmuQ(Z.now, phi.now, gamma.now, nu.now, lambda.now, eta.now, Psi.now, M, ind.scale, muX.now=muX.now, QX.now=QX.now)
  d.muQ.now <- ans.muQ.now$dens
  pi.muQ.now <- ans.muQ.now$dprior
  ldet.now <- ans.muQ.now$ldet.vec
  ldet11.now <- ans.muQ.now$ldet11.vec

  like.p <- get.like.Z(Z.now, muX.p, QX.p, ldet.p, phi.now, nb)
  like.now <- get.like.Z(Z.now, muX.now, QX.now, ldet.now, phi.now, nb)

  pi.p <- pi.muQ.p + pi.gamma.p
  pi.now <- pi.muQ.now + pi.gamma.now
  d.p <- d.gamma.p + d.muQ.p
  d.now <- d.gamma.now + d.muQ.now

  accept <- 0
  if(it==1)
    MH.ratio <- 1
  else
    MH.ratio <- exp(like.p + pi.p - d.p - like.now - pi.now + d.now)
  if(runif(1) < MH.ratio){
    gamma.now <- gamma.p
    muX.now <- muX.p
    QX.now <- QX.p
    ldet.now <- ldet.p
    ldet11.now <- ldet11.p
    accept <- 1
  }
  return(list(gamma.now=gamma.now, muX.now=muX.now, QX.now=QX.now, ldet.now=ldet.now, ldet11.now=ldet11.now, accept=accept))
}



## evaluate density of current Z
get.like.Z <- function(Z.now, mu.now, QX.now, ldet.now, phi.now, nb){

  M <- max(phi.now)
  n <- nrow(Z.now)
  
## for each missing j, draw Zij (and effectively Xij) | rest, and accept/reject
  blocks <- list()
  ind.now <- 1
  inc <- ceiling(n/nb)
  for(b in 1:nb){
    blocks[[b]] <- ind.now:min(ind.now+inc-1, n)
    ind.now <- ind.now+inc
  }

  like <- foreach(b=1:nb, .combine=sum)%dopar%{
    like.b <- 0
    for(m in 1:M){
      ind.m <- blocks[[b]][phi.now[blocks[[b]]]==m]
      if(length(ind.m)>0)
        like.b <- like.b + dmvnorm(Z.now[ind.m,], mu.now[m,], S.inv=QX.now[m,,], log.det=ldet.now[m])
    }
    like.b
  }
  return(like)
}



## Gibbs update for vee (omega)
update.vee <- function(phi.now, delta.now, M){

  if(M==1){
    vee <- 0
    omega <- 1
  }
  else{
    vee <- rep(0,M)
    for(m in 1:M){
      a.star <- sum(phi.now==m)+1
      b.star <- sum(phi.now>m)+delta.now
      vee[m] <- rbeta(1,b.star, a.star)
    }
    vee[M] <- 0
    omega <- makeprobs(vee)
  }
  return(list(vee=vee,omega=omega))
}


## Gibbs update for delta
update.delta <- function(vee.now, A.delta, B.delta){

  M <- length(vee.now)
  a.star <- M + A.delta
  b.star <- -sum(log(vee.now[-M])) + B.delta
  delta <- rgamma(1, a.star, b.star)
  return(delta)
}





## MH update for psi
update.Psi <- function(Psi.now, muX.now, QX.now, A.Psi, B.Psi, eta.now, nu.now, lambda.now, gamma.now, prop.sd.Psi, dfp, ind.scale){

  accept <- 0
  M <- dim(QX.now)[1]
  p <- dim(QX.now)[2]
  Psi.p <- Psi.now
  psi.now <- Psi.now[1,1]
  psi.p <- r.pos.proposal(psi.now, df=dfp, sigma=prop.sd.Psi)
  Psi.p <- diag(psi.p,p)
  ans.muQ.p <- try(get.pi.muQ(gamma.now, nu.now, lambda.now, eta.now, Psi.p, M, ind.scale, muX.now, QX.now), silent=FALSE)
  if(is.character(ans.muQ.p[1]))
    like.p <- -Inf
  else
    like.p <- ans.muQ.p
  ans.muQ.now <- try(get.pi.muQ(gamma.now, nu.now, lambda.now, eta.now, Psi.now, M, ind.scale, muX.now, QX.now), silent=TRUE)
  if(is.character(ans.muQ.now[1]))
    like.now <- -1E300
  else
    like.now <- ans.muQ.now

  pi.p <- dgamma(psi.p, A.Psi, B.Psi, log=TRUE)
  pi.now <- dgamma(psi.now, A.Psi, B.Psi, log=TRUE)
  dprop.p <- d.pos.proposal(psi.p, psi.now, df=dfp, sigma=prop.sd.Psi)
  dprop.now <- d.pos.proposal(psi.now, psi.p, df=dfp, sigma=prop.sd.Psi)
  MH.ratio <- exp((like.p + pi.p - dprop.p) - (like.now + pi.now - dprop.now))
  if(runif(1) < MH.ratio){
    Psi.now <- Psi.p
    accept <- 1
  }
  return(list(Psi.now=Psi.now, accept=accept))
}




## MH update for eta
update.eta <- function(muX.now, QX.now, gamma.now, nu.now, lambda.now, eta.now, Psi.now, A.eta, B.eta, prop.sd.eta, dfp, etagp, ind.scale){
 
  accept <- 0
  M <- dim(QX.now)[1]
  p <- dim(QX.now)[2]
  c.vars <- (gamma.now==1)
  nc.vars <- !c.vars
  p1 <- sum(c.vars)
  p2 <- p-p1

  eta.p <- r.pos.proposal(eta.now-p-etagp, df=dfp, sigma=prop.sd.eta) + p + etagp

  ans.muQ.p <- try(get.pi.muQ(gamma.now, nu.now, lambda.now, eta.p, Psi.now, M, ind.scale, muX.now, QX.now), silent=TRUE)
  if(is.character(ans.muQ.p[1]))
    like.p <- -Inf
  else
    like.p <- ans.muQ.p
  ans.muQ.now <- try(get.pi.muQ(gamma.now, nu.now, lambda.now, eta.now, Psi.now, M, ind.scale, muX.now, QX.now), silent=TRUE)
  if(is.character(ans.muQ.now[1]))
    like.now <- -1E300
  else
    like.now <- ans.muQ.now
  pi.p <- dgamma(eta.p-p-etagp, A.eta, B.eta, log=TRUE)
  pi.now <- dgamma(eta.now-p-etagp, A.eta, B.eta, log=TRUE)
  dprop.p <- d.pos.proposal(eta.p-p-etagp, eta.now-p-etagp, df=dfp, sigma=prop.sd.eta)
  dprop.now <- d.pos.proposal(eta.now-p-etagp, eta.p-p-etagp, df=dfp, sigma=prop.sd.eta)
  MH.ratio <- exp((like.p + pi.p - dprop.p) - (like.now + pi.now - dprop.now))
  if(runif(1) < MH.ratio){
    eta.now <- eta.p
    accept <- 1
  }
  return(list(eta.now=eta.now, accept=accept))
}




## MH update for lambda
update.lambda <- function(muX.now, QX.now, ldet11.now, gamma.now, nu.now, lambda.now, eta.now, Psi.now, A.lambda, B.lambda, prop.sd.lambda, dfp, ind.scale){
 
  accept <- 0
  M <- dim(QX.now)[1]
  p <- dim(QX.now)[2]
  c.vars <- (gamma.now==1)
  nc.vars <- !c.vars
  p1 <- sum(c.vars)
  p2 <- p-p1

  lambda.p <- r.pos.proposal(lambda.now, df=dfp, sigma=prop.sd.lambda)

  ans.muQ.p <- try(get.pi.muQ(gamma.now, nu.now, lambda.p, eta.now, Psi.now, M, ind.scale, muX.now, QX.now), silent=TRUE)
  if(is.character(ans.muQ.p[1]))
    like.p <- -Inf
  else
    like.p <- ans.muQ.p
  ans.muQ.now <- try(get.pi.muQ(gamma.now, nu.now, lambda.now, eta.now, Psi.now, M, ind.scale, muX.now, QX.now), silent=TRUE)
  if(is.character(ans.muQ.now[1]))
    like.now <- -1E300
  else
    like.now <- ans.muQ.now
  
  pi.p <- dgamma(lambda.p, A.lambda, B.lambda, log=TRUE)
  pi.now <- dgamma(lambda.now, A.lambda, B.lambda, log=TRUE)
  dprop.p <- d.pos.proposal(lambda.p, lambda.now, df=dfp, sigma=prop.sd.lambda)
  dprop.now <- d.pos.proposal(lambda.now, lambda.p, df=dfp, sigma=prop.sd.lambda)
  MH.ratio <- exp((like.p + pi.p - dprop.p) - (like.now + pi.now - dprop.now))
  if(runif(1) < MH.ratio){
    lambda.now <- lambda.p
    accept <- 1
  }
  return(list(lambda.now=lambda.now, accept=accept))
}





###############################
##### Utility Functions #######
###############################


## generate a vecotr of fold ids for cross validation
get.foldid <- function(n, nfolds=10, seed=220){

  replace.seed <- T
  if(missing(seed))
    replace.seed <- F

  if(replace.seed){
   ## set seed to specified value
    if(!any(ls(name='.GlobalEnv', all.names=T)=='.Random.seed')){
      set.seed(1)
    }
    save.seed <- .Random.seed
    set.seed(seed)  
  }

  perm <- sample(1:n, n)
  n.cv <- rep(floor(n/nfolds),nfolds)
  rem <- n - n.cv[1]*nfolds
  n.cv[index(1,rem)] <- n.cv[index(1,rem)]+1
  foldid <- rep(0,n)
  ind2 <- 0
  
  for(i in 1:nfolds){
    ind1 <- ind2+1
    ind2 <- ind2+n.cv[i]
    foldid[perm[ind1:ind2]] <- i
  }
  if(replace.seed){
   ## restore random seed to previous value
    .Random.seed <<- save.seed
  }
  return(foldid)
}




## go from vee to omega
makeprobs <- function(v){ 
   #compute the stick-breaking weights
   M <- length(v)
   if(M==1)
     return(1)
   probs <- 1-v
   probs[2:M] <- probs[2:M]*cumprod(v[2:M-1])
return(probs)}


## go from omega to vee
probs2v <- function(probs){ 
   #compute the stick-breaking weights
   M <- length(probs)
   v <- probs
   v[2:M] <- probs[2:M]/(1-cumsum(probs[2:M-1]))
   v[v>=1] <- 1-1E-16
   v[probs==0] <- .5
   return(1-v)
}



## generate truncated normal deviates
myrtruncnorm <- function(n, mu, sd, T, greater=TRUE){

  ans <- rep(0,n)
  if(length(greater)==1)
    greater <- rep(greater,n)
  less <- !greater
  if(length(mu)==1)
    mu <- rep(mu,n)
  if(length(sd)==1)
    sd <- rep(sd,n)
  if(length(T)==1)
    T <- rep(T,n)
  if(any(greater)){
    ng <- sum(greater)
    lp <- log(runif(ng))+pnorm(T[greater],mu[greater],sd[greater],lower.tail=FALSE,log.p=TRUE)
    ans[greater] <- qnorm(lp,mu[greater],sd[greater],log.p=TRUE,lower.tail=FALSE)
  }
  if(any(less)){
    nl <- sum(less)
    lp <- log(runif(nl))+pnorm(T[less],mu[less],sd[less],lower.tail=TRUE,log.p=TRUE)
    ans[less] <- qnorm(lp,mu[less],sd[less],lower.tail=TRUE,log.p=TRUE)
  }
  return(ans)
}


## index function that cretes an empty vector if trying to index from 1:0
index <- function(m,n){
  if(m<=n) return(m:n)
  else return(numeric(0))
}



## Generate a random Wishart deviate
rwish <- function(v, S){
  S <- as.matrix(S)
  if (v < nrow(S)){
      stop(message = "v is less than the dimension of S in rwish().\n")
  }
  if(nrow(S)==0)
    return(matrix(0,0,0))
  p <- nrow(S)
  CC <- chol(S)
  Z <- matrix(0, p, p)
  diag(Z) <- sqrt(rchisq(p, v:(v - p + 1)))
  if(p > 1){
    pseq <- 1:(p - 1)
    Z[rep(p*pseq, pseq) + unlist(lapply(pseq, seq))] <- rnorm(p*(p-1)/2)
  }
  return(crossprod(Z%*%CC))
}



## Generate a random inverse-Wishart deviate
riwish <- function(v, S){
  if(nrow(S)==0)
    return(matrix(0,0,0))
  Sinv <- ginv.gp(S)$inv
  W <- rwish(v, Sinv)
  Winv <- ginv.gp(W)$inv
  return(Winv)
}



## Evaluate Wishart density
dwish <- function (W, v, S, ans.W=NULL, ans.S=NULL) 
{
    if (!is.matrix(S)) 
        S <- matrix(S)
    if (!is.matrix(W)) 
        W <- matrix(W)
    if (v < nrow(S)) {
        stop(message = "v is less than the dimension of S in  dwish()\n\n")
    }
    k <- nrow(S)
    if(k==0)
      return(0)
    gammapart <- 0
    for (i in 1:k) {
        gammapart <- gammapart + lgamma((v + 1 - i)/2)
    }
    logdenom <- gammapart + (v * k/2)*log(2) + (k * (k - 1)/4)*log(pi)
    if(is.null(ans.S))
      ans.S <- ginv.gp(S)
    if(is.null(ans.W))
      ans.W <- ginv.gp(W)
    logdetS <- ans.S$log.det
    logdetW <- ans.W$log.det
    hold <- ans.S$inv%*%W
    tracehold <- sum(hold[row(hold) == col(hold)])
    lognum <- (-v/2)*logdetS + ((v - k - 1)/2)*logdetW - 1/2*tracehold
    return(lognum-logdenom)
}




## Evaluate inverse-Wishart density
diwish <- function (W, v, S, Winv=NULL, logdetW=NULL, logdetS=NULL) 
{
    if (!is.matrix(S)) 
        S <- matrix(S)
    if (!is.matrix(W)) 
        W <- matrix(W)
    if (v < nrow(S)) {
        stop("v is less than the dimension of S in  diwish().\n")
    }
    k <- nrow(S)
    gammapart <- 0
    for (i in 1:k) {
        gammapart <- gammapart + lgamma((v + 1 - i)/2)
    }
    logdenom <- gammapart + (v * k/2)*log(2) + (k * (k - 1)/4)*log(pi)
    if(is.null(logdetS)){
      ans.S <- ginv.gp(S)
      logdetS <- ans.S$log.det
    }    
    if(is.null(logdetW) || is.null(Winv)){
      ans.W <- ginv.gp(W)
      logdetW <- ans.W$log.det
      Winv <- ans.W$inv
    }
    hold <- S%*%Winv
    tracehold <- sum(hold[row(hold) == col(hold)])
    lognum <- (v/2)*logdetS - (v + k + 1)/2*logdetW -1/2*tracehold
    return(lognum - logdenom)
}



## di-gamma function (used in beta updates)
gamma2 <- function(x)
  return(digamma(x)*gamma(x))


## tri-gamma function (used in beta updates)
gamma3 <- function(x)
  return(gamma(x)*(trigamma(x)+1/(gamma(x)^2)*gamma2(x)^2))



## obtain inverse/sqrt, etc of a symmetric PD matrix via Cholesky decomp
ginv.gp <- function (X, eps = 1e-12){
    X <- as.matrix(X)
    if(any(X==Inf))
      return(list(inv = diag(nrow(X)), log.det = Inf))
    if(nrow(X)==0)
      return(list(inv=matrix(0,0,0), log.det=0, sqrt=matrix(0,0,0), sqrt.inv=matrix(0,0,0)))
    L <- chol(X)
    inv <- chol2inv(L)
    sqrt <- t(L)
    sqrt.inv <- backsolve(L,diag(1,nrow(L)))
    log.det <- as.numeric(determinant(L)$modulus)*2
    return(list(inv = inv, log.det = log.det, sqrt=sqrt, sqrt.inv=sqrt.inv))
}




## generate random MV normal deviate
rmvnorm <- function(n=1, mu=0, Sigma, S.sqrt=NULL){

  if(is.null(S.sqrt)){
    ans <- ginv.gp(Sigma)
    S.sqrt <- ans$sqrt
  }
  p <- nrow(S.sqrt)
  if(length(mu)==1) 
    mu <- rep(mu,p)
  X <- matrix(0,n,p)
  for(i in 1:n)
    X[i,] <- as.numeric(S.sqrt%*%rnorm(length(mu),0,1) + mu)
  return(X[1:n,])
}


## Evaluate MV normal density
dmvnorm <- function(x, mu, Sigma, S.inv=NULL, log.det=NULL){

  if(is.null(S.inv[1]) || is.null(log.det)){
    ans <- ginv.gp(Sigma)
    S.inv <- ans$inv
    log.det <- ans$log.det
  }
  if(is.null(dim(x))){
    n <- 1
    p <- length(x)
    x <- matrix(x,1,p)
  }
  else{
    n <- nrow(x)
    p <- ncol(x)
  }
  dens <- -n*p/2*log(2*pi)-n/2*log.det
  for(i in index(1,n))
    dens <- dens - 0.5*as.numeric(t(x[i,]-mu)%*%S.inv%*%(x[i,]-mu))
  return(dens)
}



## generate random normal-inverseWishart deviate
rniw <- function(nu,lambda,P,eta){

  Q <- rwish(eta, ginv.gp(P)$inv)
  ans.Q <- ginv.gp(Q)
  Sigma <- ans.Q$inv
  S.sqrt <-  1/sqrt(lambda)*ans.Q$sqrt.inv
  mu <- rmvnorm(1,nu,S.sqrt=1/sqrt(lambda)*ans.Q$sqrt.inv)
  return(list(mu=mu,Sigma=Sigma,Q=Q))
}


## Evlauate density of normal-inverseWishart
dniw <- function(mu, Q, nu, lambda, P, eta){

  ans.Q <- ginv.gp(Q)
  dens <- dwish(Q, eta, ginv.gp(P)$inv)
  S.inv <- lambda*Q
  p <- dim(Q)[1]
  log.det <- -p*log(lambda) - ans.Q$log.det
  dens <- dens + dmvnorm(mu, nu, S.inv=S.inv, log.det=log.det)
  return(dens)
}




## generate random matrix-normal deviate
rmatnorm <- function(M,U,V,U.sqrt=NULL,V.sqrt=NULL){

  n <- nrow(M)
  p <- ncol(M)
  if(n==0 || p==0)
    return(matrix(0,n,p))
  if(is.null(U.sqrt))
    U.sqrt <- ginv.gp(U)$sqrt
  if(is.null(V.sqrt))
    V.sqrt <- ginv.gp(V)$sqrt
  Z <- matrix(rmvnorm(1,as.numeric(M),S.sqrt=V.sqrt%x%U.sqrt),n,p)
  return(Z)
}



## Evlauate density of matrix-normal
dmatnorm <- function(Z,M,U,V,ans.U=NULL,ans.V=NULL){

  n <- nrow(M)
  p <- ncol(M)
  if(n==0 || p==0)
    return(0)
  if(is.null(ans.U))
    ans.U <- ginv.gp(U)
  if(is.null(ans.V))
    ans.V <- ginv.gp(V)
  S.inv <- ans.V$inv%x%ans.U$inv
  log.det <- n*ans.V$log.det + p*ans.U$log.det
  dens <- dmvnorm(as.numeric(Z),as.numeric(M),S.inv=S.inv, log.det=log.det)
  return(dens)
}



## Evaluate density of a normal mixture (used in phi updates)
dmvmixnorm <- function(z, omega, mu, S.inv, log.det){
  p <- length(z)
  M <- nrow(mu)
  dens <- prob <- rep(0,M)
  for(m in 1:M)
    dens[m] <- log(omega[m])-p/2*log(2*pi)-1/2*log.det[m] - 0.5*as.numeric(t(z-mu[m,])%*%S.inv[m,,]%*%(z-mu[m,]))
  for(m in 1:M){
    if(dens[m]== -Inf)
      prob[m] <- 0
    else
      prob[m] <- 1/(sum(exp(dens - dens[m])))
  }
  return(prob)
}



## evaluate BIC from a kmeans clustering
kmeansBIC <- function(fit){

  m = ncol(fit$centers)
  n = length(fit$cluster)
  k = nrow(fit$centers)
  D = fit$tot.withinss
  return(D + log(n)*m*k)
}


## select a kmeans solution with #components optimizing BIC (used in cluster initialization)
mykmeans <- function(Z, Mvec){

  BIC <- foreach(m=Mvec, .combine=c)%dopar%{
    foo.m <- kmeans(Z,m)
    kmeansBIC(foo.m)
  }
  Mopt <- Mvec[which(BIC==min(BIC))[1]]
  fopt <- kmeans(Z,Mopt)
  return(fopt)
}






############################################################################################
####### Remaining functions are for other methods used in the simulation study #############
############################################################################################

############################################################################################
########### Factorize missingness and use glmnet with missing indicators ###################
############################################################################################

convert2factors <- function(X, bin.X){

  Xf <- X
  Xd <- matrix(0,nrow(X),0)
  p <- ncol(X)
  for(j in 1:p){
    if(bin.X[j]==1){
      Xf[is.na(Xf[,j]),j] <- -1
      Xf[,j] <- Xf[,j]+1
    }
    else{
      cuts <- quantile(X[,j],prob=c(.33,.67), na.rm=TRUE)
      Xf[X[,j]<=cuts[1],j] <- 1
      Xf[X[,j]>cuts[1] & X[,j]<=cuts[2],j] <- 2
      Xf[X[,j]>cuts[2],j] <- 3
      Xf[is.na(X[,j]),j] <- 0
    }
  }
  Xf <- as.data.frame(Xf)
  return(list(Xf=Xf, Xd=Xd))
}





convertframe2factors <- function(X){

  Xf <- X
  p <- ncol(X)
  for(j in 1:p){
    if(any(is.na(Xf[,j]))){
      if(is.factor(Xf[,j])){
        levels(Xf[,j]) <- c(levels(Xf[,j]), "missing")
        Xf[is.na(Xf[,j]),j] <- "missing"
      }
      else{
        cuts <- quantile(X[,j],prob=c(.33,.67), na.rm=TRUE)
        blah <- rep("missing", n)
        blah <- as.factor(blah)
        levels(blah) <- c("1", "2", "3", "missing")
        blah[X[,j]<=cuts[1]] <- "1"
        blah[X[,j]>cuts[1] & X[,j]<=cuts[2]] <- "2"
        blah[X[,j]>cuts[2]] <- "3"
        blah[is.na(X[,j])] <- "missing"
        Xf[,j] <- blah
      }
    }
  }
  Xf <- as.data.frame(Xf)
  return(Xf)
}



add.missing.ind.cols <- function(X){

  Xf <- X
  xnames <- names(Xf)
  p <- ncol(X)
  for(j in 1:p){
    if(any(is.na(Xf[,j]))){
      if(is.factor(Xf[,j])){
        levels(Xf[,j]) <- c(levels(Xf[,j]), "missing")
        Xf[is.na(Xf[,j]),j] <- "missing"
      }
      else{
        Xf[is.na(Xf[,j]),j] <- 0
	Xf <- cbind(Xf, is.na(X[,j])*1)
        xnames <- c(xnames, paste(xnames[j],".miss",sep=""))
      }
    }
  }
  Xf <- as.data.frame(Xf)
  names(Xf) <- xnames
  return(Xf)
}






#############################################################
############# RF interpolation and GLMNET ###################
#############################################################

RFEN <- function(X, y, outcome=NULL, bin.X=NULL, maxiter=5, use=NULL, parallelize="no", family="cox", alpha=.5){

 # Fill in X with RF
  p <- ncol(X)
  if(is.null(bin.X)){
    bin.X <- rep(0,ncol(X))
    for(j in 1:ncol(X))
      bin.X[j] <- is.factor(X[,j])
  }
  Xd <- as.data.frame(X)
  for(j in 1:ncol(X)){
    if(bin.X[j]==1)
      Xd[,j] <- as.factor(Xd[,j])
  }
  if(!is.null(outcome)){
    ind.p <- (outcome==0 & y<1E-6)
    yfoo<- y
    yfoo[ind.p] <- NA
    outcomefoo <- outcome
    outcomefoo[ind.p] <- NA
    outcomefoo <- as.factor(outcomefoo)
    Xd <- cbind(Xd,yfoo,outcomefoo)
  }
  else
    Xd <- cbind(Xd,y)
  if(any(is.na(X)))
    Xh <- missForest(Xd, strata=as.list(c(1:ncol(Xd))), maxiter=maxiter, parallelize=parallelize)$ximp
  else
    Xh <- Xd
  Xh <- matrix(as.numeric(as.matrix(Xh)), nrow(Xh), ncol(Xh))[,1:ncol(X)]

 # Fit glmnet
  if(family=="cox"){
    yS <- cbind(y,outcome)
    colnames(yS) <- c("time", "status")
  }
  else
    yS <- y
  if(is.null(use))
    use <- 1:p
  XXh <- Xh[,use]
  ans.fit <- cv.glmnet(XXh, yS, family=family, alpha=alpha)
  return(list(fit=ans.fit, X=X, Xh=XXh, y=y, outcome=outcome, bin.X=bin.X))
}



predict.RFEN <- function(obj, Xp, maxiter=5, K=1, parallelize="no"){

  np <- nrow(Xp)
  if(K > np)
    K <- np
  p <- ncol(Xp)
  X <- obj$X
  n <- nrow(X)
  bin.X <- obj$bin.X
  
 # Fill in Xp using both X and Xp
  Xd <- as.data.frame(rbind(Xp,X))
  for(j in 1:ncol(X)){
    if(bin.X[j]==1)
      Xd[,j] <- as.factor(Xd[,j])
  }
  ind.list <- list()
  ind.now <- 1
  inc <- ceiling(np/K)
  for(k in 1:K){
    ind.list[[k]] <- ind.now:(ind.now+inc-1)
    ind.now <- ind.now + inc
  }
  Xhp <- matrix(0,np,p)
  for(k in 1:K){
    np.k <- length(ind.list[[k]])
    ind.k <- c(ind.list[[k]],(np+1):(np+n))
    Xd.k <- Xd[ind.k,]
    Xh.k <- missForest(Xd.k, strata=as.list(c(1:ncol(X))), maxiter=maxiter, parallelize=parallelize)$ximp
    Xhp[ind.list[[k]],] <- matrix(as.numeric(as.matrix(Xh.k[1:np.k,])), np.k, p)
  }
  lamhat.RFEN <- as.numeric(predict(obj$fit, Xhp))
  return(lamhat.RFEN)
}







##################################################################
############# Tuning by CV concordance for GBM ###################
##################################################################



cv.gbm <- function(X, y, outcome=NULL, n.trees=seq(150,500,by=10), interaction.depth, n.minobsinnode = 5, shrinkage=10^(seq(-3,-1,length=10)), bag.fraction = .5, distribution="bernoulli", foldid=NULL, nfolds=10, seed=220, verbose=FALSE){

  n <- nrow(X)
  if(is.null(foldid))
    foldid <- get.foldid(n,10, seed=seed)
  nfolds <- max(foldid)
  ns <- length(shrinkage)
  nnt <- length(n.trees)
  concord.mat <- foreach(s=1:ns, .combine=rbind)%dopar%{
    concord.s <- rep(0,nnt)
    yhat <- matrix(0,n,nnt)
    for(k in 1:nfolds){
      X.k <- X[foldid==k,,drop=F]
      y.k <- y[foldid==k]
      X.mk <- X[foldid!=k,,drop=F]
      y.mk <- y[foldid!=k]
      if(distribution=="poisson"){
        outcome.mk <- outcome[foldid!=k]
        ans.sk <- gbm.fit(X.mk, outcome.mk, offset=log(y.mk), n.trees=n.trees[nnt], interaction.depth=interaction.depth, n.minobsinnode=n.minobsinnode, shrinkage=shrinkage[s], bag.fraction=bag.fraction, distribution=distribution, verbose=verbose)
      }
      else{
        ans.sk <- try(gbm.fit(X.mk, y.mk, n.trees=n.trees[nnt], interaction.depth=interaction.depth, n.minobsinnode=n.minobsinnode, shrinkage=shrinkage[s], bag.fraction=bag.fraction, distribution=distribution, verbose=verbose), silent=TRUE)
      }
      if(distribution=="multinomial"){
        for(t in 1:nnt){
	  foo <- predict(ans.sk, X.k, n.trees=n.trees[t],type="response", verbose=verbose)
	  fuh <- foo[cbind(1:length(y.k),as.numeric(y.k),1)]
	  fuh[is.na(fuh)] <- 0
          yhat[foldid==k,t] <- fuh
	}
      }
      else{
        for(t in 1:nnt){
	  if(!is.character(ans.sk[1]))
            yhat[foldid==k,t] <- predict(ans.sk, X.k, n.trees=n.trees[t])
	  else
	    yhat[foldid==k,t] <- 0
	}
      }
    }
    for(t in 1:nnt){
      if(distribution=="coxph")
        concord.s[t] <- survConcordance(y~yhat[,t])$concord
      else if(distribution=="poisson")
        concord.s[t] <- survConcordance(Surv(y,outcome)~yhat[,t])$concord
      else if(distribution=="multinomial")
        concord.s[t] <- sum(log(yhat[,t]))
      else if(distribution=="gaussian")
        concord.s[t] <- 1-sum((y-yhat[,t])^2)/sum((y-mean(y))^2)
      else
        concord.s[t] <- get.ROC(y, yhat[,t])$auc
    }
    concord.s
  }

  inds <- which(concord.mat==max(concord.mat), arr.ind=TRUE)
  opt.shrink <- shrinkage[inds[1]]
  opt.n.trees <- n.trees[inds[2]]
  opt.concord <- concord.mat[inds[1], inds[2]]
  return(list(shrinkage=opt.shrink, n.trees=opt.n.trees, concord=opt.concord))
}














