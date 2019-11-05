
## For each case 1-8, one table
## out-of-sample concordance, log(lambda) MSE, avg model size, avg number of correct covariates, avg number of spurious covariates, and compute time.

# cd /data5/bsi/clip/s201827.clip_hpc/BSPR/main_code/simulations_12-18-16

# assuming you are in the "jasa_code" directory,
# otherwise you may have to alter some paths.
source("weibull_response_sparse_DPM.R")
registerDoMC()
options(cores=8)



##############################################################
## Read in all of the result files produced by sim_cases.R  ##
## and the necessary summaries in an object called sumarray ##
##############################################################

sumarray <- list()
inform.miss.list <- list()


for(case in 1:4){

vpcut <- NULL

pattern <- paste("CaseSeed_", case, sep="")
fnames <- dir("./Results", pattern=pattern, full.names=TRUE)

load(fnames[1])
inform.miss.list[[case]] <- inform.col


nmeth <- 11
nstats <- 11
nf <- length(fnames)
sumarray.case <- foreach(f=1:nf)%dopar%{
print(f)
  load(fnames[f])

  sumarray.f <- matrix(0,nmeth,nstats)
  ind.inform <- (beta[-1] !=0)
  if(case==6 || case==8)
    ind.inform <- (beta[-1] !=0)[-c(2,17)]
  cut <- ifelse(is.null(vpcut), rho, vpcut)

  sumarray.f[1,1] <- survConcordance(Surv(ycp, event=outcomecp)~pred.MVN.MAR)$concord
  sumarray.f[2,1] <- survConcordance(Surv(ycp, event=outcomecp)~pred.MVN.MNAR)$concord
  sumarray.f[3,1] <- survConcordance(Surv(ycp, event=outcomecp)~pred.sDPM.MAR)$concord
  sumarray.f[4,1] <- survConcordance(Surv(ycp, event=outcomecp)~pred.sDPM.MNAR)$concord
  sumarray.f[5,1] <- survConcordance(Surv(ycp, event=outcomecp)~lamhat.gbm1)$concord
  sumarray.f[6,1] <- survConcordance(Surv(ycp, event=outcomecp)~lamhat.gbm)$concord
  sumarray.f[7,1] <- survConcordance(Surv(ycp, event=outcomecp)~lamhat.RFEN1)$concord
  sumarray.f[8,1] <- survConcordance(Surv(ycp, event=outcomecp)~lamhat.RFEN2)$concord
  sumarray.f[9,1] <- survConcordance(Surv(ycp, event=outcomecp)~lamhat.fac.sw)$concord
  sumarray.f[10,1] <- survConcordance(Surv(ycp, event=outcomecp)~lamhat.cont.sw)$concord
  sumarray.f[11,1] <- survConcordance(Surv(ycp, event=outcomecp)~riskp)$concord


  sumarray.f[1,2] <- mean(((log(riskp)-mean(log(riskp)))-(pred.MVN.MAR-mean(pred.MVN.MAR)))^2)
  sumarray.f[2,2] <- mean(((log(riskp)-mean(log(riskp)))-(pred.MVN.MNAR-mean(pred.MVN.MNAR)))^2)
  sumarray.f[3,2] <- mean(((log(riskp)-mean(log(riskp)))-(pred.sDPM.MAR-mean(pred.sDPM.MAR)))^2)
  sumarray.f[4,2] <- mean(((log(riskp)-mean(log(riskp)))-(pred.sDPM.MNAR-mean(pred.sDPM.MNAR)))^2)
  sumarray.f[5,2] <- mean(((log(riskp)-mean(log(riskp)))-(lamhat.gbm-mean(lamhat.gbm1)))^2)
  sumarray.f[6,2] <- mean(((log(riskp)-mean(log(riskp)))-(lamhat.gbm-mean(lamhat.gbm)))^2)
  sumarray.f[7,2] <- mean(((log(riskp)-mean(log(riskp)))-(lamhat.RFEN1-mean(lamhat.RFEN1)))^2)
  sumarray.f[8,2] <- mean(((log(riskp)-mean(log(riskp)))-(lamhat.RFEN2-mean(lamhat.RFEN2)))^2)
  sumarray.f[9,2] <- mean(((log(riskp)-mean(log(riskp)))-(lamhat.fac.sw-mean(lamhat.fac.sw)))^2)
  sumarray.f[10,2] <- mean(((log(riskp)-mean(log(riskp)))-(lamhat.cont.sw-mean(lamhat.cont.sw)))^2)


  sumarray.f[1,3] <- cor(log(riskp),pred.MVN.MAR)^2
  sumarray.f[2,3] <- cor(log(riskp),pred.MVN.MNAR)^2
  sumarray.f[3,3] <- cor(log(riskp),pred.sDPM.MAR)^2
  sumarray.f[4,3] <- cor(log(riskp),pred.sDPM.MNAR)^2
  sumarray.f[5,3] <- cor(log(riskp),lamhat.gbm1)^2
  sumarray.f[6,3] <- cor(log(riskp),lamhat.gbm)^2
  sumarray.f[7,3] <- cor(log(riskp),lamhat.RFEN1)^2
  sumarray.f[8,3] <- cor(log(riskp),lamhat.RFEN2)^2
  sumarray.f[9,3] <- cor(log(riskp),lamhat.fac.sw)^2
  sumarray.f[10,3] <- cor(log(riskp),lamhat.cont.sw)^2
  sumarray.f[11,3] <- 1


  gamma <- (apply(ans.MVN.MAR$beta[,2:(p+1)]!=0,2,mean) > cut)*1
  sumarray.f[1,4] <- sum(gamma)
  sumarray.f[1,5] <- (sum(gamma[ind.inform]) + sum(1-gamma[!ind.inform]))/p
  
  gamma <- (apply(ans.MVN.MNAR$beta[,2:(p+1)]!=0,2,mean) > cut)*1
  sumarray.f[2,4] <- sum(gamma)
  sumarray.f[2,5] <- (sum(gamma[ind.inform]) + sum(1-gamma[!ind.inform]))/p
  
  gamma <- (apply(ans.sDPM.MAR$beta[,2:(p+1)]!=0,2,mean) > cut)*1
  sumarray.f[3,4] <- sum(gamma)
  sumarray.f[3,5] <- (sum(gamma[ind.inform]) + sum(1-gamma[!ind.inform]))/p
  
  gamma <- (apply(ans.sDPM.MNAR$beta[,2:(p+1)]!=0,2,mean) > cut)*1
  sumarray.f[4,4] <- sum(gamma)
  sumarray.f[4,5] <- (sum(gamma[ind.inform]) + sum(1-gamma[!ind.inform]))/p

  ans.gbm$var.names <- ans.gbm1$var.names <- paste("x",1:p,sep="")
  blah <- summary(ans.gbm1, plotit=FALSE)
  xord <- match(paste("x",1:p,sep=""), blah[,1])
  var.inf <- blah[xord,2]
  gamma <- (var.inf > .25)*1
  sumarray.f[5,4] <- sum(gamma)
  sumarray.f[5,5] <- (sum(gamma[ind.inform]) + sum(1-gamma[!ind.inform]))/p

  blah <- summary(ans.gbm, plotit=FALSE)
  xord <- match(paste("x",1:p,sep=""), blah[,1])
  var.inf <- blah[xord,2]
  gamma <- (var.inf > .25)*1
  sumarray.f[6,4] <- sum(gamma)
  sumarray.f[6,5] <- (sum(gamma[ind.inform]) + sum(1-gamma[!ind.inform]))/p

  ind.lam1 <- which(ans.RFEN1$fit$glmnet.fit$lambda == ans.RFEN1$fit$lambda.1se)
  gamma <- (ans.RFEN1$fit$glmnet.fit$beta[,ind.lam1]!=0)*1
  sumarray.f[7,4] <- sum(gamma)
  sumarray.f[7,5] <- (sum(gamma[ind.inform]) + sum(1-gamma[!ind.inform]))/p

  ind.lam2 <- which(ans.RFEN2$fit$glmnet.fit$lambda == ans.RFEN2$fit$lambda.1se)
  gamma <- (ans.RFEN2$fit$glmnet.fit$beta[,ind.lam1]!=0)*1
  sumarray.f[8,4] <- sum(gamma)
  sumarray.f[8,5] <- (sum(gamma[ind.inform]) + sum(1-gamma[!ind.inform]))/p

  if(case<5)
    x.inc <- match(attr(foo$terms,"term.labels"), paste("V",1:p,sep=""))
  else
    x.inc <- match(attr(foo0$terms,"term.labels"), paste("X.",1:p,sep=""))
  ind.inc <- rep(FALSE, p)
  ind.inc[x.inc] <- TRUE
  gamma <- ind.inc*1
  sumarray.f[9,4] <- sum(gamma)
  sumarray.f[9,5] <- (sum(gamma[ind.inform]) + sum(1-gamma[!ind.inform]))/p

  if(case<5)
    x.inc <- match(attr(foo0$terms,"term.labels"), paste("V",1:p,sep=""))
  else
    x.inc <- match(attr(foo0$terms,"term.labels"), paste("X.",1:p,sep=""))
  ind.inc <- rep(FALSE, p)
  ind.inc[x.inc] <- TRUE
  gamma <- ind.inc*1
  sumarray.f[10,4] <- sum(gamma)
  sumarray.f[10,5] <- (sum(gamma[ind.inform]) + sum(1-gamma[!ind.inform]))/p

  sumarray.f[11,4] <- sum(ind.inform)
  sumarray.f[11,5] <- 1

  if(case==3 || case==4){
    sumarray.f[1,6:11] <- apply(ans.MVN.MAR$beta[,inform.col+1],2,mean)
    sumarray.f[2,6:11] <- apply(ans.MVN.MNAR$beta[,inform.col+1],2,mean)
    sumarray.f[3,6:11] <- apply(ans.sDPM.MAR$beta[,inform.col+1],2,mean)
    sumarray.f[4,6:11] <- apply(ans.sDPM.MNAR$beta[,inform.col+1],2,mean)
    sumarray.f[7,6:11] <- ans.RFEN1$fit$glmnet.fit$beta[inform.col,ind.lam1]
    sumarray.f[8,6:11] <- ans.RFEN2$fit$glmnet.fit$beta[inform.col,ind.lam2]

    ind.fuh <- match( paste("V",inform.col,sep=""), attr(foo0$terms,"term.labels"))
    fuh0 <- foo0$coef[ind.fuh]
    fuh0[is.na(fuh0)] <- 0
    sumarray.f[10,6:11] <- fuh0
    
    sumarray.f[11,6:11] <- beta[inform.col+1]
  }
  sumarray.f
}

sumarray[[case]] <- array(0, c(nf,nmeth,nstats))
for(f in 1:nf)
  sumarray[[case]][f,,] <- sumarray.case[[f]]


apply(sumarray[[case]],c(2,3),mean)
}




## Compare methods via wilcox.test ##

alpha <- .01
mean.array <- sd.array <- se.array <- diff.best <- array(0,c(8,nmeth-1,nstats))
ind.inform <- (beta[-1] !=0)
for(ss in 1:8){
  nreal <- length(sumarray[[ss]][,1,1])
  for(j in 1:nstats){
    for(mm in 1:(nmeth-1)){
      foo <- sumarray[[ss]][,mm,j]
      mean.array[ss,mm,j] <- mean(foo)
      sd.array[ss,mm,j] <- sd(foo)
      se.array[ss,mm,j] <- sd.array[ss,mm,j]/length(foo)
    }
    target <- c(1,0,1,sum(ind.inform),1,beta[inform.miss.list[[3]]+1])
#    ind.best <- order((mean.array[ss,,j]-target[j])^2)[1]
    ind.best <- order((mean.array[ss,,j]-target[j])^2+sd.array[ss,,j]^2)[1]
    for(mm in (1:(nmeth-1))[-ind.best]){
      diff.best[ss,mm,j] <- (wilcox.test((sumarray[[ss]][,mm,j] - target[j])^2 +rnorm(nreal,0,1E-12) - (sumarray[[ss]][,ind.best,j] - target[j])^2, alt="greater", mu=0)$p.val < alpha/2)*1
    }
  }
}



meth.names <- c("MVN-MAR", "MVN-MNAR", "sDPM-MAR", "sDPM-MNAR", "GBM1", "GBM", "RFEN-MAR", "RFEN-MNAR", "Fac-SW", "CPH-SW", "TRUE")
#meth.ord <- c(1,2,3,4,5,6,7,8,9,10,11)
meth.ord <- c(1,2,3,4,6,7,10,11)
sum.names <- c("Concordance", "Risk Score R-sq", "Model Size", "PVC")
results <- list()
for(ss in 1:8){
  results[[ss]] <- apply(signif(apply(sumarray[[ss]],c(2,3),mean)[,c(1,3:5)], 3),2,format,digits=c(3,2,2,2))
  results[[ss]] <- cbind(results[[ss]], rbind(diff.best[ss,,c(1,3:5)],rep(NA,4)))
  rownames(results[[ss]]) <- meth.names
  colnames(results[[ss]]) <- c(sum.names,rep("SigDiff",4))
  results[[ss]] <- results[[ss]][meth.ord,c(1,5,2,6,3,7,4,8)]
}



################################################################
## Boxplots for Rsq to True log-risk (Figure 3 in manuscript) ##
################################################################

#par(mfrow=c(2,4))
height <- 2.7
width <- rep(.82*height,8)
width[c(1,4,5,8)] <- .82*height*1.423
for(ss in 1:4){
pdf(paste0("Rsq_plot_",ss,".pdf"), height=height,width=width[ss])
if(ss==1 || ss==5)
  par(mar=c(3.3,6,.5,.25), mgp=c(2.4,.8,0))
if(ss==4 || ss==8)
  par(mar=c(3.3,.25,.5,6), mgp=c(2.4,.8,0))
if(any(c(2,3,6,7)==ss))
  par(mar=c(3.3,.25,.5,.25), mgp=c(2.4,.8,0))
col1 <- rgb(0,.8,0)
col2 <- rgb(0,.2,1)
blah <- sumarray[[ss]][,meth.ord[(nrow(results[[ss]])-1):1],3]
foo <- apply(blah,2,median)
voo <- max(foo)
cols <- rep(col1,nrow(results[[ss]]))
cols[diff.best[ss,meth.ord[(nrow(results[[ss]])-1):1],2]==0] <- col2
ax.names <- row.names(results[[ss]])[(nrow(results[[ss]])-1):1]
boxplot(blah,horizontal=TRUE, yaxt='n', las=1, border=cols, xlab="", lwd=1.5, cex.axis=.9, boxwex=.7)
if(ss==1 || ss==5)
axis(2, labels=ax.names, at=1:(nrow(results[[ss]])-1), cex.axis=.9, las=1)
if(ss==4 || ss==8)
axis(4, labels=ax.names, at=1:(nrow(results[[ss]])-1), cex.axis=.9, las=1)
#title(main=paste0("Case ",ss), line=.7)
title(xlab="Risk Score R-sq", line=1.75)
abline(v=voo,col=col2,lty=1)
dev.off()
#readline("enter")
}






################################################
## Print Latex code for Table 1 in manuscript ##
################################################

case.ord <- c(1,3,2,4,5,7,6,8)
for(ss in (1:2)*2-1){
  if(ss==1)
    cat("& & \\multicolumn{4}{c}{Case 1: $\\bz\\sim$MVN and MAR} &&  \\multicolumn{4}{c}{Case 2: $\\bz\\sim$MVN and MNAR} \\\\[.03in]\n")
  if(ss==3)
    cat("& & \\multicolumn{4}{c}{Case 3: $\\bz\\sim$sDPM and MAR} &&  \\multicolumn{4}{c}{Case 4: $\\bz\\sim$sDPM and MNAR} \\\\[.03in]\n")
  if(ss==5)
    cat("& & \\multicolumn{4}{c}{Case 5: $\\bz\\sim$BPR and MAR} &&  \\multicolumn{4}{c}{Case 6: $\\bz\\sim$BPR and MNAR} \\\\[.03in]\n")
  if(ss==7)
    cat("& & \\multicolumn{4}{c}{$\\!\\!\\!$Case 7: $\\!\\bz\\sim$BPR, MAR, $x_2,x_{27}$ removed$\\!\\!\\!$} &&  \\multicolumn{4}{c}{$\\!\\!\\!$Case 8: $\\!\\bz\\sim$BPR, MNAR, $x_2,x_{27}$ removed$\\!\\!\\!$} \\\\[.03in]\n")

  cat("\\cline{1-1}\\cline{3-6}\\cline{8-11}\n")
  for(mm in 1:nrow(results[[case.ord[ss]]])){
    if(mm==nrow(results[[case.ord[ss]]]))
      cat("\\cline{1-1}\\cline{3-6}\\cline{8-11}\n")
    cat("\\!\\!",row.names(results[[case.ord[ss]]])[mm],"\\!\\!\\!&")
    for(jj in 1:4){
      if(mm<nrow(results[[case.ord[ss]]]) && results[[case.ord[ss]]][mm,2*jj]==0)
        cat(" & {\\bf ", results[[case.ord[ss]]][mm,2*jj-1], "}", sep="")
      else
        cat(" &", results[[case.ord[ss]]][mm,2*jj-1])
    }
    cat(" &")
    for(jj in 1:4){
      if(mm<nrow(results[[case.ord[ss+1]]]) && results[[case.ord[ss+1]]][mm,2*jj]==0)
        cat(" & {\\bf ", results[[case.ord[ss+1]]][mm,2*jj-1], "}", sep="")
      else
        cat(" &", results[[case.ord[ss+1]]][mm,2*jj-1])
    }
    cat(" \\\\\n")
  }
  cat("\\cline{1-1}\\cline{3-6}\\cline{8-11}\\\\[-.07in]") 
  cat("\n")
}      






################################################################
## boxplots of coefficient estimates (Figure 4 in manuscript) ##
################################################################

pdf("inform_beta_boxplots.pdf", height=5, width=8)

par(mfrow=c(2,3), mar=c(4.5, 6.0, 0.5, 0.5))
ss <- 4
for(j in 1:6){

beta.num <- c(1,2,3,26,27,28)
meth.names <- c("MVN-MAR", "MVN-MNAR", "sDPM-MAR", "sDPM-MNAR", "GBM1", "GBM", "RFEN-MAR", "RFEN-MNAR", "Fac-SW", "CPH-SW", "TRUE")
#meth.ord <- c(1,2,3,4,5,6,7,8,9,10,11)
meth.ord <- c(1,2,3,4,7,10)
beta.j <- list()
for(ll in 1:length(meth.ord)){
  mm <- meth.ord[ll]
  beta.j[[ll]] <- sumarray[[ss]][,mm,6:11][,j]
}
beta.true <- sumarray[[ss]][1,11,6:11][j]
col1 <- rgb(0,.8,0)
col2 <- rgb(0,.2,1)
cols <- rep(col1,length(meth.ord))
cols[diff.best[ss,meth.ord,j+5]==0] <- col2
#cols[4] <- col2

boxplot(beta.j[length(meth.ord):1], main='', names=meth.names[meth.ord][length(meth.ord):1], horizontal=TRUE, las=1, border=cols[length(meth.ord):1], lwd=1.5, cex.axis=.9)
title(xlab=paste0("Beta_",beta.num[j]), lin=1.8)
abline(v=beta.true, col=4)
}

dev.off()



