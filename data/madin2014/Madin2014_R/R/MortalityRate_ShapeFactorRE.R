
# last modified by SRC on 05/02/14

# Fits random effects meta-regression using REML (default) or SJ (if zero mortalities
# are used with no variances)
################################################################################

rm(list=ls())

# setwd("c:/Users/jc129774/My Documents/Sean/Lizard Corals/MortalityRate vs ShapeFactor")

library(metafor) # for rma.uni fxn
library(plotrix) # for addtable2plot fxn

load("data/sean.Rdata")

large <- sean[sean$group=="largest",]
small <- sean[sean$group=="smallest",]

large.lis <- vector("list",5)
small.lis <- vector("list",5)

large.lis[[5]] <- large[large$proportion==0.5,]
large.lis[[4]] <- large[large$proportion==0.4,]
large.lis[[3]] <- large[large$proportion==0.3,]
large.lis[[2]] <- large[large$proportion=="0.2",]
large.lis[[1]] <- large[large$proportion=="0.1",]

small.lis[[5]] <- small[small$proportion==0.5,]
small.lis[[4]] <- small[small$proportion==0.4,]
small.lis[[3]] <- small[small$proportion==0.3,]
small.lis[[2]] <- small[small$proportion=="0.2",]
small.lis[[1]] <- small[small$proportion=="0.1",]

large <- large.lis
rm(large.lis)

small <- small.lis
rm(small.lis)

  # Caution - There are zero sampling variances for GR and GP.
  # Help on rma.uni:
  # "Outcomes with non-positive sampling variances are problematic.
  # If a sampling variance is equal to zero, then its weight will be 1/0
  # for fixed-effects models when using weighted estimation.
  # Switching to unweighted estimation is a possible solution."
  # Switching to method="SJ" works without error

  # Sean: the limit on sampling variance for an all-zero sample is 1/n
  # use this and see if it makes any differences.
  
modify.var <- T
  # TRUE: plug in 1/n as the variance for all-zero samples
  # FALSE: leave variance=0 as is


quartz(width=9,height=5)
par(mfrow=c(2,5),mar=c(4,4,6,1),oma=c(1,1,0,0))
# largest
for(i in 1:5) {
  # shape factors
  meancsf <- large[[i]]$csf_mean
  # mean mortality rate estimates
  mest <- large[[i]]$mort_est
  n <- large[[i]]$n

  # sampling variance for each mest : NOTE: this is because I gave Josh the
  # wrong formula for calculating sampling variance -- should have divided by 
  # n, not n-1, so multiply by (n-1)/n to fix.
  sampvar <- large[[i]]$mort_var*(n-1)/n
  #sampvar <- large[[i]]$mort_var
  if(modify.var) sampvar[mest==0&sampvar==0] <- 1/n[mest==0&sampvar==0]
  lci <- mest-sqrt(sampvar)    # SE, not CI
  uci <- mest+sqrt(sampvar)

  fit <- rma.uni(mest,sampvar,mods=~meancsf) #,method="SJ") # need to use method=SJ
              # if including points with zero variance
              # this function does a random-effects meta-regression of mortality 
              # rate (weighted by uncertainty) against CSF 
  print(fit)
              
  newx <- seq(min(meancsf),max(meancsf),len=10)
  preds <- predict(fit,newmods=newx)
  if(!modify.var) yli <- c(0,0.5)#range(0,preds$ci.ub,uci)
  if(modify.var) yli <- c(0,0.5)
  plot(meancsf,mest,ylim=yli,xlab="",ylab="")
  title(main=paste("Largest\nFraction =",large[[i]]$prop[1]),line=3.5)
  sapply(seq_along(meancsf),function(j)lines(rep(meancsf[j],2),c(lci[j],uci[j])))
  lines(newx,preds$pred)
  lines(newx,preds$ci.lb,lty=3)
  lines(newx,preds$ci.ub,lty=3)
  tabl <- formatC(as.matrix(data.frame(est=fit$b,se=fit$se,pval=fit$pval)),dig=3)
  tabl <- cbind(c("intrcpt","slope"),tabl)
  addtable2plot(min(par()$usr[1:2]),max(par()$usr[3:4]),table=tabl,cex=0.7,
    display.rownames=F,yjust=1.1)
}
# smallest
for(i in 1:5) {
  # shape factors
  meancsf <- small[[i]]$csf_mean
  # mortality rate estimates
  mest <- small[[i]]$mort_est
  n <- small[[i]]$n
  # sampling variance for each mest
  sampvar <- small[[i]]$mort_var*(n-1)/n
  #sampvar <- small[[i]]$mort_var
  if(modify.var) sampvar[mest==0&sampvar==0] <- 1/n[mest==0&sampvar==0]
  lci <- mest-sqrt(sampvar)  # SE, not CI
  uci <- mest+sqrt(sampvar)

  fit <- rma.uni(mest,sampvar,mods=~meancsf) #,method="SJ")
  print(fit)
  newx <- seq(min(meancsf),max(meancsf),len=10)
  preds <- predict(fit,newmods=newx)
  yli <- c(0,0.8)#range(0,preds$ci.ub,uci)
  plot(meancsf,mest,ylim=yli,xlab="",ylab="")
  title(main=paste("Smallest\nFraction =",small[[i]]$prop[1]),line=3.5)
  sapply(seq_along(meancsf),function(j)lines(rep(meancsf[j],2),c(lci[j],uci[j])))
  lines(newx,preds$pred)
  lines(newx,preds$ci.lb,lty=3)
  lines(newx,preds$ci.ub,lty=3)
  tabl <- formatC(as.matrix(data.frame(est=fit$b,se=fit$se,pval=fit$pval)),dig=3)
  tabl <- cbind(c("intrcpt","slope"),tabl)
  addtable2plot(min(par()$usr[1:2]),max(par()$usr[3:4]),table=tabl,cex=0.7,
    display.rownames=F,yjust=1.1)
}
title(xlab="Mean csf",outer=T,line=-1,cex.lab=2)
title(ylab="Mean mortality rate (+/- 95% CI)",outer=T,line=-1,cex.lab=2)

plotfile <- paste("MortalityRate_ShapeFactor",if(modify.var)"_AddedVar",sep="")
#savePlot(plotfile,type="pdf")
# saved on 04/02/14 and 05/02/14


#-------------------------------------------------------------------------------
# see what happens when I omit a couple of outliers (with large mean csf)
# from the largest sets

quartz(width=9,height=2.5)
par(mfrow=c(1,5),mar=c(4,4,6,1),oma=c(1,2.5,0,0))
# largest
for(i in 1:5) {
  # shape factors
  meancsf <- large[[i]]$csf_mean
  # mean mortality rate estimates
  mest <- large[[i]]$mort_est
  n <- large[[i]]$n
  # sampling variance for each mest
  sampvar <- large[[i]]$mort_var*(n-1)/n
  #sampvar <- large[[i]]$mort_var
  if(modify.var) sampvar[mest==0&sampvar==0] <- 1/n[mest==0&sampvar==0]
  lci <- mest-sqrt(sampvar)   # SRC: modified to SE, not conf ints
  uci <- mest+sqrt(sampvar)

  # removing the biggest two meancsf values
  omitindx <- order(meancsf,decreasing=T)[1:2]
  meancsf <- meancsf[-omitindx]
  mest <- mest[-omitindx]
  n <- n[-omitindx]
  sampvar <- sampvar[-omitindx]
  lci <- lci[-omitindx]
  uci <- uci[-omitindx]
  
  fit <- rma.uni(mest,sampvar,mods=~meancsf) #,method="SJ")
  print(fit)
  newx <- seq(min(meancsf),max(meancsf),len=10)
  preds <- predict(fit,newmods=newx)
  if(!modify.var) yli <- c(0,0.5)#range(0,preds$ci.ub,uci)
  if(modify.var) yli <- c(0,0.5)
  plot(meancsf,mest,ylim=yli,xlab="",ylab="")
  title(main=paste("Largest\nFraction =",large[[i]]$prop[1]),line=3.5)
  sapply(seq_along(meancsf),function(j)lines(rep(meancsf[j],2),c(lci[j],uci[j])))
  lines(newx,preds$pred)
  lines(newx,preds$ci.lb,lty=3)
  lines(newx,preds$ci.ub,lty=3)
  tabl <- formatC(as.matrix(data.frame(est=fit$b,se=fit$se,pval=fit$pval)),dig=3)
  tabl <- cbind(c("intrcpt","slope"),tabl)
  addtable2plot(min(par()$usr[1:2]),max(par()$usr[3:4]),table=tabl,cex=0.7,
    display.rownames=F,yjust=1.1)
}
title(xlab="Mean csf",outer=T,line=-1,cex.lab=2)
title(ylab="Mean mortality rate\n(+/- 95% CI)",outer=T,line=-1,cex.lab=2)

plotfile <- paste("MortalityRate_ShapeFactor_OutliersRemoved",
            if(modify.var)"_AddedVar",sep="")
#savePlot(plotfile,type="pdf")
# saved on 04/02/14 and 05/02/14
