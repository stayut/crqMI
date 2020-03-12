library(crqMI)
set.seed(1234)
datasim=sim(n=100)
data = datasim$obs.data

tau=c(0.3,0.5)
out=crqMI(Y= data$Y, Z= data$Z, X=data$X, d= datasim$d,
          delta.v=data$delta.v, delta.t=data$delta.t,  Imp=NULL,
          M = 10, tau=tau, ctype= datasim$ctype,seed=T,
          pi.fam="logistic",se=TRUE,Boot=100,shrinkage=TRUE)


MI=out$estimate$MI  ## MI estimates
CC=out$estimate$CC  ## CC estimates
SMI=out$estimate$SMI ## shrinkage estimates
## estimation error and SE
true.par =datasim$true.par[,paste(tau)]
print(round(cbind(MI-true.par,SMI-true.par,CC-true.par),3))
print(round(cbind(out$se$MI,out$se$SMI,out$se$CC),3))

