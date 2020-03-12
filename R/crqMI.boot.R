crqMI.boot <-
function(Y, Z=NULL, X=NULL, d, delta.v, delta.t,  Imp=NULL,
           M = 10, tau=c(0.25,0.5,0.75), ctype="left",
           Boot=100,seed=T,pi.fam="logistic"){
    require(survival)
    require(quantreg)
    require(foreach)
    require(doSNOW)
    # method="linear": linear imputation model (crqMI)
    # method="nonlinear": nonlinear imputation model (scrqMI)


    Z=as.matrix(Z)
    X=as.matrix(X)
    delta.v=as.matrix(delta.v)
    vtemp=NULL
    if(seed){set.seed(528)}
    n=length(Y)
    resample.index.list=list()
    for (boot in 1:Boot){
      resample.index=sample(1:n,size = n,replace = T)
      resample.index.list[[boot]]=resample.index
    }

    cpu1<-makeCluster(4)
    #change the number of CPU cores
    registerDoSNOW(cpu1)
    if(seed){set.seed(528)}

    out.resample=foreach(boot=1:Boot,
                         .export=c("crq","Surv","crqMI.nose","pi.func",
                                   "rearrange","mycoef.crq.rearrange")) %dopar%{
                                     try(suppressWarnings(crqMI.nose(Y=Y[resample.index.list[[boot]]],
                                                                     Z=Z[resample.index.list[[boot]],],
                                                                     X=X[resample.index.list[[boot]],],
                                                                     d=d,
                                                                     delta.v=delta.v[resample.index.list[[boot]],],
                                                                     delta.t=delta.t[resample.index.list[[boot]]],
                                                                     Imp=Imp[resample.index.list[[boot]],],
                                                                     M = M, tau=tau, ctype=ctype,pi.fam=pi.fam)))}


    stopCluster(cpu1)

    error.remove= sapply(1:length(out.resample),function(x)!inherits(out.resample[[x]],"try-error"))
    out.resample=out.resample[error.remove]

    vtemp= t(sapply(1:length(out.resample), function(x){as.vector(out.resample[[x]])}))
    #vtemp=na.omit(vtemp)
    na.vtemp.pos=c(1:(max(which((sum(delta.t)/n)>tau))*nrow(out.resample[[1]])),
                   ((ncol(out.resample[[1]])/3*2)*
                      nrow(out.resample[[1]]))+1:(max(which((sum(delta.t)/n)>tau))*nrow(out.resample[[1]])))
    na.vtemp=unique(which(is.na(vtemp[,na.vtemp.pos]),arr.ind = T)[,1])

    if(length(na.vtemp)!=0){vtemp=vtemp[-na.vtemp,]}
    temp.SE= apply(vtemp,2,sd)
    temp.SE= matrix(temp.SE, byrow=F, nrow=nrow(out.resample[[1]]))
    colnames(temp.SE) = colnames(out.resample[[1]])
    return(list(SE=temp.SE,vtemp=vtemp))
  }
