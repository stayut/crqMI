sim <-
function(n,true.alpha = c(0.5,0.5),
           true.beta1=c(-0.5,0.5),
           true.eta = c(0.5,0.5),
           nc.t=c(0.2,0.2),ctype=c("left","right"), Cmax=5,
           r=0.2,
           sigma.t=0.5){

    if(!any(ls(all=TRUE)=="n")){stop("please give the sample size n")}
    if(!any(ls(all=TRUE)=="nc.t")){
      stop("please give the theoretical value for censoring rate (nc.t) in terms of univariate censored covariate")}
    if(!any(ls(all=TRUE)=="Cmax")){stop("please give the the potential max value (Cmax) for censoring in terms of event times")}

    if(length(nc.t)==1)nc.t=rep(nc.t,2)
    if(length(ctype)==1)ctype=rep(ctype,2)

    ## full observed covariates
    x1=runif(n, min = -1, max = 1)
    x2=rbinom(n, size = 1, prob = 0.5)
    X = cbind(x1,x2)
    colnames(X)=paste0("X",1:ncol(X))

    V = t(sapply(1:n,function(idx){(X%*%true.eta)[idx]+
        mvtnorm::rmvnorm(1,mean = rep(0,2), sigma = matrix(c(1,0.4,0.4,0.6),2,2))}))

    if( length(which(ctype =="left"))==2){

      d=sapply(1:2,function(idx)seq(min(V[,idx]),max(V[,idx]),length.out = (n*100)))
      d.nc= sapply(1:2,function(idx)
        sapply(1:length(d[,idx]),FUN = function(dd){sum(V[,idx]<=d[dd,idx])/n}))
      d.pos=sapply(1:2,function(idx)min(which(d.nc[,idx]>nc.t[idx])))
      d=sapply(1:2,function(idx)d[d.pos[idx],idx])
      nc.r=sapply(1:2,function(idx)d.nc[d.pos[idx],idx]) ## real value for censoring rate
      Z = sapply(1:2,function(idx)pmax(V[,idx],d[idx]))
      delta.v = (V-matrix(d,nrow = nrow(V),ncol = 2,byrow = T)>= 0)
    } else if(length(which(ctype =="right"))==2){
      d=sapply(1:2,function(idx)seq(min(V[,idx]),max(V[,idx]),length.out = (n*100)))
      d.nc= sapply(1:2,function(idx)
        sapply(1:length(d[,idx]),FUN = function(dd){sum(V[,idx]>=d[dd,idx])/n}))
      d.pos=sapply(1:2,function(idx)min(which(d.nc[,idx]<nc.t[idx])))
      d=sapply(1:2,function(idx)d[d.pos[idx],idx])
      nc.r=sapply(1:2,function(idx)d.nc[d.pos[idx],idx]) ## real value for censoring rate
      Z = sapply(1:2,function(idx)pmin(V[,idx],d[idx]))## V is observed by Z
      delta.v = (V-matrix(d,nrow = nrow(V),ncol = 2,byrow = T)<= 0)
    } else if(length(which(ctype =="left"))==1&length(which(ctype =="right"))==1){
      left.idx = which(ctype =="left")
      d.vec= seq(min(V[,left.idx]),max(V[,left.idx]),length.out = (n*100))
      d.nc= sapply(1:length(d.vec),FUN = function(dd){sum(V[,left.idx]<=d.vec[dd])/n})
      d.pos=min(which(d.nc>nc.t))
      d = rep(NA,2)
      d[left.idx]=d.vec[d.pos]
      nc.r= rep(NA,2)
      nc.r[left.idx]=d.nc[d.pos]
      Z= matrix(NA,nrow=nrow(V),ncol=ncol(V))## V is observed by Z
      Z[,left.idx]=pmax(V[,left.idx],d[left.idx])
      delta.v = matrix(NA,nrow=nrow(V),ncol=ncol(V))
      delta.v[,left.idx] = as.numeric(V[,left.idx]>= d[left.idx])

      right.idx = which(ctype =="right")
      d.vec= seq(min(V[,right.idx]),max(V[,right.idx]),length.out = (n*100))
      d.nc= sapply(1:length(d.vec),FUN = function(dd){sum(V[,right.idx]>=d.vec[dd])/n})
      d.pos=min(which(d.nc<nc.t))
      d[right.idx]=d.vec[d.pos]
      nc.r[right.idx]=d.nc[d.pos]
      Z[,right.idx]=pmin(V[,right.idx],d[right.idx])
      delta.v[,right.idx] = as.numeric(V[,right.idx]<= d[right.idx])

    }

    colnames(Z)=c("Z1","Z2")
    if(is.null(r)){stop("please input r for model 2")}
    epsilon = rnorm(n, sd = sigma.t)*(1+r*x1)
    T = drop(exp(V%*%true.alpha + X %*% true.beta1+epsilon))
    summary(T)

    ## censoring distribution is covariate-dependent
    C = rep(NA, n)
    C[x1<=0] = runif(sum(x1<=0), min=0, max=Cmax)
    C[x1>0] = runif(sum(x1>0), min=1, max=Cmax)

    # C = runif(n, min=0, max=Cmax)

    ## T is observed by Y
    Y = pmin(T,C)
    # hist(Y)
    delta.t = as.numeric(T <= C)


    names(true.alpha)= paste0("alpha_",1:length(true.alpha))
    names(true.beta1)= paste0("beta1_",1:length(true.beta1))
    names(true.eta)= paste0("eta_",1:length(true.eta))

    obs.data = list(Y=Y,Z=Z,X=X,delta.v=delta.v,delta.t=delta.t)

    true.data = list(T=T,V=V,X=X,C=C)

    true.par=c(true.alpha,true.beta1)
    true.par = matrix(true.par,nr=length(true.par),nc=100)
    true.par = rbind(qnorm((1:100)/100,mean = 0,sd=sigma.t),
                     true.par)

    true.par[4,] = true.par[4,]+r*qnorm((1:100)/100,mean = 0,sd=sigma.t)
    colnames(true.par)=(1:100)/100
    rownames(true.par)=c("Intercept","alpha1","alpha2","beta1","beta2")

    # censoring rate for covariate and outcome
    cr=c(nc.r,1-sum(delta.t)/n)
    # cat("censoring rate for covariates",nc.r,".\n")
    # cat("censoring rate for outcomes",1-sum(delta.t)/n,".\n")
    # cat("par",true.par,".\n")

    re=list(obs.data=obs.data,
            true.data=true.data,
            true.par=true.par,
            d=d,
            cr=cr,
            ctype=ctype)
    return(re)
  }
