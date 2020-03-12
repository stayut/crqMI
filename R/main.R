crqMI <-
  function(Y, Z=NULL, X=NULL, d, delta.v, delta.t, Imp=NULL,
           M = 10, tau=c(0.25,0.5,0.75), ctype="left",seed=T,
           pi.fam="logistic",se=FALSE,Boot=100,shrinkage=FALSE){
    require(survival)
    require(quantreg)
    # Y: outcomes
    # Z: covariates may be censored
    # V: potentially true covariates for Z
    # X: fully observed covariates
    # d: detection limits for censored covariates Z
    # delta.v: censoring indicators for Z: =1 if observed
    # delta.t: censoring indicators for Y: =1 if observed
    # Imp: imputation model; by default, Imp=cbind(log(Y), delta.t, X)
    # M: the times of multiple imputation
    # ctype: covariates' censoring type, left or right censored
    # pi.fam: logistic or probit regression is fitted for censoring probability for V


    out=crqMI.nose(Y=Y, Z=Z, X=X, d=d, delta.v=delta.v,
                   delta.t=delta.t,  Imp=Imp,
                   M = M, tau=tau, ctype=ctype,seed=seed,
                   pi.fam=pi.fam)
    fit=list()
    fit$estimate$MI=out[,which(colnames(out)%in%
                                 paste0("MI",".tau=",tau))]

    if(se){
      out.boot=crqMI.boot(Y=Y, Z=Z, X=X, d=d, delta.v=delta.v,
                          delta.t=delta.t,
                          Imp=Imp, M = M, tau=tau, ctype=ctype,
                          Boot=Boot,seed=seed,pi.fam=pi.fam)
      SE=out.boot$SE
      rownames(SE)=rownames(out)
      fit$se$MI=SE[,which(colnames(out)%in%
                            paste0("MI",".tau=",tau))]
    }
    if(se&shrinkage){
      vtemp=out.boot$vtemp
      vtemp.pos=which(colnames(SE)%in%
                        paste0(c("CC","MI"),".tau=",rep(tau,each=2)))
      vtemp.pos=c(sapply(1:length(vtemp.pos),function(pos){
        (nrow(SE)* (vtemp.pos[pos]-1)+1):(nrow(SE)*vtemp.pos[pos])}))

      cov.vtemp= cov(vtemp[,vtemp.pos])

      cov.dim=nrow(cov.vtemp)/2
      cov1=cov.vtemp[1:cov.dim,1:cov.dim]
      cov2=cov.vtemp[(cov.dim+1):(2*cov.dim),(cov.dim+1):(2*cov.dim)]
      cov12=cov.vtemp[1:cov.dim,(cov.dim+1):(2*cov.dim)]
      var.est.cc= matrix(diag(
        cov1+cov2-2*cov12),byrow = F,ncol = length(tau))

      var.est.cc=var.est.cc/(var.est.cc+(out[,which(colnames(out)%in%
                                                      paste0("MI",".tau=",tau))]-
                                           out[,which(colnames(out)%in%
                                                        paste0("CC",".tau=",tau))])^2)

      shrink= out[,which(colnames(out)%in%
                           paste0("CC",".tau=",tau))]+
        var.est.cc*(
          out[,which(colnames(out)%in%
                       paste0("MI",".tau=",tau))]-
            out[,which(colnames(out)%in%
                         paste0("CC",".tau=",tau))])

      SE.shrink= sqrt(diag(cbind(diag(1- c(var.est.cc)),diag(c(var.est.cc)))%*%
                             cov.vtemp%*%t(cbind(diag(1- c(var.est.cc)),diag(c(var.est.cc))))))

      shrink.se=matrix(SE.shrink, byrow=F, nrow=nrow(SE))
      colnames(shrink)=colnames(shrink.se)=paste0("SMI.tau=",tau)
      rownames(shrink.se)=rownames(shrink)
      fit$se$SMI=shrink.se
      fit$estimate$SMI=shrink
    }

    fit$estimate$CC=out[,which(colnames(out)%in%
                                 paste0("CC",".tau=",tau))]

    fit$estimate$filld=out[,which(colnames(out)%in%
                                    paste0("naive",".tau=",tau))]

    if(se){
      fit$se$CC=SE[,which(colnames(out)%in%
                            paste0("CC",".tau=",tau))]
      fit$se$filld=SE[,which(colnames(out)%in%
                               paste0("naive",".tau=",tau))]
    }

    return(fit)


  }

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

crqMI.nose <-
  function(Y, Z=NULL, X=NULL, d, delta.v, delta.t, Imp=NULL,
           M = 10, tau=c(0.25,0.5,0.75), ctype="left",seed=T,
           pi.fam="logistic"){
    require(survival)
    require(quantreg)
    # Y: outcomes
    # Z: covariates may be censored
    # V: potentially true covariates for Z
    # X: fully observed covariates
    # d: detection limits for censored covariates Z
    # delta.v: censoring indicators for Z: =1 if observed
    # delta.t: censoring indicators for Y: =1 if observed
    # Imp: imputation model; by default, Imp=cbind(log(Y), delta.t, X)
    # M: the times of multiple imputation
    # ctype: covariates' censoring type, left or right censored
    # pi.fam: logistic or probit regression is fitted for censoring probability for V


    Y=drop(Y)
    n=length(Y)
    if(is.null(X)){
      X=as.matrix(rep(1,n))
      dimX=0   # NO X, only one censored covariate V
    }else{
      X=as.matrix(X)
      dimX=1  # having X and V
    }
    Z=as.matrix(Z)
    d= drop(d)

    if(is.vector(d)){
      if(ncol(Z)!=length(d)& ncol(Z)!=length(ctype)) stop("please check the dimension of censored covariates")
      d=matrix(d,nrow = nrow(Z),ncol = ncol(Z),byrow = T)
    }
    if(is.matrix(d)){
      if(nrow(Z)!=nrow(d)&ncol(Z)!=ncol(d)& ncol(Z)!=length(ctype)) stop("please check the dimension of censored covariates")
    }

    if(!any(ls(all=TRUE)=="delta.v")){
      delta.v=c()
      for(k in 1:ncol(Z)){
        delta.v = cbind(delta.v,(Z[,k]!= d[,k]))
      }

    }
    delta.v=as.matrix(delta.v)
    sgn= rep(1,ncol(Z))

    if(any(ctype=="left")){
      Z[,which(ctype=="left")]=-Z[,which(ctype=="left")]
      d[,which(ctype=="left")]=-d[,which(ctype=="left")]
      sgn[which(ctype=="left")]=-1
    }

    if(is.null(Imp)&dimX==1){Imp=cbind(log(Y), delta.t, X)}
    if(is.null(Imp)&dimX==0){Imp=cbind(log(Y), delta.t)}
    Imp=as.matrix(Imp)

    #### step 1
    #### estimate the censoring probability
    pi0 = suppressWarnings(sapply(1:ncol(Z),function(k){pi.func(delta=as.numeric(delta.v[,k]),
                                                                A=cbind(1,Imp),fam=pi.fam)}))

    #### step 2
    #### estimate coefficients in the imputation model for V
    crqfit1 <- try(crq(Surv(Z[,1], delta.v[,1]) ~ Imp,method = "Portnoy"))


    if(inherits(crqfit1,"try-error"))stop("Errors in crq procedure")
    # ### some of the solutions are not stable (remove them, another option is to smooth them out)
    # tmp = crqfit1$sol[nrow(crqfit1$sol),]
    # idx.outlier = which((tmp-median(tmp))>5*mad(tmp))
    # if(length(idx.outlier)>0) {
    #   crqfit1$sol = crqfit1$sol[,-idx.outlier]
    #   # cat("removed", paste(length(idx.outlier)), "outliers", "\n")
    # }

    Qbar10=predict(crqfit1,stepfun=TRUE)
    Qbar1=rearrange(Qbar10) # rearranged


    #### V.im: imputed values for the censored covariate
    est = NULL
    # Vim = matrix(0, ncol=M, nrow=n)

    b=1
    if(seed){set.seed(528)}
    repeat{
      #cat("MI",b,".\n")
      V.im = Z
      V.im[delta.v==0] = NA


      # for those crossed cases
      if(length(crqfit1$Isplit)>0){
        for(i in 1:length(crqfit1$Isplit)){
          idx = crqfit1$Isplit[i]

          MI.times=1
          repeat{
            u = max(0.0001, runif(1, 1-pi0[idx,1], 1))
            if(MI.times<= 10){
              if(u<=max(crqfit1$sol[1,])+0.001){
                V.im[idx,1] = mycoef.crq.rearrange(crqfit1, Qbar=Qbar1, tau=u, pos=idx)
                break
              }else{
                MI.times=MI.times+1
              }
            }else{
              V.im[idx,1] = NA
              break
            }
          }

        }
      }


      # for those deleted censored
      if(sum(crqfit1$status==2)>0){
        for(idx in which(crqfit1$status==2))
        {
          u = runif(1, min(crqfit1$sol[1,]), 1)
          V.im[idx,1] = mycoef.crq.rearrange(crqfit1, Qbar=Qbar1, tau=u, pos=idx)

        }
      }
      V.im[which(delta.v[,1]==0 & V.im[,1]<d[,1]),1] = NA  # in case that some of imputed values are <d
      ind1 = which( is.na(V.im[,1])==FALSE)
      ind.list=list()
      ind.list[[1]]=ind1

      if(ncol(Z)>1){
        for(idx.colZ in 2:ncol(Z)){
          #### for the second column of V
          ind.last=1:n
          for(ind.len in 1:length(ind.list)){
            ind.last=ind.last[ind.list[[ind.len]]]
          }

          imp.matrix=cbind(Z[,idx.colZ], delta.v[,idx.colZ],
                           V.im[,1:(idx.colZ-1)],Imp)[ind.last,]
          crqfit2 =crq(Surv(imp.matrix[,1],imp.matrix[,2])~imp.matrix[,3:(idx.colZ+1)]+imp.matrix[,-(1:(idx.colZ+1))],
                       method = "Portnoy")

          if(inherits(crqfit2,"try-error"))stop("Errors in crq procedure")
          # ### some of the solutions are not stable (remove them, another option is to smooth them out)
          # tmp = crqfit2$sol[nrow(crqfit2$sol),]
          # idx.outlier = which((tmp-median(tmp))>5*mad(tmp))
          # if(length(idx.outlier)>0) {
          #   crqfit2$sol = crqfit2$sol[,-idx.outlier]
          #   # cat("removed", paste(length(idx.outlier)), "outliers", "\n")
          # }
          Qbar20=predict(crqfit2,stepfun=TRUE)
          Qbar2=rearrange(Qbar20) # rearranged

          if(length(crqfit2$Isplit)>0){
            for(i in 1:length(crqfit2$Isplit)){
              idx = crqfit2$Isplit[i]

              MI.times=1
              repeat{
                u = max(0.0001, runif(1, 1-pi0[ind.last,idx.colZ][idx], 1))

                if(MI.times<= 10){
                  if(u<=max(crqfit2$sol[1,])+0.001){
                    V.im[ind.last,idx.colZ][idx] = mycoef.crq.rearrange(crqfit2, Qbar=Qbar2, tau=u, pos=idx)
                    break
                  }else{
                    MI.times=MI.times+1
                  }
                }else{
                  V.im[ind.last,idx.colZ][idx] = NA

                  break
                }
              }


            }
          }


          # for those deleted censored
          if(sum(crqfit2$status==2)>0){
            for(idx in which(crqfit2$status==2))
            {
              u = runif(1, min(crqfit2$sol[1,]), 1)
              V.im[ind.last,idx.colZ][idx]=mycoef.crq.rearrange(crqfit2, Qbar=Qbar2, tau=u, pos=idx)


            }
          }
          V.im[which(delta.v[,idx.colZ]==0 & V.im[,idx.colZ]<=d[,idx.colZ]),idx.colZ] = NA  ### in case some of imputed values are <d
          ind2 = which(is.na(V.im[ind.last,idx.colZ])==FALSE)
          ind.list=c(ind.list,list(ind2))
        }
      }


      #### step 3: Find the solution of the b-th estimating equation
      data.last= cbind(log(Y), delta.t, V.im,X)
      for(ind.len in 1:length(ind.list)){
        data.last=data.last[ind.list[[ind.len]],]
      }
      coeff = c(coef(crq(Surv(data.last[,1], data.last[,2]) ~
                           data.last[,3:(ncol(Z)+2)] + data.last[,-(1:(ncol(Z)+2))],
                         method="PengHuang"),taus=tau))


      est = rbind(est, coeff)

      rm(V.im,crqfit2,ind1,ind2,ind.list)

      if(b>= M){break
      }else{b=b+1}
    }

    #### step 4: the multiple imputation estimate
    MI = apply(est, 2, mean)
    MI = matrix(MI, byrow=F, ncol=length(tau))
    colnames(MI) = paste0("MI.tau=",tau)


    #### Complete case analysis
    CC= try(coef(crq(Surv(log(Y)[apply(delta.v,1,sum)==ncol(Z)],
                          delta.t[apply(delta.v,1,sum)==ncol(Z)]) ~
                       Z[apply(delta.v,1,sum)==ncol(Z),] + X[apply(delta.v,1,sum)==ncol(Z),],
                     method="PengHuang"),taus=tau))


    naive= try(coef(crq(Surv(log(Y), delta.t) ~ Z + X,
                        method="PengHuang"),taus=tau))

    if(inherits(CC,"try-error")|inherits(naive,"try-error")){
      warning("Errors in crq procedure")
    }else{
      colnames(CC)=paste0("CC.tau=",tau)
      colnames(naive)=paste0("naive.tau=",tau)
      MI= cbind(CC,naive,MI)
    }


    MI[2:(ncol(Z)+1),]= MI[2:(ncol(Z)+1),]*sgn
    if(length(colnames(X))==0)colnames(X)=paste0("X",1:ncol(X))
    if(length(colnames(Z))==0)colnames(Z)=paste0("Z",1:ncol(Z))
    rownames(MI)=c("Intercept",colnames(Z),colnames(X))
    return(MI)
  }
mycoef.crq.rearrange <-
  function (object, Qbar, tau = 0.99, pos,...)
  {
    # Qbar is rearranged
    if (min(tau) < 0 || max(tau) > 1)
      stop("tau out of range [0,1]")
    if (length(object$tau) == 1) {
      Q <- object$fitted.values
      return(Q)
    }

    r <- object$sol[1, ]
    r = pmin(pmax(r, 0), 1)

    Qbar=Qbar[[pos]]

    if (is.unsorted(r))r <- r[-length(r)]

    if(tau >= max(r))
    {
      if((tau-max(r))<=0.001)
      {

        Q=Qbar(max(r))

      }
      else Q=NA
    } else
    {
      bin <- findInterval(tau, r)
      wgt <- (tau - r[bin])/(r[bin + 1] - r[bin])

      Q <- wgt * Qbar(r[bin]) + (1 - wgt) * Qbar(r[bin+1])
    }
    return(Q)
  }
pi.func <-
  function(delta, A,fam="logistic")
  {
    require(stats)
    # calculate pi(A) = P(delta=0|A) using logistic regression
    # A2 = cbind(1,A)
    A2=as.matrix(A)
    if(fam=="logistic"){
      #logistic regression
      theta <- glm(1-delta ~ A2-1, family=binomial("logit"))$coef
      pi = exp(A2%*%theta)/(1 +exp(A2%*%theta))
    }else if(fam=="probit"){
      #probit regression
      theta <- glm(1-delta ~ A2-1, family=binomial("probit"))$coef
      pi = pnorm(A2%*%theta)
    }
    return(pi)
  }
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

