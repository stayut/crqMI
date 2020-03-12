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
