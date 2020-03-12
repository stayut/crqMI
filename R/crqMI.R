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
