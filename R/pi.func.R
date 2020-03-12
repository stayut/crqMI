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
