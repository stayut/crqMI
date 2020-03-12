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
