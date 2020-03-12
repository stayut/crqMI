\name{crqMI}
\alias{crqMI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Functions to fit censored quantile regression models with covariates subject to detection limits
}
\description{
The method is developed using some techniques in Peng and Huang (2008) and Wang and Feng (2012).
}
\usage{
crqMI(Y, Z, X, d, delta.v, delta.t, Imp = NULL, M = 10,
tau=c(0.25,0.5,0.75), ctype="left", seed=TRUE,
pi.fam="logistic",se=FALSE, Boot=100, shrinkage=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{ survival outcomes.
%%     ~~Describe \code{Y} here~~
}
  \item{Z}{ observations of censored covariates.
%%     ~~Describe \code{Z} here~~
}
  \item{X}{ fully observed covariates.
%%     ~~Describe \code{X} here~~
}
  \item{d}{ detection limits for censored covariates Z.
%%     ~~Describe \code{d} here~~
}
  \item{delta.v}{ censoring indicators for Z.  1 indicates observed otherwise censored.
%%     ~~Describe \code{delta.v} here~~
}
  \item{delta.t}{ censoring indicators for Y. 1 indicates observed otherwise censored.
%%     ~~Describe \code{delta.t} here~~
}
  \item{Imp}{ imputation matrix. Default is \code{Imp=NULL}, regarding (X,log(Y),delta.c) as the impuation model.
%%     ~~Describe \code{Imp} here~~
}
  \item{M}{ the number of imputation times
%%     ~~Describe \code{M} here~~
}
  \item{tau}{the quantile(s) at which the model is to be estimated.
%%     ~~Describe \code{tau} here~~
}
  \item{ctype}{ censoring type(s) of censored covariate(s), left or right censored. It can be a character or a vector of characters. By default \code{ctype="left"}, regarding all censored covariates as left-censored.
%%     ~~Describe \code{ctype} here~~
}
  \item{seed}{ If \code{seed=TRUE}, see seed in the computation.
%%     ~~Describe \code{seed} here~~
}
  \item{pi.fam}{Default is \code{pi.fam="logistic"}, using logistic regression to estimate the conditional censoring probability Pr(delta.v=0|imp). Another choice is \code{pi.fam="probit"}
%%     ~~Describe \code{pi.fam} here~~
}
\item{se}{ If \code{se=TRUE}, returning the estimates of standard error (SE) based on the Bootstrap method. Default is \code{se=FALSE}
}
 \item{Boot}{ the number of bootstrap resamples to estimate SE. Default is \code{Boot=100}
}
\item{shrinkage}{If \code{shrinkage=TRUE}, returning the shrinkage estimate. Default is \code{shrinkage=FALSE}
}

}
\details{
%The proposed method(MI and SMI)
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
 \item{estimate }{a list of estimates obtained from complete case analysis (CC), naive method ignoring censoring of covariates (filld), the proposed method (MI) and its shrinkage version (SMI) }
 \item{se}{a list of estimated stardard errors obtained from complete case analysis (CC), naive method ignoring censoring of covariates (filld), the proposed method (MI) and its shrinkage version (SMI)}

}
\references{
Peng, L. and Huang, Y. (2008). Survival analysis with quantile regression models, J. Am. Stat. Assoc., 103, 637-649.

Wang, H.J. and Feng, X. (2012) Multiple imputation for m-regression with censored covariates, J. Am. Stat. Assoc., 104, 1117-1128.
}
\author{
Yu Tonghui, Xiang Liming and Wang Huixia Judy.
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{crq}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
library(crqMI)
set.seed(1234)
datasim=sim(n=100)
data = datasim$obs.data

tau=c(0.3,0.5)

############  if SEs are not required #####
out=crqMI(Y= data$Y, Z= data$Z, X=data$X, d= datasim$d,
          delta.v=data$delta.v, delta.t=data$delta.t,  Imp=NULL,
          M = 10, tau=tau, ctype= datasim$ctype,seed=T,
          pi.fam="logistic",se=FALSE,shrinkage=TRUE)

MI=out$estimate$MI  ## MI estimates
CC=out$estimate$CC  ## CC estimates
true.par =datasim$true.par[,paste(tau)]
print(round(cbind(MI-true.par,CC-true.par),3)) # the estimation errors

############  need to compute SEs and shrinkage estimates#####
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


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
