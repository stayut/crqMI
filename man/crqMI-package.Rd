\name{crqMI-package}
\alias{crqMI-package}
%\alias{crqMI}
\docType{package}
\title{
\packageTitle{crqMI}
}
\description{
\packageDescription{crqMI}
}
\details{
%The DESCRIPTION file:
%\packageDESCRIPTION{crqMI}
%\packageIndices{crqMI}
%~~ An overview of how to use the package, including the most important ~~
%~~ functions ~~
}
\author{
\packageAuthor{crqMI}

Maintainer: \packageMaintainer{crqMI}
}
\references{
%~~ Literature or other references for background information ~~
}
~~ Optionally other standard keywords, one per line, from file KEYWORDS ~~
~~ in the R documentation directory ~~
\keyword{ package }
\seealso{
%#~~ Optional links to other man pages, e.g. ~~
%$~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{
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
}
