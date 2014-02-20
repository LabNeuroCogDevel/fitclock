\docType{package}
\name{fitclock}
\alias{fitclock}
\alias{fitclock-package}
\title{R package to fit data from variants of Michael Frank clock task}
\description{
Object-oriented implementation of the time clock algorithm,
a multiple regression-style prediction of reaction time on
each trial of the clock task. Implementation fits the
standard model including strategic exploration (adapting
toward uncertainty) from Frank et al. 2009 Nat Neuro, and
also allows for parameters to vary by contingency or run
condition (e.g., emotion). Functions are also provided to
convolve trial-wise regressors with a hemodynamic response
function to generate an fMRI design matrix.
}
\details{
\tabular{ll}{ Package: \tab fitclock\cr Type: \tab
Package\cr Version: \tab 0.7\cr Date: \tab 2014-01-27\cr
License: \tab GPL\cr Depends: \tab methods\cr }
}
\author{
Michael Hallquist, code adapted from Michael Frank

Maintainer: Michael Hallquist <michael.hallquist@gmail.com>
}
\references{
Frank, M. J., Doll, B. B., Oas-Terpstra, J., & Moreno, F.
(2009). Prefrontal and striatal dopaminergic genes predict
individual differences in exploration and exploitation.
Nature Neuroscience, 12(8), 1062-1068. doi:10.1038/nn.2342
}

