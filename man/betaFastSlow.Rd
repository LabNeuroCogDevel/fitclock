\name{betaFastSlow}
\alias{betaFastSlow}
\title{beta distribution worker class (handled as environment shared by epsilon and rho)}
\usage{
betaFastSlow(w = NULL, alpha_slow = 1.01, beta_slow = 1.01,
  alpha_fast = 1.01, beta_fast = 1.01, decay = 1, lastUpdateTrial = 1L,
  local_RT_learning_rate = 0.1)
}
\arguments{
  \item{w=NULL}{Shared workspace environment}

  \item{alpha_slow}{alpha shape parameter for slow (above
  average) RTs}

  \item{beta_slow}{beta shape parameter for slow RTs}

  \item{alpha_fast}{alpha parameter for fast RTs}

  \item{beta_fast}{beta parameter for fast RTs}

  \item{decay}{decay parameter (knowledge of fast and slow
  reward dists decays with time)}

  \item{lastUpdateTrial}{what is the latest trial number
  for which beta updates were computed (necessary because
  betas are shared by rho and epsilon -- avoid double
  update)}

  \item{local_RT_learning_rate}{learning rate modulating
  running calculation of recent mean RT.}
}
\description{
beta distribution worker class (handled as environment
shared by epsilon and rho)
}
\keyword{internal}

