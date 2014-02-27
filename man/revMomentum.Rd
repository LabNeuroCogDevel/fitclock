\name{revMomentum}
\alias{revMomentum}
\title{gamma: reverse momentum parameter, adjust RTs incrementally in a direction before switching}
\usage{
revMomentum(min_value = 0, max_value = 1, init_value = 0.5,
  cur_value = init_value, par_scale = 0.1, by = NULL)
}
\arguments{
  \item{min_value}{Lower bound for gamma parameter.}

  \item{max_value}{Upper bound for gamma parameter.}

  \item{init_value}{Initial value for gamma parameter.}

  \item{cur_value}{Current value for gamma parameter.}

  \item{par_scale}{Expected parameter scale.}

  \item{by}{character vector defining one or more run-level
  fields over which this parameter varies.}
}
\description{
gamma: reverse momentum parameter, adjust RTs incrementally
in a direction before switching
}

