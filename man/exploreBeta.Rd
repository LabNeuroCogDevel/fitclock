\name{exploreBeta}
\alias{exploreBeta}
\title{epsilon: adapt to fast/slow responses in proportion to uncertainty about other dist}
\usage{
exploreBeta(min_value = 0, max_value = 1e+05, init_value = 2000,
  cur_value = init_value, par_scale = 1000, by = NULL)
}
\arguments{
  \item{min_value}{Lower bound for epsilon parameter.}

  \item{max_value}{Upper bound for epsilon parameter.}

  \item{init_value}{Initial value for epsilon parameter.}

  \item{cur_value}{Current value for epsilon parameter.}

  \item{par_scale}{Expected parameter scale.}

  \item{by}{character vector defining one or more run-level
  fields over which this parameter varies.}
}
\description{
epsilon: adapt to fast/slow responses in proportion to
uncertainty about other dist
}

