\name{meanRT}
\alias{meanRT}
\title{K: baseline response speed parameter.}
\usage{
meanRT(min_value = 100, max_value = 5000, init_value = 1000,
  cur_value = init_value, par_scale = 1000, by = NULL)
}
\arguments{
  \item{min_value}{Lower bound for K parameter.}

  \item{max_value}{Upper bound for K parameter.}

  \item{init_value}{Initial value for K parameter.}

  \item{cur_value}{Current value for K parameter.}

  \item{par_scale}{Expected parameter scale (log scale)}

  \item{by}{character vector defining one or more run-level
  fields over which this parameter varies.}
}
\description{
Essentially the RT intercept in the multiple regression
prediction equation.
}

