\docType{class}
\name{clockdata_run-class}
\alias{clockdata_run-class}
\title{dataset object for run-level data (multiple trials)}
\description{
Note that the workspace used during the fitting process
resides at the level of the clockdata_run object since
prediction equation is principally intended to fit multiple
trials in a given run. Optimal parameters across runs or
subjects and runs are derived by summing SSEs for each run,
essentially finding a single set of parameter values that
fit all runs/subjects reasonably well.
}
\section{Fields}{
  \describe{ \item{\code{w}:}{ shared workspace environment
  used during fit. } \item{\code{run_number}:}{ number of
  this run in a multi-run sequence. Important when some
  values (e.g., V) carry across runs. } \item{\code{SSE}:}{
  sum of squared errors for this run: predicted versus
  observed RTs } \item{\code{RTobs}:}{ vector of observed
  RTs } \item{\code{Reward}:}{ vector of obtained rewards.
  } \item{\code{avg_RT}:}{ average reaction time for this
  run (maybe overridden by global RT) }
  \item{\code{global_trial_number}:}{ vector of trial
  numbers in the overall experiment (1..runs x trials/run)
  } \item{\code{rew_function}:}{ string denoting run reward
  contingency (e.g., "IEV") } \item{\code{run_condition}:}{
  string denoting some other aspect of run (e.g., emotion,
  such as "happy") } \item{\code{by_lookup}:}{ at fit-time,
  clock_model copies in a named character vector that is
  the union of all relevant fields for run definition
  (usually rew_function + run_condition) }
  \item{\code{orig_data_frame}:}{ optional data.frame from
  original experiment run containing full saved data (in
  case there are additional variables of interest) } }
}

\section{Methods}{
  \describe{ \item{\code{initialize_workspace(prior_w)}:}{
  sets up shared environment for fitting, \code{w}, based
  on \code{prior_w}, workspace of prior run. }
  \item{\code{plot_RTs()}:}{ plot observed RTs (partially
  implemented) } }
}

