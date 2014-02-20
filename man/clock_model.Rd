\docType{class}
\name{clock_model}
\alias{clock_model}
\alias{clock_model-class}
\title{Object-oriented implementation of the Frank clock algorithm.}
\description{
The clock_model class defines a predictive model to fit RT
and reward data using data from the clock task. The
algorithm can fit group-level data (optimize theta over all
subjects and runs), subject-level multi-run data (optimize
theta over multiple runs within subject), or single-run
data.
}
\details{
At initialization, data to be fit can be passed using the
\code{clock_data} parameter, and should be of class
clockdata_group (for group data), clockdata_subject (for
multi-run subject data), or clockdata_run (for single-run
data).

After initializing the object, one or more parameter
objects are added to the prediction equation using the
$add_params(...) method. Bounds on the parameters can be
specified at initialization.
}
\section{Fields}{
  \describe{ \item{\code{params}:}{ \code{list} of param
  objects used to predict RTs. Note: parameters are
  evaluated in order. } \item{\code{noiseWt}:}{ A numeric
  scalar (in milliseconds) specifying how much random
  uniform noise to add to predicted RTs. Default: 0 }
  \item{\code{SSE}:}{ Total sum of squared errors for
  fitted theta values. } \item{\code{AIC}:}{ Akaike
  Information Criterion for fitted theta values. }
  \item{\code{clock_data}:}{ Clock data to be fit. Should
  be of class clockdata_group, clockdata_subject, or
  clockdata_run. DO NOT ASSIGN YOURSELF. }
  \item{\code{use_global_avg_RT}:}{ Whether to use the
  global reaction time across all runs within subject as
  reference for scale and K. }
  \item{\code{global_avg_RT}:}{ Average reaction time
  across all runs within subject. } \item{\code{all_by}:}{
  character vector that is the union of all run-level
  fields that define vary-by-run parameters }
  \item{\code{clock_fit}:}{ clock_fit object storing all
  information about result of algorithm fit (e.g., expected
  value) } }
}

\section{Methods}{
  \describe{ \item{\code{add_params(...)}:}{ Method to add
  parameter objects to prediction equation. Should be a
  list of params or a single param. }
  \item{\code{fit(toFit=NULL, random_starts=NULL,
  profile=TRUE)}:}{ Fit \code{clock_data} using the list of
  \code{params}. Data to fit can be set at call to
  \code{$fit()} using the \code{toFit} parameter. }
  \item{\code{incremental_fit=function(toFit=NULL,
  plot=TRUE, njobs=1L)}:}{ Based on the \code{params} list,
  test model subsets building from 1:p parameters. Useful
  for checking AICs to verify utility of each parameter. }
  \item{\code{list_params()}:}{ Returns a data.frame
  describing the parameters, current values, and bounds of
  this model. } \item{\code{params_current()}:}{ Named
  vector of current values for all parameters. }
  \item{\code{params_minimum()}:}{ Named vector of lower
  bounds for all parameters. }
  \item{\code{params_maximum()}:}{ Named vector of upper
  bounds for all parameters. }
  \item{\code{params_initial()}:}{ Named vector of initial
  values for all parameters. }
  \item{\code{params_par_scale()}:}{ Named vector of
  expected log parameter scale for all parameters. }
  \item{\code{reorder_params(neworder)}:}{ Re-orders the
  parameter list according to a named character vector,
  \code{neworder}. Useful because params evaluated in
  order. } \item{\code{set_data(cdata)}:}{ Sets the clock
  dataset to be fit. cdata should be a clockdata object.
  Use this method if you wish to alter the clock_data
  field. }

  }
}
\examples{
#use example fMRI dataset distrbuted with package
exSubj <- clockdata_subject(subject_ID="008", dataset=clocksubject_fMRI_008jh)

exModel <- clock_model()
exModel$add_params(
   meanRT(max_value=4000),
   autocorrPrevRT(),
   goForGold(),
   go(),
   noGo(),
   meanSlowFast(),
   exploreBeta()
)

exModel$set_data(exSubj)
incrFit <- exModel$incremental_fit(njobs=6)
fit <- exModel$fit(random_starts=5)
}

