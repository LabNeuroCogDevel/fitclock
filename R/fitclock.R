#' @title R package to fit data from variants of Michael Frank clock task
#' @name fitclock
#' @docType package
#' @description Object-oriented implementation of the time clock algorithm, a
#'     multiple regression-style prediction of reaction time on each trial of the
#'     clock task. Implementation fits the standard model including strategic
#'     exploration (adapting toward uncertainty) from Frank et al. 2009 Nat Neuro,
#'     and also allows for parameters to vary by contingency or run condition
#'     (e.g., emotion). Functions are also provided to convolve trial-wise
#'     regressors with a hemodynamic response function to generate an fMRI design
#'     matrix.
#' 
#' @details \tabular{ll}{
#' Package: \tab fitclock\cr
#' Type: \tab Package\cr
#' Version: \tab 0.7\cr
#' Date: \tab 2014-01-27\cr
#' License: \tab GPL\cr
#' Depends: \tab methods\cr
#' }
#' @author Michael Hallquist, code adapted from Michael Frank
#' 
#' Maintainer: Michael Hallquist <michael.hallquist@@gmail.com>
#' @references Frank, M. J., Doll, B. B., Oas-Terpstra, J., & Moreno,
#' F. (2009). Prefrontal and striatal dopaminergic genes predict individual
#' differences in exploration and exploitation. Nature Neuroscience, 12(8),
#' 1062-1068. doi:10.1038/nn.2342
NULL

#' Clock behavioral data for 400 trials across 8 fMRI runs from a single subject (008)
#' 
#' Dataset is a .csv file consisting of the following fields 
#' 
#' \itemize{
#'   \item run. fMRI run (1:8)  
#'   \item trial. Global trial number (1:400) 
#'   \item rewFunc. Reward contingency (DEV, IEV, CEV, CEVR)
#'   \item emotion. Face emotion of central stimulus (scram, fear, happy) 
#'   \item magnitude. Magnitude of expected reward given RT  
#'   \item probability. Probability of expected reward given RT
#'   \item score. Obtained reward (probabilistic receipt of payoff)
#'   \item ev. Expected value of response given RT (magnitude*probability)
#'   \item rt. Reaction time (ms)
#'   \item clock_onset. Run onset time of clock stimulus (sec) 
#'   \item isi_onset. Run onset time of 50ms ISI (sec)
#'   \item feedback_onset. Run onset time of 850ms reward feedback (sec)
#'   \item iti_onset. Run onset time of ITI (sec)
#'   \item iti_ideal. Desired ITI duration (sec) based on fMRI design optimization.
#'   \item image. Image file displayed on screen.
#' }
#' 
#' @docType data
#' @keywords datasets
#' @name clocksubject_fMRI_008jh
#' @usage data(clocksubject_fMRI_008jh)
#' @format A data.frame with 400 rows and 15 variables
NULL

#' @name clocksubject_1000
#' @title clocksubject_1000
#' @description This dataset is an example of multi-run single subject data
#'       from the clock task.
#' @docType data
#' @keywords datasets
#' @usage data(clocksubject_1000)
#' @format a data.frame containing 9 runs of 40 trials, representing a 3x3 design:
#'       rewFunc: IEV, DEV, CEVR; and emotion: fear, happy, scrambled.
#' @source Laboratory for Developmental Personality Neuroscience, Penn State University
#' @author Michael Hallquist, 2013-01-27
NULL
