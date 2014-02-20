#' fitclock
#' @name fitclock
#' @docType package
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
#' @source Laboratory for Neurocognitive Development, University of Pittsburgh
#' @author Michael Hallquist, 2013-01-27
NULL
