fitclock
===============

The `fitclock` package implements related algorithms to fit
reaction times and rewards from the clock task developed by
Michael Frank and colleagues. 

Installation
------------

To install the package from github, run the following in R:
 
```r
if (!require(devtools)) { install.packages("devtools"); require(devtools) }

install_github("fitclock", "LabNeuroCogDevel", args="--byte-compile")

#load package
library(fitclock)
```

Background and Architecture
--------

The package includes an example dataset of a single subject who completed
eight runs of the clock task (50 trials per run). This dataset is called
`clocksubject_fMRI_008jh`, and you can see the basic structure:

```r
str(clocksubject_fMRI_008jh)
head(clocksubject_fMRI_008jh)
```

The fitclock package can fit a vector of parameters for a group (optimize
parameters of multiple subjects with multiple runs), a subject (many runs),
or a single run. The data storage structure for the package defines three object:
`clockdata_group`, `clockdata_subject`, and `clockdata_run`. Not surprisingly,
`clockdata_group` objects are composed of multiple `clockdata_subject` objects,
which are in turn composed of multiple `clockdata_run` objects.

The basic time-clock algorithm is designed to fit a series of trials with a
consistent contingency (i.e., a run). Thus, observed reaction times and rewards
are stored at the `clockdata_run` level, whereas the `clockdata_subject` and
`clockdata_group` objects are more storage lists for multiple runs or subjects,
respectively.

At this time, the `clockdata_group` functionality is not fully implemented,
but scaffolding is in place. Fitting subjects or runs, however, is complete,
and examples for each are provided below.

Importing data into fitclock
----------------------------

To import behavioral data into the fitclock package, write the behavioral data
for all runs as a .csv file, perhaps using the ClockToCSV.m script. This csv file should
have approximately the following form:

```r
'data.frame':	400 obs. of  15 variables:
 $ run           : int  1 1 1 1 1 1 1 1 1 1 ...
 $ trial         : int  1 2 3 4 5 6 7 8 9 10 ...
 $ rewFunc       : Factor w/ 4 levels "CEV","CEVR","DEV",..: 3 3 3 3 3 3 3 3 3 3 ...
 $ emotion       : Factor w/ 3 levels "fear","happy",..: 3 3 3 3 3 3 3 3 3 3 ...
 $ magnitude     : int  77 81 67 79 67 73 66 64 73 74 ...
 $ probability   : num  0.66 0.359 0.79 0.295 0.738 ...
 $ score         : int  0 0 67 0 67 73 66 0 73 0 ...
 $ ev            : num  50.8 29.2 53 23.3 49.3 ...
 $ rt            : int  934 3001 401 3818 605 384 201 334 934 501 ...
 $ clock_onset   : num  8.11 15.12 20.19 23.67 33.56 ...
 $ isi_onset     : num  9.06 18.14 20.61 27.51 34.17 ...
 $ feedback_onset: num  9.11 18.19 20.66 27.56 34.22 ...
 $ iti_onset     : num  10 19.1 21.6 28.5 35.1 ...
 $ iti_ideal     : num  5.1 1.1 2.1 5.1 0.1 0.1 1.1 0.1 3.1 3.1 ...
 $ image         : Factor w/ 75 levels "fear_1.png","fear_10.png",..: 62 73 63 54 61 70 57 67 71 75 ...
```

where `run` is the run number (from first to last), `trial` is the global trial number (across all runs),
`rewFunc` is the reward contingency function (so far: CEV, CEVR, DEV, IEV), `emotion` is the emotion
condition for the run (if relevant), `magnitude` is the magnitude of probabilistic reward, `probability`
is the probability of receiving the reward, `score` is the actual score obtained by the player, and `ev` is
the long-run average expected value for the RT given the contingency (i.e., magnitude*probability).

In addition, for the fMRI implementation, run timing information, including `clock_onset` (time of the onset of
clock stimulus within the run, in seconds), `isi_onset` (onset of 50ms isi between clock and feedback),
`feedback_onset` (onset of reward outcome), and `iti_onset` (time of fixation cross between runs) can be included
in the .csv file. These fields are not relevant for the behavioral fits, but are preserved for use in the creation of
a design matrix for fMRI analysis (see `clock_fit$build_design_matrix`).

Finally, other fields, such as the ideal ITI length or the image displayed can be included, but are not used at this time.

Setting up an analysis model
----------------------------

Brief worked examples

```r
library(fitclock)

#setup subject
jh <- clockdata_subject(subject_ID="008_jh", dataset=clocksubject_fMRI_008jh)

#setup model to fit RT differences
expDiff_model <- clock_model(fit_RT_diffs=TRUE)
expDiff_model$add_params(
    meanRT(max_value=4000),
    autocorrPrevRT(),
    goForGold(),
    go(),
    noGo(),
    meanSlowFast(),
    exploreBeta()
)

#tell model which dataset to use
expDiff_model$set_data(jh)

#test the incremental contribution of each parameter to AIC (fit)
incr_fit <- expDiff_model$incremental_fit(njobs=6)

#vector of AIC values
sapply(incr_fit$incremental_fits, "[[", "AIC")

#fit full model, using 5 random starts and choosing the best fit
f <- expDiff_model$fit(random_starts=5)

#design matrix matching Badre et al. 2012 Neuron
#writes 3-column timing files (onset, duration, parametric value) for each run to the "run_timing"
#directory within the current working directory
d <- f$build_design_matrix(regressors=c("mean_uncertainty", "rel_uncertainty", "rpe_pos", "rpe_neg", "rt"), 
    event_onsets=c("clock_onset", "clock_onset", "feedback_onset", "feedback_onset", "feedback_onset"), 
    durations=c("rt", "rt", "feedback_duration", "feedback_duration", 0), baselineCoefOrder=2, writeTimingFiles=TRUE)

#fit moving average smoothed RTs
expDiff_model$smooth <- 5 #window size of 5 trials
expDiff_model$set_data(jh) #need to set dataset again to trigger use of smoothed RTs
fsmooth <- expDiff_model$fit()

###
#model variant allowing for negative epsilon (exploit)
#sticky choice variant (see Badre et al. 2012 Neuron)

#test sticky choice
negEps <- clock_model()
negEps$add_params(
    K=meanRT(max_value=4000),
    stickyChoice=stickyChoice(),
    alphaG=go(),
    alphaN=noGo(),
    rho=meanSlowFast(),
    epsilonBeta=exploreBeta(min_value=-100000) #allow negative epsilon
)

negEps$set_data(jh)
fnegEps <- negEps$fit()

```

------------------------------------------
Example of a simple delta_v learning model

The basic delta_v model tries to find the learning rate that minimizes reward prediction errors
(i.e., the discrepancy between expected value and obtained rewards). This model is divorced from
fitting/predicting RTs per se, and instead tries to understand whether the agent is tracking value.

```r

library(fitclock)

jh <- clockdata_subject(subject_ID="008_jh", dataset=clocksubject_fMRI_008jh)

vm <- deltavalue_model(clock_data=jh, alphaV=0.1)
vm$predict() #SSE of predicted - observed rewards at learning rate of 0.1
V_0p1 <- vm$V #v matrix for learning rate of 0.1
f <- vm$fit() #estimate learning rate as a free parameter
V_free <- vm$V #v matrix for free parameter


#design matrix: clock onset, feedback_onset, EV, PE-, PE+
d <- f$build_design_matrix(regressors=c("clock", "feedback", "ev", "rpe_neg", "rpe_pos"), 
    event_onsets=c("clock_onset", "feedback_onset", "feedback_onset", "feedback_onset", "feedback_onset"), 
    durations=c(0, 0, "feedback_duration", "feedback_duration", "feedback_duration"), baselineCoefOrder=2, writeTimingFiles="AFNI",
    runVolumes=c(223,273,280,244,324,228,282,310))

```