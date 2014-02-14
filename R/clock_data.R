#dataset object is supposed to represent a group of data over which a single set of parameters (theta) are to be found
#so in the group case, get best-fitting parameters for all subjects and runs
#for a single subject, best fitting-parameters across subject's runs
#for a single run, best-fitting parameters for this run

#have predict in alg detect the class of the dataset to determine how to proceed

#' dataset object for group-level data (multiple subjects with multiple runs)
#' 
#' @section Fields:
#'    \describe{
#'      \item{\code{subjects}:}{ \code{list} of clockdata_subject objects defining the group. }
#'      \item{\code{dataset_name}:}{ optional character vector identifying the dataset. }
#'    }
#' 
#' @section Methods:
#'    \describe{
#'      \item{\code{add_subjects(...)}:}{ add one or more clockdata_subject objects to the dataset. }
#'      \item{\code{delete_subjects(subjects_to_delete=NULL)}:}{ delete subjects with IDs specified by character vector \code{subjects_to_delete}. }  
#'    }
#' 
#' @importFrom methods setRefClass
#' @export clockdata_group
#' @exportClass clockdata_group
clockdata_group <- setRefClass(
    Class="clockdata_group",
    fields=list(
        subjects="list",
        dataset_name="character"
    ),
    methods=list(
        initialize=function(dataset_name=NULL, subjects=NULL) {
          if (!is.null(dataset_name)) { dataset_name <<- dataset_name }
          if (!is.null(subjects)) {
            if (!is.list(subjects) && class(subjects) == "clockdata_subject") {
              subjects <<- list(subjects) #single subject
            } else if (is.list(subjects) && all(sapply(subjects, class) == "clockdata_subject")) {
              subjects <<- subjects
            } else { stop ("subjects should be a list of clockdata_subject objects or a single clockdata_subject object")}
          }
        },
        add_subjects=function(...) {
          if (!is.null(subjects)) {
            existingSubjectIDs <- sapply(subjects, function(s) { s$subject_ID } )
          } else { existingSubjectNums <- c() }
          
          s_objs <- as.list(match.call())[-1L] #first element of match.call is a class for refMethodDef for the alg object
          
          #verify that all arguments are clockdata_subject objects
          if (!all(sapply(subjects, class) == "clockdata_subject")) { stop("All arguments to add_subjects must be clockdata_subject objects.") }
          
          newSubjectIDs <- sapply(s_objs, function(s) { s$subject_ID } )
          
          if (any(dupes <- newSubjectIDs %in% existingSubjectIDs)) {
            cat("Cannot add a subject with the same ID as an existing subject\n")
            stop("The following subjects conflict: ", paste(newSubjectIDs[dupes], collapse=", "))
          }
          
          subjects <<- c(subjects, s_objs) #concatenate
          #subjects <<- subjects[sort(sapply(subjects, function(s) { s$subject_ID }))] #sort in ascending order
          
        },
        delete_subjects=function(subjects_to_delete=NULL) {
          if (is.null(subjects_to_delete)) { return(invisible()) }
          if (!is.character(subjects_to_delete)) { stop("delete_subjects expects a character vector of subject IDs to delete") }
          
          subjects_to_keep <- which(! sapply(subjects, function(s) { s$subject_ID } ) %in% subjects_to_delete)
          unmatched <- which(! subjects_to_delete %in% subjects)
          
          subjects <<- subjects[subjects_to_keep]
          
          if (length(unmatched) > 0L) {
            warning("Unable to find the following subjects for deletion: ", paste0(unmatched, collapse=", "))
          }
        }
    
    
    )
)

#' dataset object for subject-level data (multiple runs within a single subject)
#' 
#' If the \code{csv_file} field is passed at instantiation, data will be imported from the file specified.
#' 
#' @section Fields:
#'    \describe{
#'      \item{\code{subject_ID}:}{ unique identifier string for subject. }
#'      \item{\code{runs}:}{ \code{list} of clockdata_run objects defining the subject. }
#'      \item{\code{csv_file}:}{ comma-separated file containing behavior for subject (one or more runs, each with many trials) }
#'    }
#' 
#' @section Methods:
#'    \describe{
#'      \item{\code{add_runs(r_objs)}:}{ add a list of clockdata_run objects, \code{r_objs} to the subject dataset. }
#'      \item{\code{delete_subjects(subjects_to_delete=NULL)}:}{ delete subjects with IDs specified by character vector \code{subjects_to_delete}. }
#'      \item{\code{delete_runs(runs_to_delete=NULL)}:}{ a numeric vector of run numbers to delete. }
#'      \item{\code{plot_runs()}:}{ plot observed RTs and reward for all runs. (not implemented yet }  
#'    }
#' 
#' @importFrom ggplot2 ggplot
#' @importFrom methods setRefClass
#' @importFrom fmri fmri.stimulus
#' @importFrom fmri fmri.design
#' @export clockdata_subject
#' @exportClass clockdata_subject
clockdata_subject <- setRefClass(
    Class="clockdata_subject",
    fields=list(
        subject_ID="character",
        runs="list",
        csv_file="character"
    ),
    methods=list(
        initialize=function(subject_ID=NULL, runs=NULL, csv_file=NULL, ...) {
          if (is.null(subject_ID)) { 
            stop("clockdata_subject requires subject_ID at initialization.") 
          } else {
            subject_ID <<- subject_ID
          }
          if (!is.null(runs)) {
            if (!is.list(runs) && class(runs) == "clockdata_run") {
              runs <<- list(runs) #single run
            } else if (is.list(runs) && all(sapply(runs, class) == "clockdata_run")) {
              runs <<- runs
            } else { stop ("runs should be a list of clockdata_run objects or a single clockdata_run object")}
          }
          if (!is.null(csv_file)) {
            csv_file <<- csv_file
            message("Importing data from csv file: ", csv_file)
            .self$import_runs_from_csv()
          }
          callSuper(...)
        },
        import_runs_from_csv=function(fname=NULL) {          
          if (is.null(fname)) { fname <- csv_file } #use subject data file
          stopifnot(file.exists(fname))
          #for now, assuming a relatively fixed format for these CSV files
          
          sdata <- read.csv(fname, header=TRUE)
          d_split <- split(sdata, sdata$run)
          
          for (d in d_split) {
            r <- clockdata_run(
                run_number=d$run[1L],
                RTobs=d$rt, 
                Reward=d$score,
                global_trial_number=d$trial,
                rew_function=as.character(d$rewFunc[1L]),
                run_condition=as.character(d$emotion[1L]),
                orig_data_frame=d)
            
            if ("clock_onset" %in% names (d)) { r$clock_onset <- d$clock_onset }
            if ("feedback_onset" %in% names (d)) { r$feedback_onset <- d$feedback_onset }
            if ("iti_onset" %in% names (d)) { r$iti_onset <- d$iti_onset }
            
            .self$add_runs(r)
          }
          
        },
        #add_runs=function(...) {
        add_runs=function(r_objs) {
          if (!is.null(runs)) {
            existingRunNums <- sapply(runs, function(r) { r$run_number } )
          } else { existingRunNums <- c() }
          
          #this is leading to all kinds of problems here -- objects are staying as language
          #then, when evaluated, they expect the (local) data from the initialization
          #core problem is that match.call seems to always get a language expression, whereas I need an object
          #makes me think that the add_params, although working, is actually double-instantiating all objects
          #once as parameters to the function
          #once at the eval step
          
          #r_objs <- as.list(match.call())[-1L] #first element of match.call is a class for refMethodDef for the alg object
          #r_objs <- lapply(r_objs, eval) #force initialization of run objects (otherwise stays as a language expression)
          
          #verify that all arguments are clockdata_run objects
          
          if (!is.list(r_objs)) {
            if (class(r_objs) != "clockdata_run") { stop("add_runs expects a clockdata_run object or a list of such objects.") }
            r_objs <- list(r_objs)
          } else if (is.list(r_objs) && !all(sapply(r_objs, class) == "clockdata_run")) { 
            stop("add_runs expects a clockdata_run object or a list of such objects.") 
          }
          
          newRunNums <- sapply(r_objs, function(r) { r$run_number } )
          
          if (any(dupes <- newRunNums %in% existingRunNums)) {
            cat("Cannot add a run with the same number as an existing run\n")
            stop("The following runs conflict: ", paste(newRunNums[dupes], collapse=", "))
          }
          
          runs <<- c(runs, r_objs) #concatenate
          runs <<- runs[sort(sapply(runs, function(r) { r$run_number }))] #sort in ascending order
          
        },
        delete_runs=function(runs_to_delete=NULL) {
          if (is.null(runs_to_delete)) { return(invisible()) }
          if (!is.numeric(runs_to_delete)) { stop("delete_runs expects a numeric vector of run numbers to delete") }
          
          runs_to_keep <- which(! sapply(runs, function(r) { r$run_number} ) %in% runs_to_delete)
          unmatched <- which(! runs_to_delete %in% runs)
          
          runs <<- runs[runs_to_keep]
          
          ##TODO: Need to sequentially re-number remaining runs
          
          if (length(unmatched) > 0L) {
            warning("Unable to find the following runs for deletion: ", paste0(unmatched, collapse=", "))
          }
        },
        plot_runs=function() {
          
        }
    )
)

#' an object to store the results of a model being fit to a dataset
#' 
#' Returned by the $fit() method of clock_model. 
#' 
#' @section Fields:
#'    \describe{
#'      \item{\code{RTobs}:}{ matrix of observed reaction times. Can be subjects x runs x trials, runs x trials, or just 1 x trials. }
#'      \item{\code{RTpred}:}{ matrix of predicted reaction times. }
#'      \item{\code{Reward}:}{ matrix of obtained rewards. }
#'      \item{\code{total_SSE}:}{ total sum of squared errors across all fitted data. }
#'      \item{\code{SSE}:}{ vector or matrix of SSEs for each subject and run }
#'      \item{\code{AIC}:}{ vector or matrix of AICs for each subject and run }
#'      \item{\code{elapsed_time}:}{ numeric vector of elapsed time for optimizing the parameters }
#'      \item{\code{profile_data}:}{ \code{list} of profile timing information from Rprof, if requested during $fit }
#'      \item{\code{opt_data}:}{ \code{list} of output from optimizer, nlminb. }
#'    }
#' 
#' @importFrom methods setRefClass
#' @importFrom fmri fmri.design
#' @importFrom Hmisc Lag
#' @export clock_fit
#' @exportClass clock_fit
clock_fit <- setRefClass(
    Class="clock_fit",
    fields=list(
        RTobs="matrix", #trial vector (run), run x trial matrix (subject), or subject x run x trial matrix (group)
        RTpred="matrix",
        Reward="matrix",
        total_SSE="numeric", #scalar of total sum of squared errors
        theta="matrix", #named matrix of parameters and bounds
        SSE="numeric", #vector or matrix of SSEs over subjects and runs
        AIC="numeric", 
        elapsed_time="numeric",
        profile_data="list",
        opt_data="list", #list of results from optimizer
        pred_contrib="list",
        clock_onset="matrix",
        feedback_onset="matrix",
        iti_onset="matrix",
        bfs_var_fast="matrix", #trialwise variance of beta distribution for fast responses
        bfs_var_slow="matrix",
        bfs_mean_fast="matrix",
        bfs_mean_slow="matrix",
        ev="matrix", #expected value
        rpe="matrix" #reward prediction error
    ),
    methods=list(
        initialize=function(...) {
          callSuper(...) #default assignment of fields
        },
        #extended from fmri package to support parametric regressor
        fmri.stimulus=function(scans=1, onsets=c(1), durations=c(1), values=c(1),
            rt=3, times=NULL, mean=TRUE, a1 = 6, a2 = 12, b1 = 0.9, b2 = 0.9, cc = 0.35) {
          
          mygamma <- function(x, a1, a2, b1, b2, c) {
            d1 <- a1 * b1
            d2 <- a2 * b2
            c1 <- ( x/d1 )^a1
            c2 <- c * ( x/d2 )^a2
            res <- c1 * exp(-(x-d1)/b1) - c2 * exp(-(x-d2)/b2)
            res
          }
          
          if (is.null(times)) {
            scale <- 1
          } else {
            #upsample time grid by a factor of 100 to get best estimate of hrf at each volume 
            scale <- 100
            onsets <- times/rt*scale
            durations <- durations/rt*scale
            rt <- rt/scale
            scans <- scans*scale
          }
          numberofonsets <- length(onsets)
          
          if (length(durations) == 1) {
            durations <- rep(durations,numberofonsets)
          } else if (length(durations) != numberofonsets)  {
            stop("Length of durations vector does not match the number of onsets!")
          }
          
          if (length(values) == 1) {
            #use the same regressor height (usually 1.0) for all onsets
            values <- rep(values, numberofonsets)
          } else if (length(values) != numberofonsets) {
            stop("Length of values vector does not match the number of onsets!")
          }
          
          stimulus <- rep(0, scans)
          
          for (i in 1:numberofonsets) {
            for (j in onsets[i]:(onsets[i]+durations[i]-1)) {
              stimulus[j] <- values[i]
            }
          }
          stimulus <- c(rep(0,20*scale),stimulus,rep(0,20*scale))
          
          #  zero pad stimulus vector to avoid bounding/edge effects in convolve
          hrf <- convolve(stimulus,mygamma(((40*scale)+scans):1, a1, a2, b1/rt, b2/rt, cc))/scale
          hrf <- hrf[-(1:(20*scale))][1:scans]
          hrf <- hrf[unique((scale:scans)%/%scale)*scale]
          
          dim(hrf) <- c(scans/scale,1)
          
          if (mean) {
            hrf - mean(hrf)
          } else {
            hrf
          }
        },
        build_design_matrix=function(
            regressors=NULL,
            event_onsets=NULL,
            durations=NULL,
            convolve=TRUE) {
          
#          if (is.null(regressors)) {
#            regressors <- c("buttonPress", dimnames(f$pred_contrib[[1L]])[[1L]]) #use all param names from first run
#          }
          
          browser()
          if (length(unique(sapply(list(regressors, event_onsets, durations), length))) != 1L) { stop("regressors, event_onsets, and durations must have the same length") }
          #require(fmri)
          
          #determine the last fMRI volume to be analyzed
          last_fmri_volume <- apply(.self$iti_onset, 1, function(itis) {
                ceiling(itis[length(itis)] + 12.0 ) #fixed 12-second ITI after every run
              })
          
#          if ("buttonPress" %in% regressors) {
#            #browser()
#            press_time <- .self$clock_onset + .self$RTobs/1000
#            press_time.hrf <- lapply(1:nrow(press_time), function(run) { 
#                  fmri.stimulus(scans=last_fmri_volume[run], times=press_time[run,], durations=0, rt=1.0) #hard-coded 1.0s TR for now
#                })
#          }
#          #drop button press
#          regressors <- regressors[-which(regressors=="buttonPress")]
          
          designCols <- c()
          
          #build design matrix in the order specified
          dmat <- sapply(1:length(regressors), function(r) {
                #Note: all of the regressor computations here should result in nruns x ntrials matrices
                #regressor values
                if (regressors[r]=="rel_uncertainty") {
                  reg <- abs(bfs_var_fast - bfs_var_slow)
                } else if (regressors[r]=="mean_uncertainty") {
                  reg <- abs(bfs_var_fast + bfs_var_slow)/2 #trialwise average uncertainty for fast and slow responses
                } else if (regressors[r]=="button_press") {
                  reg <- matrix(1.0, dim(RTobs)) #standard task indicator regressor
                } else if (regressors[r]=="rpe_pos") {
                  reg <- apply(rpe, c(1,2), function(x) { if (x > 0) x else 0 })  
                } else if (regressors[r]=="rep_neg") {
                  reg <- apply(rpe, c(1,2), function(x) { if (x < 0) abs(x) else 0 }) #abs so that greater activation scales with response to negative PE
                }
                
                #replace missing values with 0 for clarity in convolution
                reg[which(is.na(reg))] <- 0
                
                #determine onset times for events
                if (event_onsets[r] == "clock_onset") {
                  times <- clock_onset
                } else if (event_onsets[r] == "feedback_onset") {
                  times <- feedback_onset
                } else if (event_onsets[r] == "iti_onset") {
                  times <- iti_onset
                } else if (event_onsets[r] == "rt") {
                  times <- clock_onset + RTobs/1000
                }
                
                if (durations[r] == "rt") {
                  durations <- RTobs/1000
                } else if (durations[r] == "feedback_duration") {
                  durations <- iti_onset - feedback_onset
                } else if (durations[r] == "iti_duration") {
                  #need to lag clock matrix: iti duration is: clock_onset(t+1) - iti_onset(t)
                  #lag_clock <- apply(clock_onset)
                } else if (!is.na(suppressWarnings(as.numeric(durations[r])))) {
                  #user-specified scalar duration (not currently supporting a user-specified per-regressor vector of durations) 
                  durations <- matrix(as.numeric(durations[r]), dim(RTobs))
                }
                
                conv.reg <- lapply(1:nrow(reg), function(run) { 
                      fmri.stimulus(scans=last_fmri_volume[run], values=reg[run,], times=times[run,], durations=durations[run,], rt=1.0) #hard-coded 1.0s TR for now
                    })
                
                conv.reg
              })
              
              #returns a 2-d list of runs x regressors. Needs to stay as list since runs vary in length, so aggregate is not rectangular
              #bind by run
              runMats <- sapply(1:dim(dmat)[1L], function(r) {
                    m <- do.call(cbind, dmat[r,])
                    dm <- fmri.design(m, order=2)
                    dm
                  })
              
              browser()
        }
    )

)


#' dataset object for run-level data (multiple trials)
#' 
#' Note that the workspace used during the fitting process resides at the level of the
#' clockdata_run object since prediction equation is principally intended to fit
#' multiple trials in a given run. Optimal parameters across runs or subjects and runs
#' are derived by summing SSEs for each run, essentially finding a single set of parameter
#' values that fit all runs/subjects reasonably well.
#' 
#' @section Fields:
#'    \describe{
#'      \item{\code{w}:}{ shared workspace environment used during fit. }
#'      \item{\code{run_number}:}{ number of this run in a multi-run sequence. Important when some values (e.g., V) carry across runs. }
#'      \item{\code{SSE}:}{ sum of squared errors for this run: predicted versus observed RTs }
#'      \item{\code{RTobs}:}{ vector of observed RTs }
#'      \item{\code{Reward}:}{ vector of obtained rewards. }
#'      \item{\code{avg_RT}:}{ average reaction time for this run (maybe overridden by global RT) }
#'      \item{\code{global_trial_number}:}{ vector of trial numbers in the overall experiment (1..runs x trials/run) }
#'      \item{\code{rew_function}:}{ string denoting run reward contingency (e.g., "IEV") }
#'      \item{\code{run_condition}:}{ string denoting some other aspect of run (e.g., emotion, such as "happy") }
#'      \item{\code{by_lookup}:}{ at fit-time, clock_model copies in a named character vector that is the union of all relevant fields for run definition (usually rew_function + run_condition) }
#'      \item{\code{orig_data_frame}:}{ optional data.frame from original experiment run containing full saved data (in case there are additional variables of interest) }
#'    }
#' 
#' @section Methods:
#'    \describe{
#'      \item{\code{initialize_workspace(prior_w)}:}{ sets up shared environment for fitting, \code{w}, based on \code{prior_w}, workspace of prior run. }
#'      \item{\code{plot_RTs()}:}{ plot observed RTs (partially implemented) }  
#'    }
#' 
#' @importFrom ggplot2 ggplot
#' @importFrom methods setRefClass
#' @export clockdata_run
#' @exportClass clockdata_run
clockdata_run <- setRefClass(
    Class="clockdata_run",
    fields=list(
        w="environment",
        SSE="numeric",
        run_number="numeric", #number of this run within a multi-run session
        RTobs="numeric", #vector of observed RTs
        Reward="numeric", #vector of obtained rewards
        avg_RT="numeric", #average reaction time, used for some parameter fits
        global_trial_number="numeric", #vector of trial numbers in the overall experiment (1..runs x trials/run)
        rew_function="character",
        run_condition="character", #optional string specifying the conditions for this run (e.g., fear faces)
        by_lookup="character", #at fit-time, alg copies in a named character vector that is the union of all relevant fields for run definition (usually rew_function + run_condition) 
        orig_data_frame="data.frame", #optional data.frame from original experiment run containing full saved data (in case there are additional variables of interest)
        clock_onset="numeric", #vector of clock stimulus onset times (in seconds)
        feedback_onset="numeric", #vector of feedback stimulus onset times (in seconds)
        iti_onset="numeric" #vector of fixation iti onset times (in seconds)
    ),
    methods=list(
        initialize=function(run_number=NA_integer_, RTobs=NA_integer_, Reward=NA_integer_, global_trial_number=NA_integer_,
            rew_function=NULL, run_condition=NA_character_, ...) {
          
          if (is.na(run_number[1L])) { stop("At this point, clockdata_run must have a run number indicating temporal order.") }
          if (is.na(RTobs[1L])) { stop("construction of clockdata_run requires observed reaction times (RTobs)") }
          if (is.na(Reward[1L])) { stop("construction of clockdata_run requires observed rewards (Reward)") }
          
          run_number <<- run_number
          RTobs <<- RTobs
          avg_RT <<- mean(RTobs, na.rm=TRUE)
          Reward <<- Reward
          
          if (is.na(global_trial_number[1L])) { 
            warning("global_trial_number not provided. Defaulting to 1..t")
            global_trial_number <<- 1:length(RTobs)
          } else {
            global_trial_number <<- global_trial_number
          }
          
          if (length(unique(sapply(list(RTobs, Reward, global_trial_number), length))) > 1L) {
            stop("RTobs, Reward, and global_trial_number must have the same length")
          }
          
          if (is.null(rew_function)) { stop("Construction of clockdata_run requires rew_function") }
          rew_function <<- rew_function
          run_condition <<- run_condition
          
          initialize_workspace() #initial setup of workspace
          
          callSuper(...) #assign unmatched args
          
        },
        initialize_workspace=function(prior_w) {
          #initialize shared workspace
          #clear out old values if refitting
          w <<- new.env(parent = baseenv()) #make sure we do not accidentally inherit objects from .Globalenv as parent 
          
          w$ntrials <<- length(Reward)
          w$RTobs   <<- RTobs #initialize fields
          w$Reward  <<- Reward
          w$RTpred  <<- rep(NA_real_, w$ntrials) #initialize empty predicted RT
          w$RTpred[1L] <<- w$RTobs[1L] #cannot predict first trial behavior per se, so use observed RT so that t=1 doesn't contribute to SSE
          w$rpe     <<- rep(NA_real_, w$ntrials) #vector of reward prediction errors
          w$V       <<- rep(NA_real_, w$ntrials) #expected value vector
          
          if (!missing(prior_w) && !is.null(prior_w) && is.environment(prior_w)) {
            #cat("using priorV: ", prior_w$V[length(prior_w$V)], "\n") #carry forward expected value from last trial
            w$V[1L]   <<- prior_w$V[length(prior_w$V)] #carry forward expected value from last trial
          } else {
            #cat("using 0 V\n")
            w$V[1L]   <<- 0 #no assumption on prior expected value
          }
          
          w$avg_RT  <<- avg_RT
          w$alphaV  <<- 0.1 # learning rate for critic (V).  Just set this to avoid degeneracy
        },
        plot_RTs=function() {
          stopifnot(exists(w$RTobs))
          if (exists(w$RTpred)) {
            rtDf <- data.frame(run_condition=run_condition, rew_function=rew_function,
                rt=c(w$RTobs, w$RTpred), rt_type=rep(c("observed", "predicted"), each=length(w$RTobs))
            )
            rtPlot <- ggplot(rtDf, aes(x=trial, y=rt, color=rt_type))
          } else {
            message("Predicted RTs not found. Plotting only observed RTs.")
            rtDf <- data.frame(run_condition=run_condition, rew_function=rew_function, trial=1:length(w$RTobs),
                rt=w$RTobs, rt_type=rep("observed", length(w$RTobs))
            )
            rtPlot <- ggplot(rtDf, aes(x=trial, y=rt))
          }
          
          rtPlot <- rtPlot + geom_line(size=1.1) + ggtitle(paste(rtDf$rew_function, rtDf$run_condition, sep=", ")) +
              theme_bw(base_size=14)
          
          print(rtPlot)
          return(invisible(rtPlot))
          
        }
    )
)

