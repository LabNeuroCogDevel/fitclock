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
#' If the \code{dataset} field is passed at instantiation, data will be imported from the file specified.
#' 
#' @section Fields:
#'    \describe{
#'      \item{\code{subject_ID}:}{ unique identifier string for subject. }
#'      \item{\code{runs}:}{ \code{list} of clockdata_run objects defining the subject. }
#'      \item{\code{dataset}:}{ comma-separated file containing behavior for subject (one or more runs, each with many trials) }
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
#' @importFrom tools file_ext
#' @export clockdata_subject
#' @exportClass clockdata_subject
clockdata_subject <- setRefClass(
    Class="clockdata_subject",
    fields=list(
        subject_ID="character",
        runs="list",
        dataset="character"
    ),
    methods=list(
        initialize=function(subject_ID=NULL, runs=NULL, dataset=NULL, ...) {
          if (is.null(subject_ID)) { 
            warning("clockdata_subject requires subject_ID at initialization.") #converted to warning so that $copy works
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
          if (!is.null(dataset)) {
            .self$import_runs(dataset)
          }
          callSuper(...)
        },
        import_runs=function(df=NULL) {
          if (is.null(df)) { df <- dataset } #use subject data file, if available
          
          if (inherits(df, "data.frame")) {
            sdata <- df
          } else if (inherits(df, "character")) {
            stopifnot(file.exists(df))
            dataset <<- df
            if (tolower(file_ext(df) == "csv")) {
              message("Importing data from csv file: ", dataset)
              sdata <- read.csv(df, header=TRUE)
            }
          }
          
          #for now, assuming a relatively fixed format for these CSV files
          d_split <- split(sdata, sdata$run)
          
          for (d in d_split) {
            r <- clockdata_run(
                run_number=d$run[1L],
                RTraw=d$rt, 
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
#' @importFrom car vif
#' @importFrom Hmisc Lag
#' @export clock_fit
#' @exportClass clock_fit
clock_fit <- setRefClass(
    Class="clock_fit",
    fields=list(
        ntrials="numeric", #scalar of total trials
        RTobs="matrix", #trial vector (run), run x trial matrix (subject), or subject x run x trial matrix (group)
        RTraw="matrix", #original RTs. RTobs is the RTs that were fitted (e.g., could be raw, smoothed, differenced, etc.)
        RTpred="matrix",
        Reward="matrix",
        total_SSE="numeric", #scalar of total sum of squared errors
        theta="matrix", #named matrix of parameters and bounds
        SSE="numeric", #vector or matrix of SSEs over subjects and runs
        AIC="numeric", 
        nparams="numeric",
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
        rpe="matrix", #reward prediction error
        run_condition="character", #vector of run conditions
        rew_function="character" #vector of reward contingencies        
    ),
    methods=list(
        initialize=function(...) {
          callSuper(...) #default assignment of fields
        },
        populate_fit=function(clock_data=NULL) {
          if (is.null(clock_data)) { return(invisible(NULL)) }
          #based on the clock_data object, populate the fit object with various relevant fields
          #such as the iti_onset times, observed and predicted reaction times, etc.
          #only populate field if it is not already initialized
          
          if (class(clock_data) == "clockdata_run") {
            #convert to list object so thta the lapply calls below work for both single and multiple runs
            clock_data <- list(runs=clock_data)
          }
          
          if (length(ntrials)==0L) { ntrials <<- sum(unlist(lapply(clock_data$runs, function(r) { r$w$ntrials } ))) }
          if (length(RTobs)==0L) { RTobs <<- do.call(rbind, lapply(clock_data$runs, function(r) { r$RTobs })) } #FIXME
          if (length(RTraw)==0L) { RTraw <<- do.call(rbind, lapply(clock_data$runs, function(r) { r$RTraw })) } #FIXME
          if (length(RTpred)==0L) { RTpred <<- do.call(rbind, lapply(clock_data$runs, function(r) { r$w$RTpred })) }
          if (length(Reward)==0L) { Reward <<- do.call(rbind, lapply(clock_data$runs, function(r) { r$Reward }))  }
          if (length(clock_onset)==0L) { clock_onset <<- do.call(rbind, lapply(clock_data$runs, function(r) { if(length(r$clock_onset) == 0 ) NA else r$clock_onset })) }
          if (length(feedback_onset)==0L) { feedback_onset <<- do.call(rbind, lapply(clock_data$runs, function(r) { if(length(r$feedback_onset) == 0 ) NA else r$feedback_onset })) }
          if (length(iti_onset)==0L) { iti_onset <<- do.call(rbind, lapply(clock_data$runs, function(r) { if(length(r$iti_onset) == 0 ) NA else r$iti_onset })) }
          if (length(rew_function)==0L) { rew_function <<- do.call(c, lapply(clock_data$runs, function(r) { r$rew_function })) }
          if (length(run_condition)==0L) { run_condition <<- do.call(c, lapply(clock_data$runs, function(r) { r$run_condition })) }
          if (length(ev)==0L) { ev <<- do.call(rbind, lapply(clock_data$runs, function(r) { r$w$V })) } #expected value
          if (length(rpe)==0L) { rpe <<- Reward - ev } #better or worse than expected?
            
          #get a list prediction contributions of each parameter per run: each element is a params x trials matrix
          if (length(pred_contrib)==0L) { 
            pred_contrib <<- lapply(clock_data$runs, function(r) { 
                  if (!is.null(r$w$pred_contrib)) { do.call(rbind, r$w$pred_contrib) } else { NULL } 
                }) 
          }
          #flatten this? 
          #arr_pred_contrib <- do.call(abind, list(along=0, pred_contrib))
            
          if (exists("betaFastSlow", envir=clock_data$runs[[1L]]$w, inherits=FALSE)) {
            hasBeta <- TRUE
            bfs_var_fast <<- do.call(rbind, lapply(clock_data$runs, function(r) { r$w$betaFastSlow$var_fast })) #trialwise estimate of uncertainty re: fast responses
            bfs_var_slow <<- do.call(rbind, lapply(clock_data$runs, function(r) { r$w$betaFastSlow$var_slow })) #trialwise estimate of uncertainty re: slow responses
            bfs_mean_fast <<- do.call(rbind, lapply(clock_data$runs, function(r) { r$w$betaFastSlow$mean_fast })) #trialwise estimate of mean value for fast responses
            bfs_mean_slow <<- do.call(rbind, lapply(clock_data$runs, function(r) { r$w$betaFastSlow$mean_slow })) #trialwise estimate of mean value for slow responses
          } else { hasBeta <- FALSE }
        },
        build_design_matrix=function(
            regressors=NULL,
            event_onsets=NULL,
            durations=NULL,
            baselineCoefOrder=-1L,
            runVolumes=NULL, #vector of total fMRI volumes for each run (used for convolved regressors)
            plot=TRUE,
            writeTimingFiles=NULL,
            output_directory="run_timing") {
          
          if (is.null(regressors)) {
            stop("regressor names, event onsets, and event durations must be specified.")
          }
          
          if (length(unique(sapply(list(regressors, event_onsets, durations), length))) != 1L) { 
            stop("regressors, event_onsets, and durations must have the same length") 
          }
          
          #require(fmri)
          
          if (is.null(runVolumes)) {
            #determine the last fMRI volume to be analyzed
            last_fmri_volume <- apply(.self$iti_onset, 1, function(itis) {
                  ceiling(itis[length(itis)] + 12.0 ) #fixed 12-second ITI after every run
                })
            message("Assuming that last fMRI volume was 12 seconds after the onset of the last ITI.")
            message(paste0("Resulting lengths: ", paste(last_fmri_volume, collapse=", ")))
          } else {
            if (length(runVolumes) != nrow(RTobs)) { stop("Length of runVolumes is ", length(runVolumes), ", but number of runs in fit object is ", nrow(RTobs)) }
            last_fmri_volume <- runVolumes
          }
          
          #build design matrix in the order specified
          dmat <- sapply(1:length(regressors), function(r) {
                #Note: all of the regressor computations here should result in nruns x ntrials matrices
                #regressor values
                if (regressors[r] == "rel_uncertainty") {
                  #trialwise relative uncertainty about fast versus slow responses: subtract variances 
                  reg <- abs(bfs_var_fast - bfs_var_slow)
                } else if (regressors[r] == "mean_uncertainty") {
                  reg <- abs(bfs_var_fast + bfs_var_slow)/2 #trialwise average uncertainty for fast and slow responses
                } else if (regressors[r] == "ev") {
                  reg <- ev #expected value
                } else if (regressors[r] == "rpe_pos") {
                  reg <- apply(rpe, c(1,2), function(x) { if (x > 0) x else 0 }) #positive prediction error
                } else if (regressors[r] == "rpe_neg") {
                  reg <- apply(rpe, c(1,2), function(x) { if (x < 0) abs(x) else 0 }) #abs so that greater activation scales with response to negative PE
                } else if (regressors[r] == "rt") {
                  reg <- RTobs #parametric regressor for overall reaction time (from Badre)
                } else {
                  message("Assuming that ", regressors[r], " is a task indicator function.")
                  reg <- array(1.0, dim=dim(RTobs)) #standard task indicator regressor
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
                
                if (durations[r] %in% c("rt", "clock_duration")) { #time of clock on screen
                  durmat <- RTobs/1000
                } else if (durations[r] == "feedback_duration") {
                  durmat <- iti_onset - feedback_onset
                } else if (durations[r] == "iti_duration") {
                  #need to lag clock matrix: iti duration is: clock_onset(t+1) - iti_onset(t)
                  #lag_clock <- apply(clock_onset)
                } else if (!is.na(suppressWarnings(as.numeric(durations[r])))) {
                  #user-specified scalar duration (not currently supporting a user-specified per-regressor vector of durations) 
                  durmat <- array(as.numeric(durations[r]), dim=dim(RTobs))
                } else {
                  stop("Unknown duration keyword:", durations[r])
                }
                
                #fsl-style 3-column format: onset duration value
                output <- lapply(1:nrow(reg), function(run) {
                      cbind(onset=times[run,], duration=durmat[run,], value=reg[run,])
                    })
                
                names(output) <- paste0("run", 1:nrow(reg))
                
                return(output)
              })
          
          dimnames(dmat)[[2L]] <- regressors
          
          
          #returns a 2-d list of runs x regressors. Needs to stay as list since runs vary in length, so aggregate is not rectangular
          #each element in the 2-d list is a 2-d matrix: trials x (onset, duration, value) 
          
          #create an HRF-convolved version of the list
          dmat.convolve <- lapply(1:dim(dmat)[1L], function(run) {
                run.convolve <- lapply(dmat[run,], function(reg) {
                      fmri.stimulus(scans=last_fmri_volume[run], values=reg[,"value"], times=reg[,"onset"], durations=reg[,"duration"], rt=1.0) #hard-coded 1.0s TR for now      
                    })
                do.call(data.frame, run.convolve) #pull into a data.frame with ntrials rows and nregressors cols (convolved)
              })
          
          #dmat.convolve should now be a 1-d runs list where each element is a data.frame of convolved regressors.
          names(dmat.convolve) <- paste0("run", 1:length(dmat.convolve))
          
          #Write timing files to disk for analysis by AFNI, FSL, etc.
          if (!is.null(writeTimingFiles)) {
            dir.create(output_directory, recursive=TRUE, showWarnings=FALSE)
            if ("FSL" %in% writeTimingFiles) {
              for (run in 1:dim(dmat)[1L]) {
                for (reg in 1:dim(dmat)[2L]) {
                  fname <- paste0("run", run, "_", dimnames(dmat)[[2L]][reg], "_FSL3col.txt")
                  write.table(dmat[[run,reg]], file=file.path(output_directory, fname), sep="\t", eol="\n", col.names=FALSE, row.names=FALSE)
                }
              }
            }

            if ("AFNI" %in% writeTimingFiles) {
              #use dmBLOCK-style regressors: time*modulation:duration
              for (reg in 1:dim(dmat)[2L]) {
                regMat <- c()
                for (run in 1:dim(dmat)[1L]) {
                  dmStr <- apply(dmat[[run,reg]], 1, function(row) {
                        row <- plyr::round_any(row, .000001) #keep precision reasonable for output format 
                        paste0(row["onset"], "*", row["value"], ":", row["duration"])
                      })
                  regMat <- rbind(regMat, dmStr)                  
                }
                fname <- paste0(dimnames(dmat)[[2L]][reg], "_dmBLOCK.txt")
                write.table(regMat, file=file.path(output_directory, fname), sep="\t", eol="\n", quote=FALSE, col.names=FALSE, row.names=FALSE)
              }
              
              #write convolved regressors
              #AFNI amplitude modulation forces a mean and deviation from the mean regressor for each effect
              #as a result, if two parametric influences occur at a given time, it leads to perfect collinearity.
              conv_concat <- list()
              lapply(1:length(dmat.convolve), function(r) {
                    lapply(1:length(dmat.convolve[[r]]), function(v) {
                          regName <- names(dmat.convolve[[r]])[v]
                          fname <- paste0(names(dmat.convolve)[r], "_", regName, ".1D")
                          toWrite <- plyr::round_any(dmat.convolve[[r]][[v]], .000001)
                          conv_concat[[regName]] <<- c(conv_concat[[regName]], toWrite) #add for concatenated 1D file
                          write.table(toWrite, file=file.path(output_directory, fname), sep="\n", eol="\n", quote=FALSE, col.names=FALSE, row.names=FALSE)
                        })              
                  })
              
              #write run-concatenated convolved regressors (for use in AFNI)
              lapply(1:length(conv_concat), function(v) {
                    fname <- paste0(names(conv_concat)[v], "_concat.1D")
                    write.table(conv_concat[[v]], file=file.path(output_directory, fname), sep="\n", eol="\n", quote=FALSE, col.names=FALSE, row.names=FALSE)
                  })
            }
          }
          
          collinearityDiag.raw <- apply(dmat, 1, function(run) {
                #check correlations among regressors for trial-wise estimates
                cmat <- do.call(data.frame, lapply(run, function(regressor) {
                          regressor[,"value"]
                        }))
                
                corvals <- cor(cmat, use="pairwise.complete.obs")
                vifMat <- data.frame(cbind(const=rep(1,nrow(cmat)), cmat)) #add dummy constant for vif
                vifForm <- as.formula(paste("const ~ 1 +", paste(names(cmat), collapse=" + ")))
                
                varInfl <- tryCatch(car::vif(lm(vifForm, data=vifMat)), error=function(e) { NA }) #return NA if failure
                list(r=corvals, vif=varInfl)
              })
          
          collinearityDiag.convolve <- lapply(dmat.convolve, function(run) { #apply(dmat.convolve, 1, function(run) {                    
                corvals <- cor(run, use="pairwise.complete.obs")
                vifMat <- data.frame(cbind(const=rep(1,nrow(run)), run)) #add dummy constant for vif 
                vifForm <- as.formula(paste("const ~ 1 +", paste(names(run), collapse=" + ")))
                
                varInfl <- tryCatch(car::vif(lm(vifForm, data=vifMat)), error=function(e) { NA }) #return NA if failure
                list(r=corvals, vif=varInfl)
              })
          
          #add baseline terms to convolved design matrices
          if (baselineCoefOrder > -1L) {
            dmat.convolve <- lapply(dmat.convolve, function(r) {
                  n <- names(r)
                  d <- data.frame(fmri.design(r, order=2))
                  names(d) <- c(n, paste0("base", 0:baselineCoefOrder))
                  d
                })
            
          }
          
          return(list(design=dmat, design.convolve=dmat.convolve, collin.raw=collinearityDiag.raw, collin.convolve=collinearityDiag.convolve))
          
        },
        plotRTs=function() {
          rtDf <- data.frame(
              run=rep(1:nrow(RTobs), each=length(RTobs)),
              run_condition=rep(run_condition, each=length(RTobs)), 
              rew_function=rep(rew_function, each=length(RTobs)),
              reward=Reward,
              rt=c(RTobs, RTpred),
              rt_type=rep(c("observed", "predicted"), each=length(RTobs))
          )
          
          rtPlot <- ggplot(rtDf, aes(x=trial, y=rt, color=rt_type))
          
          rtPlot <- rtPlot + geom_line(size=1.1) + ggtitle(paste(rtDf$rew_function, rtDf$run_condition, sep=", ")) +
              theme_bw(base_size=14)
          
          print(rtPlot)
          return(invisible(rtPlot))
          
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
        RTobs="numeric", #vector of RTs to be fit (raw, differenced, smoothed, etc.)
        RTraw="numeric", #vector of raw RTs
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
        initialize=function(run_number=NA_integer_, RTraw=NA_integer_, Reward=NA_integer_, global_trial_number=NA_integer_,
            rew_function=NULL, run_condition=NA_character_, ...) {
          
          if (is.na(run_number[1L])) { stop("At this point, clockdata_run must have a run number indicating temporal order.") }
          if (is.na(RTraw[1L])) { stop("construction of clockdata_run requires observed reaction times (RTraw)") }
          if (is.na(Reward[1L])) { stop("construction of clockdata_run requires observed rewards (Reward)") }
          
          run_number <<- run_number
          RTraw <<- RTraw
          avg_RT <<- mean(RTraw, na.rm=TRUE)
          Reward <<- Reward
          
          if (is.na(global_trial_number[1L])) { 
            warning("global_trial_number not provided. Defaulting to 1..t")
            global_trial_number <<- 1:length(RTraw)
          } else {
            global_trial_number <<- global_trial_number
          }
          
          if (length(unique(sapply(list(RTraw, Reward, global_trial_number), length))) > 1L) {
            stop("RTraw, Reward, and global_trial_number must have the same length")
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

