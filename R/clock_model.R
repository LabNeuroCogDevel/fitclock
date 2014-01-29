#' Object-oriented implementation of the Frank clock algorithm.
#' 
#' The clock_model class defines a predictive model to fit RT and reward data using data from the clock
#' task. The algorithm can fit group-level data (optimize theta over all subjects and runs), subject-level 
#' multi-run data (optimize theta over multiple runs within subject), or single-run data. 
#' 
#' At initialization, data to be fit can be passed using the \code{clock_data} parameter, and should be of class 
#' clockdata_group (for group data), clockdata_subject (for multi-run subject data), or clockdata_run (for single-run data).
#' 
#' After initializing the object, one or more parameter objects are added to the prediction equation
#' using the $add_params(...) method. Bounds on the parameters can be specified at initialization.
#'  
#' @section Fields:
#'    \describe{
#'      \item{\code{params}:}{ \code{list} of param objects used to predict RTs. Note: parameters are evaluated in order. }
#'      \item{\code{noiseWt}:}{ A numeric scalar (in milliseconds) specifying how much random uniform noise to add to predicted RTs. Default: 0 }
#'      \item{\code{SSE}:}{ Total sum of squared errors for fitted theta values. }
#'      \item{\code{AIC}:}{ Akaike Information Criterion for fitted theta values. }
#'      \item{\code{clock_data}:}{ Clock data to be fit. Should be of class clockdata_group, clockdata_subject, or clockdata_run. DO NOT ASSIGN YOURSELF. }
#'      \item{\code{use_global_avg_RT}:}{ Whether to use the global reaction time across all runs within subject as reference for scale and K. }
#'      \item{\code{global_avg_RT}:}{ Average reaction time across all runs within subject. }
#'      \item{\code{all_by}:}{ #character vector that is the union of all run-level fields that define vary-by-run parameters }
#'    }
#' 
#' @section Methods:
#'    \describe{
#'      \item{\code{add_params(...)}:}{ Method to add parameter objects to prediction equation. Should be a list of params or a single param. }
#'      \item{\code{fit(toFit=NULL, random_starts=NULL, profile=TRUE)}:}{ Fit \code{clock_data} using the list of \code{params}. Data to fit can be set at call to \code{$fit()} using the \code{toFit} parameter. } 
#'      \item{\code{list_params()}:}{ Returns a data.frame describing the parameters, current values, and bounds of this model. }
#'      \item{\code{params_current()}:}{ Named vector of current values for all parameters. }
#'      \item{\code{params_minimum()}:}{ Named vector of lower bounds for all parameters. }
#'      \item{\code{params_maximum()}:}{ Named vector of upper bounds for all parameters. }
#'      \item{\code{params_initial()}:}{ Named vector of initial values for all parameters. }
#'      \item{\code{params_par_scale()}:}{ Named vector of expected log parameter scale for all parameters. }
#'      \item{\code{reorder_params(neworder)}:}{ Re-orders the parameter list according to a named character vector, \code{neworder}. Useful because params evaluated in order. }
#'      \item{\code{set_data(cdata)}:}{ Sets the clock dataset to be fit. cdata should be a clockdata object. Use this method if you wish to alter the clock_data field. }
#' 
#'    }
#'
#' @importFrom methods setRefClass
#' @export

#add a reset method here to set values back to initial (clear for re-using the same alg object with a new dataset)
#actually, couldn't we just note when the dataset is updated and zero out parameters?
clock_model <- setRefClass(
    Class="clock_model",
    fields=list(
        #updateEquation="expression",
        params="list",
        noiseWt="numeric",
        SSE="numeric",
        AIC="numeric",
        clock_data="ANY", #allow this to be dataset, subject, or run,
        use_global_avg_RT="logical",
        global_avg_RT="numeric",
        all_by="character" #character vector that is the union of all run-level fields defining vary-by-run parameters
    ),
    methods=list(
        initialize=function(clock_data=NULL, ...) {
          cat("Initializing clock_model\n")
          
          params <<- list() #initialize empty list of model parameters
          noiseWt <<- 0 #do not add noise to RT prediction
          use_global_avg_RT <<- TRUE #whether to use average RT across all blocks in fit (e.g., "go for gold" scales wrt avg_RT). 
          
          if (!is.null(clock_data)) { set_data(clock_data) }
          callSuper(...) #for classes inheriting from this, pass through unmatched iniitalization values
        },

        #Compute the global reaction time, used by some parameters to scale relative differences 
        #  For subject-level data, global RT is average across all runs and trials
        #  For run-level data, global RT is average across trials.
        #  TODO: Figure out what to do in case of group data (not currently handled). 
        set_global_avg_RT=function() {
          if (use_global_avg_RT && !is.null(clock_data)) {
            if (class(clock_data) == "clockdata_subject") {
              #compute average reaction time across all runs
              global_avg_RT <<- mean(do.call("c", lapply(clock_data$runs, "[[", "RTobs")), na.rm=TRUE)
            } else if (class(clock_data) == "clockdata_run") {
              global_avg_RT <<- mean(clock_data$RTobs, na.rm=TRUE)
            }
          } 
        },
        
        #function to append elements from one list to another, warning about overwriting existing elements
        lappend=function(target, append) {
          if (!is.list(target)) stop("target is not a list.")
          if (!is.list(append)) stop("append is not a list.")
          
          #force names of list to match corresponding parameter names
          #this discards any list names provided by user
          names(append) <- sapply(append, function(p) { return(p$name) })
          
          for (elementName in names(append)) {
            if (!is.null(target[[elementName]])) warning("Element is already present in target list: ", elementName)
            target[[elementName]] <- append[[elementName]]
          }
          
          return(target)
        },
        
        #method to add one or more parameter objects to the model. 
        add_params=function(...) {
          #accept a variable number of parameter objects
          p_objs <- as.list(match.call())[-1L] #first element of match.call is a class for refMethodDef for the clock_model object
          p_objs <- lapply(p_objs, eval) #force initialization of param objects (otherwise stays as a language expression)
          
          #pass along address of shared workspace to sub-objects
          #now diverting workspace issues until just before fit
          #allow clock_model object to be completely setup independent of data
          #p_objs <- lapply(p_objs, function(p) { p$w <- .self$w; p$reset_workspace(); return(p) })
          
          params <<- lappend(params, p_objs)
        },
        
        #Set the dataset to be fit. Important to expose method because after setting clock_data field,
        #we need to compute the global RT and setup any vary-by-run parameters given the details of the data to be fit.
        set_data=function(cdata) {
          if (!class(cdata) %in% c("clockdata_group", "clockdata_subject", "clockdata_run")) {
            stop("data for clock_model object must be one of the following types: clockdata_group, clockdata_subject, or clockdata_run")
          }
          
          clock_data <<- cdata
          set_global_avg_RT()
          
          setup_param_by()
        },
        
        #expand parameters to accommodate per-condition variation, depending on clock_data
        setup_param_by=function() {
          if (inherits(clock_data, "uninitializedField")) { return(invisible(NULL)) } #do nothing
          
          #determine by settings for all parameters
          if (!is.null( uniq_by <- unique(unlist(lapply(params, function(p) { p$by }))) )) { 
            all_by <<- uniq_by 
          } else {
            return(invisible(NULL)) #no use of by-condition parameters
          }
          
          #develop a lookup matrix of unique conditions combinations
          if (class(clock_data) == "clockdata_subject") {
            conMat <- lapply(clock_data$runs, function(r) {
                  #set by_lookup (union of all by fields) for each run object
                  #and return a list of conditions for each run to be expanded for each parameter below, as needed
                  r$by_lookup <- sapply(all_by, function(b) { r[[b]]} )
                })
            conMat_df <- data.frame(do.call(rbind, conMat))
            for (n in 1:ncol(conMat_df)) {
              conMat_df[,n] <- paste(names(conMat_df)[n], conMat_df[,n], sep=":")
            }
            params <<- lapply(params, function(p) {
                  p$base_name <- p$name[1L] #just basic parameter name without condition (used for getTheta lookup) 
                  if (!is.null(p$by)) {
                    #make sure that by is sorted alpha to correspond with getTheta
                    p$by <- sort(p$by)
                    #format parameter names into <param>/condition:value/condition:value/.../ for passing to optimizer (parse within predict)
                    #for now, leave off param name, just keep condition names
                    #when calling params_initial, etc., lapply will include param name in vector
                    p$name <- paste(p$name, unique(apply(conMat_df[,p$by, drop=FALSE], 1, paste, collapse="/")), sep="/") #alpha sort is important so that getTheta can assume that order of theta condition string is alphabetical
                    #p$name <- unique(apply(conMat_df[,p$by, drop=FALSE], 1, paste, collapse="/"))
                    #I suppose using the first element as the value to replicate could backfire if an clock_model object is recycled across datasets.
                    p$init_value <- rep_len(p$init_value[1L], length(p$name))
                    p$cur_value <- rep_len(p$cur_value[1L], length(p$name))
                    p$max_value <- rep_len(p$max_value[1L], length(p$name))
                    p$min_value <- rep_len(p$min_value[1L], length(p$name))
                    p$par_scale <- rep_len(p$par_scale[1L], length(p$name))
                    names(p$init_value) <- names(p$cur_value) <- names(p$max_value) <- names(p$min_value) <- names(p$par_scale) <- p$name
                  }
                  p
                })
          }
          
        },
        
        #obtain a data.frame describing the parameters, values, and bounds of this model.
        list_params=function() {
          #cat("Current parameters in algorithm\n\n")
          df <- cbind(min_value=params_minimum(),
              max_value=params_maximum(),
              init_value=params_initial(),
              cur_value=params_current(),
              par_scale=params_par_scale()
          )
          rownames(df) <- unlist(lapply(params, function(p) { p$name }))
          #cat(names(params), sep="\n")
          #print(format(df,digits=3,scientific=FALSE))
          return(df)
        },
        
        reorder_params=function(neworder) {
          #reorders the parameter collection
          #potentially useful if there are dependencies among parameters or workspace values: params are evaluated in sequential order 
          left_out <- names(params)[! names(params) %in% neworder]
          if (length(left_out) > 0L) {
            warning("reorder_params call did not include position for: ", paste(left_out, collapse=", "), ". Defaulting to end")
            neworder <- c(neworder, left_out)
          }
          if (length(unmatched <- neworder[! neworder %in% names(params)]) > 0L) {
            stop("reorder_params unable to find parameters in clock_model: ", paste(unmatched, collapse=", "))
          }
          params <<- params[neworder]
        },
        
        #these functions get values, bounds, and scales for all parameters in the current model 
        params_current=function() {
          unlist(unname(lapply(params, function(p) { p$cur_value }))) #need to unname lapply vector to avoid appending names(params) to each element
        },        
        params_minimum=function() {
          unlist(unname(lapply(params, function(p) { p$min_value })))          					
        },        
        params_maximum=function() {
          unlist(unname(lapply(params, function(p) { p$max_value })))
        },        
        params_initial=function() {
          unlist(unname(lapply(params, function(p) { p$init_value })))
        },        
        params_par_scale=function() {
          unlist(unname(lapply(params, function(p) { p$par_scale })))
        },
        
        fit=function(toFit=NULL, random_starts=NULL, profile=TRUE) {
          
          #can pass in new dataset at $fit call
          #if not passed in, use the current clock_data field
          if (!is.null(toFit)) {
            set_data(toFit)
          }
          
          initialValues <- params_initial()
          lower <- params_minimum()
          upper <- params_maximum()
          
          if (inherits(clock_data, "uninitializedField")) { stop("clock dataset must be provided prior to $fit(). Set using $set_data()") }
          
          #require(optimx)
          
          #call predict using optimizer
          if (length(initialValues) == 1L) { method="Brent"
          } else { method="L-BFGS-B" }
          
          #lapply(params, function(p) { p$value_history <- numeric(0) }) #reset value history for all parameters prior to fit REFCLASS VERSION
          params <<- lapply(params, function(p) { p$value_history <- numeric(0); return(p) }) #reset value history for all parameters prior to fit S3 VERSION
          
          #initialValues <- initialValues / params_par_scale() #scale prior to optimization
          
          #system.time(optResult <- optim(par=initialValues, fn=.self$predict, method=method,
          #        lower=lower, upper=upper,
          #        control=list(parscale=params_par_scale()),
          #        updateFields=FALSE, trackHistory=TRUE))
          
          if (!is.null(random_starts) && is.numeric(random_starts)) {
            require(foreach)
            require(doMC)
            njobs <- min(parallel::detectCores, random_starts)
            registerDoMC(njobs)
            
          }
          
          #track optimization
          if (profile) { Rprof("tc_prof.out") }
          
          #this is the most sensible, and corresponds to optim above (and is somewhat faster)
          elapsed_time <- system.time(optResult <- nlminb(start=initialValues, objective=.self$predict, 
                  lower=lower, upper=upper, scale=1/params_par_scale(),
                  updateFields=FALSE, trackHistory=TRUE))
          
          if (profile) {
            Rprof(NULL) #stop profiling
            
            prof <- summaryRprof("tc_prof.out")#, lines = "both")
            unlink("tc_prof.out") #remove profile file
          } else {
            prof <- list()
          }
          
          #this is roughly what is suggested in the PORT documentation (1/max scaling). But it's slower than 1/magnitude above.
          #system.time(optResult <- nlminb(start=initialValues, objective=.self$predict, lower=lower, upper=upper, scale=1/upper, #scale=c(1, 100, 100, 100, 100),
          #        updateFields=FALSE, trackHistory=TRUE))
          
          #even though this is faster, I have no idea why one would use a scale of 100 for all...
          #system.time(optResult <- nlminb(start=initialValues, objective=.self$predict, lower=lower, upper=upper, scale=100L, #scale=c(1, 100, 100, 100, 100),
          #        updateFields=FALSE, trackHistory=TRUE))
          
          #system.time(optResult <- nlminb(start=initialValues, objective=.self$predict, lower=lower, upper=upper, scale=10L, #scale=c(1, 100, 100, 100, 100),
          #        updateFields=FALSE, trackHistory=TRUE))
          
          #system.time(optResult <- nlminb(start=initialValues, objective=.self$predict, lower=lower, upper=upper,
          #        updateFields=FALSE, trackHistory=TRUE))
          
          #scalecheck()
          #optimx does not appear to be suited to pass through scale to nlminb...
          #elapsed_time <- system.time(optResult <- optimx(par=initialValues, fn=.self$predict, method="nlminb",
          #        lower=lower, upper=upper, 
          #        control=list(
          #            scale=1/params_par_scale(), #parameter scaling since values are quite different in magnitude 
          #            usenumDeriv=FALSE, 
          #            all.methods=FALSE
          #            ),
          #        updateFields=FALSE, trackHistory=TRUE))
          
          #time <- system.time(p <- optimx(par=initialValues, fn=.self$predict, lower=lower, upper=upper, method=method, control=list(parscale=c(1e3, 1e0, 1e-1, 1e-1))))
          #time <- system.time(p <- optimx(par=initialValues, fn=.self$predictNoLookup, lower=lower, upper=upper, method=method, updateFields=FALSE))
          
          if (optResult$convergence == 0) {
            #success
            #lapply(params, function(p) { p$cur_value <- optResult$par[p$name] }) #set current parameter values based on optimization REFCLASS
            #need to unname the parameter value to avoid getting a named vector, which lapply concatenates with list names (e.g., lambda.lambda)
            params <<- lapply(params, function(p) { p$cur_value <- unname(optResult$par[p$name]); return(p) }) #set current parameter values based on optimization S3
            
            #compute trialwise predictions with optimized parameters
            .self$predict(updateFields=TRUE, trackHistory=FALSE)
            
            #set SSE for fit
            SSE <<- optResult$objective
            
            ##TODO: Would be nice if $fit returned a fit object that had parameter values, the constituent environments that were fit,
            # and key workspace variables as data.frames (e.g., predicted RTs, observed RTs, trial-by-trial param contributions, etc.
            
            ##only works for subject and run fits, not group
            if (class(clock_data) == "clockdata_subject") {
              ntrials <- sum(unlist(lapply(clock_data$runs, function(r) { r$w$ntrials } )))
              RTobs <- do.call(rbind, lapply(clock_data$runs, function(r) { r$RTobs }))
              RTpred <- do.call(rbind, lapply(clock_data$runs, function(r) { r$w$RTpred }))
              Reward <- do.call(rbind, lapply(clock_data$runs, function(r) { r$Reward }))
            } else if (class(clock_data) == "clockdata_run") {
              ntrials <- clock_data$w$ntrials
              RTobs <- clock_data$RTobs
              RTpred <- clock_data$w$RTpred
              Reward <- clock_data$Reward
            }
            
            nparams <- length(optResult$par) #number of free parameters
            
            AIC <<- ntrials*(log(2*pi*(SSE/ntrials))+1) + 2*nparams
            
            fitResult <- clock_fit(RTobs=RTobs, RTpred=RTpred, Reward=Reward, total_SSE=SSE, AIC=AIC, 
                elapsed_time=unclass(elapsed_time), opt_data=optResult, profile_data=prof, theta=as.matrix(list_params()))
          } else {
            warning("Optimization failed.")
            fitResult <- NULL
          }
          
          return(fitResult)
          
        },
        
        predict=function(theta=params_current(), updateFields=FALSE, trackHistory=FALSE) {
          ##TODO: clock_model object shoud allow for some sort of symbolic specification of how
          #values from prior runs are carried forward.
          #this is an clock_model-level decision, not subject, run, parameter, etc.
          
          totalSSE <- 0
          
          #track value of parameters over the course of optimization
          if (trackHistory) {
            params <<- lapply(params, function(p) { p$value_history <- c(p$value_history, theta[p$name]); return(p) }) #set current parameter values based on optimization S3 version
            #lapply(params, function(p) { p$value_history <- c(p$value_history, theta[p$name]) }) #set current parameter values based on optimization REFCLASS VERSION
          }
          
          if (class(clock_data)=="clockdata_group") {
            #loop over subjects and runs (group fit)
            for (s in clock_data$subjects) { #make this parallel
              for (r in s$runs) {
                
              }
            }
          } else if (class(clock_data)=="clockdata_subject") {            
            for (r in 1:length(clock_data$runs)) {
              if (r > 1L) {
                prior_w <- clock_data$runs[[r - 1]]$w #pass along environment from previous run
              } else { prior_w <- NULL }
              
              #pass forward theta (likely from optim)
              run_SSE <- fit_run(theta=theta, to_fit=clock_data$runs[[r]], prior_w=prior_w, updateFields=updateFields)
              clock_data$runs[[r]]$SSE <- run_SSE #copy run
              totalSSE <- totalSSE + run_SSE
            }
            
          } else if (class(clock_data)=="clockdata_run") {
            #pass forward theta (likely from optim)
            totalSSE <- fit_run(theta=theta, to_fit=clock_data)
          }
          
          return(totalSSE)
        },
        
        fit_run=function(to_fit=NULL, prior_w=NULL, theta=params_current(), updateFields=FALSE) {
          if (is.null(to_fit) || !class(to_fit) == "clockdata_run") { stop("fit_run requires a clockdata_run object") }
          
          #setup basics of workspace -- shouldn't this be in clock_model conceptually?
          to_fit$initialize_workspace(prior_w)
          
          #copy environment of run to local variable for looping over trials
          w <- to_fit$w
          
          #if using global RT average, set environment here
          #N.B. Need to set here before calling reset in each param because this is where bestRT[1] is set...
          #would be better to make this more robust
          if (use_global_avg_RT) { w$avg_RT <- global_avg_RT }
          
          #lapply(.self$params, function(p) { p$w <- .self$w; p$reset_workspace() }) #setup workspace variables for each parameter REFCLASS VERSION
          #For speed, theta obj + condition lookup should occur here, outside of trial loop
          params <<- lapply(.self$params, function(p) { 
                p$w <- w
                reset_workspace(p)
                p$theta_lookup <- getTheta(p, to_fit$by_lookup) #define which element of theta is relevant for each parameter (based on $by and this run's condition)
                return(p)
              }) #setup workspace variables for each parameter (NEED TO UPDATE PARAMS SINCE THEY ARE PASS BY VALUE)
          
          #optim using Brent, or optimize fail to pass names of parameters into predict.
          #names are critical for the params to lookup properly.
          if (is.null(names(theta))) {
            #warning("Unnamed theta vector passed to predict. Looking up from parameter list.")
            names(theta) <- names(params) #sapply(params, function(p) { return(p$name) })
          }
          
          #reset workspace before proceeding!
          
          w$RT_new  <- w$RTobs[1L] #use observed RT as predicted RT for t=1
          w$RT_last <- w$RTobs[1L]
          w$RT_last2 <- w$RTobs[1L]
          
          ##TODO: Risk that we initialize some vals above (like initial V), but predict is called iteratively for optimization, so values will not be reset here if they were set in initialize
          
          for (t in 2:w$ntrials) {
            w$cur_trial <- t
            w$lastTrial <- t - 1
            
            eval(
                quote({
                      Rew_last <- Reward[lastTrial]
                      RT_last  <- RTobs[lastTrial]
                      if (cur_trial > 2) { RT_last2 <- RTobs[cur_trial - 2] }
                      
                      V_last <- V[lastTrial]
                      V_new = V_last + alphaV*(Rew_last - V_last) # update critic expected value
                      V[cur_trial] <- V_new
                    }), 
                w)
            
            #refclass version
            #w$RT_new <<- sum(sapply(params, function(p) { p$getRTUpdate(theta, updateFields=updateFields) } )) + noiseWt*(runif(1,-0.5,0.5)) #add or subtract noise according to noiseWt (0 for now)
            #w$RT_new <- sum(sapply(params, function(p) { getRTUpdate(p, theta, updateFields=updateFields) } )) + noiseWt*(runif(1,-0.5,0.5)) #add or subtract noise according to noiseWt (0 for now)
            
            #need to pass along run condition to getRTUpdate functions
            #sapply(all_by, function(b) { to_fit[[b]]} )
            
            noiseContrib <- if (noiseWt > 0.0) { noiseWt*(runif(1,-0.5,0.5)) } else { 0.0 }
            w$RT_new <- sum(pred_contrib <- unlist(lapply(params, function(p) { getRTUpdate(p, theta) } ))) + noiseContrib
            
            if (updateFields) { #track trial-by-trial contribution of each param to RT prediction
              pnames <- names(params)
              for ( p in 1:length(pnames) ) {
                w$pred_contrib[[ pnames[p] ]][w$cur_trial] <- pred_contrib[p]
              }
            }
            
            w$RTpred[w$cur_trial] <- w$RT_new
            
            #TODO: incorporate this code from matlab. Affects beta dist local_RT and autocorrelation parameters 
            #if RT_last==0, RT_last = RT_last2; end; %% if last trial there
            #was no response, use trial before that for updating RT avg and autocorrelation effects (otherwise counted as 0)
            
            w$rpe[t-1] <- w$Rew_last - w$V_last
            
          }
          
          #cat("SSE: ", SSE, "\n")
          return(SSE)
        })

)
