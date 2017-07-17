#dataset object is supposed to represent a group of data over which a single set of parameters (theta) are to be found
#so in the group case, get best-fitting parameters for all subjects and runs
#for a single subject, best fitting-parameters across subject's runs
#for a single run, best-fitting parameters for this run

#have predict in alg detect the class of the dataset to determine how to proceed

#abstract build_design_matrix function from clock_fit object to allow for more rapid tweaks without re-fitting data.
#' Create fmri design matrix for AFNI, FSL, or internal convolved design.
#'
#' @importFrom data.table as.data.table
#' @importFrom plyr round_any 
#' @importFrom orthopolynom legendre.polynomials polynomial.values
#' @export
build_design_matrix=function(
    fitobj=NULL,
    regressors=NULL,
    event_onsets=NULL,
    durations=NULL,
    normalizations=NULL, #normalization of HRF
    center_values=FALSE, #whether to center parametric regressors prior to convolution
    baselineCoefOrder=-1L,
    baselineParameterization="Legendre",
    runVolumes=NULL, #vector of total fMRI volumes for each run (used for convolved regressors)
    dropVolumes=0L, #vector of how many volumes to drop from the beginning of a given run
    runsToOutput=NULL,
    plot=TRUE,
    writeTimingFiles=NULL,
    output_directory="run_timing",
    tr=1.0, #TR of scan in seconds
    convolve=TRUE, #whether to convolve the regressors with the HRF
    convolve_wi_run=TRUE, #whether to mean center parametric regressors within runs before convolution
    add_derivs=FALSE, #whether to add temporal derivatives for each column after convolution
    parmax1=FALSE, #whether to rescale a convolved regressor to max=1.0 after convolution (normalize scaling across runs and subjects)
    high_pass=NULL #whether to apply a high-pass filter to the design matrix (e.g., to match fmri preprocessing)
) {
  
  if (is.null(regressors)) {
    stop("regressor names, event onsets, and event durations must be specified.")
  }
  
  if (length(unique(sapply(list(regressors, event_onsets, durations), length))) != 1L) { 
    stop("regressors, event_onsets, and durations must have the same length") 
  }
  
  if (is.null(normalizations)) {
    normalizations <- rep("none", length(regressors))
  }
  
  if (is.null(runVolumes)) {
    #determine the last fMRI volume to be analyzed
    runVolumes <- sapply(fitobj$iti_onset, function(itis) {
          ceiling((itis[length(itis)] + 12.0)/tr ) #fixed 12-second ITI after every run
        })
    message("Assuming that last fMRI volume was 12 seconds after the onset of the last ITI.")
    message(paste0("Resulting lengths: ", paste(runVolumes, collapse=", ")))
  } else {
    if (length(runVolumes) != length(fitobj$RTraw)) { warning("Length of runVolumes is ", length(runVolumes), ", but number of runs in fit object is ", nrow(fitobj$RTraw)) }
  }
  
  if (length(dropVolumes) < length(runVolumes)) {
    message("Using first element of dropVolumes for all runs: ", dropVolumes[1L])
    dropVolumes <- rep(dropVolumes[1L], length(runVolumes)) 
  }
  
  #determine which run fits should be output for fmri analysis
  if (is.null(runsToOutput)) {
    message("Assuming that all runs should be fit and run numbers are sequential ascending")
    runsToOutput <- 1:length(runVolumes)
  }
  
  timeOffset <- tr*dropVolumes
  runVolumes <- runVolumes - dropVolumes #elementwise subtraction of dropped volumes from full lengths
    
  #build design matrix in the order specified
  #use an explicit cbind of the lapply to create a runs x regressors list
  dmat <- do.call(cbind, lapply(1:length(regressors), function(r) {
            if (grepl("^custom_", regressors[r])) {

              #custom regressor
              #expect that the field in fitobj will be a list of nruns
              #each element is a matrix-like object (can be data.frame) 
              #containing columns 'onset', 'duration', and 'value' fields
              regname <- sub("^custom_", "", regressors[r], perl=TRUE)
              reg <- fitobj$sceptic[[regname]] #yoked to sceptic field for now
              
              #create a run-level concatenated data.frame for each run
              output <-  lapply(1:nrow(reg), function(run) {
                    do.call(rbind, lapply(reg[run,], function(tlist) { as.matrix(tlist[,c("onset", "duration", "value")]) })) #bind all trial-level matrix objects together
                  })
              
              names(output) <- paste0("run", 1:nrow(reg))
              #stopifnot(all(c("onset", "duration", "value") %in% names(output[[1]]))) #check that field names match expectation
            } else {
              #Note: all of the regressor computations here should result in nruns x ntrials matrices
              #regressor values
              if (regressors[r] == "rel_uncertainty") {
                #trialwise relative uncertainty about fast versus slow responses: subtract variances
                reg <- lapply(1:length(fitobj$bfs_var_fast), function(run) {
                      abs(fitobj$bfs_var_fast[[run]] - fitobj$bfs_var_slow[[run]])
                    })
              } else if (regressors[r] == "mean_uncertainty") {
                reg <- lapply(1:length(fitobj$bfs_var_fast), function(run) {
                      abs(fitobj$bfs_var_fast[[run]] + fitobj$bfs_var_slow[[run]])/2
                    })
              } else if (regressors[r] == "ev") {
                reg <- fitobj$ev #expected value
              } else if (regressors[r] == "rpe") {
                reg <- fitobj$rpe #rpe maintaining sign
              } else if (regressors[r] == "rpe_pos") {
                reg <- lapply(fitobj$rpe, function(run) { sapply(run, function(x) { if (is.na(x) || x > 0) x else 0 }) }) #positive prediction error
              } else if (regressors[r] == "rpe_pos_sqrt") {
                reg <- lapply(fitobj$rpe, function(run) { sapply(run, function(x) { if (is.na(x) || x > 0) sqrt(x) else 0 }) }) #positive prediction error
              } else if (regressors[r] == "rpe_pos_bin") {
                reg <- lapply(fitobj$rpe, function(run) { sapply(run, function(x) { if (is.na(x)) NA else if (x > 0) 1 else 0 }) }) #1/0 binarized PPE
              } else if (regressors[r] == "rpe_neg") {
                reg <- lapply(fitobj$rpe, function(run) { sapply(run, function(x) { if (is.na(x) || x < 0) abs(x) else 0 }) }) #abs so that greater activation scales with response to negative PE
              } else if (regressors[r] == "rpe_neg_sqrt") {
                reg <- lapply(fitobj$rpe, function(run) { sapply(run, function(x) { if (is.na(x) || x < 0) sqrt(abs(x)) else 0 }) }) #abs so that greater activation scales with response to negative PE
              } else if (regressors[r] == "rpe_neg_bin") {
                reg <- lapply(fitobj$rpe, function(run) { sapply(run, function(x) { if (is.na(x)) NA else if (x < 0) 1 else 0 }) }) #1/0 binarized NPE
              } else if (regressors[r] == "rt") {
                reg <- fitobj$RTraw #parametric regressor for overall reaction time (from Badre)
              } else if (grepl("^sceptic_", regressors[r])) {
                #all sceptic parameters just get copied across directly
                reg <- fitobj$sceptic[[sub("^sceptic_", "", regressors[r], perl=TRUE)]] #name must match 
              } else {
                message("Assuming that ", regressors[r], " is a task indicator function.")
                reg <- lapply(fitobj$RTraw, function(run) { rep(1.0, length(run)) }) #standard task indicator regressor
              }
              
              #replace missing values with 0 for clarity in convolution
              reg <- lapply(reg, function(run) { run[which(is.na(run))] <- 0; return(run) })
              
              #determine onset times for events
              if (event_onsets[r] == "clock_onset") {
                times <- fitobj$clock_onset
              } else if (event_onsets[r] == "feedback_onset") {
                times <- fitobj$feedback_onset
              } else if (event_onsets[r] == "iti_onset") {
                times <- fitobj$iti_onset
              } else if (event_onsets[r] == "rt") {
                times <- lapply(1:fitobj$nruns, function(run) { fitobj$clock_onset[[run]] + fitobj$RTraw[[run]]/1000 })
              }
              
              if (durations[r] %in% c("rt", "clock_duration")) { #time of clock on screen
                durmat <- lapply(fitobj$RTraw, function(run) { run/1000 })
              } else if (durations[r] == "feedback_duration") {
                durmat <- lapply(1:fitobj$nruns, function(run) { fitobj$iti_onset[[run]] - fitobj$feedback_onset[[run]] })
              } else if (durations[r] == "iti_duration") {
                #need to lag clock matrix: iti duration is: clock_onset(t+1) - iti_onset(t)
                #lag_clock <- apply(clock_onset)
              } else if (!is.na(suppressWarnings(as.numeric(durations[r])))) {
                #user-specified scalar duration (not currently supporting a user-specified per-regressor vector of durations)
                durmat <- lapply(fitobj$RTraw, function(run) { rep(as.numeric(durations[r]), length(run)) })
              } else {
                stop("Unknown duration keyword:", durations[r])
              }
              
              #fsl-style 3-column format: onset duration value
              output <- lapply(1:length(reg), function(run) {
                    cbind(onset=times[[run]] - timeOffset[run], duration=durmat[[run]], value=reg[[run]])
                  })
              
              names(output) <- paste0("run", 1:length(reg))
            }
                        
            return(output)
          }))
  
  dimnames(dmat)[[2L]] <- regressors
  
  #only retain runs to be analyzed
  dmat <- dmat[runsToOutput,,drop=FALSE] 
  
  #returns a 2-d list of runs x regressors. Needs to stay as list since runs vary in length, so aggregate is not rectangular
  #each element in the 2-d list is a 2-d matrix: trials x (onset, duration, value) 
  
  #subfunction used by the summed runs and separate runs convolution steps
  convolve_regressor <- function(reg, vols, normalization="none", high_pass=NULL) {
    #check for the possibility that the onset + event duration exceeds the number of good volumes in the run (e.g., if truncated for high movement) 
    if (any(whichHigh <- (reg[,"onset"] + reg[,"duration"]) > vols)) {
      reg <- reg[!whichHigh,]
    }
    
    #see hrf_convolve_normalize for implementation details
    if (normalization == "evtmax_1.0") {
      #each event is 1.0-normalized
      x <- hrf_convolve_normalize(scans=vols, values=reg[,"value"], times=reg[,"onset"], durations=reg[,"duration"], rt=tr, normeach=TRUE, center_values=center_values, parmax1=parmax1)
    } else if (normalization == "durmax_1.0") {
      #peak amplitude of hrf is 1.0 (before multiplying by parametric value) based on stimulus duration (stims > 10s approach 1.0)
      x <- hrf_convolve_normalize(scans=vols, values=reg[,"value"], times=reg[,"onset"], durations=reg[,"duration"], rt=tr, normeach=FALSE, center_values=center_values, parmax1=parmax1)
    } else {
      #no normalization
      x <- fmri.stimulus(scans=vols, values=reg[,"value"], times=reg[,"onset"], durations=reg[,"duration"], rt=tr, center_values=center_values, convolve=convolve, parmax1=parmax1)
    }
    
    #apply high-pass filter after convolution if requested
    if (!is.null(high_pass)) {
      x <- fir1Bandpass(x, TR=tr, low=high_pass, high=1/tr/2, plotFilter=FALSE, forward_reverse=TRUE, padx=1, detrend=1)
    }
    
    return(x)
    
  }
  
  #concatenate regressors across runs by adding timing from MR files.
  runtiming <- cumsum(runVolumes)*tr #timing in seconds of the start of successive runs
  
  #for visualization of events, return concatenated onsets (some code redundancy below)
  #note that we want to add the run timing from the r-1 run to timing for run r.
  #example: run 1 is 300 volumes with a TR of 2.0
  # thus, run 1 has timing 0s .. 298s, and the first volume of run 2 would be 300s
  concat_onsets <- lapply(1:dim(dmat)[2L], function(reg) {
        thisreg <- dmat[,reg]
        concat_reg <- do.call(c, lapply(1:length(thisreg), function(run) {
                  timing <- thisreg[[run]]
                  timing[,"onset"] <- timing[,"onset"] + ifelse(run > 1, runtiming[run-1], 0)
                  timing[timing[,"value"] != 0, "onset"] #remove exact 0 values from "events"
                }))
        
        concat_reg
      })
  names(concat_onsets) <- dimnames(dmat)[[2L]]
  
  if (convolve_wi_run) {
    #create an HRF-convolved version of the list
    dmat.convolve <- lapply(1:dim(dmat)[1L], function(i) {
          run.convolve <- lapply(1:dim(dmat)[2L], function(j) {
                reg <- dmat[[i,j]] #regressor j for a given run i
                convolve_regressor(reg, runVolumes[i], normalizations[j], high_pass=high_pass)
              })
          
          df <- do.call(data.frame, run.convolve) #pull into a data.frame with ntrials rows and nregressors cols (convolved)
          names(df) <- dimnames(dmat)[[2L]]
          return(df)
        })    
  } else {
    #issue with convolution of each run separately is that the normalization and mean centering are applied within-run
    #in the case of EV, for example, this will always scale the regressor in terms of relative differences in value within run, but
    #will fail to capture relative differences across run (e.g., if value tends to be higher in run 8 than run 1)
    
    dmat_sumruns <- lapply(1:dim(dmat)[2L], function(reg) {
          thisreg <- dmat[,reg]
          concattiming <- do.call(rbind, lapply(1:length(thisreg), function(run) {
                    timing <- thisreg[[run]]
                    timing[,"onset"] <- timing[,"onset"] + ifelse(run > 1, runtiming[run-1], 0)
                    timing
                  }))
          
          #convolve concatenated events with hrf
          all.convolve <- convolve_regressor(concattiming, sum(runVolumes), normalizations[reg], high_pass=high_pass)
          
          #now, to be consistent with code below (and elsewhere), split back into runs
          splitreg <- split(all.convolve, do.call(c, sapply(1:length(runVolumes), function(x) { rep(x, runVolumes[x]) })))
          
          return(splitreg)
          
        })
    
    #now have a list of length(regressors) with length(runs) elements.
    #need to reshape into a list of data.frames where each data.fram is a run with regressors
    dmat.convolve <- lapply(1:length(dmat_sumruns[[1L]]), function(run) {
          df <- data.frame(lapply(dmat_sumruns, "[[", run))
          names(df) <- dimnames(dmat)[[2L]]
          return(df)
        })
    
  }
  
  #handle the addition of temporal dervitives
  if (add_derivs) {
    message("Adding temporal derivatives of each substantive regressor orthogonalized against design matrix")
    
    dmat.convolve <- lapply(dmat.convolve, function(run) {
          dmat <- as.matrix(run) #need as a matrix for lm call
          
          dmat_derivatives <- do.call(cbind, lapply(1:ncol(dmat), function(col) {
                    dx <- c(0, diff(dmat[,col]))
                    
                    ##orthogonalize wrt design
                    return(residuals(lm(dx ~ dmat)))                
                  }))
          
          colnames(dmat_derivatives) <- paste0("d_", colnames(dmat))
          
          cbind(run, dmat_derivatives) #return design matrix with derivatives added
        })
  }
  
  
#  quick code to test that the design matrix regressors are sensible in their heights and (de)correlation with each other after convolution. 
#  In particular, verify that zero is not being treated as an "event" by the parametric convolution after mean centering     
#  browser()
#  dm1 <- dmat[1,]
#  dcon1 <- dmat.convolve[[1]]
#  dcon1$time <- 1:nrow(dcon1)
#  cor(dcon1)
#  m <- melt(dcon1, id.vars="time")
#  ggplot(m, aes(x=time, y=value)) + geom_line() + facet_grid(variable~., scales="free_y")
#  lapply(dcon1, range)
#  lapply(dm1, function(reg) range(reg[,"value"]))
#  vtest <- dm1$rpe_pos[,"value"]
#  vtest <- vtest[vtest != 0.0]
#  vtest - mean(vtest)
  
  #dmat.convolve should now be a 1-d runs list where each element is a data.frame of convolved regressors.
  names(dmat.convolve) <- paste0("run", runsToOutput)
  
  #Write timing files to disk for analysis by AFNI, FSL, etc.
  if (!is.null(writeTimingFiles)) {
    dir.create(output_directory, recursive=TRUE, showWarnings=FALSE)
    
    if ("convolved" %in% writeTimingFiles) {
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
    
    if ("FSL" %in% writeTimingFiles) {
      for (i in 1:dim(dmat)[1L]) {
        for (reg in 1:dim(dmat)[2L]) {
          regout <- dmat[[i,reg]]
          if (center_values && !all(regout[,"value"] == 0.0)) {
            #remove zero-value events from the regressor
            regout <- regout[regout[,"value"] != 0, ]
            
            #now mean center values (unless there is no variation, such as a task indicator function)
            if (sd(regout[,"value"]) > 0) { regout[,"value"] <- regout[,"value"] - mean(regout[,"value"]) }            
          } 
          
          fname <- paste0("run", runsToOutput[i], "_", dimnames(dmat)[[2L]][reg], "_FSL3col.txt")
          write.table(regout, file=file.path(output_directory, fname), sep="\t", eol="\n", col.names=FALSE, row.names=FALSE)
        }
      }
    }
    
    if ("AFNI" %in% writeTimingFiles) {
      #use dmBLOCK-style regressors: time*modulation:duration. One line per run
      
      #need to unify multiple values that share onset time
      #time*modulation1,modulation2:duration
      regonsets <- apply(dmat, 2, function(reg) {
            unname(do.call(c, lapply(reg, function(run) { run[,"onset"]})))
          })
      
      regdurations <- apply(dmat, 2, function(reg) {
            unname(do.call(c, lapply(reg, function(run) { run[,"duration"]})))
          })
      
      #my code (longer than data.table)
#      combcols <- list()
#      colstomatch <- 1:ncol(regonsets)
#      while(length(colstomatch) > 0L) {
#        matches <- colstomatch[1L] #always match self
#        if (length(colstomatch) > 1L) {
#          for (coltest in colstomatch[2:length(colstomatch)]) {
#            if (identical(regonsets[,matches[1L]], regonsets[,coltest])) {
#              matches <- c(matches, coltest)
#            }
#          }
#        }
#        combcols[[paste(matches, collapse=".")]] <- matches
#        colstomatch <- colstomatch[-1*which(colstomatch %in% matches)]
#      }
      
      #magical code from here: http://stackoverflow.com/questions/22993637/efficient-r-code-for-finding-indices-associated-with-unique-values-in-vector
      first_onset_duration <- paste(regonsets[1,], regdurations[1,]) #use combination of onset and duration to determine unique dmBLOCK regressors
      dt = as.data.table(first_onset_duration)[, list(comb=list(.I)), by=first_onset_duration] #just matches on first row (should work in general)
      
      lapply(dt$comb, function(comb) {
            combmat <- dmat[,comb, drop=F]
            #onsets and durations are constant, amplitudes vary
            runvec <- c() #character vector of AFNI runs (one row per run)
            for (i in 1:dim(combmat)[1L]) {
              runonsets <- combmat[[i,1]][,"onset"] #just use first regressor of combo to get vector onsets and durations (since combinations, by definition, share these)
              rundurations <- combmat[[i,1]][,"duration"]
              runvalues <- do.call(cbind, lapply(combmat[i,], function(reg) { reg[,"value"] }))
              
              #AFNI doesn't like us if we pass in the boxcar ourselves in the dmBLOCK format (since it creates this internally). Filter out.
              indicatorFunc <- apply(runvalues, 2, function(col) { all(col == 1.0)} )
              if (any(indicatorFunc)) { runvalues <- runvalues[,-1*which(indicatorFunc), drop=FALSE] }
              
              #if the indicator regressor was the only thing present, revert to the notation TIME:DURATION notation for dmBLOCK (not TIME*PARAMETER:DURATION)
              if (ncol(runvalues) == 0L) {
                runvec[i] <- paste(sapply(1:length(runonsets), function(j) { 
                          paste0(plyr::round_any(runonsets[j], .000001), ":", plyr::round_any(rundurations[j], .000001))
                        }), collapse=" ")
              } else {
                runvec[i] <- paste(sapply(1:length(runonsets), function(j) { 
                          paste0(plyr::round_any(runonsets[j], .000001), "*", paste(plyr::round_any(runvalues[j,], .000001), collapse=","), ":", plyr::round_any(rundurations[j], .000001))
                        }), collapse=" ")
              }              
            }
            
            writeLines(runvec, file.path(output_directory, paste0(paste(dimnames(combmat)[[2L]], collapse="_"), "_dmBLOCK.txt")))
          })
      
    }
  }
  
  collinearityDiag.raw <- apply(dmat, 1, function(run) {
        #custom regressors are not guaranteed to have ntrials as the dimension. Remove them from raw diagnostics since they're not on the same trial grid
        #rlengths <- sapply(run, length)
        custom_reg <- grep("^custom_", names(run))
        if (length(custom_reg) > 0L) { run <- run[-1*custom_reg]}
            
        #check correlations among regressors for trial-wise estimates
        cmat <- do.call(data.frame, lapply(run, function(regressor) {
                  regressor[,"value"]
                }))
        
        #remove any constant columns (e.g., task indicator regressor for clock) so that VIF computation is sensible
        zerosd <- sapply(cmat, sd)
        cmat_noconst <- cmat[,zerosd != 0.0, drop=FALSE]

        if (ncol(cmat_noconst) == 0L) {
          #if only task indicator regressors are included, then they cannot be collinear before convolution since there is no variability
          return(list(r=NA_real_, vif=NA_real_))
        } else {
          corvals <- suppressWarnings(cor(cmat, use="pairwise.complete.obs")) #drop warnings about SD of 0 since we want it to populate NAs there.
          vifMat <- data.frame(cbind(dummy=rnorm(nrow(cmat_noconst)), cmat_noconst)) #add dummy constant for vif
          vifForm <- as.formula(paste("dummy ~ 1 +", paste(names(cmat_noconst), collapse=" + ")))
          
          varInfl <- tryCatch(car::vif(lm(vifForm, data=vifMat)), error=function(e) { rep(NA, ncol(cmat_noconst)) }) #return NA if failure
          return(list(r=corvals, vif=varInfl))  
        }
      })
  
  collinearityDiag.convolve <- lapply(dmat.convolve, function(run) { #apply(dmat.convolve, 1, function(run) {                    
        corvals <- cor(run, use="pairwise.complete.obs")
        vifMat <- data.frame(cbind(dummy=rnorm(nrow(run)), run)) #add dummy constant for vif 
        vifForm <- as.formula(paste("dummy ~ 1 +", paste(names(run), collapse=" + ")))
        
        varInfl <- tryCatch(car::vif(lm(vifForm, data=vifMat)), error=function(e) { NA }) #return NA if failure
        list(r=corvals, vif=varInfl)
      })
  
  #add baseline terms to convolved design matrices
  if (baselineCoefOrder > -1L) {
    dmat.convolve <- lapply(dmat.convolve, function(r) {
          n <- names(r)
          if (baselineParameterization == "Legendre") {
            #consistent with AFNI approach, use Legendre polynomials, which are centered at zero to allow for appropriate baseline 
            unnormalized.p.list <- legendre.polynomials( baselineCoefOrder, normalized=FALSE )
            #evaluate polynomial between -1 and 1, which is where its desirable centered behavior exists.
            #although the functions are centered at 0, evaluating them on a given grid may lead to slight deviations from mean=0. Thus, center perfectly.
            baseline <- polynomial.values(polynomials=unnormalized.p.list, x=seq(-1,1, length.out=nrow(r)))
            baseline[2:length(baseline)] <- lapply(baseline[2:length(baseline)], function(v) { v - mean(v) }) #don't center constant!
            
            names(baseline) <- paste0("base", 0:baselineCoefOrder)
            d <- cbind(r, baseline)
          } else {
            #compute polynomials that are orthogonal to design effects
            #this may be somewhat dubious. for example, a linear trend in a task regressor would be preserved because it is given preference over baseline 
            d <- data.frame(fmri.design(r, order=baselineCoefOrder))
            names(d) <- c(n, paste0("base", 0:baselineCoefOrder))
          }          
          d
        })
    
  }
  
  return(list(design=dmat, design.convolve=dmat.convolve, collin.raw=collinearityDiag.raw, collin.convolve=collinearityDiag.convolve, concat_onsets=concat_onsets, runVolumes=runVolumes))
  
}

#' Concatenate design matrices for each run to form a single design with unique baselines per run (ala AFNI)
#'
#' @importFrom plyr rbind.fill
#' @export
concatDesignRuns <- function(d) {
  
  d_allruns <- do.call(rbind.fill, lapply(1:length(d$design.convolve), function(r) {
            thisrun <- d$design.convolve[[r]]
            basecols <- grepl("base", names(thisrun))
            ##note that this will rename repeated names into something like ev and ev.1, which is good
            names(thisrun) <- gsub("base", paste0("run", r, "base"), names(thisrun))
            thisrun
          }))
  
  d_allruns[which(is.na(d_allruns), arr.ind=TRUE)] <- 0
  d_allruns <- as.matrix(d_allruns) #needed for lm
  
  d_allruns
}


#' Visualize design matrix, including event onset times and run boundaries
#'
#' @importFrom ggplot2 ggplot aes geom_line theme_bw facet_grid geom_vline ggsave
#' @importFrom reshape2 melt
#' @export
visualizeDesignMatrix <- function(d, outfile=NULL, runboundaries=NULL, events=NULL, includeBaseline=TRUE) {
  
  if (!includeBaseline) {
    d <- d[,!grepl("(run[0-9]+)*base", colnames(d))]
  }
  
  print(round(cor(d), 3))
  d <- as.data.frame(d)
  d$volume <- 1:nrow(d)
  d.m <- melt(d, id.vars="volume")
  g <- ggplot(d.m, aes(x=volume, y=value)) + geom_line(size=1.2) + theme_bw(base_size=15) + facet_grid(variable ~ ., scales="free_y")
  
  colors <- c("black", "blue", "red", "orange") #just a hack for color scheme right now
  
  if (!is.null(runboundaries)) {
    g <- g + geom_vline(xintercept=runboundaries, color=colors[1L])
  }
  
  if (!is.null(events)) {
    for (i in 1:length(events)) {
      g <- g + geom_vline(xintercept=events[[i]], color=colors[i+1])
    }
  }
  
  if (!is.null(outfile)) {
    ggsave(filename=outfile, plot=g, width=21, height=9)
  }
  
  return(invisible(g))
}


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
                orig_data_frame=d)
            
            if ("emotion" %in% names(d)) { r$run_condition <- as.character(d$emotion[1L]) }            
            if ("clock_onset" %in% names(d)) { r$clock_onset <- d$clock_onset }
            if ("feedback_onset" %in% names(d)) { r$feedback_onset <- d$feedback_onset }
            if ("iti_onset" %in% names(d)) { r$iti_onset <- d$iti_onset }
            
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
#' @importFrom ggplot2 geom_line ggtitle theme_bw
#' @export clock_fit
#' @exportClass clock_fit
clock_fit <- setRefClass(
    Class="clock_fit",
    fields=list(
        nruns="numeric", #scalar number of total runs
        ntrials="numeric", #scalar of total trials
        RTobs="list", #trial vector (run), run x trial matrix (subject), or subject x run x trial matrix (group)
        RTraw="list", #original RTs. RTobs is the RTs that were fitted (e.g., could be raw, smoothed, differenced, etc.)
        RTpred="list",
        Reward="list",
        total_SSE="numeric", #scalar of total sum of squared errors
        theta="matrix", #named matrix of parameters and bounds
        SSE="numeric", #vector or matrix of SSEs over subjects and runs
        AIC="numeric", 
        nparams="numeric",
        elapsed_time="numeric",
        profile_data="list",
        opt_data="list", #list of results from optimizer
        pred_contrib="list",
        clock_onset="list",
        feedback_onset="list",
        iti_onset="list",
        bfs_var_fast="list", #trialwise variance of beta distribution for fast responses
        bfs_var_slow="list",
        bfs_mean_fast="list",
        bfs_mean_slow="list",
        ev="list",
        rpe="list",
        sceptic="list",
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
            #convert to list of list objects so that the lapply calls below work for both single and multiple runs
            clock_data <- list(runs=list(clock_data))
          }
          
          #Nov 2015: transitioned to list storage for multi-run attributes (e.g., RTs) to allow for unequal numbers of trials across runs.
          if (length(nruns)==0L) { nruns <<- length(clock_data$runs) }
          if (length(ntrials)==0L) { ntrials <<- sum(unlist(lapply(clock_data$runs, function(r) { r$w$ntrials } ))) }
          if (length(RTobs)==0L) { RTobs <<- lapply(clock_data$runs, function(r) { r$RTobs }) } #FIXME
          if (length(RTraw)==0L) { RTraw <<- lapply(clock_data$runs, function(r) { r$RTraw }) } #FIXME
          if (length(RTpred)==0L) { RTpred <<- lapply(clock_data$runs, function(r) { r$w$RTpred }) }
          if (length(Reward)==0L) { Reward <<- lapply(clock_data$runs, function(r) { r$Reward })  }
          if (length(clock_onset)==0L) { clock_onset <<- lapply(clock_data$runs, function(r) { if(length(r$clock_onset) == 0 ) NA else r$clock_onset }) }
          if (length(feedback_onset)==0L) { feedback_onset <<- lapply(clock_data$runs, function(r) { if(length(r$feedback_onset) == 0 ) NA else r$feedback_onset }) }
          if (length(iti_onset)==0L) { iti_onset <<- lapply(clock_data$runs, function(r) { if(length(r$iti_onset) == 0 ) NA else r$iti_onset }) }
          if (length(rew_function)==0L) { rew_function <<- do.call(c, lapply(clock_data$runs, function(r) { r$rew_function })) }
          if (length(run_condition)==0L) { run_condition <<- do.call(c, lapply(clock_data$runs, function(r) { r$run_condition })) }
          if (length(ev)==0L) { ev <<- lapply(clock_data$runs, function(r) { r$w$V }) } #expected value
          if (length(rpe)==0L) { rpe <<- lapply(1:length(Reward), function(r) { Reward[[r]] - ev[[r]] }) } #better or worse than expected?
          
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
            bfs_var_fast <<- lapply(clock_data$runs, function(r) { r$w$betaFastSlow$var_fast }) #trialwise estimate of uncertainty re: fast responses
            bfs_var_slow <<- lapply(clock_data$runs, function(r) { r$w$betaFastSlow$var_slow }) #trialwise estimate of uncertainty re: slow responses
            bfs_mean_fast <<- lapply(clock_data$runs, function(r) { r$w$betaFastSlow$mean_fast }) #trialwise estimate of mean value for fast responses
            bfs_mean_slow <<- lapply(clock_data$runs, function(r) { r$w$betaFastSlow$mean_slow }) #trialwise estimate of mean value for slow responses
          } else { hasBeta <- FALSE }
        },
        plotRTs=function() {
          rtDf <- data.frame(
              trial=rep(1:ncol(RTobs), nrow(RTobs)),
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

