#' convolve a regressor with a normalized HRF with peak amplitude = 1.0
#' 
#' extends fmri.stimulus by allowing for two normalization approaches (building on AFNI dmUBLOCK):
#'   1) pre-convolution HRF max=1.0 normalization of each stimulus regardless of duration: identical to dmUBLOCK(1)
#'   2) pre-convolution HRF max=1.0 normalization for long events (15+ sec) -- height of HRF is modulated by duration of event: identical to dmUBLOCK(0)
#' @export
hrf_convolve_normalize <- function(scans, times, durations, values, rt=1.0, normeach=FALSE, rmzeros=TRUE,
    demean_events=TRUE, demean_convolved=FALSE, a1 = 6, a2 = 12, b1 = 0.9, b2 = 0.9, cc = 0.35) {
  
  #this is my hacky way to figure out the peak value 
  #obtain an estimate of the peak HRF height of a long event convolved with the HRF at these settings of a1, b1, etc.
  hrf_boxcar <- fmri.stimulus(scans=300, values=1.0, times=100, durations=100, rt=1.0, mean=FALSE,  #don't mean center for computing max height 
      a1=a1, a2=a2, b1=b1, b2=b2, cc=cc)
  hrf_max <- max(hrf_boxcar)
  
  #remove zeros from events vector to avoid these becoming non-zero values in the hrf convolution after grand mean centering
  #for example, for negative RPEs, if one does not occur, it will have a value of zero. But after grand mean centering below, this
  #will now be treated as an "event" by the HRF convolution since the amplitude is no longer zero.
  if (rmzeros) {
    zeros <- which(values==0.0)
    if (length(zeros) > 0L) {
      times <- times[-1.0*zeros]
      values <- values[-1.0*zeros]
      durations <- durations[-1.0*zeros]
    }
  }
  
  #demean parametric regressor (good idea to remove collinearity)
  if (demean_events && !all(values==1.0)) {
    values <- values - mean(values)
  }
  
  #for each event, convolve it with hrf, normalize, then sum convolved events to get full timecourse
  normedEvents <- sapply(1:length(times), function(i) {
        #obtain unit-convolved duration-modulated regressor to define HRF prior to modulation by parametric regressor
        stim_conv <- fmri.stimulus(scans=scans, values=1.0, times=times[i], durations=durations[i], rt=rt, mean=FALSE, a1=a1, a2=a2, b1=b1, b2=b2, cc=cc)
        if (normeach) {
          stim_conv <- stim_conv/max(stim_conv) #rescale HRF to a max of 1.0 for each event, regardless of duration -- EQUIVALENT TO dmUBLOCK(1)
        } else {
          stim_conv <- stim_conv/hrf_max #rescale HRF to a max of 1.0 for long event -- EQUIVALENT TO dmUBLOCK(0)
        }
        
        stim_conv <- stim_conv*values[i] #for each event, multiply by parametric regressor value
      })
  
  tc_conv <- apply(normedEvents, 1, sum)
  
  #grand mean center convolved regressor
  if (demean_convolved) { tc_conv <- tc_conv - mean(tc_conv) }
  
  tc_conv
}

#' convolve a regressor with a hemodynamic response function for fMRI analysis.
#' 
#' extended from fmri package to allow for continuous-valued regressor,
#' which is passed using the values parameter.
#' @export
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
}

#' compute a moving average smooth over a time series (here, a vector of RTs)
#'  
#' used to fit smoothed RTs (\code{clock_model} object).
#' Function is an adapted version of \code{runmean} from the \code{caTools} package.
#' @keywords internal
runmean <- function(x, k=5) {
  n <- length(x)
  k2 = k%/%2
  k1 = k - k2 - 1
  y = c(sum(x[1:k]), diff(x, k))
  y = cumsum(y)/k
  y = c(rep(0, k1), y, rep(0, k2))
  
  idx1 = 1:k1
  idx2 = (n - k2 + 1):n
  
  #use mean of available data at tails
  for (i in idx1) y[i] = mean(x[1:(i + k2)])
  for (i in idx2) y[i] = mean(x[(i - k1):n])
  
  return(y)
}
