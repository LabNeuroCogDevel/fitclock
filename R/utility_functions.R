#' convolve a regressor with a normalized HRF with peak amplitude = 1.0
#' 
#' extends fmri.stimulus by allowing for two normalization approaches (building on AFNI dmUBLOCK):
#'   1) pre-convolution HRF max=1.0 normalization of each stimulus regardless of duration: identical to dmUBLOCK(1)
#'   2) pre-convolution HRF max=1.0 normalization for long events (15+ sec) -- height of HRF is modulated by duration of event: identical to dmUBLOCK(0)
#' @export
hrf_convolve_normalize <- function(scans, times, durations, values, rt=1.0, normeach=FALSE, rm_zeros=TRUE,
    center_values=TRUE, parmax1=FALSE, demean_convolved=FALSE, a1 = 6, a2 = 12, b1 = 0.9, b2 = 0.9, cc = 0.35) {
  
  #this is my hacky way to figure out the peak value 
  #obtain an estimate of the peak HRF height of a long event convolved with the HRF at these settings of a1, b1, etc.
  hrf_boxcar <- fmri.stimulus(scans=300, values=1.0, times=100, durations=100, rt=rt, demean=FALSE,  #don't mean center for computing max height 
      a1=a1, a2=a2, b1=b1, b2=b2, cc=cc)
  hrf_max <- max(hrf_boxcar)
  
  #N.B. this code is duplicated below in fmri.stimulus to allow the core convolution function to mean-center parametric regressors
  #without any form of height normalization, as implemented here. may be able to reduce this code here, but because each even receives
  #a unique convolution in the code below, would need to maintain centering at least.  
  
  #remove zeros from events vector to avoid these becoming non-zero values in the hrf convolution after grand mean centering
  #for example, for negative RPEs, if one does not occur, it will have a value of zero. But after grand mean centering below, this
  #will now be treated as an "event" by the HRF convolution since the amplitude is no longer zero.
  if (rm_zeros) {
    zeros <- which(values==0.0)
    if (length(zeros) > 0L) {
      times <- times[-1.0*zeros]
      values <- values[-1.0*zeros]
      durations <- durations[-1.0*zeros]
    }
  }
  
  if (length(times) == 0L) {
    warning("No non-zero events for regressor to be convolved. Returning all-zero result for fMRI GLM.")
    return(rep(0, scans))
  }
  
  #handle mean centering of parametric values prior to convolution
  #this is useful when one wishes to dissociate variance due to parametric modulation versus stimulus occurrence
  if (center_values && !all(values==1.0)) {
    values <- values - mean(values)
  }
      
  #for each event, convolve it with hrf, normalize, then sum convolved events to get full timecourse
  normedEvents <- sapply(1:length(times), function(i) {
        #obtain unit-convolved duration-modulated regressor to define HRF prior to modulation by parametric regressor
        stim_conv <- fmri.stimulus(scans=scans, values=1.0, times=times[i], durations=durations[i], rt=rt, demean=FALSE, center_values=FALSE, a1=a1, a2=a2, b1=b1, b2=b2, cc=cc)
        if (normeach) {
          stim_conv <- stim_conv/max(stim_conv) #rescale HRF to a max of 1.0 for each event, regardless of duration -- EQUIVALENT TO dmUBLOCK(1)
        } else {
          stim_conv <- stim_conv/hrf_max #rescale HRF to a max of 1.0 for long event -- EQUIVALENT TO dmUBLOCK(0)
        }
        
        stim_conv <- stim_conv*values[i] #for each event, multiply by parametric regressor value
      })
  
  tc_conv <- apply(normedEvents, 1, sum)
  
  #allow for the convolved regressor to be rescaled to have a maximum height of 1.0. This normalizes the range of the regressor across runs and subjects
  #such that the betas for the regressor are on similar scales, just with a change of gain... experimenting with this currently given SCEPTIC DAUC regressors
  #which tend to be highly skewed, potentially driven by behavioral parameter scaling challenges. Note that this interacts with the normeach setting such that
  #for normeach=FALSE (durmax_1.0), the height of 1.0 will reflect a combination of the parameter and the duration. The interpretation is somewhat cleaner under
  #normeach=TRUE (evtmax_1.0) since the HRF has height of 1.0 for each stimulus, and then rescaling to max 1.0 after convolution with the parametric value
  #will change the relative heights of the stimuli solely according to the parameter.
  if (parmax1) {
    tc_conv <- tc_conv/max(tc_conv)
  }
  
  #grand mean center convolved regressor
  if (demean_convolved) { tc_conv <- tc_conv - mean(tc_conv) }
  
  tc_conv
}

#' convolve a regressor with a hemodynamic response function for fMRI analysis.
#' 
#' extended from fmri package to allow for continuous-valued regressor,
#' which is passed using the values parameter.
#' also support mean centering of parametric regressor prior to convolution 
#' to dissociate it from stimulus occurrence (when event regressor also in model)
#' 
#' @export
fmri.stimulus=function(scans=1, onsets=c(1), durations=c(1), values=c(1), center_values=FALSE, rm_zeros=TRUE, convolve=TRUE,
    rt=3, times=NULL, demean=TRUE, a1 = 6, a2 = 12, b1 = 0.9, b2 = 0.9, cc = 0.35, parmax1=FALSE) {
  
  mygamma <- function(x, a1, a2, b1, b2, c) {
    d1 <- a1 * b1
    d2 <- a2 * b2
    c1 <- ( x/d1 )^a1
    c2 <- c * ( x/d2 )^a2
    res <- c1 * exp(-(x-d1)/b1) - c2 * exp(-(x-d2)/b2)
    res
  }
  
  #remove zeros from events vector to avoid these becoming non-zero values in the hrf convolution after grand mean centering
  #for example, for negative RPEs, if one does not occur, it will have a value of zero. But after grand mean centering below, this
  #will now be treated as an "event" by the HRF convolution since the amplitude is no longer zero.
  if (rm_zeros) {
    zeros <- which(values==0.0)
    if (length(zeros) > 0L) {
      times <- times[-1.0*zeros]
      values <- values[-1.0*zeros]
      durations <- durations[-1.0*zeros]
    }
  }
  
  if (length(times) == 0L) {
    warning("No non-zero events for regressor to be convolved. Returning all-zero result for fMRI GLM.")
    return(rep(0, scans))
  }

  #handle mean centering of parametric values prior to convolution
  #this is useful when one wishes to dissociate variance due to parametric modulation versus stimulus occurrence
  if (center_values && !all(values==1.0)) {
    values <- values - mean(values)
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
  
  #rescale regressor to maximum height of 1.0 for scaling similarity across instances (see hrf_convolve_normalize for details)
  #only applicable to convolve regressors 
  if (parmax1 && convolve) { hrf <- hrf/max(hrf) }
  
  if (!convolve) {
    #just return the box car without convolving by HRF
    stimulus <- stimulus[-(1:(20*scale))][1:scans] #remove zero padding
    stimulus <- stimulus[unique((scale:scans)%/%scale)*scale] #subset the elements of the upsampled grid back onto the observed TRs
    return(stimulus)
  } else if (demean) {
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

#' detrend time series up to quadratic trend. Used by fir1Bandpass prior to filtering
#'  
#' @keywords internal
detrendts <- function(x, order=0) {
  #order 0=demean; order 1=linear; order 2=quadratic
  lin <- 1:length(x)
  quad <- lin^2
  
  if (order == 0)
    residuals(lm(x~1))
  else if (order == 1)
    residuals(lm(x ~ 1 + lin))
  else if (order == 2)
    residuals(lm(x ~ 1 + lin + quad))
  else
    stop("order not supported:", order)
}

#' function to FIR-1 bandpass filter a signal. Can low- or high-pass filter by specifying 0 for low or >= Nyquist for high.
#'  
#' @keywords internal
#' @importFrom signal fir1 filtfilt freqz
#' @export
fir1Bandpass <- function(x, TR=2.0, low=.009, high=.08, n=250, plotFilter=FALSE, forward_reverse=TRUE, padx=0, detrend=1) {
  
  #check for all NA
  if (all(is.na(x))) return(x)
  
  #n refers to filter order. 250 does quite well with typical signals
  Fs <- 1/TR
  nyq <- Fs/2
  
  #enforce filter upper bound at 1.0 (nyquist)
  if (high/nyq > 1.0) { high <- nyq }
  
  #coefficients are specified in the normalized 0-1 range.
  fir1Coef <- fir1(n, c(low/nyq, high/nyq), type="pass")
  
  if (plotFilter) print(freqz(fir1Coef, Fs=Fs))
  
  origLen <- length(x)
  
  #handle detrending (almost always a good idea to demean, if not detrend, for fourier series to be valid!)
  if (!is.null(detrend) && detrend >= 0)
    x <- detrendts(x, order=detrend)
  
  #zero-pad data, if requested
  x <- c(x, rep(0*x, padx))
  
  #as the order of the filter exceeds the length of the time series,
  #some sort of phase distortion is introduced.
  #forward+reverse filtering cleans it up
  if (forward_reverse) xfilt <- filtfilt(fir1Coef, x) 
  else xfilt <- filter(fir1Coef, x)
  
  return(xfilt[1:origLen])
}
