#' convolve a regressor with a hemodynamic response function for fMRI analysis.
#' 
#' extended from fmri package to allow for continuous-valued regressor,
#' which is passed using the values parameter.
#' @keywords internal
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
