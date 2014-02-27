\name{fmri.stimulus}
\alias{fmri.stimulus}
\title{convolve a regressor with a hemodynamic response function for fMRI analysis.}
\usage{
fmri.stimulus(scans = 1, onsets = c(1), durations = c(1), values = c(1),
  rt = 3, times = NULL, mean = TRUE, a1 = 6, a2 = 12, b1 = 0.9,
  b2 = 0.9, cc = 0.35)
}
\description{
extended from fmri package to allow for continuous-valued
regressor, which is passed using the values parameter.
}
\keyword{internal}

