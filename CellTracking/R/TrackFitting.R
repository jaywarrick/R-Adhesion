#' @return the sum square error between the TrackList data and a sweep with the provided amplitude and phaseShift
sseBulk <- function(trackList, trackMatrix, amplitude, phaseShift, timeScalingFactor, ti)
{
     # Cols are frames, rows are track ids
     # we pass the 'trackMatrix' so we don't have to obtain it each iteration

     # get the sweep (for this we need the 'trackList')
     sweep <- getSweep(amplitude=amplitude, phaseShift=phaseShift, offset=pi, sin=trackList$sin, ti, fi=trackList$fi, ff=trackList$ff, sweepDuration=timeScalingFactor*trackList$sweepDuration, t=trackList$tAll, guess=NULL)
     # For each index in tAll (i.e., for each frame)
     sse <- sum((t(trackMatrix)-sweep$v)^2, na.rm=TRUE) # Do the transpose because the subtract function typically makes the subtracted vector vertical
     cat("(", amplitude, ",", timeScalingFactor, ",", ti, ") = ", sse, "\n", sep="")
     return(sse)
}

#' @return The best fit of a sweep fitted to all the data simultaneously to obtain an ensemble prediction of the phase shift to be used for all the Tracks
getBulkPhaseShift <- function(trackList, tiGuess=0)
{
     phaseShift <- pi
     require(stats)
     trackMatrix <- trackList$getMatrix()
     amplitude <- max(as.numeric(trackList$getProp(fun=function(x){r <- x$range('x', rel=TRUE); r <- (r[2]-r[1])/5; return(r)})))
     guess <- c(timeScalingFactor=1, ti=tiGuess)
     #guess <- c(amplitude=amplitude, phaseShift=1.2*pi)
     #amplitudeLimits=c(0.1*amplitude, amplitude)
     #phaseShiftLimits=c(1*pi,1.4*pi)
     timeScalingFactorLimits=c(0.99, 1.01)
     tiLimits=c(-0.5/trackList$fi, 0.5/trackList$fi)
     #sweepDurationLimits=c(0.98*trackList$sweepDuration, 1.02*trackList$sweepDuration)
     #fiLimits=c(0.98*trackList$fi, 1.02*trackList$fi)
     #ffLimits=c(0.98*trackList$ff, 1.02*trackList$ff)
     bestFit <- optim(par=guess,
                      function(par, trackList, trackMatrix){sseBulk(trackList=trackList, trackMatrix=trackMatrix, amplitude=amplitude, phaseShift=phaseShift, timeScalingFactor=par['timeScalingFactor'], ti=par['ti'])},
                      method='L-BFGS-B',
                      #method='SANN',
                      lower=c(min(timeScalingFactorLimits), min(tiLimits)),
                      upper=c(max(timeScalingFactorLimits), max(tiLimits)),
                      control=list(trace=0),
                      trackList=trackList,
                      trackMatrix=trackMatrix)
     return(list(par=c(phaseShift=phaseShift, amplitude=amplitude, timeScalingFactor=as.numeric(bestFit$par['timeScalingFactor']), ti=as.numeric(bestFit$par['ti']), offset=as.numeric(bestFit$par['offset'])), fit=bestFit))
}

#'Fist SSE BULK

#' @return The best fit of a sweep fitted to all the data simultaneously to obtain an ensemble prediction of the phase shift to be used for all the Tracks
getBulkPhaseShift2 <- function(trackList, resolution=0.05)
{
     require(stats)
     trackMatrix <- trackList$getMatrix()
     amplitude <- max(as.numeric(trackList$getProp(fun=function(x){r <- x$range('x', rel=TRUE); r <- (r[2]-r[1])/2; return(r)})))
     #guess <- c(amplitude=amplitude, phaseShift=0)
     guess <- c(amplitude=amplitude)
     amplitudeLimits=c(0.5, 50*amplitude)
     phaseShiftLimits=c(1*pi,1.4*pi)
     phaseShifts <- seq(-pi,pi,resolution*2*pi)
     sseList <- c()
     amplitudeList <- c()
     for(phaseShift in phaseShifts)
     {
          bestFit <- optim(par=guess,
                           function(par, trackList, trackMatrix){sseBulk(trackList=trackList, trackMatrix=trackMatrix, amplitude=par['amplitude'], phaseShift=phaseShift)},
                           method='L-BFGS-B',
                           lower=c(min(amplitudeLimits)),
                           upper=c(max(amplitudeLimits)),
                           control=list(trace=0),
                           trackList=trackList,
                           trackMatrix=trackMatrix)
          amplitudeList <- c(amplitudeList, bestFit$par['amplitude'])
          sseList <- c(sseList, bestFit$value)
     }

     i <- which.min(sseList);
     return(list(par=c(phaseShift=phaseShifts[i], amplitudeList[i]), sseList=sseList))
}
