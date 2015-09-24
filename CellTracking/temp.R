setValidFrames = function(trackList, amplitude=1, phaseShift=-pi, validStart=0.01, validEnd=0.99)
{
     "
     #' Takes the fit and determines where the cells switch directions
     #'
     #' 'validStart' and 'validEnd' represent precentages. In other words, after the switch in direction,
     #' valid points start at 'validStart' % of the way to the next switch in direction while 'validEnd'
     #' occurs 'validEnd' % of the way to the that same next switch in direction.
     #'
     #' This filter is good for removing inaccurate values of the velocity when using a triangle waveform because
     #' the estimate of velocity will be artificially be lower due to sample aliasing of particle motion.
     #' Possible values are 0, 1, 2, and 3 for inflectionPtsToFilter. Multiple at once can be provided.
     #' 0 Represents the upward zero-crossing of the position and max positive velocity
     #' 1 Represents the max position and downward zero-crossing of the velocity
     #' 2 Represents the downward zero-crossing of the position and min (i.e most negative) velocity
     #' 4 Represents the min (i.e., most negative) position and upward zero-crossing of the velocity
     #'
     #' 'inflectionPtsToFilter' refers to the first point of the numbered sections (i.e., 1 refers to the max point)
     #' Thus if 1 is included in 'inflectionPtsToFilter', points adjacent to this point will be considered for wobble and filtered.
     #'
     #' We use 'inflectionPtsToFilter' of c(1,3) to mark the points at which the flow switches direction
     "

     # Helper function: get data points adjacent to specified inflection timepoint
     getNearests <- function(t, inflectionPoint)
     {
          upper <- which.max((t-inflectionPoint) >= 0)
          if(is.na(upper) || upper == 1)
          {
               return(NA)
          }
          lower <- upper - 1
          return(c(lower, upper))
     }

     sweep <- getSweep(amplitude=amplitude, phaseShift=phaseShift, offset=0, sin=trackList$sin, fi=trackList$fi, ff=trackList$ff, sweepDuration=trackList$sweepDuration, t=trackList$tAll, guess=NULL)
     inflectionsToAddress <- sweep$inflectionNums %in% c(1,3) # These are times at which flow switches directions
     indicesToRemove <- numeric(0)
     for(i in which(inflectionsToAddress))
     {
          nearests <- getNearests(sweep$t, sweep$inflections[i])
          if(!is.na(nearests)[1])
          {
               indicesToRemove <- c(indicesToRemove, nearests)
          }
     }
     validFrames0 <- trackList$allFrames
     validFrames0 <- validFrames0[-indicesToRemove]

     validFrames1 <- numeric(0)
     for(tIndex in 1:base::length(trackList$tAll))
     {
          # Get the nearest inflection at or beyond this time

          infIndex <- which.max((sweep$inflections >= trackList$tAll[tIndex]) & inflectionsToAddress)
          if(is.na(infIndex)) next

          # Get the bounding inflections that represent changes in fluid direction
          infT2 <- sweep$inflections[infIndex] # take the inflection we found
          if((infIndex-2) < 1)
          {
               infT1 <- 0
          }
          else
          {
               infT1 <- sweep$inflections[infIndex-2] # also take two inflections prior because each inflection represents pi/2 and we want to go back to the last change in direction which is pi ago.
          }
          dInfT <- infT2-infT1 # define the time interval size between these two inflections

          # Within the if statement calculate the fractional location of this time index in the interval between the two inflections.
          if( (trackList$tAll[tIndex] >= (infT1 + validStart*dInfT)) & (trackList$tAll[tIndex] <= (infT1 + validEnd*dInfT)) )
          {
               # If it is within the startValid and endValid bounds, add it to the list of the valid frames
               validFrames1 <- c(validFrames1, trackList$allFrames[tIndex]) # (tIndex-1) = frame because frames are indices that start at 0
          }
     }

     validFrames <<- sort(intersect(validFrames1, validFrames0))

     for(track in trackList$tracks)
     {
          track$setValidFrames(validFrames)
     }
}
