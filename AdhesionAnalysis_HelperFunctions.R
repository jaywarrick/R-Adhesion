
library(signal)
library(zoo)
library(pracma)
# source("http://bioconductor.org/biocLite.R")
# biocLite('GraphAlignment')
library('GraphAlignment')

##### Analysis Methods #####

interpolateDerivative <- function(f0, f1, f2, x0, x1, x2, xj)
{
     # This is a three point interpolation of the derivative where the interpolated point
     # is the middle of the 3 points. This simplifies to the the three-point midpoint formula
     # when the time steps are equal but can handle when timepoints are missing (i.e., the 
     # time-step on either side of the 3 points is not equal)
     term1 <- f0*((2*xj-x1-x2)/((x0-x1)*(x0-x2)))
     term2 <- f1*((2*xj-x0-x2)/((x1-x0)*(x1-x2)))
     term3 <- f2*((2*xj-x0-x1)/((x2-x0)*(x2-x1)))
     v <- term1 + term2 + term3
}

first <- function(x)
{
     return(x[1])
}

last <- function(x)
{
     return(x[numel(x)])
}

localDerivative <- function(x, t, i)
{
     if(i == 1)
     {
          #return((x[i+1]-x[i])/(t[i+1]-t[i]))
          return(interpolateDerivative(x[i], x[i+1], x[i+2], t[i], t[i+1], t[i+2], t[i]))
     }
     else if(i == length(x))
     {
          #return((x[i]-x[i-1])/(t[i]-t[i-1]))
          return(interpolateDerivative(x[i-2], x[i-1], x[i], t[i-2], t[i-1], t[i], t[i]))
     }
     else
     {
          return(interpolateDerivative(x[i-1], x[i], x[i+1], t[i-1], t[i], t[i+1], t[i]))
     }
}

getDerivative <- function(x, t)
{
     v <- numeric(0)
     for(i in 1:length(x))
     {
          v <- c(v, localDerivative(x, t, i))
     }
     return(v)
}

normalize <- function(data, abs=FALSE)
{
     if(abs)
     {
          data2 <- abs(data)
          return((data)/(max(abs(data))))
     }
     else
     {
          return((data-min(data))/(max(data)-min(data)))
     }
}

trackLengthFilter <- function(track, min=0, max=1000000)
{
     l <- length(track)
     if(l >= min & l <=max)
     {
          return(TRUE)
     }
     else
     {
          return(FALSE)
     }
}

trackFrameFilter <- function(track, startMin=0, startMax=1000000, endMin=0, endMax=1000000)
{
     # All start and end points are inclusive
     startFrame <- first(track@frames)
     endFrame <- last(track@frames)
     
     if(startFrame >= startMin & startFrame <= startMax & endFrame >= endMin & endFrame <= endMax)
     {
          return(TRUE)
     }
     else
     {
          return(FALSE)
     }
}

filterTrackList <- function(trackList, filterFunction, ...)
{
     trackList@tracks <- Filter(function(x){filterFunction(x,...)}, trackList@tracks)
     return(trackList)
}

sortTrackList <- function(object, trackPropFun=function(object){return(length(object@frames))}, decreasing = TRUE) 
{
     ret <- getProp(object, trackPropFun=trackPropFun)
     sorted <- sort(unlist(ret), index.return=T, decreasing=decreasing)
     newTracks <- object@tracks[sorted$ix]
     object@tracks <- newTracks
     return(object)
}

applyTrackList <- function(object, trackFun=function(track){return(calculateVelocities(track))}) 
{
     objectret <- getProp(object, trackPropFun=trackPropFun)
     sorted <- sort(unlist(ret), index.return=T, decreasing=decreasing)
     newTracks <- object@tracks[sorted$ix]
     object@tracks <- newTracks
     return(object)
}

setTrackListValidTimeIndices <- function(trackList, inflectionPtsToFilter=c(1,3), applyWobbleFilter=TRUE, wobbleStart=0.25, wobbleEnd=0.99)
{
     for(i in 1:length(trackList@tracks))
     {
          trackList@tracks[[i]] <- setTrackValidTimeIndices(trackList@tracks[[i]], inflectionPtsToFilter=c(1,3), applyWobbleFilter=applyWobbleFilter, wobbleStart=wobbleStart, wobbleEnd=wobbleEnd)
     }
     return(trackList)
}

getNearests <- function(t, inflectionPoint)
{
     diff <- t-inflectionPoint
     upper <- which(diff >= 0)[1]
     if(is.na(upper) || upper == 1)
     {
          return(NA)
     }
     lower <- upper - 1
     return(c(lower, upper))
}

setTrackValidTimeIndices <- function(track, inflectionPtsToFilter=c(1,3), applyWobbleFilter=TRUE, wobbleStart=0.25, wobbleEnd=0.99)
{
     #' Takes the track and determines which points are adjacent to the types of points listed in 
     #' the 'inflectionPtsToFilter' vector of integers and removes them. This filter is good for
     #' removing inaccurate values of the velocity when using a triangle waveform because the estimate
     #' of velocity will be artificially be lower due to sample aliasing of particle motion.
     #' Possible values are 0, 1, 2, and 3 for inflectionPtsToFilter. Multiple at once can be provided.
     #' 0 Represents the upward zero-crossing of the position and max positive velocity
     #' 1 Represents the max position and downward zero-crossing of the velocity
     #' 2 Represents the downward zero-crossing of the position and min (i.e most negative) velocity
     #' 4 Represents the min (i.e., most negative) position and upward zero-crossing of the velocity
     
     sweep <- getTrackSweep(track)
     indicesToRemove <- c()
     inflectionsToAddress <- sweep$inflectionNums %in% inflectionPtsToFilter # These are times at which flow switches directions
     for(i in which(inflectionsToAddress))
     {
          nearests <- getNearests(sweep$t, sweep$inflections[i])
          if(!is.na(nearests)[1])
          {
               indicesToRemove <- c(indicesToRemove, nearests)
          }
     }
     #print(indicesToRemove)
     validTimeIndices <- c(1:length(track@t))[-indicesToRemove]
     
     
     ### Account for the "wobble" of cells when flow switches directions
     if(applyWobbleFilter)
     {
          otherValidTimes <- numeric(0) # This 
          for(tIndex in 1:length(track@t))
          {
               #browser()
               #print(paste('tIndex:',tIndex))
               infIndex <- which((sweep$inflections >= track@t[tIndex]) & inflectionsToAddress)[1]
               #print(paste('t2Index:',t2Index))
               if(is.na(infIndex)) next
               infT2 <- sweep$inflections[infIndex]
               infT1 <- sweep$inflections[infIndex-2]
               dInfT <- infT2-infT1
               if( (track@t[tIndex] >= (infT1 + wobbleStart*dInfT)) & (track@t[tIndex] <= (infT1 + wobbleEnd*dInfT)) )
               {
                    otherValidTimes <- c(otherValidTimes, tIndex)
               }
               #browser()
          }
     }
     if(applyWobbleFilter)
     {
          #browser()
          validTimeIndices <- sort(intersect(validTimeIndices, otherValidTimes))
     }
     track@validTimeIndices <- validTimeIndices
     return(track)
}

##### Sweep Functions #####

getSweep <- function(sin=TRUE, amplitude=1, phaseShift=0, offset=0, fi=0.2, ff=0.01, tAll=seq(1,500,1), frames=-1, guess=NULL)
{
     if(frames[1] < 0 | frames[1] > length(tAll)-1)
     {
          frames <- (1:length(t)) - 1 #Subtract 1 to be consistent with frames from JEX, which start at 0 instead of 1
     }
     t <- tAll[frames + 1]
     ti <- first(t)
     tf <- last(t)
     N <- log(ff/fi)/log(2)
     R <- N / (tf-ti)
     
     if(is.null(guess))
     {
          A <- amplitude
          phi <- phaseShift
          b <- offset
     }
     else
     {
          A <- guess['amplitude']
          phi <- guess['phaseShift']
          b <- guess['offset']
     }
     
     ### Create "saw-tooth" version ###
     # line essentially means position. With a triangular waveform, the velocity equals 4*A*f(t) because in a single cycle, the position of a particle travels a total distance of 4*A.
     # If we integrate the velocity, we get the position (i.e., line), we shift this line by a phase shift described below.
     # Whenever the position travels a distance of A, we have essentially passed through pi/2 of a cycle.
     # So we do mod division of the position by A to determine the zero crossings and peaks of the triangular wave
     f <- fi*(ff/fi)^(tAll/tf)
     if(ff==fi | A == 0)
     {
          line <- 4*A*f*tAll - A*(phi/(pi/2)) # the last term shifts the line up and down, essentially dictating where the "sections" of the line are placed (i.e., the initial phase of the triangle wave)
          line1 <- 4*1*f*tAll - 1*(phi/(pi/2)) # WITH DUMMY AMPLITUDE OF 1 in case A is 0 then we have something to calculate sections with that isn't a simple flat line
     }
     else
     {
          line <- 4*A*fi*(  ((ff/fi)^(tAll/tf)*tf)/log(ff/fi)-tf/log(ff/fi)  ) - A*(phi/(pi/2)) # the last term shifts the line up and down, essentially dictating where the "sections" of the line are placed (i.e., the initial phase of the triangle wave)
          line1 <- 4*1*fi*(  ((ff/fi)^(tAll/tf)*tf)/log(ff/fi)-tf/log(ff/fi)  ) - 1*(phi/(pi/2)) # WITH DUMMY AMPLITUDE OF 1 in case A is 0 then we have something to calculate sections with that isn't a simple flat line
     }
     
     if(A == 0)
     {
          x <- b + line # Flatline
     }
     else
     {
          x <- b + (line %% A) #*sign(fSweep$x)
     }
     # Fix each pi/2 section of x
     sections <- (line1 %/% 1) %% 4
     x[sections==1] <- -1*(x[sections==1]-1)
     x[sections==2] <- -1*(x[sections==2])
     x[sections==3] <- (x[sections==3]-1)
     
     # Find pi/2 intervals on frequency sweep
     
     offsetInflections <- (  (-1*(phi/(pi/2))) %/% 1  ) 
     ni <- seq(1,1000,1) + offsetInflections
     suppressWarnings(inflections <- (log((log(ff/fi)*phi)/(2*fi*pi*tf)+(log(ff/fi)*ni)/(4*fi*tf)+1)*tf)/log(ff/fi))
     inflections <- inflections[!is.nan(inflections)]
     inflectionNums <- (  seq(0,length(inflections)-1,1) + (  offsetInflections %% 4  )  ) %% 4
     
     startTimeI <- which(inflections >= t[1])[1]
     endTimeI <- which(inflections > t[length(t)])[1] - 1
     inflections <- inflections[startTimeI:endTimeI]
     inflectionNums <- inflectionNums[startTimeI:endTimeI]
     
     predicted <- A*sin(2*pi*((fi*(-1+2^(R*tAll)))/(R*log(2))) - phi) + b
     
     #browser()
     #print((-A*(phi/(pi/2))) %/%A %% 4)
     if(sin)
     {
          v=getDerivative(x=predicted, t=tAll)[frames+1]
          return(list(A=A, t=t, x=predicted[frames+1], line=line, v=v, sections=sections, inflections=inflections, inflectionNums=inflectionNums))
     }
     else
     {
          v=getDerivative(x=x, t=tAll)[frames+1]
          return(list(A=A, t=t, x=x[frames+1], line=line, v=v, sections=sections, inflections=inflections, inflectionNums=inflectionNums))
     }
}

getTrackSweep <- function(track, slot='vx', amplitude=1, offset=mean(slot(track,slot)))
{
     return(getSweep(sin=track@parent@sin, amplitude=amplitude, phaseShift=track@phaseShift, offset=offset, fi=track@parent@fi, ff=track@parent@ff, tAll=track@parent@tAll, frames=track@frames, guess=NULL))
}

# 
# ##### TrackList #####
# 
# setClass('trackList', 
#          representation(tracks='list', tracksTable='data.frame', sin='logical', fi='numeric', ff='numeric', tAll='numeric', phaseShift='numeric'),
#          prototype(tracks=list(), tracksTable=c(), sin=c(), fi=c(), ff=c(), tAll=c(), phaseShift=c()))
# 
# setMethod('initialize', 'trackList', function(.Object, sin, fi, ff, tAll)
# {
#      .Object@sin <- sin
#      .Object@fi <- fi
#      .Object@ff <- ff
#      .Object@tAll <- tAll
#      return(.Object)
# })
# 
# trackListFromFile <- function(object, file=NULL, sin=FALSE, fi, ff, tAll)
# {
#      library(foreign)
#      tracksFile <- read.arff(file)
#      tracksFile2 <- reorganizeTable(tracksFile, nameCol='Metadata')
#      
#      duh <- list()
#      for(row in 1:length(tracksFile2[,1]))
#      {
#           id <- tracksFile2[row,'Track']
#           start <- tracksFile2[row,'polygonPts']
#           pattern <- tracksFile2[row,'patternPts']
#           newTrack <- trackFromTrackROI(id=id, sin=sin, fi=fi, ff=ff, tAll=tAll, start=start, pattern=pattern)
#           duh[as.character(tracksFile2[row,'Track'])] <- newTrack
#      }
#      
#      object <- new('trackList', sin=sin, fi=fi, ff=ff, tAll=tAll)
#      object@tracks <- duh
#      object@tracksTable <- tracksFile2
#      return(object)
# }
# 
# setMethod('length', 'trackList', function(x)
# {
#      length(x@tracks)
# }
# )
# 
# setMethod('plot', 'trackList', function(x, slot='vx', rel=FALSE, ...)
# {
#      
#      xRanges <- getProp(x, function(track){range(track@t)})
#      xRanges <- matrix(unlist(xRanges), ncol=2, byrow=TRUE)
#      yRanges <- getProp(x, function(track){range(slot(track, slot))})
#      yRanges <- matrix(unlist(yRanges), ncol=2, byrow=TRUE)
#      
#      Xmin <- min(xRanges[,1], na.rm=TRUE)
#      Xmax <- max(xRanges[,2], na.rm=TRUE)
#      Ymin <- min(yRanges[,1], na.rm=TRUE)
#      Ymax <- max(yRanges[,2], na.rm=TRUE)
#      
#      xlim <- c(Xmin, Xmax)
#      ylim <- c(Ymin, Ymax)
#      
#      first <- TRUE
#      for(track in x@tracks)
#      {
#           #print(track@id)
#           #print(xlim)
#           #print(ylim)
#           if(first)
#           {
#                plot(track, slotY=slot, relY=rel, add=F, xlim=xlim, ylim=ylim, withTitle=FALSE, main='All Tracks', ...)
#                first = FALSE
#           }
#           else
#           {
#                plot(track, slotY=slot, relY=rel, add=T, ...)
#           }
#      }
# }
# ) 
# 
# setGeneric('getTrack', function(object, id) standardGeneric('getTrack')) 
# setMethod('getTrack', 'trackList', function(object, id)
# {
#      return(object@tracks[[as.character(id)]])
# })  
# 
# setGeneric('setTrack', function(object, track) standardGeneric('setTrack')) 
# setMethod('setTrack', 'trackList', function(object, track)
# {
#      track@parent <- object
#      object@tracks[[as.character(track@id)]] <- track
#      return(object)
# })  
# 
# setGeneric('removeTrack', function(object, id) standardGeneric('removeTrack')) 
# setMethod('removeTrack', 'trackList', function(object, id)
# {
#      object@tracks[[as.character(id)]] <- NULL
#      return(object)
# })  
# 
# setGeneric('getProp', function(object, trackPropFun=function(track){return(length(track))}) standardGeneric('getProp')) 
# setMethod('getProp', 'trackList', function(object, trackPropFun)
# {
#      ret <- list()
#      for(track in object@tracks)
#      {
#           ret[[as.character(track@id)]] <- trackPropFun(track)
#      }
#      return(ret)
# })
# 
# setGeneric('addTrackPoint', function(object, id, x, y, t, frame) standardGeneric('addTrackPoint')) 
# setMethod('addTrackPoint', 'trackList', function(object, id, x, y, t, frame)
# {
#      track <- getTrack(object, id)
#      if(is.null(track))
#      {
#           track <- new('track', id=id)
#      }
#      track <- addPoint(track, x=x, y=y, t=t, frame=frame)
#      track@parent <- object
#      object@tracks[[as.character(id)]] <- track
#      return(object)
# })
# 
# ##### Track #####
# 
# setClass('track',
#          representation(id='numeric', x='numeric', y='numeric', t='numeric', frames='numeric', vx='numeric', vy='numeric', positionFit='list', velocityFit='list', validFrames='numeric', parent='trackList'),
#          prototype(id=c(), x=c(), y=c(), t=c(), frames=c(), vx=c(), vy=c(), positionFit=list(), velocityFit=list(), validFrames=c(), parent=c()))
# 
# setMethod('initialize', 'track', function(.Object, id)
# {
#      .Object@id <- id
#      return(.Object)
# })
# 
# trackFromTrackROI <- function(id, start, pattern, parent)
# {
#      pairs <- strsplit(pattern,';')[[1]]
#      x <- numeric(0)
#      x0 <- numeric(0)
#      y <- numeric(0)
#      y0 <- numeric(0)
#      frames <- numeric(0)
#      first <- TRUE
#      for(pair in pairs)
#      {
#           if(first)
#           {
#                nums <- strsplit(start,',')[[1]]
#                x0 <- as.numeric(nums[1])
#                y0 <- as.numeric(nums[2])
#                i0 <- as.numeric(nums[3])
#                first <- FALSE
#                x <- append(x,x0)
#                y <- append(y,y0)
#                frames <- append(frames,i0)
#           }
#           else
#           {
#                nums <- strsplit(pair,',')[[1]]
#                x <- append(x, x0+as.numeric(nums[1]))
#                y <- append(y, y0+as.numeric(nums[2]))
#                frames <- append(frames,as.numeric(nums[3]))
#                # print(nums)
#           }
#      }
#      
# {
#      #      if(!is.na(samplingFreq))
#      #      {
#      #           t <- seq(ti, tf, 1/samplingFreq)[index + 1] # Index from JEX starts at 0, whereas first index in R starts at 1
#      #           if(length(which(is.na(t)) > 0))
#      #           {
#      #                stop("The sampling frequency, ti, and tf do not match the roi data. Check to see if they are correct. Last frame in track (", max(index), ") may exceed last frame suggested by sampling = (", length(seq(ti, tf, 1/samplingFreq)), ").")
#      #           }
#      #      }
#      #      else
#      #      {
#      #           if(max(index)[1] > length.out)
#      #           {
#      #                stop("The number of frames in this track (", max(index), ") exceeds the suggested number of frames of the dataset, as indicated by length.out (", length.out, ").")
#      #           }
#      #           t <- seq(ti, tf, length.out=length.out)[index + 1]
#      #           samplingFreq <- 1 / (t[2]-t[1])
#      #      }
# }
# 
# object <- new('track', id=as.numeric(as.character(id)), parent=parent)
# object@x <- x
# object@y <- y
# object@t <- tAll[frames + 1] # frames start at 0 whereas indices start at 1 in R
# object@frames <- frames
# object@positionFit <- list()
# object@velocityFit <- list()
# object@validFrames <- frames
# object <- calculateVelocities(object)
# return(object)
# }
# 
# setMethod('show', 'track', function(object)
# {
#      for(name in slotNames(object))
#      {
#           if(!(name %in% c()))#c('x','y','t','vx','vy')))
#           {
#                cat("\n", name, " = ", paste(slot(object, name), collapse=","), sep="")
#           }
#           else
#           {
#                cat("\n", name, " = ...", sep="")
#           }
#      }
# })
# 
# setMethod('plot', 'track', function(x, slotX='t', slotY='x', relX=FALSE, relY=TRUE, add=FALSE, withTitle=TRUE, col='black', lwd=1, lty=1, xlab=slotX, ylab=slotY, ...)
# {
#      xData <- getSlot(x, slotX, relX)
#      yData <- getSlot(x, slotY, relY)
#      if(relX)
#      {
#           xlab <- paste(xlab, ' - mean(', xlab, ')', sep='')
#      }
#      if(relY)
#      {
#           ylab <- paste(ylab, ' - mean(', ylab, ')', sep='')
#      }
#      
#      if(is.na(xData) || is.na(yData))
#      {
#           print(paste('Nothing to plot for track ', x@id, '.', sep=''))
#      }
#      else
#      {
#           if(add)
#           {
#                
#                lines(xData, yData, col=col, lwd=lwd, lty=lty, ...)
#           }
#           else
#           {
#                if(withTitle)
#                {
#                     plot(xData, yData, type='l', col=col, lwd=lwd, lty=lty, xlab=xlab, ylab=ylab, main=as.character(x@id), ...)
#                }
#                else
#                {
#                     plot(xData, yData, type='l', col=col, lwd=lwd, lty=lty, xlab=xlab, ylab=ylab, ...)
#                }
#           }
#      }
# })
# 
# setMethod('length', 'track', function(x)
# {
#      numel(x@x)
# })
# 
# setMethod('range', 'track', function(x, slot, rel=FALSE)
# {
#      range(getSlot(x, slot, rel=rel))
# })
# 
# setGeneric('getSlot', function(object, slot, rel=FALSE, validOnly=FALSE) standardGeneric('getSlot')) 
# setMethod('getSlot', 'track', function(object, slot, rel, validOnly)
# {
#      if(rel)
#      {
#           temp <- slot(object, slot) - mean(slot(object, slot))
#           if(validOnly)
#           {
#                temp <- temp[object@validIndices]
#           }
#      }
#      else
#      {
#           temp <- slot(object, slot)
#           if(validOnly)
#           {
#                temp <- temp[object@validIndices]
#           }
#      }
#      return(temp)
# })
# 
# setGeneric('addPoint', function(object, x, y, t, frame) standardGeneric('addPoint')) 
# setMethod('addPoint', 'track', function(object, x, y, t, frame)
# {
#      object@x <- c(object@x, x)
#      object@y <- c(object@y, y)
#      object@t <- c(object@t, t)
#      object@frames <- c(object@frames, frame)
#      return(object)
# })
# 
# setGeneric('calculateVelocities', function(object) standardGeneric('calculateVelocities')) 
# setMethod('calculateVelocities', 'track', function(object)
# {
#      object@vx <- getDerivative(object@x, object@t)
#      object@vy <- getDerivative(object@y, object@t)
#      return(object)
# })
# 

##### Track Fitting #####

# sse <- function(sin=TRUE, amplitude=1, phaseShift=0, offset=0, fi=0.2, ff=0.01, ti=0, tf=500, samplingFreq=0.5, data, indices)
# {
#      predicted <- getSweep(sin=sin, amplitude=amplitude, phaseShift=phaseShift, offset=offset, fi=fi, ff=ff, ti=ti, tf=tf, samplingFreq=samplingFreq, indices=indices)
#      return(sum((data-predicted$v)^2))
# }

sseTrack <- function(object, slot='vx', amplitude=50, phaseShift=0, offset=0, frames=-1)
{
     if(frames[1] < 0)
     {
          frames <- slot(object, 'frames')
     }
     predicted <- getSweep(sin=object@parent@sin, phaseShift=phaseShift, amplitude=amplitude, offset=offset, fi=object@parent@fi, ff=object@parent@ff, tAll=object@parent@tAll, frames=object@frames)
     data <- slot(object, slot)
     indicesToGet <- which(object@frames == frames)
     result <- sum((data[indicesToGet]-predicted$v)^2)
     #print(paste(result,first(frames), ' to ', last(frames)))
     if(is.na(result))
     {
          browser()
     }
     return(result)
     #sse(sin=object@sin, phaseShift=phaseShift, amplitude=amplitude, offset=offset, fi=object@fi, ff=object@ff, ti=object@ti, tf=object@tf, samplingFreq=object@samplingFreq, data=slot(object, slot), indices=indices)
}

getTrackFitFixedPhase <- function(object, slot='vx', initialAmplitude=40, amplitudeLimits=c(0,500), phaseShift=0, initialOffset=0, offsetLimits=c(0, 10000), guess=NULL)
{
     if(is.null(guess))
     {
          bestFit <- optim(par=c(amplitude=initialAmplitude, offset=initialOffset), 
                           function(par, object, slot, phaseShift, offset){sseTrack(object=object, slot=slot, amplitude=par['amplitude'], phaseShift=phaseShift, offset=par['offset'])}, 
                           method='L-BFGS-B', 
                           lower=c(min(amplitudeLimits), min(offsetLimits)), 
                           upper=c(max(amplitudeLimits), max(offsetLimits)), 
                           control=list(trace=0), 
                           object=object,
                           slot=slot,
                           phaseShift=phaseShift)
     }else
     {
          bestFit <- optim(par=c(amplitude=as.numeric(guess['amplitude']), offset=as.numeric(guess['offset'])), 
                           function(par, object, slot, phaseShift, offset){sseTrack(object=object, slot=slot, amplitude=par['amplitude'], phaseShift=phaseShift, offset=par['offset'])}, 
                           method='L-BFGS-B', 
                           lower=c(min(amplitudeLimits), min(offsetLimits)), 
                           upper=c(max(amplitudeLimits), max(offsetLimits)), 
                           control=list(trace=0), 
                           object=object,
                           slot=slot,
                           phaseShift=phaseShift)
     } 
     return(list(par=c(phaseShift=as.numeric(phaseShift), amplitude=as.numeric(bestFit$par['amplitude']), offset=as.numeric(bestFit$par['offset'])), fit=bestFit))
}

getWindowIndices <- function(sweep, window, i)
{
     # Is startTimeI a vector index or a frame number?
     startLine <- sweep$line[i]
     endI <- which(sweep$line > startLine + window*as.numeric(sweep$A))[1]-1 # window*sweep$A represents window number of pi/2 sections of signal
     if(endI > i)
     {
          endI <- i
     }
     if(endI > length(sweep$t))
     {
          return(NULL)
     }
     return(i:endI)
     
     #      startInfTime <- which(inflections >= startTime)[1]
     #      endInfTime <- startInfTime + window - 1
     #      #      validInflections <- inflections[i:(i + window - 1)]
     #      #      startInfTime <- validInflections[1]
     #      #      endInfTime <- validInflections[length(validInflections)]
     #      #      startTimeI <- which(times >= startInfTime)[1]
     #      endTimeI <- which(times > endInfTime)[1] - 1
     #      #      if(is.na(startTimeI) | is.na(endTimeI))
     #      #      {
     #      #           browser()
     #      #      }
     #      return(startTimeI:endTimeI)
}

getTrackRollFit <- function(object, window=7, slot='vx', amplitudeLimits=c(0.25,500), offset=0, guess=NULL)
{
     if(is.null(guess))
     {
          guess <- getTrackParamGuess(object)
     }
     
     amplitude <- guess['amplitude']
     indices <- object@index # Index means frames in JEX
     fullSweep <- getTrackSweep(track=object, amplitude=amplitude)
     inflections <- fullSweep$inflections
     times <- fullSweep$t
     line <- fullSweep$line
     
     start <- 1
     end <- length(times)
     
     currentFit <- NULL
     currentAmplitude <- amplitude
     results <- list()
     browser()
     for(i in start:end)
     {
          if(!is.null(currentFit))
          {
               currentAmplitude <- currentFit$par['amplitude']
          }
          windowI <- getWindowIndices(fullSweep, window=window, i=i)
          if(is.null(windowI))
          {
               browser()
               break
          }
          print(paste('Starting next timepoint', start, end, i))
          currentFit <- getTrackFitAmplitudeOnly(object, slot=slot, initialAmplitude=currentAmplitude, amplitudeLimits=amplitudeLimits, phaseShift=object@phaseShift, offset=offset, indices=object@index[windowI], guess=FALSE)
          results[as.character(indices[i + (windowI%/%2) + 1])] <- list(currentFit)
     }
     return(results)
}

getTrackFitAmplitudeOnly <- function(object, slot='vx', initialAmplitude=40, amplitudeLimits=c(0,500), phaseShift=0, offset=0, indices=-1, guess=FALSE)
{
     if(indices[1] < 0)
     {
          indices <- object@index
     }
     if(guess[1]==TRUE)
     {
          initialAmplitude <- getTrackParamGuess(object)['amplitude']
     }
     
     ### Out this into calc of bestFit so we can specify indices
     #sse(sin=object@sin, phaseShift=phaseShift, amplitude=amplitude, offset=offset, fi=object@fi, ff=object@ff, ti=object@ti, tf=object@tf, samplingFreq=object@samplingFreq, data=slot(object, slot), indices=indices)
     
     bestFit <- optim(par=c(amplitude=as.numeric(initialAmplitude)), 
                      function(par, object, slot, phaseShift, offset, indices){sseTrack(object=object, slot=slot, phaseShift=phaseShift, amplitude=par['amplitude'], offset=offset, indices=indices)}, 
                      method='L-BFGS-B', 
                      lower=c(min(amplitudeLimits)), 
                      upper=c(max(amplitudeLimits)), 
                      control=list(trace=0), 
                      object=object,
                      slot=slot,
                      phaseShift=phaseShift,
                      offset=offset,
                      indices=indices)
     return(bestFit)
}

getTrackFitAll <- function(object, slot='vx', initialAmplitude=40, amplitudeLimits=c(0,500), initialPhaseShift=0, phaseShiftLimits=c(-pi, pi), initialOffset=0, offsetLimits=c(0,10000), guess=TRUE)
{
     if(guess[1]==TRUE)
     {
          trackGuess <- getTrackParamGuess(object)
          bestFit <- optim(par=c(phaseShift=as.numeric(trackGuess['phaseShift']), amplitude=as.numeric(trackGuess['amplitude']), offset=as.numeric(trackGuess['offset'])), 
                           function(par, object, slot){sseTrack(object=object, slot=slot, phaseShift=par['phaseShift'], amplitude=par['amplitude'], offset=par['offset'])}, 
                           method='L-BFGS-B', 
                           lower=c(min(phaseShiftLimits), min(amplitudeLimits), min(offsetLimits)), 
                           upper=c(max(phaseShiftLimits), max(amplitudeLimits), max(offsetLimits)), 
                           control=list(trace=0), 
                           object=object,
                           slot=slot)
          return(bestFit)
     }else
     {
          bestFit <- optim(par=c(phaseShift=initialPhaseShift, amplitude=initialAmplitude, offset=initialOffset), 
                           function(par, object, slot){sseTrack(object=object, slot=slot, phaseShift=par['phaseShift'], amplitude=par['amplitude'], offset=par['offset'])}, 
                           method='L-BFGS-B', 
                           lower=c(min(phaseShiftLimits), min(amplitudeLimits), min(offsetLimits)), 
                           upper=c(max(phaseShiftLimits), max(amplitudeLimits), max(offsetLimits)), 
                           control=list(trace=0), 
                           object=object,
                           slot=slot)
          return(bestFit)
     }
     return(list(par=c(phaseShift=as.numeric(bestFit$par['phaseShift']), amplitude=as.numeric(bestFit$par['amplitude']), offset=as.numeric(bestFit$par['offset'])), fit=bestFit))
}

getTrackParamGuess <- function(object, slot='x', amplitudeLimits=c(0,500), phaseShiftLimits=c(-pi, pi), offsetLimits=c(-10000,10000))
{
     offset <- mean(slot(object,slot))
     amplitude <- mean(abs(getSlot(object,slot,rel=T)))
     #     amplitude <- (amplitude[2]-amplitude[1])/2
     phaseShift <- 0
     return(c(phaseShift=phaseShift, amplitude=amplitude, offset=offset))
}

plotFit <- function(track, fit=NULL, ...)
{
     if(is.null(fit))
     {
          print('Generating fit to plot...')
          fit <- getTrackFitAll(track)
     }
     plot(track, slotY='vx', relY=FALSE, add=FALSE, ...)
     guessPar <- getTrackParamGuess(track)
     guess <- getSweep(sin=track@parent@sin, amplitude=guessPar['amplitude'], phaseShift=guessPar['phaseShift'], offset=guessPar['offset'], fi=track@parent@fi, ff=track@parent@ff, tAll=track@parent@tAll, frames=track@frames)
     #lines(guess$t, guess$v, col='blue')
     predicted <- getSweep(sin=track@parent@sin, amplitude=fit$par['amplitude'], phaseShift=fit$par['phaseShift'], offset=fit$par['offset'], fi=track@parent@fi, ff=track@parent@ff, tAll=track@parent@tAll, frames=track@frames)
     lines(predicted$t, predicted$v, col='red')
     #lines(predicted$t, predicted$x-mean(predicted$x), col='green')
     points(predicted$inflections, rep(0,length(predicted$inflections)))
     points(track@t[track@validFrames], track@vx[track@validFrames], col='blue')
     #print("Guess:")
     #print(guessPar)
     print("Fit:")
     print(fit$par)
}

##### TrackList Fitting #####

sseTrackList <- function(object, slot, phaseShift, initialAmplitude=50, amplitudeLimits=c(0,500), initialOffset=0, offsetLimits=c(0,10000), guess=TRUE)
{
     print(paste("Optimizing phaseshift: ", phaseShift, sep=''))
     sseRet <- 0
     count <- 1
     total <- length(object@tracks)
     if(guess==TRUE)
     {
          for(track in object@tracks)
          {
               bestFit <- getTrackFitFixedPhase(object=track, slot=slot, initialAmplitude=initialAmplitude, amplitudeLimits=amplitudeLimits, phaseShift=phaseShift, initialOffset=initialOffset, offsetLimits=offsetLimits, guess=getTrackParamGuess(track))
               sseRet <- sseRet + bestFit$fit$value #bestFit$value gives the value of the sseTrack at the best fit param
               print(paste(count, " of ", total, " - ID:", track@id, ", phaseShift:", bestFit$par['phaseShift'], ", amplitude:", bestFit$par['amplitude'], ", offset:", bestFit$par['offset']))
               count <- count + 1
          }
     }
     else
     {
          for(track in object@tracks)
          {
               bestFit <- getTrackFitFixedPhase(object=track, initialAmplitude=initialAmplitude, amplitudeLimits=amplitudeLimits, phaseShift=phaseShift, initialOffset=initialOffset, offsetLimits=offsetLimits) # Don't use the getTrackParamGuessFunction
               sseRet <- sseRet + bestFit$fit$value #bestFit$value gives the value of the sseTrack at the best fit param
               print(paste(count, " of ", total, " - ID:", track@id, ", phaseShift:", bestFit$par['phaseShift'], ", amplitude:", bestFit$par['amplitude'], ", offset:", bestFit$par['offset']))
               count <- count + 1
          }
     }
     return(sseRet)
}

getTrackListPhaseShiftGuess <- function(object, slot='x')
{
     ranges <- numeric(0)
     ids <- numeric(0)
     for(track in object@tracks)
     {
          range <- range(track, 'x')
          ranges <- c(ranges, range[2]-range[1])
          ids <- c(ids, track@id)
     }
     testerIndex <- which(ranges == max(ranges))[1]
     testerId <- ids[testerIndex]
     testerTrack <- getTrack(object, testerId)
     fit <- getTrackFitAll(testerTrack)
     plotFit(testerTrack, fit)
     return(as.numeric(fit$par['phaseShift']))
}

setTrackListPhaseShift <- function(object, slot='vx', initialPhaseShift=0, phaseShiftLimits=c(-pi, pi), guess=TRUE)
{
     if(guess==TRUE)
     {
          initialPhaseShift <- getTrackListPhaseShiftGuess(object, slot)
     }
     bestFit <- optim(par=c(phaseShift=initialPhaseShift),
                      function(par, object, slot) {sseTrackList(object=object, slot=slot, phaseShift=par['phaseShift'], guess=TRUE)},
                      method='L-BFGS-B',
                      lower=min(phaseShiftLimits),
                      upper=max(phaseShiftLimits),
                      control=list(trace=3),
                      object=object,
                      slot=slot)
     object@phaseShift <- bestFit$par['phaseShift']
     for(i in 1:length(object@tracks))
     {
          object@tracks[[i]]@phaseShift <- bestFit$par['phaseShift']
     }
     return(object)
}


##### Butterworth Filters #####

getFilter <- function(samplingFreq=0.5, passFreq=0.01, stopFreq=0.005, Rp=0.5, Rs=30)
{
     return(butter(buttord(Wp=passFreq/(samplingFreq/2), Ws=stopFreq/(samplingFreq/2), Rp=Rp, Rs=Rs)))
}

compare <- function(stopFreq=0.005, passFreq=0.01, samplingFreq=0.5, atFreq=0.01)
{
     t <- seq(0,5*1/passFreq,1/samplingFreq)
     bf <- getFilter(stopFreq=stopFreq, passFreq=passFreq, samplingFreq=samplingFreq)
     raw <- sin(2*pi*atFreq*t)
     filtered <- as.numeric(filter(bf, raw))
     plot(t, raw, type='l', xlab='t', ylab='amplitude')
     lines(t, filtered, col='red')
}

plotFreqResponse <- function(filter, samplingFreq=0.5, passFreq=0.01, stopFreq=0.005, Rp=0.5, Rc=30)
{
     plot(c(0, passFreq, passFreq, 0, 0), c(0, 0, -Rp, -Rp, 0),
          type = "l", xlab = "Frequency (Hz)", ylab = "Attenuation (dB)",
          col = "green", ylim = c(-Rc-Rp,0), xlim = c(0,4*passFreq))
     lines(c(0, stopFreq, stopFreq, 0, 0), c(-Rc-Rp, -Rc-Rp, -Rc, -Rc, -Rc-Rp),
           col = "red", ylim = c(-Rc-Rp,0), xlim = c(0,2000))
     hf <- freqz(filter, n=5000, Fs = samplingFreq)
     lines(hf$f, 20*log10(abs(hf$h)))
}

filterWandering <- function(object, slot='x', stopFreq=0.005, passFreq=0.01)
{
     bwf = getFilter(samplingFreq=object@samplingFreq, stopFreq=stopFreq, passFreq=passFreq)
     ret <- filter(bwf, slot(object, slot))
     newTrack <- object
     slot(newTrack, slot) <- as.numeric(ret)
     return(newTrack)
}