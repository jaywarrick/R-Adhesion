
library(signal)
library(zoo)
source("http://bioconductor.org/biocLite.R")
biocLite('GraphAlignment')
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

trackFrameFilter <- function(track, startBy=0, endBy=1000000)
{
    startFrame <- track@index[1]
    endFrame <- track@index[length(track)]
    if(startFrame >= startBy & endFrame <= endBy)
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

getSweep <- function(sin=TRUE, amplitude=1, phaseShift=0, offset=0, fi=0.2, ff=0.01, ti=0, tf=500, samplingFreq=0.5, indices=-1, guess=NULL)
{
    t <- seq(ti, tf, 1/samplingFreq)
    if(indices[1] < 0 || indices[1] > length(t))
    {
        indices <- (1:length(t)) - 1 #Subtract 1 to be consistent with indices from JEX, which start at 0 instead of 1
    }
    t <- t[indices + 1]
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
    f <- fi*(ff/fi)^(t/tf)
    line <- 4*A*fi*(  ((ff/fi)^(t/tf)*tf)/log(ff/fi)-tf/log(ff/fi)  ) - A*(phi/(pi/2)) # the last term shifts the line up and down, essentially dictating where the "sections" of the line are placed (i.e., the initial phase of the triangle wave)
    line1 <- 4*1*fi*(  ((ff/fi)^(t/tf)*tf)/log(ff/fi)-tf/log(ff/fi)  ) - 1*(phi/(pi/2)) # WITH DUMMY AMPLITUDE OF 1 in case A is 0 then we have something to calculate sections with that isn't a simple flat line
    x <- b + (line %% A) #*sign(fSweep$x)
    # Fix each pi/2 section of x
    sections <- (line1 %/% 1) %% 4
    x[sections==1] <- -1*(x[sections==1]-1)
    x[sections==2] <- -1*(x[sections==2])
    x[sections==3] <- (x[sections==3]-1)
    
    # Find pi/2 intervals on frequency sweep
    n <- seq(0,1000,1)
    suppressWarnings(inflections <- (log((log(ff/fi)*phi)/(2*fi*pi*tf)+(log(ff/fi)*n)/(4*fi*tf)+1)*tf)/log(ff/fi))
    inflections <- inflections[!is.nan(inflections)]
    inflectionNums <- (  seq(0,length(inflections)-1,1) + (  (-A*(phi/(pi/2))) %/% A %% 4  )  ) %% 4
    
    predicted <- A*sin(2*pi*((fi*(-1+2^(R*t)))/(R*log(2))) - phi) + b
    
    #browser()
    #print((-A*(phi/(pi/2))) %/%A %% 4)
    if(sin)
    {
        return(list(t=t, x=predicted, v=getDerivative(x=predicted, t=t), sections=sections, inflections=inflections, inflectionNums=inflectionNums))
    }
    else
    {
        return(list(t=t, x=x, v=getDerivative(x=x, t=t), sections=sections, inflections=inflections, inflectionNums=inflectionNums))
    }
}

getTrackSweep <- function(track, slot='vx', sin=TRUE, amplitude=1, offset=mean(slot(track,slot)))
{
    return(getSweep(sin=sin, amplitude=amplitude, phaseShift=track@phaseShift, offset=offset, fi=track@fi, ff=track@ff, ti=track@ti, tf=track@tf, samplingFreq=track@samplingFreq, indices=track@index, guess=NULL))
}

##### Track #####

setClass('track', representation(id='numeric', x='numeric', y='numeric', t='numeric', index='numeric', vx='numeric', vy='numeric', fi='numeric', ff='numeric', ti='numeric', tf='numeric', samplingFreq='numeric', phaseShift='numeric', validTimeIndices='numeric', amplitude='numeric'), prototype(id=0, x=0, y=0, t=0, index=0, vx=0, vy=0, fi=1, ff=0.1, ti=0, tf=1000, samplingFreq=1, phaseShift=0, validTimeIndices=0))

library(foreign)
tracksFile <- read.arff(path)
tracksFile2 <- reorganizeTable(tracksFile, nameCol='Metadata')
row=6076
duh <- new('track',tracksFile2[row,'Track'], tracksFile2[row,'polygonPts'], tracksFile2[row,'patternPts'], fi=0.2, ff=0.01, ti=0, tf=500, length.out=2361)

setMethod('initialize', 'track', function(.Object, id='0', start='0,0,0', pattern='0,0,0', fi, ff, ti, tf, samplingFreq=NA, length.out=1)
{
    .Object@id <- as.numeric(as.character(id))
    pairs <- strsplit(pattern,';')[[1]]
    x <- numeric(0)
    x0 <- numeric(0)
    y <- numeric(0)
    y0 <- numeric(0)
    index <- numeric(0)
    first <- TRUE
    for(pair in pairs)
    {
        if(first)
        {
            nums <- strsplit(start,',')[[1]]
            x0 <- as.numeric(nums[1])
            y0 <- as.numeric(nums[2])
            i0 <- as.numeric(nums[3])
            first <- FALSE
            x <- append(x,x0)
            y <- append(y,y0)
            index <- append(index,i0)
        }
        else
        {
            nums <- strsplit(pair,',')[[1]]
            x <- append(x, x0+as.numeric(nums[1]))
            y <- append(y, y0+as.numeric(nums[2]))
            index <- append(index,as.numeric(nums[3]))
            # print(nums)
        }
    }
    
    if(!is.na(samplingFreq))
    {
        t <- seq(ti, tf, 1/samplingFreq)[index + 1] # Index from JEX starts at 0, whereas first index in R starts at 1
        if(length(which(is.na(t)) > 0))
        {
            stop("The sampling frequency, ti, and tf do not match the roi data. Check to see if they are correct. Last frame in track (", max(index), ") may exceed last frame suggested by sampling = (", length(seq(ti, tf, 1/samplingFreq)), ").")
        }
    }
    else
    {
        if(max(index)[1] > length.out)
        {
            stop("The number of frames in this track (", max(index), ") exceeds the suggested number of frames of the dataset, as indicated by length.out (", length.out, ").")
        }
        t <- seq(ti, tf, length.out=length.out)[index + 1]
        samplingFreq <- 1 / (t[2]-t[1])
    }
    
    
    #     print(x)
    #     print(y)
    #     print(index)
    #     print(t)
    .Object@x <- x
    .Object@vx <- getDerivative(x, t)
    .Object@y <- y
    .Object@vy <- getDerivative(y, t)
    .Object@t <- t
    .Object@index <- index
    .Object@fi <- fi
    .Object@ff <- ff
    .Object@ti <- ti
    .Object@tf <- tf
    .Object@samplingFreq <- samplingFreq
    .Object
}
)

setMethod('show', 'track', function(object)
{
    print(paste('Track:', object@id))
    print(paste('[fi=', object@fi, ', ff=', object@ff, ', ti=', object@ti, ', tf=', object@tf, ', samplingFreq=', object@samplingFreq, ', phaseShift=', object@phaseShift, sep=''))
    for(i in 1:length(object@x))
    {
        print(paste('[',object@x[i],',',object@y[i],',',object@t[i],']',sep=''))
    }
    print(object@validTimeIndices)
}
)

setMethod('plot', 'track', function(x, slotX='t', slotY='x', relX=FALSE, relY=TRUE, add=FALSE, withTitle=TRUE, col='black', lwd=1, lty=1, xlab=slotX, ylab=slotY, ...)
{
    xData <- getSlot(x, slotX, relX)
    yData <- getSlot(x, slotY, relY)
    if(relX)
    {
        xlab <- paste(xlab, ' - mean(', xlab, ')', sep='')
    }
    if(relY)
    {
        ylab <- paste(ylab, ' - mean(', ylab, ')', sep='')
    }
    
    if(is.na(xData) || is.na(yData))
    {
        print(paste('Nothing to plot for track ', x@id, '.', sep=''))
    }
    else
    {
        if(add)
        {
            
            lines(xData, yData, col=col, lwd=lwd, lty=lty, ...)
        }
        else
        {
            if(withTitle)
            {
                plot(xData, yData, type='l', col=col, lwd=lwd, lty=lty, xlab=xlab, ylab=ylab, main=as.character(x@id), ...)
            }
            else
            {
                plot(xData, yData, type='l', col=col, lwd=lwd, lty=lty, xlab=xlab, ylab=ylab, ...)
            }
        }
    }
}
)

setMethod('length', 'track', function(x)
{
    length(x@x)
}
)

setMethod('range', 'track', function(x, slot, rel=FALSE)
{
    range(getSlot(x, slot, rel=rel))
}
)

setGeneric('getSlot', function(object, slot, rel=FALSE, validOnly=FALSE) standardGeneric('getSlot')) 
setMethod('getSlot', 'track', function(object, slot, rel, validOnly)
{
    if(rel)
    {
        temp <- slot(object, slot) - mean(slot(object, slot))
        if(validOnly)
        {
            temp <- temp[object@validIndices]
        }
    }
    else
    {
        temp <- slot(object, slot)
        if(validOnly)
        {
            temp <- temp[object@validIndices]
        }
    }
    return(temp)
}
)

##### TrackList #####

setClass('trackList', representation(tracks='list', fi='numeric', ff='numeric', ti='numeric', tf='numeric', samplingFreq='numeric', phaseShift='numeric'))

setMethod('initialize', 'trackList', function(.Object, template=NULL, tracksFilePath, fi, ff, ti, tf, samplingFreq=NA, length.out)
{
    # Determine samplingFreq if necessary
    if(is.na(samplingFreq))
    {
        temp <- seq(ti, tf, length.out=length.out)
        samplingFreq <- 1 / (temp[2]-temp[1])
    }
    
    # Initialize this way if given a trackList as a template
    if(is.null(template))
    {
        library(foreign)
        tracksFile <- read.arff(path)
        tracksFile2 <- reorganizeTable(tracksFile, nameCol='Metadata')
        
        duh <- list()
        for(row in 1:length(tracksFile2[,1]))
        {
            duh[as.character(tracksFile2[row,'Track'])] <- new('track',tracksFile2[row,'Track'], tracksFile2[row,'polygonPts'], tracksFile2[row,'patternPts'], fi=fi, ff=ff, ti=ti, tf=tf, samplingFreq=samplingFreq, length.out=length.out)
        }
        
        .Object@tracks <- duh
        .Object@fi <- fi
        .Object@ff <- ff
        .Object@ti <- ti
        .Object@tf <- tf
        .Object@samplingFreq <- samplingFreq
        .Object@phaseShift <- 0
        .Object
    }
    else # initialize this way if using the other possible parameters
    {
        .Object@fi <- template@fi
        .Object@ff <- template@ff
        .Object@ti <- template@ti
        .Object@tf <- template@tf
        .Object@samplingFreq <- template@samplingFreq
        .Object@phaseShift <- template@phaseShift
        .Object
    }
}
)

setMethod('length', 'trackList', function(x)
{
    length(x@tracks)
}
)

setMethod('plot', 'trackList', function(x, slot='vx', rel=FALSE)
{

    xRanges <- getProp(x, function(track){range(track@t)})
    xRanges <- matrix(unlist(xRanges), ncol=2, byrow=TRUE)
    yRanges <- getProp(x, function(track){range(slot(track, slot))})
    yRanges <- matrix(unlist(yRanges), ncol=2, byrow=TRUE)

    Xmin <- min(xRanges[,1], na.rm=TRUE)
    Xmax <- max(xRanges[,2], na.rm=TRUE)
    Ymin <- min(yRanges[,1], na.rm=TRUE)
    Ymax <- max(yRanges[,2], na.rm=TRUE)
    
    xlim <- c(Xmin, Xmax)
    ylim <- c(Ymin, Ymax)
    
    first <- TRUE
    for(track in x@tracks)
    {
        print(track@id)
        print(xlim)
        print(ylim)
        if(first)
        {
            plot(track, slotY=slot, relY=rel, add=F, xlim=xlim, ylim=ylim, withTitle=FALSE, main='All Tracks')
            first = FALSE
        }
        else
        {
            plot(track, slotY=slot, relY=rel, add=T)
        }
    }
}
)

setGeneric('getTrack', function(object, id) standardGeneric('getTrack')) 
setMethod('getTrack', 'trackList', function(object, id)
{
    return(object@tracks[[as.character(id)]])
}
)  

setGeneric('addTrack', function(object, id) standardGeneric('addTrack')) 
setMethod('addTrack', 'trackList', function(object, id)
{
    for(track in object@tracks)
    {
        if(track@id == id)
        {
            return(track)
        }
    }
    return(object)
}
)  

setGeneric('removeTrack', function(object, id) standardGeneric('removeTrack')) 
setMethod('removeTrack', 'trackList', function(object, id)
{
    for(track in object@tracks)
    {
        if(track@id == id)
        {
            return(track)
        }
    }
    return(object)
}
)  

setGeneric('getProp', function(object, trackPropFun=function(track){return(length(track))}) standardGeneric('getProp')) 
setMethod('getProp', 'trackList', function(object, trackPropFun)
{
    ret <- list()
    for(track in object@tracks)
    {
        ret[[as.character(track@id)]] <- trackPropFun(track)
    }
    return(ret)
}
)

setGeneric('sortTrackList', function(object, trackPropFun=function(object){return(length(object@tracks))}, decreasing = TRUE) standardGeneric('sortTrackList')) 
setMethod('sortTrackList', 'trackList', function(object, trackPropFun, decreasing)
{
    ret <- getProp(object, trackPropFun=trackPropFun)
    sorted <- sort(unlist(ret), index.return=T, decreasing=decreasing)
    newTracks <- object@tracks[sorted$ix]
    object@tracks <- newTracks
    return(object)
}
)

##### Track Fitting #####

sse <- function(amplitude=1, phaseShift=0, offset=0, fi=0.2, ff=0.01, ti=0, tf=500, samplingFreq=0.5, data, indices)
{
    predicted <- getSweep(sin=sin, amplitude=amplitude, phaseShift=phaseShift, offset=offset, fi=fi, ff=ff, ti=ti, tf=tf, samplingFreq=samplingFreq, indices=indices)
    return(sum((data-predicted$v)^2))
}

sseTrack <- function(object, sin=FALSE, slot='vx', amplitude=50, phaseShift=0, offset=0, indices=-1)
{
    if(indices[1] < 0)
    {
        indices <- slot(object, 'index')
    }
    sse(sin=sin, phaseShift=phaseShift, amplitude=amplitude, offset=offset, fi=object@fi, ff=object@ff, ti=object@ti, tf=object@tf, samplingFreq=object@samplingFreq, data=slot(object, slot), indices=indices)
}

getTrackFitFixedPhase <- function(object, sin=FALSE, slot='vx', initialAmplitude=40, amplitudeLimits=c(0,500), phaseShift=0, initialOffset=0, offsetLimits=c(0, 10000), guess=NULL)
{
    if(is.null(guess))
    {
        bestFit <- optim(par=c(amplitude=initialAmplitude, offset=initialOffset), 
                         function(par, object, sin, slot, phaseShift, offset){sseTrack(object=object, sin=sin, slot=slot, amplitude=par['amplitude'], phaseShift=phaseShift, offset=par['offset'])}, 
                         method='L-BFGS-B', 
                         lower=c(min(amplitudeLimits), min(offsetLimits)), 
                         upper=c(max(amplitudeLimits), max(offsetLimits)), 
                         control=list(trace=0), 
                         object=object, 
                         sin=sin,
                         slot=slot,
                         phaseShift=phaseShift)
    }else
    {
        bestFit <- optim(par=c(amplitude=as.numeric(guess['amplitude']), offset=as.numeric(guess['offset'])), 
                         function(par, object, sin, slot, phaseShift, offset){sseTrack(object=object, sin=sin, slot=slot, amplitude=par['amplitude'], phaseShift=phaseShift, offset=par['offset'])}, 
                         method='L-BFGS-B', 
                         lower=c(min(amplitudeLimits), min(offsetLimits)), 
                         upper=c(max(amplitudeLimits), max(offsetLimits)), 
                         control=list(trace=0), 
                         object=object, 
                         sin=sin,
                         slot=slot,
                         phaseShift=phaseShift)
    } 
    return(list(par=c(phaseShift=as.numeric(phaseShift), amplitude=as.numeric(bestFit$par['amplitude']), offset=as.numeric(bestFit$par['offset'])), fit=bestFit))
}

getTrackRollFit <- function(object, window=5, slot='vx', sin=FALSE, amplitudeLimits=c(0,500), offset=0, guess=NULL)
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
    
    getWindowIndicies <- function(times, indicies, inflections, window, i)
    {
        validInflections <- inflections[i:(i + window - 1)]
        startInfTime <- validInflections[1]
        endInfTime <- validInflections[length(validInflections)]
        startTimeI <- which(times >= startInfTime)[1]
        endTimeI <- which(times >= endInfTime)[1]
        return(indicies[startTimeI:endTimeI])
    }
    
    start <- 1
    end <- length(indices)-(window-1)
    
    currentFit <- NULL
    currentAmplitude <- amplitude
    results <- list()
    for(i in start:end)
    {
        if(!is.null(currentFit))
        {
            currentAmplitude <- currentFit$par['amplitude']
        }
        currentFit <- getTrackFitAmplitudeOnly(object, sin=sin, slot=slot, initialAmplitude=currentAmplitude, amplitudeLimits=amplitudeLimits, phaseShift=object@phaseShift, offset=offset, indices=-1, guess=FALSE)
        results[as.character(i)] <- currentFit
    }
}

getTrackFitAmplitudeOnly <- function(object, sin=FALSE, slot='vx', initialAmplitude=40, amplitudeLimits=c(0,500), phaseShift=0, offset=0, guess=FALSE)
{
    if(indices < 0)
    {
        indices <- object@index
    }
    if(guess[1]==TRUE)
    {
        amplitude <- getTrackParamGuess(object)['amplitude']
    }
    
    ### Out this into calc of bestFit so we can specify indices
    sse(sin=sin, phaseShift=phaseShift, amplitude=amplitude, offset=offset, fi=object@fi, ff=object@ff, ti=object@ti, tf=object@tf, samplingFreq=object@samplingFreq, data=slot(object, slot), indices=indices)
    
    bestFit <- optim(par=c(amplitude=initialAmplitude), 
                     function(par, object, sin, slot, phaseShift, offset){sseTrack(object=object, sin=sin, slot=slot, phaseShift=phaseShift, amplitude=par['amplitude'], offset=offset)}, 
                     method='L-BFGS-B', 
                     lower=c(min(phaseShiftLimits), min(amplitudeLimits), min(offsetLimits)), 
                     upper=c(max(phaseShiftLimits), max(amplitudeLimits), max(offsetLimits)), 
                     control=list(trace=0), 
                     object=object, 
                     sin=sin,
                     slot=slot,
                     phaseShift=phaseShift,
                     offset=offset)
    return(bestFit)
}

getTrackFitAll <- function(object, sin=FALSE, slot='vx', initialAmplitude=40, amplitudeLimits=c(0,500), initialPhaseShift=0, phaseShiftLimits=c(-pi, pi), initialOffset=0, offsetLimits=c(0,10000), guess=TRUE)
{
    if(guess[1]==TRUE)
    {
        trackGuess <- getTrackParamGuess(object)
        bestFit <- optim(par=c(phaseShift=as.numeric(trackGuess['phaseShift']), amplitude=as.numeric(trackGuess['amplitude']), offset=as.numeric(trackGuess['offset'])), 
                         function(par, object, sin, slot){sseTrack(object=object, sin=sin, slot=slot, phaseShift=par['phaseShift'], amplitude=par['amplitude'], offset=par['offset'])}, 
                         method='L-BFGS-B', 
                         lower=c(min(phaseShiftLimits), min(amplitudeLimits), min(offsetLimits)), 
                         upper=c(max(phaseShiftLimits), max(amplitudeLimits), max(offsetLimits)), 
                         control=list(trace=0), 
                         object=object, 
                         sin=sin,
                         slot=slot)
        return(bestFit)
    }else
    {
        bestFit <- optim(par=c(phaseShift=initialPhaseShift, amplitude=initialAmplitude, offset=initialOffset), 
                         function(par, object, sin, slot){sseTrack(object=object, sin=sin, slot=slot, phaseShift=par['phaseShift'], amplitude=par['amplitude'], offset=par['offset'])}, 
                         method='L-BFGS-B', 
                         lower=c(min(phaseShiftLimits), min(amplitudeLimits), min(offsetLimits)), 
                         upper=c(max(phaseShiftLimits), max(amplitudeLimits), max(offsetLimits)), 
                         control=list(trace=0), 
                         object=object, 
                         sin=sin,
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

plotFit <- function(track, fit=NULL)
{
    if(is.null(fit))
    {
        print('Generating fit to plot...')
        fit <- getTrackSinFitAll(track)
    }
    plot(track, slotY='vx', relY=FALSE, add=FALSE)
    guessPar <- getTrackParamGuess(track)
    guess <- getSweep(sin=TRUE, amplitude=guessPar['amplitude'], phaseShift=guessPar['phaseShift'], offset=guessPar['offset'], fi=track@fi, ff=track@ff, ti=track@ti, tf=track@tf, samplingFreq=track@samplingFreq, indices=track@index)
    #lines(guess$t, guess$v, col='blue')
    predicted <- getSweep(sin=TRUE, amplitude=fit$par['amplitude'], phaseShift=fit$par['phaseShift'], offset=fit$par['offset'], fi=track@fi, ff=track@ff, ti=track@ti, tf=track@tf, samplingFreq=track@samplingFreq, indices=track@index)
    lines(predicted$t, predicted$v, col='red')
    #lines(predicted$t, predicted$x-mean(predicted$x), col='green')
    points(predicted$inflections, rep(0,length(predicted$inflections)))
    points(track@t[track@validTimeIndices], track@vx[track@validTimeIndices], col='blue')
    #print("Guess:")
    #print(guessPar)
    print("Fit:")
    print(fit$par)
}

##### TrackList Fitting #####

sseTrackListSin <- function(object, slot, phaseShift, initialAmplitude=50, amplitudeLimits=c(0,500), initialOffset=0, offsetLimits=c(0,10000), guess=TRUE)
{
    print(paste("Optimizing phaseshift: ", phaseShift, sep=''))
    sseRet <- 0
    if(guess==TRUE)
    {
        for(track in object@tracks)
        {
            bestFit <- getTrackSinFitFixedPhase(object=track, slot=slot, initialAmplitude=initialAmplitude, amplitudeLimits=amplitudeLimits, phaseShift=phaseShift, initialOffset=initialOffset, offsetLimits=offsetLimits, guess=getTrackParamGuess(track))
            sseRet <- sseRet + bestFit$fit$value #bestFit$value gives the value of the sseTrack at the best fit param
            #print(paste("ID:", track@id, ", phaseShift:", bestFit$par['phaseShift'], ", amplitude:", bestFit$par['amplitude'], ", offset:", bestFit$par['offset']))
        }
    }
    else
    {
        for(track in object@tracks)
        {
            bestFit <- getTrackSinFitFixedPhase(object=track, initialAmplitude=initialAmplitude, amplitudeLimits=amplitudeLimits, phaseShift=phaseShift, initialOffset=initialOffset, offsetLimits=offsetLimits) # Don't use the getTrackParamGuessFunction
            sseRet <- sseRet + bestFit$fit$value #bestFit$value gives the value of the sseTrack at the best fit param
            #print(paste("ID:", track@id, ", phaseShift:", bestFit$par['phaseShift'], ", amplitude:", bestFit$par['amplitude'], ", offset:", bestFit$par['offset']))
        }
    }
    return(sseRet)
}

sseTrackListTri <- function(object, slot, phaseShift, initialAmplitude=50, amplitudeLimits=c(0,500), initialOffset=0, offsetLimits=c(0,10000), guess=TRUE)
{
    print(paste("Optimizing phaseshift: ", phaseShift, sep=''))
    sseRet <- 0
    if(guess==TRUE)
    {
        for(track in object@tracks)
        {
            bestFit <- getTrackTriFitFixedPhase(object=track, slot=slot, initialAmplitude=initialAmplitude, amplitudeLimits=amplitudeLimits, phaseShift=phaseShift, initialOffset=initialOffset, offsetLimits=offsetLimits, guess=getTrackParamGuess(track))
            sseRet <- sseRet + bestFit$fit$value #bestFit$value gives the value of the sseTrack at the best fit param
            #print(paste("ID:", track@id, ", phaseShift:", bestFit$par['phaseShift'], ", amplitude:", bestFit$par['amplitude'], ", offset:", bestFit$par['offset']))
        }
    }
    else
    {
        for(track in object@tracks)
        {
            bestFit <- getTrackTriFitFixedPhase(object=track, initialAmplitude=initialAmplitude, amplitudeLimits=amplitudeLimits, phaseShift=phaseShift, initialOffset=initialOffset, offsetLimits=offsetLimits) # Don't use the getTrackParamGuessFunction
            sseRet <- sseRet + bestFit$fit$value #bestFit$value gives the value of the sseTrack at the best fit param
            #print(paste("ID:", track@id, ", phaseShift:", bestFit$par['phaseShift'], ", amplitude:", bestFit$par['amplitude'], ", offset:", bestFit$par['offset']))
        }
    }
    return(sseRet)
}

getTrackListSinPhaseShiftGuess <- function(object, slot='x')
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
    fit <- getTrackSinFitAll(testerTrack)
    plotFit(testerTrack, fit)
    return(as.numeric(fit$par['phaseShift']))
}

getTrackListTriPhaseShiftGuess <- function(object, slot='x')
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
    fit <- getTrackTriFitAll(testerTrack)
    plotFit(testerTrack, fit)
    return(as.numeric(fit$par['phaseShift']))
}

setTrackListSinPhaseShift <- function(object, slot='vx', initialPhaseShift=0, phaseShiftLimits=c(-pi, pi), guess=TRUE)
{
    if(guess==TRUE)
    {
        initialPhaseShift <- getTrackListSinPhaseShiftGuess(object, slot)
    }
    bestFit <- optim(par=c(phaseShift=initialPhaseShift), function(par, object, slot) {sseTrackListSin(object=object, slot=slot, phaseShift=par['phaseShift'], guess=TRUE)}, method='L-BFGS-B', lower=min(phaseShiftLimits), upper=max(phaseShiftLimits), control=list(trace=3), object=object, slot=slot)
    object@phaseShift <- bestFit$par['phaseShift']
    for(i in 1:length(object@tracks))
    {
        object@tracks[[i]]@phaseShift <- bestFit$par['phaseShift']
    }
    return(object)
}

setTrackListTriPhaseShift <- function(object, slot='vx', initialPhaseShift=0, phaseShiftLimits=c(-pi, pi), guess=TRUE)
{
    if(guess==TRUE)
    {
        initialPhaseShift <- getTrackListTriPhaseShiftGuess(object, slot)
    }
    bestFit <- optim(par=c(phaseShift=initialPhaseShift), function(par, object, slot) {sseTrackListTri(object=object, slot=slot, phaseShift=par['phaseShift'], guess=TRUE)}, method='L-BFGS-B', lower=min(phaseShiftLimits), upper=max(phaseShiftLimits), control=list(trace=3), object=object, slot=slot)
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