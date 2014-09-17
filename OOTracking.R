library(foreign)
library(pracma)
library(plyr)
library(methods)
library(stats)

##### Maxima #####

# This is a test and a change by jaywarrick

Maxima <- setRefClass("Maxima",
                      fields = list(frame="numeric", points="data.frame"),
                      methods = list(
                           initializeWithROI = function(frame, polygon=NULL)
                           {
                                frame <<- as.numeric(as.character(frame))
                                if(!is.null(polygon))
                                {
                                     pairs <- strsplit(polygon,';')[[1]]
                                     x <- numeric(0)
                                     x0 <- numeric(0)
                                     y <- numeric(0)
                                     y0 <- numeric(0)
                                     index <- numeric(0)
                                     first <- TRUE
                                     for(pair in pairs)
                                     {
                                          nums <- strsplit(pair,',')[[1]]
                                          x <- append(x, as.numeric(nums[1]))
                                          y <- append(y, as.numeric(nums[2]))
                                          index <- append(index,as.numeric(nums[3]))
                                     }
                                     
                                     points <<- data.frame(id=index, x=x, y=y)
                                }
                                else
                                {
                                     points <<- data.frame(id=numeric(0), x=numeric(0), y=numeric(0))
                                }
                           },
                           length = function()
                           {
                                return(nrow(points))
                           },
                           addPoint = function(id, x, y)
                           {
                                points <<- rbind(points, data.frame(id=id, x=x, y=y))
                           },
                           plotMaxima = function(IDs=TRUE, pch=20, xlab='X [pixel]', ylab='Y [pixel]', cex=0.5, ...)
                           {
                                plot(points$x, points$y, pch=20, xlab=xlab, ylab=ylab, cex=cex, ...)
                                if(IDs)
                                {
                                     text(x=points$x, y=points$y, labels=as.character(points$id), cex=0.5, adj=c(0,0))
                                }
                           },
                           getXYZ = function()
                           {
                                data <- cbind(points, data.frame(z=0))
                                return(data[,2:4])
                           },
                           getPoint = function(index=NULL, id=NULL)
                           {
                                if(is.null(index) & is.null(id))
                                {
                                     return(NULL)
                                }
                                else if(!is.null(index))
                                {
                                     return(points[index,])
                                }
                                else
                                {
                                     return(points[points$id==id,])
                                }
                           }
                      )
)

##### MaximaList #####

MaximaList <- setRefClass('MaximaList',
                          fields = list(maxima='list', maxID='numeric', trackBackStart='numeric', trackBackEnd='numeric'),
                          methods = list(
                               initializeWithFile = function(path, timeDimName='Time')
                               {
                                    require(foreign)
                                    if(!is.null(path))
                                    {
                                         maximaFile <- read.arff(path)
                                         maximaFile2 <- reorganizeTable(maximaFile, nameCol='Metadata')
                                         maximaFile2[,timeDimName] <- as.numeric(as.character(maximaFile2[,timeDimName]))
                                         temp <- list()
                                         for(r in 1:nrow(maximaFile2))
                                         {
                                              newMaxima <- new('Maxima')
                                              newMaxima$initializeWithROI(frame=maximaFile2[r,timeDimName], polygon=maximaFile2[r,'polygonPts'])
                                              temp[[as.character(maximaFile2[r,timeDimName])]] <- newMaxima
                                         }
                                         temp <- temp[order(as.numeric(names(temp)), decreasing=FALSE)]
                                         maxima <<- temp
                                    }
                                    else
                                    {
                                         maxima <<- list()
                                    }
                                    
                                    maxID <<- -1 # Used for appending maxima
                               },
                               getMaxima = function(frame)
                               {
                                    return(maxima[[as.character(frame)]])
                               },
                               setMaxima = function(newMaxima)
                               {
                                    maxima[[as.character(newMaxima$frame)]] <<- newMaxima
                               },
                               length = function()
                               {
                                    return(base::length(maxima))
                               },
                               trackBack = function(startFrame, endFrame, maxDist=150, direction=c(1,0,0), directionality=10, uniformityDistThresh=3, digits=1)
                               {
                                    if(startFrame >= length() | startFrame < 0 | endFrame >= length() | endFrame >= startFrame | endFrame < 0)
                                    {
                                         stop("Start frame or end frame out of bounds. Start frame must be >=0 and > (yes >) end frame because we work backward. End frame must be < last frame of data set.")
                                    }
                                    trackBackStart <<- startFrame
                                    trackBackEnd <<- endFrame
                                    currentMaxima <- getMaxima(startFrame)
                                    previousMaxima <- getMaxima(startFrame-1)
                                    maxID <<- max(currentMaxima$points$id)
                                    while(!is.null(previousMaxima) && previousMaxima$frame >= endFrame)
                                    {
                                         cat("Linking frame: ", currentMaxima$frame, "\n", sep="")
                                         t1 <- currentMaxima$getXYZ()
                                         t0 <- previousMaxima$getXYZ() # maxima to prepend
                                         results <- DirectionalLinearAssignment(points=list(t0=t0, t1=t1), maxDist=maxDist, direction=direction, directionality=directionality, uniformityDistThresh=uniformityDistThresh, digits=digits)
                                         px <- results$px
                                         newMaxima <- new('Maxima')
                                         newMaxima$frame <- previousMaxima$frame
                                         #### Can potentially vectorize this process to make it faster
                                         for(i in 1:base::length(px))
                                         {
                                              if(px[i] <= nrow(t0) & i <= nrow(t1))
                                              {
                                                   #print("A")
                                                   # Then we have a link between t0 and t1 and the points should share ids
                                                   p0 <- previousMaxima$getPoint(index=px[i])
                                                   p1 <- currentMaxima$getPoint(index=i)
                                                   newMaxima$addPoint(id=p1$id, x=p0$x, y=p0$y)
                                              }
                                              else if(px[i] <= nrow(t0) & i > nrow(t1))
                                              {
                                                   #print("B")
                                                   # This point in t0 represents the start of a new track
                                                   # increment the maxID and add the point to the newMaxima
                                                   maxID <<- maxID + 1
                                                   p0 <- previousMaxima$getPoint(index=px[i])
                                                   newMaxima$addPoint(id=maxID, x=p0$x, y=p0$y)
                                              }
                                              else
                                              {
                                                   #print("C")
                                                   # This point in t1 represents the end of a track or an auxillary match to enable other types of matches
                                                   # We don't need to do anything
                                              }
                                         }
                                         setMaxima(newMaxima)
                                         currentMaxima <- newMaxima
                                         previousMaxima <- getMaxima(newMaxima$frame - 1)
                                    }
                                    cat("Done tracking. Trimming Data to start and end frame")
                                    allFrameNames <- as.numeric(names(maxima))
                                    indicesToRemove <- which(!(allFrameNames %in% endFrame:startFrame))
                                    maxima[indicesToRemove] <<- NULL
                               },
                               getTrackList = function(sin=FALSE, fi, ff, tAll, frameRange=c(0,-1))
                               {
                                    if(frameRange[1] < 0 | frameRange[2] > trackBackStart | frameRange[2] < frameRange[1])
                                    {
                                         exportFrames <- trackBackEnd:trackBackStart
                                    }
                                    else
                                    {
                                         exportFrames <- frameRange[1]:frameRange[2]
                                    }
                                    trackList <- new('TrackList', sin=sin, fi=fi, ff=ff, tAll=tAll)
                                    count = 0
                                    total = base::length(exportFrames)
                                    for(f in exportFrames)
                                    {
                                         .maxima <- getMaxima(f)
                                         for(i in 1:.maxima$length())
                                         {    
                                              data <- .maxima$points[i,]
                                              trackList$addTrackPoint(id=data$id, x=data$x, y=data$y, t=tAll[.maxima$frame+1], frame=.maxima$frame)
                                         }
                                         count <- count + 1
                                         cat("Generating TrackList: ", round(100*count/total, digits=2), "%\n", sep="")
                                    }
                                    trackList$calculateVelocities()
                                    return(trackList)
                               },
                               getProp = function(fun=function(.maxima){return(base::range(.maxima$points$x))})
                               {
                                    ret <- list()
                                    for(.maxima in maxima)
                                    {
                                         ret[[as.character(.maxima$frame)]] <- fun(.maxima)
                                    }
                                    return(ret)
                               },
                               generateMaximaPlots = function(path=NULL, baseName='Frame_')
                               {
                                    dir.create(path, recursive=TRUE)
                                    wd <- getwd()
                                    setwd(path)
                                    xlim <- base::range(unlist(getProp(fun=function(.maxima){return(base::range(.maxima$points$x))})))
                                    ylim <- base::range(unlist(getProp(fun=function(.maxima){return(base::range(.maxima$points$y))})))
                                    
                                    makeThePlots <- function()
                                    {
                                         for(.maxima in maxima)
                                         {
                                              cat("Plotting Frame = ", .maxima$frame, "\n", sep="")
                                              pdf(file=paste(path, '/', baseName, .maxima$frame, '.pdf', sep=''), width=9, height=4)
                                              .maxima$plotMaxima(IDs=TRUE, xlim=xlim, ylim=ylim, main=as.character(.maxima$frame))
                                              dev.off()
                                         }
                                    }
                                    
                                    tryCatch(expr=makeThePlots(), finally=setwd(wd))
                               },
                               offsetFrames = function(offset=0)
                               {
                                    names(maxima) <<- as.character(as.numeric(names(maxima))+offset)
                                    for(.maxima in maxima)
                                    {
                                         .maxima$frame <-.maxima$frame+offset
                                    }
                                    trackBackStart <<- trackBackStart + offset
                                    trackBackEnd <<- trackBackEnd + offset
                               }
                          )
)

##### TrackList #####

TrackList <- setRefClass('TrackList',
                         fields = list(tracks='list', sin='logical', fi='numeric', ff='numeric', tAll='numeric', phaseShift='numeric', validFrames='numeric'),
                         methods = list(
                              initializeWithFile = function(file=NULL, sin=FALSE, fi, ff, tAll)
                              {
                                   require(foreign)
                                   tracksFile <- read.arff(file)
                                   tracksFile2 <- reorganizeTable(tracksFile, nameCol='Metadata')
                                   
                                   for(row in 1:nrow(tracksFile2))
                                   {
                                        id <- tracksFile2[row,'Track']
                                        start <- tracksFile2[row,'polygonPts']
                                        pattern <- tracksFile2[row,'patternPts']
                                        newTrack <- new('Track')
                                        newTrack$initializeWithTrackROI(id=id, start=start, pattern=pattern, tAll=tAll)
                                        setTrack(newTrack)
                                   }
                                   sin <<- sin
                                   fi <<- fi
                                   ff <<- ff
                                   tAll <<- tAll
                              },
                              length = function()
                              {
                                   return(base::length(tracks))
                              },
                              plotTrackList = function(x, slot='vx', fun=NULL, rel=FALSE, ...)
                              {
                                   xRanges <- getProp(fun=function(track){track$range('t')})
                                   xRanges <- matrix(unlist(xRanges), ncol=2, byrow=TRUE)
                                   yRanges <- getProp(fun=function(track){track$range(slot)})
                                   yRanges <- matrix(unlist(yRanges), ncol=2, byrow=TRUE)
                                   
                                   Xmin <- min(xRanges[,1], na.rm=TRUE)
                                   Xmax <- max(xRanges[,2], na.rm=TRUE)
                                   Ymin <- min(yRanges[,1], na.rm=TRUE)
                                   Ymax <- max(yRanges[,2], na.rm=TRUE)
                                   
                                   args <- list(...)
                                   if(is.null(args$xlim))
                                   {
                                        xlim <- c(Xmin, Xmax)
                                   }
                                   else
                                   {
                                        xlim <- args$xlim
                                        args$xlim <- NULL
                                   }
                                   if(is.null(args$ylim))
                                   {
                                        ylim <- c(Ymin, Ymax)
                                   }
                                   else
                                   {
                                        ylim <- args$ylim
                                        args$ylim <- NULL
                                   }
                                   print(args)
                                   first <- TRUE
                                   for(track in tracks)
                                   {
                                        if(first)
                                        {
                                             do.call(track$plotTrack, c(list(slotY=slot, funY=fun, relY=rel, add=F, xlim=xlim, ylim=ylim, withTitle=FALSE, main='All Tracks'), args))
                                             first = FALSE
                                        }
                                        else
                                        {
                                             do.call(track$plotTrack, c(list(slotY=slot, funY=fun, relY=rel, add=T), args))
                                        }
                                   }
                              },
                              getTrack = function(id)
                              {
                                   return(tracks[[as.character(id)]])
                              },
                              setTrack = function(newTrack)
                              {
                                   newTrack$.parent <- .self
                                   tracks[[as.character(newTrack$id)]] <<- newTrack
                              },
                              removeTrack = function(id)
                              {
                                   tracks[[as.character(id)]] <<- NULL
                              },
                              getProp = function(fun=function(x){return(x$length())})
                              {
                                   ret <- list()
                                   for(.track in tracks)
                                   {
                                        ret[[as.character(.track$id)]] <- fun(.track)
                                   }
                                   return(ret)
                              },
                              addTrackPoint = function(id, x, y, t, frame)
                              {
                                   track <- getTrack(id)
                                   if(is.null(track))
                                   {
                                        track <- new('Track')
                                        track$id <- id
                                        # track$.parent <- .self # setTrack resets .parent
                                   }
                                   track$addPoint(x=x, y=y, t=t, frame=frame)
                                   setTrack(track)
                              },
                              calculateVelocities = function()
                              {
                                   callOnTracks('calculateVelocities')
                              },
                              smoothVelocities = function(fit, dist, maxWidth)
                              {
                                   widths <- getWindowWidths(fit=bestFit, trackList=trackList, dist=10, maxWidth=25)
                                   callOnTracks('smoothVelocities', widths=widths)
                              },
                              applyToTracks = function(fun, ...)
                              {
                                   for(track in tracks)
                                   {
                                        track <- fun(track, ...)
                                   }
                              },
                              callOnTracks = function(funName, ...)
                              {
                                   tot <- length()
                                   count <- 0
                                   for(track in tracks)
                                   {
                                        count <- count + 1
                                        cat("Calling", funName, "on track", count, "of", tot, "\n")
                                        # Have to do eval(parse()) because track[[funName]] is NULL while track$parsedFunName is not NULL, don't know why
                                        # Now that the function is loaded we can call it using the [[]] method
                                        theCall <- paste("track$'", funName, "'", sep="")
                                        theFunc <- eval(parse(text=theCall))
                                        if(is.null(theFunc))
                                        {
                                             stop(cat("Couldn't find function with name",funName))
                                        }
                                        do.call(theFunc, list(...))
                                   }
                              },
                              filterTracks = function(fun, ...)
                              {
                                   tracks <<- Filter(function(x){fun(x,...)}, tracks)
                                   if(base::length(tracks) == 0)
                                   {
                                        message("No tracks fit filter, resulting tracklist is of length 0!")
                                   }
                              },
                              sortTracks = function(fun=function(x){return(x$length())}, decreasing=TRUE, ...)
                              {
                                   ret <- getProp(fun=fun, ...)
                                   sorted <- sort(unlist(ret), index.return=T, decreasing=decreasing)
                                   tracks <<- tracks[sorted$ix]
                              },
                              getMatrix = function(slot='vx', validOnly=FALSE)
                              {
                                   #' Get a matrix of the tracklist data. Track id's are rows while 
                                   #' frame numbers are columns. Use the matrix form of TrackList to do 
                                   #' bulk operations such as determining the sse of all the data to a 
                                   #' single curve.
                                   frames <- (1:base::length(tAll))-1
                                   ids <- names(tracks)
                                   data <- matrix(NA, base::length(ids), base::length(frames), dimnames=list(id=names(tracks), frame=as.character(frames)))
                                   if(validOnly)
                                   {
                                        for(track in tracks)
                                        {
                                             data[as.character(track$id), (track$validFrames + 1)] <- track$getSlot(slot=slot, rel=FALSE, validOnly=validOnly)
                                        }
                                   }
                                   else
                                   {
                                        for(track in tracks)
                                        {
                                             data[as.character(track$id), (track$points$frame + 1)] <- track$points[,slot]
                                        }
                                   }
                                   
                                   return(data)
                              },
                              setValidFrames = function(fit, validStart=0.1, validEnd=0.99)
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
                                        upper <- which((t-inflectionPoint) >= 0)[1]
                                        if(is.na(upper) || upper == 1)
                                        {
                                             return(NA)
                                        }
                                        lower <- upper - 1
                                        return(c(lower, upper))
                                   }
                                   
                                   sweep <- getSweep(amplitude=fit$par['amplitude'], phaseShift=fit$par['phaseShift'], offset=0, sin=sin, fi=fi, ff=ff, tAll=tAll, frames=-1, guess=NULL)
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
                                   validFrames0 <- c(1:base::length(tAll)) - 1
                                   validFrames0 <- validFrames0[-indicesToRemove]
                                   
                                   validFrames1 <- numeric(0)
                                   for(tIndex in 1:base::length(tAll))
                                   {
                                        # Get the nearest inflection at or beyond this time
                                        
                                        temp <- which((sweep$inflections >= tAll[tIndex]) & inflectionsToAddress)
                                        infIndex <- temp[1]
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
                                        if( (tAll[tIndex] >= (infT1 + validStart*dInfT)) & (tAll[tIndex] <= (infT1 + validEnd*dInfT)) )
                                        {
                                             # If it is within the startValid and endValid bounds, add it to the list of the valid frames
                                             validFrames1 <- c(validFrames1, tIndex-1) # (tIndex-1) = frame because frames are indicies that start at 0
                                        }
                                   }
                                   
                                   validFrames <<- sort(intersect(validFrames1, validFrames0))
                                   
                                   for(track in tracks)
                                   {
                                        track$setValidFrames(validFrames)
                                   }
                              },
                              getPercentAdhered = function(velocityThreshold=3)
                              {
                                   trackMatrix <- getMatrix(slot='vxs', validOnly=TRUE)
                                   ret <- list()
                                   for(frame in colnames(trackMatrix))
                                   {
                                        velocities <- trackMatrix[,frame]
                                        velocities <- abs(velocities[!is.na(velocities)])
                                        if(!isempty(velocities))
                                        {
                                             adhered <- sum(velocities < velocityThreshold)/base::length(velocities)
                                             ret[[frame]] <- adhered
                                        }
                                   }
                                   percents <- 100*as.numeric(ret)
                                   times <- trackList$tAll[as.numeric(names(ret))+1]
                                   return(data.frame(time=times, percentAdhered=percents))
                              }
                         )
) 

##### Track #####

Track <- setRefClass('Track',
                     fields = list(id='numeric', points='data.frame', positionFit='list', velocityFit='list', validFrames='numeric', .parent='TrackList'),
                     methods = list(
                          initializeWithTrackROI = function(id, start, pattern, tAll)
                          {
                               "
                               #' Using the value in the 'PolygonPts' value and the 'PatternPts'
                               #' value from a JEX ROI that represents a track, create an R 'Track'
                               #' object. Assign the track the id given as 'id'.
                               "
                               pairs <- strsplit(pattern,';')[[1]]
                               x <- numeric(0)
                               x0 <- numeric(0)
                               y <- numeric(0)
                               y0 <- numeric(0)
                               frames <- numeric(0)
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
                                         frames <- append(frames,i0)
                                    }
                                    else
                                    {
                                         nums <- strsplit(pair,',')[[1]]
                                         x <- append(x, x0+as.numeric(nums[1]))
                                         y <- append(y, y0+as.numeric(nums[2]))
                                         frames <- append(frames,as.numeric(nums[3]))
                                         # print(nums)
                                    }
                               }                               
                               id <<- as.numeric(as.character(id))
                               if(base::length(frames) > base::length(tAll))
                               {
                                    stop(cat("Length of tAll seems to short for the data that is being used to initialize the track. frames length = ", base::length(frames), " tAll length = ", base::length(tAll)))
                               }
                               points <<- data.frame(frame=frames, x=x, y=y, t=tAll[frames+1])
                               calculateVelocities()
                          },
                          length = function()
                          {
                               "
                               #' Return the number of frames in this track.
                               "
                               return(nrow(points))
                          },
                          getSlot = function(slot, rel=FALSE, validOnly=FALSE)
                          {
                               "
                               #' Get a vector of the values indicated by 'slot' for this point.
                               #' This accessor method is useful because the data is stored within
                               #' a data.frame field called points. 
                               #' 'validOnly' indicates whether to return values for the validFrames only or all frames
                               "
                               
                               validIndices <- points$frame %in% validFrames
                               
                               if(rel)
                               {
                                    temp <- points[,slot] - mean(points[,slot])
                                    if(validOnly & !is.null(validFrames))
                                    {
                                         temp <- temp[validIndices] # frames start at 0 while R indices start at 1
                                    }
                               }
                               else
                               {
                                    temp <- points[,slot]
                                    if(validOnly & !isempty(validFrames))
                                    {
                                         temp <- temp[validIndices]
                                    }
                               }
                               temp <- temp[!is.na(temp)]
                               return(temp)
                          },
                          range = function(slot, rel=FALSE)
                          {
                               "
                               #' For the indicated 'slot' (i.e., column within the points field),
                               #' return the range that this value takes.
                               "
                               ret <- base::range(getSlot(slot, rel=rel))
                               return(ret)
                          },
                          calculateVelocities = function()
                          {
                               "
                               #' Internally populate/update the values 'vx' and 'vy' within
                               #' the points vield using 'getDerivative(x, y)'
                               "
                               points$vx <<- getDerivative(points$x, points$t)
                               points$vy <<- getDerivative(points$y, points$t)
                          },
                          smoothVelocities = function(widths)
                          {
                               points$vxs <<- sapply(seq_along(points$frame), 
                                                     getAverage,
                                                     frames=points$frame,
                                                     widths=widths,
                                                     data=points$vx
                               )
                               points$vys <<- sapply(seq_along(points$frame), 
                                                     getAverage,
                                                     frames=points$frame,
                                                     widths=widths,
                                                     data=points$vy
                               )
                          },
                          plotTrack = function(slotX='t', slotY='x', funX=NULL, funY=NULL, relX=FALSE, relY=TRUE, validOnly=FALSE, add=FALSE, withTitle=TRUE, col='black', lwd=1, lty=1, xlab=slotX, ylab=slotY, type='l', ...)
                          {
                               "
                               #' Plot the track using 'slotX' as the x variable and
                               #' 'slotY' as the y variable, using 'relX' and 'relY'
                               #' to indicate whether to plot the value 'relative' to
                               #' its mean value. Additional args are passed to the 
                               #' 'plot' method that is used internally.
                               #' 'validOnly' indicates whether to plot validFrames only or all frames
                               "
                               xData <- getSlot(slotX, relX, validOnly=validOnly)
                               if(!is.null(funX))
                               {
                                    xData <- funX(xData)
                               }
                               yData <- getSlot(slotY, relY, validOnly=validOnly)
                               if(!is.null(funY))
                               {
                                    yData <- funY(yData)
                               }
                               if(base::length(xData) < 0 || base::length(yData) < 0)
                               {
                                    message(paste('No valid data to plot for track ', id, '.', sep=''))
                               }
                               if(relX)
                               {
                                    xlab <- paste(xlab, ' - mean(', xlab, ')', sep='')
                               }
                               if(relY)
                               {
                                    ylab <- paste(ylab, ' - mean(', ylab, ')', sep='')
                               }
                               if(base::length(list(...)$validOnly)>0)
                               {
                                    print(list(...))
                                    stop("what?")
                               }
                               if(is.na(xData) || is.na(yData))
                               {
                                    message(paste('Nothing to plot for track ', id, '.', sep=''))
                               } else if(base::length(xData)==0 || base::length(yData) == 0)
                               {
                                    message(paste('No valid data to plot for track ', id, '.', sep=''))
                               } else
                               {
                                    if(add)
                                    {
                                         graphics::lines(xData, yData, col=col, lwd=lwd, lty=lty, type=type, ...)
                                    }
                                    else
                                    {
                                         if(withTitle)
                                         {
                                              graphics::plot(xData, yData, col=col, lwd=lwd, lty=lty, xlab=xlab, ylab=ylab, main=as.character(id), type=type, ...)
                                         }
                                         else
                                         {
                                              graphics::plot(xData, yData, col=col, lwd=lwd, lty=lty, xlab=xlab, ylab=ylab, type=type, ...)
                                         }
                                    }
                               }
                          },
                          addPoint = function(x, y, t, frame)
                          {
                               "
                              #' Add a point to the 'points' field of this track. It is
                              #' assumed you are adding things in an appropriate order.
                              #' Typically the points are listed in order of frame.
                              "
                               points <<- rbind(points, data.frame(x=x, y=y, t=t, frame=frame))
                          },
                          show = function()
                          {
                               "
                              #' Prints the information of this object to the command line
                              #' output. This overrides the basic 'show' method as a track
                              #' carries a reference to its parent trackList (if it has been
                              #' added to a trackList using the 'setTrack' method) which
                              #' results in recursive printing of TrackList information
                              #' because a 'show' for 'TrackList' calls 'show' on each 'Track'.
                              #' This method eliminates this issue. 
                              "
                               cat("Reference class object of class 'Track'\n")
                               for(name in names(Track$fields()))
                               {
                                    if(name %in% c('.parent'))#c('x','y','t','vx','vy')))
                                    {
                                         cat("Field '", name, "':\n", sep='')
                                         cat("(can't print... recursive)\n")
                                    }
                                    else
                                    {
                                         cat("Field '", name, "':\n", sep='')
                                         methods::show(.self[[name]])
                                    }
                               }
                               cat("\n")
                          },
                          setValidFrames = function(frames)
                          {
                               validFrames <<- points$frame[points$frame %in% frames]
                          },
                          getTrackSweep = function(amplitude=1, offset=mean(getSlot(slot='x')), validFramesOnly=FALSE, guess=NULL)
                          {
                               "
                               #' This method is provided as a convenience. It calls 'getSweep' using
	 				      #' parameters that exist within the 'Track' and parent 'TrackList'
	 				      #' when available or the 'getSweep' defualts.
	 				      "
                               if(validFramesOnly)
                               {
                                    frames <- validFrames
                               }
                               else
                               {
                                    frames <- points$frame
                               }
                               args <- list(sin=.parent$sin, fi=.parent$fi, ff=.parent$ff, tAll=.parent$tAll, phaseShift=.parent$phaseShift, amplitude=amplitude, offset=offset, frames=frames, guess=guess)
                               args <- args[!isempty(args)]
                               return(do.call(getSweep, args))
                          },
                          sseTrack = function(amplitude=50, phaseShift=0, offset=0, validFramesOnly=FALSE)
                          {
                               "
      				      #' Calculates the sum square error (i.e., sse) between this track
	 				      #' and a sweep funciton with the given 'amplitude', 'phaseShift',
	 				      #' and 'offset'. The parameter 'validFramesOnly' will limit the
	 				      #' calculation to just the 'validFrames' listed in this Track.
	 				      #' The parameter of 'sin' = T/F, and tAll, fi, and ff are passed
	 				      #' from the track's parent 'TrackList'.
	 				      "
                               if(validFramesOnly)
                               {
                                    frames <- validFrames
                               }
                               else
                               {
                                    frames <- points$frame
                               }
                               args <- list(sin=.parent$sin, fi=.parent$fi, ff=.parent$ff, tAll=.parent$tAll, amplitude=amplitude, offset=offset, frames=frames)
                               args <- args[!isempty(args)]
                               predicted <- do.call(getSweep, args)
                               data <- object$points$vx                ### Explicitly fitting vx ###
                               indicesToGet <- which(points$frame %in% frames)
                               result <- sum((data[indicesToGet]-predicted$v)^2)
                               return(result)
                          }
                     )
)

##### General #####

#' Take an arff file and reorganize it into a more standard 'table' format. Specifically this is used to
#' import an arff file from JEX as JEX uses a column called 'Measurement' to define the type of measurment
#' or property being stored and 'Value', the value of that property.
#' 
#' @param data An object that is the result of using foreign::read.arff(file) on an arff file
#' @param baseName An optional basename to add to whatever label is in the \code{nameCol} portion of each row entry
#' @param convertToNumeric An option to convert the columns of information within \code{data} leading up to 
#' \code{nameCol} and \code{valueCol} to numeric or to leave as text. Default is to convert to numeric (i.e., TRUE)
#' @param nameCol The name of the column that describes the nature of the value in the \code{valueCol}
#' @param valueCol The name of the column with the values of the properties listed in the \code{nameCol}
reorganizeTable <- function(data, baseName=NA, convertToNumeric=TRUE, nameCol='Measurement', valueCol='Value')
{
     require(plyr)
     idCols <- names(data)
     idCols <- idCols[-which(idCols %in% c(nameCol,valueCol))]
     newData <- data.frame(stringsAsFactors=FALSE)
     measurements <- unique(data[,nameCol])
     #     m <- measurements[1]
     #     browser()
     for(m in measurements)
     {
          if(is.na(baseName))
          {
               newColName <- m
               newColName <- gsub(' ','.',newColName, fixed=TRUE) # Get rid of extraneous spaces
          }else
          {
               newColName <- paste(baseName,'.',m, sep='')
               newColName <- gsub(' ','.',newColName, fixed=TRUE) # Get rid of extraneous spaces
          }
          
          temp <- data[data[,nameCol]==m,]
          temp2 <- temp[,c(idCols,valueCol)]
          names(temp2)[names(temp2)==valueCol] <- newColName
          if(nrow(newData) == 0)
          {
               newData <- temp2
          }else
          {
               newData <- merge(newData, temp2, by=idCols)
          }
     }
     
     if(convertToNumeric)
     {
          for(n in idCols)
          {
               newData[,n] <- as.numeric(as.character(newData[,n]))
          }
     }
     
     return(newData)
}

first <- function(x)
{
     return(x[1])
}

last <- function(x)
{
     require(pracma)
     return(x[numel(x)])
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

interpolateDerivative <- function(f0, f1, f2, x0, x1, x2, xj)
{
     # This is a three point interpolation of the derivative where the interpolated point
     # is the middle of the 3 points. This simplifies to the the three-point midpoint formula
     # when the time steps are equal but can handle when timesteps are unequal (i.e., the 
     # time-step on either side of the 3 points is not equal)
     term1 <- f0*((2*xj-x1-x2)/((x0-x1)*(x0-x2)))
     term2 <- f1*((2*xj-x0-x2)/((x1-x0)*(x1-x2)))
     term3 <- f2*((2*xj-x0-x1)/((x2-x0)*(x2-x1)))
     return(term1 + term2 + term3)
}

getWindowWidths <- function(fit, trackList, dist, maxWidth)
{
     dt <- trackList$tAll[2]-trackList$tAll[1]
     v <- 4*fit$par[['amplitude']]*(trackList$fi*(trackList$ff/trackList$fi)^(trackList$tAll/last(trackList$tAll)))
     widths <- ceiling(dist/(v*dt))
     widths[widths > maxWidth] <- maxWidth
     return(widths)
}

#' 'x' is the index within 'frames' for which we will perform the calculation
#' 'frames' are the frames in this track
#' 'widths' are the velocity dependent widths of the averaging window and has the same length as 'tAll'
#' 'data' is the vector of data for which we will calculate the windowed averages
getAverage <- function(x, frames, widths, data)
{
     # Subtract 1 to represent the number of intervals instead of number of points to average
     width <- widths[frames[x]+1] - 1
     newX <- x - floor(width/2)
     if(newX < 1)
     {
          newX <- 1
     }
     if((newX+width) > length(frames))
     {
          width <- length(frames)-newX
     }
     mean(data[newX:(newX+width)])
}

##### Track Filtering #####

#' Returns true if the track length within the indicated length bounds.
trackLengthFilter <- function(track, min=0, max=1000000)
{
     l <- track$length()
     if(l >= min & l <=max)
     {
          return(TRUE)
     }
     else
     {
          return(FALSE)
     }
}

#' Returns true if the track is defined within the indicated frame bounds.
trackFrameFilter <- function(track, startMin=0, startMax=1000000, endMin=0, endMax=1000000)
{
     startFrame <- first(track$points$frame)
     endFrame <- last(track$points$frame)
     
     if(startFrame >= startMin & startFrame <= startMax & endFrame >= endMin & endFrame <= endMax)
     {
          return(TRUE)
     }
     else
     {
          return(FALSE)
     }
}

##### Fitting #####

getSweep <- function(amplitude=1, phaseShift=0, offset=0, sin=FALSE, fi=2, ff=0.1, tAll=seq(0,500,1), frames=-1, guess=NULL)
{
     if(frames[1] < 0 | max(frames) > (length(tAll)-1))
     {
          frames <- (1:length(tAll)) - 1 #Subtract 1 to be consistent with frames from JEX, which start at 0 instead of 1
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
     if(ff==fi)
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
     x[sections==1] <- -1*(x[sections==1]-A)
     x[sections==2] <- -1*(x[sections==2])
     x[sections==3] <- (x[sections==3]-A)
     
     # Find pi/2 intervals on frequency sweep
     
     offsetInflections <- (  (-1*(phi/(pi/2))) %/% 1  ) 
     ni <- seq(0,1000,1) + offsetInflections
     suppressWarnings(inflections <- (log((log(ff/fi)*phi)/(2*fi*pi*tf)+(log(ff/fi)*ni)/(4*fi*tf)+1)*tf)/log(ff/fi))
     inflections <- inflections[!is.nan(inflections)]
     inflectionNums <- (  seq(0,length(inflections)-1,1) + (  offsetInflections %% 4  )  ) %% 4
     startTimeI <- which(inflections >= t[1])[1]
     endTimeI <- which(inflections <= t[length(t)])
     endTimeI <- endTimeI[length(endTimeI)]
     if(is.na(endTimeI))
     {
          stop(paste(t[length(t)], "is greater than last time of last inflection point. Therefore, there must be an error in fi/ff or tAll as there are inflections guaranteed in during times ti to tf given fi and ff."))
     }
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

sseBulk <- function(trackList, trackMatrix, amplitude=50, phaseShift=0)
{
     # Cols are frames, rows are track ids
     # we pass the 'trackMatrix' so we don't have to obtain it each iteration
     
     # get the sweep (for this we need the 'trackList')
     sweep <- getSweep(amplitude=amplitude, phaseShift=phaseShift, offset=0, sin=trackList$sin, fi=trackList$fi, ff=trackList$ff, tAll=trackList$tAll, frames=-1, guess=NULL)
     
     # For each index in tAll (i.e., for each frame)
     sse <- sum((t(trackMatrix)-sweep$v)^2, na.rm=TRUE) # Do the transpose because the subtract function typically makes the subtracted vector vertical
     cat("(", amplitude, ",", phaseShift, ") = ", sse, "\n", sep="")
     return(sse)
}

getBulkPhaseShift <- function(trackList)
{
     require(stats)
     trackMatrix <- trackList$getMatrix()
     amplitude <- max(as.numeric(trackList$getProp(fun=function(x){r <- x$range('x', rel=TRUE); r <- (r[2]-r[1])/2; return(r)})))
     guess <- c(amplitude=amplitude, phaseShift=0)
     amplitudeLimits=c(0.5, 10*amplitude)
     phaseShiftLimits=c(-pi, pi)
     offsetLimits=c(0,10000)
     bestFit <- optim(par=guess, 
                      function(par, trackList, trackMatrix){sseBulk(trackList=trackList, trackMatrix=trackMatrix, amplitude=par['amplitude'], phaseShift=par['phaseShift'])}, 
                      method='L-BFGS-B', 
                      lower=c(min(amplitudeLimits), min(phaseShiftLimits)), 
                      upper=c(max(amplitudeLimits), max(phaseShiftLimits)), 
                      control=list(trace=0),
                      trackList=trackList,
                      trackMatrix=trackMatrix)
     return(list(par=c(phaseShift=as.numeric(bestFit$par['phaseShift']), amplitude=as.numeric(bestFit$par['amplitude']), offset=as.numeric(bestFit$par['offset'])), fit=bestFit))
}
