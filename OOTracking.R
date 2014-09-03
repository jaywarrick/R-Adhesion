rm(list=ls())

source('~/.Rprofile')
source('~/Public/DropBox/GitHub/R-Adhesion/Tracking.R')
source('~/Public/DropBox/GitHub/R-Adhesion/ADhesionAnalysis_HelperFunctions.R')
library(foreign)
path1 <- '/Users/jaywarrick/Documents/JEX/LocalTest/temp/JEXData0000000001.arff'
path2 <- '/Users/jaywarrick/Documents/JEX/LocalTest/PC3 vs LNCaP/Cell_x0_y0/Roi-Maxima Upper/x0_y0.jxd'
fileTable <- read.arff(path1)


##### Maxima #####

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
                                    library(foreign)
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
                               trackBack = function(startFrame, endFrame, maxDist=150, direction=c(1,0,0), directionality=10, uniformityDistThresh=15, digits=1)
                               {
                                    if(startFrame >= length() | startFrame < 0 | endFrame >= length() | endFrame >= startFrame | endFrame < 0)
                                    {
                                         stop("Start frame or end frame out of bounds. Start frame must be >=0 and > (yes >) end frame because we work backward. End frame must be < last frame of data set.")
                                    }
                                    trackBackStart <<- startFrame
                                    trackBackEnd <<- endFrame
                                    currentMaxima <- .self$getMaxima(startFrame)
                                    previousMaxima <- .self$getMaxima(startFrame-1)
                                    maxID <<- max(currentMaxima$points$id)
                                    while(!is.null(previousMaxima) & previousMaxima$frame >= endFrame)
                                    {
                                         cat("Current frame: ", currentMaxima$frame, "\n", sep="")
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
                                         .self$setMaxima(newMaxima)
                                         currentMaxima <- newMaxima
                                         previousMaxima <- .self$getMaxima(newMaxima$frame - 1)
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
                                    exportIndices <- exportFrames + 1
                                    trackList <- new('TrackList', sin=sin, fi=fi, ff=ff, tAll=tAll[exportIndices])
                                    count = 0
                                    total = base::length(exportFrames)
                                    for(maxima in maximaList$maxima[exportIndices])
                                    {
                                         for(i in 1:maxima$length())
                                         {
                                              data <- maxima$points[i,]
                                              trackList$addTrackPoint(id=data$id, x=data$x, y=data$y, t=tAll[maxima$frame+1], frame=maxima$frame)
                                         }
                                         count <- count + 1
                                         cat("Percent complete: ", round(100*count/total, digits=2), "%\n", sep="")
                                    }
                                    trackList$calculateVelocities()
                                    return(trackList)
                               }
                          )
)

##### TrackList #####

TrackList <- setRefClass('TrackList',
                         fields = list(tracks='list', sin='logical', fi='numeric', ff='numeric', tAll='numeric', phaseShift='numeric'),
                         methods = list(
                              initializeWithFile = function(object, file=NULL, sin=FALSE, fi, ff, tAll)
                              {
                                   library(foreign)
                                   tracksFile <- read.arff(file)
                                   tracksFile2 <- reorganizeTable(tracksFile, nameCol='Metadata')
                                   
                                   duh <- list()
                                   for(row in 1:nrow(tracksFile2))
                                   {
                                        id <- tracksFile2[row,'Track']
                                        start <- tracksFile2[row,'polygonPts']
                                        pattern <- tracksFile2[row,'patternPts']
                                        newTrack <- trackFromTrackROI(id=id, sin=sin, fi=fi, ff=ff, tAll=tAll, start=start, pattern=pattern)
                                        duh[as.character(tracksFile2[row,'Track'])] <- newTrack
                                   }
                                   
                                   sin <<- sin
                                   fi <<- fi
                                   ff <<- ff
                                   tAll <<- tAll
                                   tracks <<- duh
                              },
                              length = function()
                              {
                                   return(base::length(tracks))
                              },
                              plotTrackList = function(x, slot='vx', rel=FALSE, ...)
                              {
                                   xRanges <- getProp(fun=function(track){track$range('t')})
                                   xRanges <- matrix(unlist(xRanges), ncol=2, byrow=TRUE)
                                   yRanges <- getProp(fun=function(track){track$range(slot)})
                                   yRanges <- matrix(unlist(yRanges), ncol=2, byrow=TRUE)
                                   
                                   Xmin <- min(xRanges[,1], na.rm=TRUE)
                                   Xmax <- max(xRanges[,2], na.rm=TRUE)
                                   Ymin <- min(yRanges[,1], na.rm=TRUE)
                                   Ymax <- max(yRanges[,2], na.rm=TRUE)
                                   
                                   xlim <- c(Xmin, Xmax)
                                   ylim <- c(Ymin, Ymax)
                                   
                                   first <- TRUE
                                   for(track in tracks)
                                   {
                                        if(first)
                                        {
                                             track$plotTrack(slotY=slot, relY=rel, add=F, xlim=xlim, ylim=ylim, withTitle=FALSE, main='All Tracks', ...)
                                             first = FALSE
                                        }
                                        else
                                        {
                                             track$plotTrack(slotY=slot, relY=rel, add=T, ...)
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
                                   track <- .self$getTrack(id)
                                   if(is.null(track))
                                   {
                                        track <- new('Track')
                                        track$id <- id
                                        # track$.parent <- .self # setTrack resets .parent
                                   }
                                   track$addPoint(x=x, y=y, t=t, frame=frame)
                                   .self$setTrack(track)
                              },
                              calculateVelocities = function()
                              {
                                   for(track in tracks)
                                   {
                                        track$calculateVelocities()
                                   }
                              }
                         )
) 

##### Track #####

Track <- setRefClass('Track',
                     fields = list(id='numeric', points='data.frame', positionFit='list', velocityFit='list', validFrames='numeric', .parent='TrackList'),
                     methods = list(
                          trackFromTrackROI = function(id, start, pattern)
                          {
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
                               calculateVelocities()
                               
                               id <<- as.numeric(as.character(id))
                               points <<- data.frame(frame=frames, x=x, y=y, vx=vx, vy=vy)
                          },
                          length = function()
                          {
                               return(nrow(points))
                          },
                          getSlot = function(slot, rel=FALSE, validOnly=FALSE)
                          {
                               if(rel)
                               {
                                    temp <- points[,slot] - mean(points[,slot])
                                    if(validOnly & !is.null(validFrames))
                                    {
                                         temp <- temp[validFrames + 1] # frames start at 0 while R indices start at 1
                                    }
                               }
                               else
                               {
                                    temp <- points[,slot]
                                    if(validOnly & !isempty(validFrames))
                                    {
                                         temp <- temp[validFrames + 1]
                                    }
                               }
                               return(temp)
                          },
                          range = function(slot, rel=FALSE)
                          {
                               ret <- base::range(getSlot(slot, rel=rel))
                               return(ret)
                          },
                          calculateVelocities = function()
                          {
                               points$vx <<- getDerivative(points$x, points$t)
                               points$vy <<- getDerivative(points$y, points$t)
                          },
                          plotTrack = function(slotX='t', slotY='x', relX=FALSE, relY=TRUE, add=FALSE, withTitle=TRUE, col='black', lwd=1, lty=1, xlab=slotX, ylab=slotY, ...)
                          {
                               xData <- getSlot(slotX, relX)
                               yData <- getSlot(slotY, relY)
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
                                    message(paste('Nothing to plot for track ', id, '.', sep=''))
                               }
                               else
                               {
                                    if(add)
                                    {
                                         
                                         graphics::lines(xData, yData, col=col, lwd=lwd, lty=lty, ...)
                                    }
                                    else
                                    {
                                         if(withTitle)
                                         {
                                              graphics::plot(xData, yData, type='l', col=col, lwd=lwd, lty=lty, xlab=xlab, ylab=ylab, main=as.character(id), ...)
                                         }
                                         else
                                         {
                                              graphics::plot(xData, yData, type='l', col=col, lwd=lwd, lty=lty, xlab=xlab, ylab=ylab, ...)
                                         }
                                    }
                               }
                          },
                          addPoint = function(x, y, t, frame)
                          {
                               points <<- rbind(points, data.frame(frame=frame, t=t, x=x, y=y))
                          },
                          show = function()
                          {
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
                          }
                     )
)


##### Testing 1 #####

# maxima1 <- new('Maxima')
# maxima2 <- maxima1$copy()
# maxima1$initializeWithROI(frame=0, polygon='1,1,1;2,2,2;3,3,3;4,4,4;5,5,5')
# maxima2$initializeWithROI(frame=1, polygon='1,2,5;2,3,4;3,4,3;4,5,2;500,600,1')
# 
# maximaList <- new('MaximaList')
# maximaList$setMaxima(maxima1)
# maximaList$setMaxima(maxima2)
# maximaList$trackBack(startFrame=1)
# 
# trackList <- maximaList$getTrackList(sin=FALSE, fi=2, ff=0.1, tAll=0:1)

##### Testing 2 #####

# maximaList <- new('MaximaList')
# maximaList$initializeWithFile(path=path2)
# mListCopy <- maximaList$copy()
# mListCopy$trackBack(startFrame=501)
# # trackList <- mListCopy$getTrackList(sin=FALSE, fi=2, ff=0.1, tAll=0:515)

##### Testing 3 #####
# 50 ms exposure. 2361 images in 500 s = 4.722 frames / sec
path3 <- '/Volumes/BeebeBig/Jay/JEX Databases/Adhesion FACS/RPMI P-Sel 5Hz-100mHz/Cell_x0_y0/Roi-Tracks Roi/x0_y0.jxd'
path4 <- '/Volumes/BeebeBig/Jay/JEX Databases/Adhesion FACS/RPMI P-Sel 5Hz-100mHz/Cell_x0_y0/Roi-Maxima/x0_y0.jxd'
maximaList <- new('MaximaList')
maximaList$initializeWithFile(path=path4)
mListCopy <- maximaList$copy()
mListCopy$trackBack(startFrame=2360, endFrame=2340)
mListCopy <- MaximaList$new(mListCopy)
mListCopy$trackBackEnd <- 2340
trackList <- mListCopy$getTrackList(sin=TRUE, fi=5, ff=0.1, tAll=seq(0,500,length.out=2361))
track <- trackList$getTrack(1)
track$plotTrack()
trackList$plotTrackList()


# trackListImport <- trackListFromFile(file=path3, sin=TRUE, fi=1, ff=0.1, t=seq(0, 515, 1))
# trackList <- filterTrackList(trackListImport, trackLengthFilter, min=10)
# trackList <- filterTrackList(trackList, trackFrameFilter, endMin=length(trackList@tAll)-1)
# trackList <- sortTrackList(trackList)
# plot(trackList, slot='vx')
# trackList <- setTrackListPhaseShift(object=trackList)
# aTrack <- getTrack(trackList, 2319)
# plot(aTrack, slotY='vx', relY=FALSE)
# plotFit(aTrack, xlab='Time [s]', ylab='Velocity [pixel/s]')
# trackListF <- filterTrackList(trackList, trackLengthFilter, min=4)
