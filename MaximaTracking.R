##### Maxima #####

setClass('maxima', representation(frame='numeric', points='data.frame'), prototype(frame=0, points=data.frame()))

setMethod('initialize', 'maxima', function(.Object, frame, polygon=NULL)
{
     .Object@frame <- as.numeric(as.character(frame))
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
          
          .Object@points <- data.frame(id=index, x=x, y=y)
     }
     else
     {
          .Object@points <- data.frame(id=numeric(0), x=numeric(0), y=numeric(0))
     }
     
     .Object
})

setMethod('plot', 'maxima', function(x, IDs=TRUE, pch=20, xlab='X [pixel]', ylab='Y [pixel]', cex=0.5, ...)
{
     plot(x@points$x, x@points$y, pch=20, xlab=xlab, ylab=ylab, cex=cex, ...)
     if(IDs)
     {
          text(x=x@points$x, y=x@points$y, labels=as.character(x@points$id), cex=0.5, adj=c(0,0))
     }
})

setMethod('length', 'maxima', function(x)
{
     return(length(x@points$frame))
})

setGeneric('getPoints', function(object) standardGeneric('getPoints')) 
setMethod('getPoints', 'maxima', function(object)
{
     return(object@points)
})

setGeneric('getPoint', function(object, index=NULL, id=NULL) standardGeneric('getPoint')) 
setMethod('getPoint', 'maxima', function(object, index, id)
{
     if(is.null(index) & is.null(id))
     {
          return(NULL)
     }
     else if(!is.null(index))
     {
          return(object@points[index,])
     }
     else
     {
          return(object@points[object@points$id==id,])
     }
})

setGeneric('addPoint', function(object, id, x, y) standardGeneric('addPoint')) 
setMethod('addPoint', 'maxima', function(object, id, x, y)
{
     object@points <- rbind(object@points, data.frame(id=id, x=x, y=y))
     return(object)
})

setGeneric('getXYZ', function(object) standardGeneric('getXYZ')) 
setMethod('getXYZ', 'maxima', function(object)
{
     data <- cbind(object@points, data.frame(z=0))
     return(data[,2:4])
})

##### Maxima List #####

setClass('maximaList', representation(maxima='list', maximaTable='data.frame', maxID='numeric'))

setMethod('initialize', 'maximaList', function(.Object, maximaFilePath=NULL, timeDimName='Time')
{
     library(foreign)
     if(!is.null(maximaFilePath))
     {
          maximaFile <- read.arff(maximaFilePath)
          maximaFile2 <- reorganizeTable(maximaFile, nameCol='Metadata')
          maximaFile2[,timeDimName] <- as.numeric(as.character(maximaFile2[,timeDimName]))
          .Object@maximaTable <- maximaFile2
          maxima <- list()
          for(r in 1:nrow(maximaFile2))
          {
               newMaxima <- new('maxima', frame=maximaFile2[r,timeDimName], polygon=maximaFile2[r,'polygonPts'])
               maxima[[as.character(maximaFile2[r,timeDimName])]] <- newMaxima
          }
          maxima <- maxima[order(as.numeric(names(maxima)), decreasing=FALSE)]
          .Object@maxima <- maxima
     }
     else
     {
          .Object@maxima <- list()
     }
     
     
     .Object@maxID <- -1 # Used for appending maxima
     
     .Object
})

setGeneric('getMaxima', function(object, frame=0) standardGeneric('getMaxima')) 
setMethod('getMaxima', 'maximaList', function(object, frame)
{
     return(object@maxima[[as.character(frame)]])
})

setGeneric('getLastMaxima', function(object) standardGeneric('getLastMaxima')) 
setMethod('getLastMaxima', 'maximaList', function(object)
{
     n <- max(as.numeric(names(object@maxima)))
     return(object@maxima[[as.character(n)]])
})

setGeneric('getFirstMaxima', function(object) standardGeneric('getFirstMaxima')) 
setMethod('getFirstMaxima', 'maximaList', function(object)
{
     n <- min(as.numeric(names(object@maxima)))
     return(object@maxima[[as.character(n)]])
})

setGeneric('setMaxima', function(object, maxima) standardGeneric('setMaxima')) 
setMethod('setMaxima', 'maximaList', function(object, maxima)
{
     object@maxima[[as.character(maxima@frame)]] <- maxima
     return(object)
})

setGeneric('getNextLinkedMaxima', function(newObject, initialObject, maxDist=150, direction=c(1,0,0), directionality=10, uniformityDistThresh=15, digits=1) standardGeneric('getNextLinkedMaxima')) 
setMethod('getNextLinkedMaxima', 'maximaList', function(newObject, initialObject, maxDist, direction, directionality, uniformityDistThresh, digits)
{
     currentMaxima <- getFirstMaxima(newObject)
     maxID <- newObject@maxID
     maxima <- getMaxima(initialObject, frame=(currentMaxima@frame-1))
     if(is.null(maxima))
     {
          return(NULL)
     }
     t1 <- getXYZ(currentMaxima)
     t0 <- getXYZ(maxima) # maxima to prepend
     results <- DirectionalLinearAssignment(points=list(t0=t0, t1=t1), maxDist=maxDist, direction=direction, directionality=directionality, uniformityDistThresh=uniformityDistThresh, digits=digits)
     px <- results$px
     newMaxima <- new('maxima', frame=maxima@frame)
     for(i in 1:length(px))
     {
          if(px[i] <= nrow(t0) & i <= nrow(t1))
          {
               # Then we have a link between t0 and t1 and the points should share ids
               p0 <- getPoint(maxima, index=px[i])
               p1 <- getPoint(currentMaxima, index=i)
               newMaxima <- addPoint(newMaxima, id=p1$id, x=p0$x, y=p0$y)
          }
          else if(px[i] <= nrow(t0) & i > nrow(t1))
          {
               # This point in t0 represents the start of a new track
               # increment the maxID and add the point to the newMaxima
               maxID <- maxID + 1
               p0 <- getPoint(maxima, index=px[i])
               newMaxima <- addPoint(newMaxima, id=maxID, x=p0$x, y=p0$y)
          }
          else
          {
               # This point in t1 represents the end of a track or an auxillary match to enable other types of matches
               # We don't need to do anything
          }
     }
     return(list(newMaxima=newMaxima, maxID=maxID))
})


