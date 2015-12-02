graphics.off()
rm(list=ls())
library(foreign)
library(pracma)
library(methods)
library(data.table)
source('/Users/jaywarrick/Public/DropBox/GitHub/R-Adhesion/CellTracking/R/LAPJV.R')
source('/Users/jaywarrick/Public/DropBox/GitHub/R-Adhesion/CellTracking/R/Maxima.R')
source('/Users/jaywarrick/Public/DropBox/GitHub/R-Adhesion/CellTracking/R/MaximaList.R')
source('/Users/jaywarrick/Public/DropBox/GitHub/R-Adhesion/CellTracking/R/PackageFunctions.R')
source('/Users/jaywarrick/Public/DropBox/GitHub/R-Adhesion/CellTracking/R/TrackList.R')
source('/Users/jaywarrick/Public/DropBox/GitHub/R-Adhesion/CellTracking/R/Track.R')
source('/Users/jaywarrick/Public/DropBox/GitHub/R-Adhesion/CellTracking/R/TrackFilters.R')
source('/Users/jaywarrick/Public/DropBox/GitHub/R-Adhesion/CellTracking/R/TrackFitting.R')
source('/Users/jaywarrick/Public/DropBox/GitHub/R-Adhesion/CellTracking/R/Tracking.R')
source('~/.Rprofile')
source('~/Desktop/A Sandbox/Test.R')

getData <- function(db, ds, x, y, type, name)
{
     ret <- list()
     ret$db <- db
     ret$tmp <- '/Users/jaywarrick/Desktop/A Sandbox/TrackingOutput/' #file.path(db,'temp','RScriptTempFolder')
     ret$ds <- ds
     ret$x <- x
     ret$y <- y
     ret$type <- type
     ret$name <- name
     ret$value <- read.arff(file.path(db, ds, paste0('Cell_x',x,'_y',y), paste0(type,'-',name), paste0('x',x,'_y',y,'.jxd')))
     return(ret)
}

getRangeX <- function(x)
{
     temp <- range(x$getSlot(slot='x', rel=TRUE))
     return(temp[2]-temp[1])
}
getRangeY <- function(x)
{
     temp <- range(x$getSlot(slot='y', rel=TRUE))
     return(temp[2]-temp[1])
}

analyze <- function(data)
{
     mList <- new('MaximaList')
     mList$initializeWithROIDataFrame(roiTable=data$value, timeDimName='Wash')
     mList$trackBack(startFrame=4, endFrame=0, maxDist=4, directionality=1, uniformityDistThresh=-1, digits=1)
     tList <- mList$getStandardTrackList(timePerFrame=1)
     tList2 <- tList$copy()
     tList2$filterTracks(fun=trackLengthFilter, min=5, max=5)

     rangeX <- mean(unlist(tList2$getProp(fun=getRangeX)))
     rangeY <- mean(unlist(tList2$getProp(fun=getRangeY)))
     rangeSDX <- sd(unlist(tList2$getProp(fun=getRangeX)))
     rangeSDY <- sd(unlist(tList2$getProp(fun=getRangeY)))
     quant99X <- quantile(unlist(tList2$getProp(fun=getRangeX)), .99)
     quant99Y <- quantile(unlist(tList2$getProp(fun=getRangeX)), .99)

     ret <- data.frame(ROI=as.factor(data$name), Measurement=c('rangeX','rangeY','rangeSDX','rangeSDY','quant99X','quant99Y'), Value=c(rangeX, rangeY, rangeSDX, rangeSDY, quant99X, quant99Y))
     path <- file.path(data$tmp, paste0('x',data$x,'_y',data$y,'_',data$name,'.arff'))
     write.arff(ret, file=path)
}

for(i in 0:2)
{
     for(dName in c('Maxima (center)', 'Maxima (side 1)', 'Maxima (side 2)'))
     {
          data <- getData(db='/Volumes/Data/JEX Databases/MM Registration Expt', ds='20150922', x=0, y=i, type='Roi', name=dName)
          analyze(data)
     }
}

results <- NULL
for(i in 0:2)
{
     for(dName in c('Maxima (center)', 'Maxima (side 1)', 'Maxima (side 2)'))
     {
          path <- file.path('/Users/jaywarrick/Desktop/A Sandbox/TrackingOutput', paste0('x0_y',i,'_',dName,'.arff'))
          print(path)
          if(is.null(results))
          {
               results <- read.arff(path)
               results$well <- i
          }
          else
          {
               temp <- read.arff(path)
               temp$well <- i
               results <- rbind(results, temp)
          }
     }
}
results <- data.table(results)
results2 <- data.table(reorganizeTable(data.frame(results), valueCol=c('Value'), convertToNumeric = FALSE))
summaryResults <- results2[,list(meanX=mean(rangeX), meanY=mean(rangeY), sdX=mean(rangeSDX), sdY=mean(rangeSDY)),by=c("well")]
write.table(x=results, file=file.path('/Users/jaywarrick/Desktop/A Sandbox/TrackingOutput','results.txt'))
write.table(x=summaryResults, file=file.path('/Users/jaywarrick/Desktop/A Sandbox/TrackingOutput','summaryResults'))

