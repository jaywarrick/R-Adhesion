graphics.off()
rm(list=ls())
library('foreign')
jexTempRFolder <- '/Volumes/Data/JEX Databases/MM Registration Expt/temp/RScriptTempFolder'
jexDBFolder <- '/Volumes/Data/JEX Databases/MM Registration Expt'
data1 <- list()
data1$type <- 'Roi'
data1$name <- 'Maxima (center)'
data1$value <- read.arff('/Volumes/Data/JEX Databases/MM Registration Expt/20150922/Cell_x0_y0/Roi-Maxima (center)/x0_y0.jxd')

library(foreign)
library(pracma)
library(methods)
source('/Users/jaywarrick/Public/DropBox/GitHub/R-Adhesion/CellTracking/R/LAPJV.R')
source('/Users/jaywarrick/Public/DropBox/GitHub/R-Adhesion/CellTracking/R/Maxima.R')
source('/Users/jaywarrick/Public/DropBox/GitHub/R-Adhesion/CellTracking/R/MaximaList.R')
source('/Users/jaywarrick/Public/DropBox/GitHub/R-Adhesion/CellTracking/R/PackageFunctions.R')
source('/Users/jaywarrick/Public/DropBox/GitHub/R-Adhesion/CellTracking/R/TrackList.R')
source('/Users/jaywarrick/Public/DropBox/GitHub/R-Adhesion/CellTracking/R/Track.R')
source('/Users/jaywarrick/Public/DropBox/GitHub/R-Adhesion/CellTracking/R/TrackFilters.R')
source('/Users/jaywarrick/Public/DropBox/GitHub/R-Adhesion/CellTracking/R/TrackFitting.R')
source('/Users/jaywarrick/Public/DropBox/GitHub/R-Adhesion/CellTracking/R/Tracking.R')

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

# data1 <- list()
# data1$type <- 'Roi'
# data1$name <- 'Maxima (center)'
# data1$value <- read.arff('/Volumes/Data/JEX Databases/MM Registration Expt/20150922/Cell_x0_y0/Roi-Maxima (center)/x0_y0.jxd')
mList <- new('MaximaList')
mList$initializeWithROIDataFrame(roiTable=data1$value, timeDimName='Wash')
mList$trackBack(startFrame=4, endFrame=0, maxDist=4, directionality=1, uniformityDistThresh=-1, digits=1)
tList <- mList$getStandardTrackList(timePerFrame=1)
tList2 <- tList$copy()
tList2$filterTracks(fun=trackLengthFilter, min=5, max=5)

rangeX <- mean(unlist(tList2$getProp(fun=getRangeX)))
rangeY <- mean(unlist(tList2$getProp(fun=getRangeY)))
rangeSDX <- sd(unlist(tList2$getProp(fun=getRangeX)))
rangeSDY <- sd(unlist(tList2$getProp(fun=getRangeY)))

ret <- data.frame(ROI=as.factor(data1$name), Measurement=c('rangeX','rangeY','rangeSDX','rangeSDY'), Value=c(rangeX, rangeY, rangeSDX, rangeSDY))
path <- file.path(jexTempRFolder, 'myOutFile.arff')
jexTempRFolder <- '/Users/jaywarrick/Desktop/A Sandbox'
write.arff(ret, file=path)
fileList1 <- c(path)
