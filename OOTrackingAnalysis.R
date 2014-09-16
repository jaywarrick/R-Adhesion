rm(list=ls())
source('~/.Rprofile')
source('~/Public/DropBox/GitHub/R-Adhesion/Tracking.R')
source('~/Public/DropBox/GitHub/R-Adhesion/OOTracking.R')
setwd('/Users/jaywarrick/Documents/MMB/Projects/Adhesion/R/Testing')
load(file="20140911.Rdata")



trackList <- TrackList$new(trackList)
for(track in trackList$tracks)
{
     trackList$setTrack(Track$new(track))
}
maximaList <- MaximaList$new(maximaList)
for(.maxima in maximaList$maxima)
{
     maximaList$setMaxima(Maxima$new(.maxima))
}

trackList$setValidFrames(fit=bestFit, validStart=0.15, validEnd=0.9)
trackList$smoothVelocities(fit=bestFit, dist=10, maxWidth=25)
trackList$plotTrackList(validOnly=TRUE, slot='vxs', fun=abs, ylim=c(0,50), xlim=c(400,500))
trackMatrix <- trackList$getMatrix(slot='vxs', validOnly=TRUE)
sum(!is.na(trackMatrix))
ret <- list()
for(frame in colnames(trackMatrix))
{
     velocities <- trackMatrix[,frame]
     velocities <- abs(velocities[!is.na(velocities)])
     if(!isempty(velocities))
     {
          adhered <- sum(velocities < 3)/length(velocities)
          if(adhered == 1)
          {
               browser()
          }
          else
          {
               ret[[frame]] <- adhered
          }
     }
     else
     {
          print(frame)
     }
}
times <- trackList$tAll[as.numeric(names(ret))+1]
plot(times, as.numeric(ret), xlab='Time [s]', ylab='Percent Adhered [%]')

aTrack <- Track$new(aTrack)
trackList$setTrack(aTrack)
aTrack$plotTrack(slotY='vx', validOnly=FALSE, type='l', col='blue', pch=20, cex=0.25)
aTrack$plotTrack(slotY='vx', validOnly=FALSE, type='p', col='black', add=TRUE, pch=20, cex=0.5)
aTrack$plotTrack(slotY='vx', validOnly=TRUE, type='p', col='red', add=TRUE, pch=20, cex=0.5)
# length(aTrack$getSlot(slot='x', validOnly=T))
# length(aTrack$getSlot(slot='x', validOnly=F))
# aTrack$length()
# sum(aTrack$getSlot(slot='x', validOnly=T) %in% aTrack$getSlot(slot='x', validOnly=F))

# fit <- bestFit
# validStart = 0.25
# validEnd = 0.75
# 
# sin <- trackList$sin
# fi <- trackList$fi
# ff <- trackList$ff
# tAll <- trackList$tAll
# 
# sweep <- getSweep(amplitude=fit$par['amplitude'], phaseShift=fit$par['phaseShift'], offset=0, sin=sin, fi=fi, ff=ff, tAll=tAll, frames=-1, guess=NULL)
# 
# sin <- FALSE
# fi <- 0.1
# ff <- 0.05
# tAll <- seq(0,1000,1)/10
# validStart = 0.1
# validEnd = 0.95
# sweep <- getSweep(amplitude=1, phaseShift=0, offset=0, sin=sin, fi=fi, ff=ff, tAll=tAll, frames=-1, guess=NULL)
# inflectionsToAddress <- sweep$inflectionNums %in% c(1,3) # These are times at which flow switches directions
# validFrames <- numeric(0)
# for(tIndex in 1:base::length(tAll))
# {
#      # Get the nearest inflection at or beyond this time
#      
#      temp <- which((sweep$inflections >= tAll[tIndex]) & inflectionsToAddress)
#      infIndex <- temp[1]
#      if(is.na(infIndex)) next
#      
#      # Get the bounding inflections that represent changes in fluid direction
#      infT2 <- sweep$inflections[infIndex] # take the inflection we found
#      if((infIndex-2) < 1)
#      {
#           infT1 <- 0
#      }
#      else
#      {
#           infT1 <- sweep$inflections[infIndex-2] # also take two inflections prior because each inflection represents pi/2 and we want to go back to the last change in direction which is pi ago.
#      }
#      dInfT <- infT2-infT1 # define the time interval size between these two inflections
#      
#      # Within the if statement calculate the fractional location of this time index in the interval between the two inflections.
#      if( (tAll[tIndex] >= (infT1 + validStart*dInfT)) & (tAll[tIndex] <= (infT1 + validEnd*dInfT)) )
#      {
#           # If it is within the startValid and endValid bounds, add it to the list of the valid frames
#           validFrames <- c(validFrames, tIndex-1) # (tIndex-1) = frame because frames are indicies that start at 0
#      }
# }
# plot(sweep$t, sweep$v, type='l')
# points(sweep$t[validFrames+1], sweep$v[validFrames + 1], type='p', pch=20, cex=1, col='red')
# sweep$inflectionNums
# 
# ni <- 0
# ti <- first(tAll)
# tf <- last(tAll)
# phi <- 0
# inflections <- (log((log(ff/fi)*phi)/(2*fi*pi*tf)+(log(ff/fi)*ni)/(4*fi*tf)+1)*tf)/log(ff/fi)



##### Testing 1 #####

# maxima1 <- new('Maxima')
# maxima2 <- maxima1$copy()
# maxima3 <- maxima1$copy()
# maxima1$initializeWithROI(frame=0, polygon='1,1,1;2,2,2;3,3,3;4,4,4;51,61,5')
# maxima2$initializeWithROI(frame=1, polygon='2,1,5;3,2,4;4,3,3;5,4,2;50,60,1')
# maxima3$initializeWithROI(frame=2, polygon='1,1,1;2,2,2;3,3,3;4,4,4;5,5,5')
# 
# maximaList <- new('MaximaList')
# maximaList$setMaxima(maxima1)
# maximaList$setMaxima(maxima2)
# maximaList$setMaxima(maxima3)
# maximaList$trackBack(startFrame=2, endFrame=0)
# 
# trackList <- maximaList$getTrackList(sin=FALSE, fi=2, ff=0.1, tAll=0:1)
# maximaList$generateMaximaPlots(path='~/Documents/MMB/Projects/Adhesion/R/Testing/Plots1')

##### Testing 2 #####

# maximaList <- new('MaximaList')
# maximaList$initializeWithFile(path=path2)
# mListCopy <- maximaList$copy()
# mListCopy$trackBack(startFrame=501)
# # trackList <- mListCopy$getTrackList(sin=FALSE, fi=2, ff=0.1, tAll=0:515)

# # ##### Testing 3 #####
# # 50 ms exposure. 2361 images in 500 s = 4.722 frames / sec
# path3 <- '/Volumes/BeebeBig/Jay/JEX Databases/Adhesion FACS/RPMI P-Sel 5Hz-100mHz/Cell_x0_y0/Roi-Tracks Roi/x0_y0.jxd'
# path4 <- '/Volumes/BeebeBig/Jay/JEX Databases/Adhesion FACS/RPMI P-Sel 5Hz-100mHz/Cell_x0_y0/Roi-Maxima/x0_y0.jxd'
# path5 <- '~/Documents/MMB/Projects/Adhesion/R/Testing/SparseMaxima.txt'
# if(!('maximaList' %in% ls()))
# {
#      maximaList <- new('MaximaList')
#      maximaList$initializeWithFile(path=path5)
# } else
# {
#      maximaList <- MaximaList$new(maximaList)
# }
# mListCopy <- maximaList$copy()
# mListCopy$trackBack(startFrame=9994, endFrame=0, maxDist=150, direction=c(1,0,0), directionality=10, uniformityDistThresh=2, digits=1)
# mListCopy$generateMaximaPlots(path='~/Documents/MMB/Projects/Adhesion/R/Testing/Plots1')
# trackList <- mListCopy$getTrackList(sin=FALSE, fi=2, ff=0.01, tAll=seq(0,500,length.out=maximaList$length()))
# trackList$filterTracks(fun = trackLengthFilter, min=500, max=1000000)
# # trackList$filterTracks(fun = trackFrameFilter, startMin=0, startMax=1000000, endMin=maximaList$length()-1, endMax=1000000)
# trackList$plotTrackList()
# bestFit <- getBulkPhaseShift(trackList)
# duh <- getSweep(amplitude=bestFit$par['amplitude'], phaseShift=bestFit$par['phaseShift'], offset=0, sin=trackList$sin, fi=trackList$fi, ff=trackList$ff, tAll=trackList$tAll, frames=-1, guess=NULL)
# lines(duh$t, duh$v, col='blue')
# # aTrack <- trackList$getTrack(0)
# widths <- getWindowWidths(fit=bestFit, trackList=trackList, dist=10, maxWidth=1000)
# # aTrack$smoothVelocities(widths)
# 
# trackList$smoothVelocities(fit=bestFit, dist=10, maxWidth=25)
# trackList$plotTrackList(slot='vxs', xlim=c(50,100))
# lines(duh$t, duh$v, col='blue')
# 
# setwd('/Users/jaywarrick/Documents/MMB/Projects/Adhesion/R/Testing')
# save(list = c('maximaList','trackList', 'bestFit'), file="20140911.Rdata")

##### Testing Sweep #####
# duh <- getSweep(amplitude=100, phaseShift=0, offset=0, sin=FALSE, fi=2, ff=0.01, tAll=seq(0,500,length.out=10001), frames=-1, guess=NULL)
# plot(duh$t, duh$x, col='red', type='l', xlim=c(230, 250))

# ##### Testing 4 #####
# 
# path3 <- '/Volumes/BeebeBig/Jay/JEX Databases/Adhesion FACS/RPMI P-Sel 5Hz-100mHz/Cell_x0_y0/Roi-Tracks Roi/x0_y0.jxd'
# trackList <- new('TrackList')
# trackList$initializeWithFile(file=path3, sin=TRUE, fi=5, ff=0.1, tAll=seq(0,500,length.out=2361))
# trackList$filterTracks(fun = trackLengthFilter, min=50, max=1000000)
# trackList$filterTracks(fun = trackFrameFilter, startMin=0, startMax=1000000, endMin=2360, endMax=1000000)
# bestFit <- getBulkPhaseShift(trackList)
# 
# ##### Testing 5 #####
# 
# path <- "/Volumes/BeebeBig/Jay/JEX Databases/Adhesion FACS/Test/Cell_x0_y0/Roi-Tracks Roi/x0_y0.jxd"
# path2 <- '/Users/jaywarrick/Documents/JEX/Raw Data/LNCaP.arff'
# path3 <- '/Users/jaywarrick/Documents/JEX/LocalTest/PC3 vs LNCaP/Cell_x0_y0/Roi-Tracks Upper/x0_y0.jxd'
# trackList <- new('TrackList')
# trackList$initializeWithFile(file=path3, sin=TRUE, fi=1, ff=0.1, tAll=seq(0, 515, 1))
# trackList$filterTracks(fun = trackLengthFilter, min=50, max=1000000)
# trackList$filterTracks(fun = trackFrameFilter, startMin=0, startMax=1000000, endMin=515, endMax=1000000)
# bestFit <- getBulkPhaseShift(trackList)

