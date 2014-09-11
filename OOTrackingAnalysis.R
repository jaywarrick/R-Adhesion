rm(list=ls())
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
