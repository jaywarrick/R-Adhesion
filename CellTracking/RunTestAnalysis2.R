# Clear the workspace of any data, variables, and functions that are currently loaded.
#rm(list=ls())

#source("http://bioconductor.org/biocLite.R")
#biocLite('GraphAlignment')
# Load the functions that we'll be using.
#source('C:/Users/biomems/Documents/GitHub/R-Adhesion/Tracking.R')
#source('C:/Users/biomems/Documents/GitHub/R-Adhesion/OOTracking.R')
source('~/.Rprofile')
source('~/Public/DropBox/GitHub/R-Adhesion/CellTracking/R/LAPJV.R')
source('~/Public/DropBox/GitHub/R-Adhesion/CellTracking/R/Tracking.R')
source('~/Public/DropBox/GitHub/R-Adhesion/CellTracking/R/TrackList.R')
source('~/Public/DropBox/GitHub/R-Adhesion/CellTracking/R/Track.R')
source('~/Public/DropBox/GitHub/R-Adhesion/CellTracking/R/TrackList.R')
source('~/Public/DropBox/GitHub/R-Adhesion/CellTracking/R/Track.R')
source('~/Public/DropBox/GitHub/R-Adhesion/CellTracking/R/Maxima.R')
source('~/Public/DropBox/GitHub/R-Adhesion/CellTracking/R/MaximaList.R')
source('~/Public/DropBox/GitHub/R-Adhesion/CellTracking/R/TrackFilters.R')
source('~/Public/DropBox/GitHub/R-Adhesion/CellTracking/R/TrackFitting.R')
source('~/Public/DropBox/GitHub/R-Adhesion/CellTracking/R/PackageFunctions.R')

setwd('~/Public/DropBox/GitHub/R-Adhesion/CellTracking/R/Test')

#well <- 'x0_y8'

#path <- '/Users/jaywarrick/Documents/MMB/Projects/Adhesion/MM AMD/'

# Create a MaximaList object to help us with a bunch of stuff
maximaList <- new('MaximaList')

# The MaximaList knows how to read a JEX ROI file, so let's read it in to initialize the dataset
# maximaList$initializeWithFile(path='~/Documents/MMB/Projects/Adhesion/R/Testing/SparseMaxima.txt', timeDimName='Time')
maximaList$initializeWithJEXROIFile(path=paste0(path, well, '.jxd'), timeDimName='Time')

# Copy it and do stuff, keeping the original around in case we want to redo stuff without having to read in the dataset again.
mListCopy <- maximaList$copy()

# mListCopy$offsetFrames(-1) # This makes frame 0 = time 0 give import method started index at 1

# Track the cells, backward in time starting at the frame called "startFrame" (latest frame) and ending at "endFrame" (the earliest frame)
# FYI Frame numbers START AT 0 while in R indices START AT 1.
# mListCopy$trackBack(startFrame=5000, endFrame=4600, maxDist=150, direction=c(1,0,0), directionality=10, uniformityDistThresh=2, digits=1)
mListCopy$trackBack(startFrame=max(as.numeric(names(maximaList$maxima))), endFrame=1, maxDist=150, direction=c(1,0,0), directionality=10, uniformityDistThresh=2, digits=1)

# Plot all the maxima to pdf plots so that you can flip through pdfs and see where cells are when (good for trouble shooting but not necessary)
# mListCopy$generateMaximaPlots(path='~/Documents/MMB/Projects/Adhesion/R/Testing/Plots1')
mListCopy$generateMaximaPlots(path='~/Public/DropBox/GitHub/R-Adhesion/CellTracking/Plots')

# Offset frames to actual start frame of the imported dataset.
# Don't need to for this set...
# mListCopy$offsetFrames(offset=3991)

# Now that the cells are tracked, we want to have the data reorganized into a list of tracks (i.e., a TrackList)
trackList <- mListCopy$getTrackList(sin=FALSE, fi=1, ff=0.01, sweepDuration=300, t0_Frame=0, timePerFrame=0.035)

# Get rid of short tracks (i.e., cells that are too hard to follow for long periods of time like ones that go on and off screen)
trackList$filterTracks(fun = trackLengthFilter, min=20, max=1000000)

# We could get rid of other tracks based on when they start or stop with this filter but for now we'll include all tracks no matter when they start or stop by skipping this step.
# trackList$filterTracks(fun = trackFrameFilter, startMin=0, startMax=1000000, endMin=maximaList$length()-1, endMax=1000000)

# Load from a saved thingy
#load('/Users/jaywarrick/Public/DropBox/GitHub/R-Adhesion/CellTracking/CellTracking.RData')

# Fit all the data points with a single curve to determine the phaseShift of the cells
#bestFit <- getBulkPhaseShift2(trackList, tiGuess=0)
bestFit <- getBulkPhaseShift(trackList, tiGuess=0)

# Smooth the velocities (i.e., average over multiple frames) to get a more accurate measure, especially for slow moving cells
# We need the bestFit information to estimate the speed of cells over time
# 'dist' refers to the distance in pixels that we want to see the cell travel before quantifying its velocity
# 'maxWidth' refers to the maximum number of frames to average together to estimate velocity (i.e., it's great to try and use more frames but we don't want to use too many)
trackList$smoothVelocities(fit=bestFit, dist=10, maxWidth=15)

# Given the flow switches direction and attached cells can 'wobble', we only want to gage whether as cell is adhered
# after it is done wobbling. To do this we set which frames are 'valid' for determining whether the cell is adhered or not.
# We need the 'bestFit' to estimate where the chages in flow direction are.
# We then want to wait before quantifying because the cell could be wobbling so 'validStart' dictates when
# after the change in flow direction we can see if the cell is adhered or not. 'validStart' is expressed as a
# fraction of the interval between changed in flow direction (i.e., validStart=0.1 says to wait 10% of the time between direction changes to mark time frames as valid)
# 'validEnd' is also a fraction of the same interval. Frames between 'validStart' and 'validEnd' within each interval bounded by changes in flow direction
# are recorded as being valid.
trackList$setValidFrames(fit=bestFit, validStart=0.1, validEnd=0.9)

# update the trackList upon editing a method
#trackList <- TrackList$new( trackList )

######### Figure out why the percent adhered isn't working.... The time column is messed up...

# Plot the tracks to get a sense of things
#trackList$getTrack(1)$plotTrack()
trackList$plotTrackList(slot='x', rel=TRUE, ylim=c(-250,250), validOnly=FALSE, xlim=c(0,300), ylab='Position [pixels]', xlab='t [s]')
trackList$plotTrackList(slot='vx', rel=TRUE, ylim=c(-500,500), validOnly=FALSE, xlim=c(0,300), ylab='Velocity [pixels/sec]', xlab='t [s]')

# f <- trackList$fi * (trackList$ff/trackList$fi)^(trackList$tAll/trackList$sweepDuration)
# A = 150*(6.45/4)*1e-6
# mu <- 0.00078
# h <- 200e-6
# tau = 8*pi*A*f*mu/h
# plot(trackList$tAll, 10*tau, xlab='t [s]', ylab='Shear Stress [dynes/cm^2]')

fitCurveData <- getSweep(amplitude=bestFit$par[['amplitude']], phaseShift=bestFit$par[['phaseShift']], sweepDuration=(1/bestFit$par[['timeScalingFactor']])*trackList$sweepDuration, offset=0, sin=trackList$sin, ti=bestFit$par[['ti']], fi=trackList$fi, ff=trackList$ff, t=trackList$tAll, guess=NULL)
trackList$plotTrackList(slot='vx', rel=TRUE, ylim=c(-500,500), validOnly=FALSE, xlim=c(0,5))
lines(fitCurveData$t, fitCurveData$v, col='red')
trackList$plotTrackList(slot='vx', rel=TRUE, ylim=c(-50,50), validOnly=FALSE, xlim=c(125,300))
lines(fitCurveData$t, fitCurveData$v, col='red')
trackList$plotTrackList(slot='vx', rel=TRUE, ylim=c(-500,500), validOnly=FALSE, xlim=c(0,300))
lines(fitCurveData$t, fitCurveData$v, col='red')
#fitCurveData <- getSweep(amplitude=150, phaseShift=bestFit$par[['phaseShift']], offset=0, sin=trackList$sin, fi=trackList$fi, ff=trackList$ff, sweepDuration=trackList$sweepDuration, tAll=trackList$tAll, frames=-1, guess=NULL)

# Get and plot the percent of cells adhered over time
results = trackList$getPercentAdhered(velocityThreshold=5)
#pdf(file=paste0('/Users/jaywarrick/Documents/MMB/Projects/Adhesion/MM AMD/',well,'_PercentAdhered.pdf'), width=6, height=4)
plot(results$time, results$percentAdhered, xlab='Time [s]', ylab='Percent Adhered [%]', pch=20, cex=0.75, ylim=c(0,70))
#dev.off()

# Now just grab all tracks that can be tracked from the last frame backward continuously
contTrackList <- trackList$copy()
tempLastFrame <- max(trackList$allFrames)
# Get rid of short tracks (i.e., cells that are too hard to follow for long periods of time like ones that go on and off screen)
contTrackList$filterTracks(fun = trackFrameFilter, endMin=tempLastFrame)
results = contTrackList$getPercentAdhered(velocityThreshold=5)
#pdf(file=paste0('/Users/jaywarrick/Documents/MMB/Projects/Adhesion/MM AMD/',well,'_PercentAdhered.pdf'), width=6, height=4)
plot(results$time, results$percentAdhered, xlab='Time [s]', ylab='Percent Adhered [%]', pch=20, cex=0.75, ylim=c(0,70))
#dev.off()
trackList <- trackList$copy()
contTrackList <- contTrackList$copy()
trackList$save(objectName=well, file=paste0(path,well,'.RData'))
contTrackList$save(objectName=well, file=paste0(path,well,'_c','.RData'))
