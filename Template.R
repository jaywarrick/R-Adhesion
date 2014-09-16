# Clear the workspace of any data, variables, and functions that are currently loaded.
rm(list=ls())
<<<<<<< HEAD

#source("http://bioconductor.org/biocLite.R")
#biocLite('GraphAlignment')
# Load the functions that we'll be using.
<<<<<<< HEAD
# source('~/Public/DropBox/GitHub/R-Adhesion/Tracking.R')
# source('~/Public/DropBox/GitHub/R-Adhesion/OOTracking.R')
=======
#MMMMMINE
#source("http://bioconductor.org/biocLite.R")
#biocLite('GraphAlignment')
# Load the functions that we'll be using.
>>>>>>> FETCH_HEAD
=======
source('~/Public/DropBox/GitHub/R-Adhesion/Tracking.R')
source('~/Public/DropBox/GitHub/R-Adhesion/OOTracking.R')
>>>>>>> parent of 0cfcc93... Create comments to switch between Liz and Jay
source('C:/Users/Elizabeth/Documents/GitHub/R-Adhesion/Tracking.R')
source('C:/Users/Elizabeth/Documents/GitHub/R-Adhesion/OOTracking.R')

# Create a MaximaList object to help us with a bunch of stuff
maximaList <- new('MaximaList')

# The MaximaList knows how to read a JEX ROI file, so let's read it in to initialize the dataset
<<<<<<< HEAD
<<<<<<< HEAD
# maximaList$initializeWithFile(path='~/Documents/MMB/Projects/Adhesion/R/Testing/SparseMaxima.txt', timeDimName='Time')
=======
>>>>>>> FETCH_HEAD
=======
maximaList$initializeWithFile(path='~/Documents/MMB/Projects/Adhesion/R/Testing/SparseMaxima.txt')
>>>>>>> parent of 0cfcc93... Create comments to switch between Liz and Jay
maximaList$initializeWithFile(path='C:/Users/Elizabeth/Desktop/x0_y0 points.txt', timeDimName='T')

# Copy it and do stuff, keeping the original around in case we wnat to redo stuff without having to read in the dataset again.
mListCopy <- maximaList$copy()

# Track the cells, backward in time starting at the frame called "startFrame" (latest frame) and ending at "endFrame" (the earliest frame)
# FYI Frame numbers START AT 0 while in R indices START AT 1.
<<<<<<< HEAD
<<<<<<< HEAD
# mListCopy$trackBack(startFrame=5000, endFrame=4600, maxDist=150, direction=c(1,0,0), directionality=10, uniformityDistThresh=2, digits=1)
mListCopy$trackBack(startFrame=479, endFrame=0, maxDist=150, direction=c(1,0,0), directionality=10, uniformityDistThresh=2, digits=1)

# Plot all the maxima to pdf plots so that you can flip through pdfs and see where cells are when (good for trouble shooting but not necessary)
# mListCopy$generateMaximaPlots(path='~/Documents/MMB/Projects/Adhesion/R/Testing/Plots1')
mListCopy$generateMaximaPlots(path='C:/Users/Elizabeth/Documents/MMB/Projects/Adhesion/R/Testing/Plots1')
=======
mListCopy$trackBack(startFrame=479, endFrame=0, maxDist=150, direction=c(1,0,0), directionality=10, uniformityDistThresh=2, digits=1)

# Plot all the maxima to pdf plots so that you can flip through pdfs and see where cells are when (good for trouble shooting but not necessary)
#mListCopy$generateMaximaPlots(path='C:/Users/Elizabeth/Documents/MMB/Projects/Adhesion/R/Testing/Plots1')
>>>>>>> FETCH_HEAD
=======
mListCopy$trackBack(startFrame=5000, endFrame=4600, maxDist=150, direction=c(1,0,0), directionality=10, uniformityDistThresh=2, digits=1)
mListCopy$trackBack(startFrame=479, endFrame=0, maxDist=150, direction=c(1,0,0), directionality=10, uniformityDistThresh=2, digits=1)

# Plot all the maxima to pdf plots so that you can flip through pdfs and see where cells are when (good for trouble shooting but not necessary)
mListCopy$generateMaximaPlots(path='~/Documents/MMB/Projects/Adhesion/R/Testing/Plots1')
#mListCopy$generateMaximaPlots(path='C:/Users/Elizabeth/Documents/MMB/Projects/Adhesion/R/Testing/Plots1')
>>>>>>> parent of 0cfcc93... Create comments to switch between Liz and Jay

# Offset frames to actual start frame of the imported dataset.
mListCopy$offsetFrames(offset=3991)

# Now that the cells are tracked, we want to have the data reorganized into a list of tracks (i.e., a TrackList)
<<<<<<< HEAD
<<<<<<< HEAD
trackList <- mListCopy$getTrackList(sin=FALSE, fi=2, ff=0.01, tAll=seq(0,500,length.out=10001))
=======
trackList <- mListCopy$getTrackList(sin=FALSE, fi=2, ff=0.01, tAll=seq(0,500,length.out=max(as.numeric(names(mListCopy$maxima)) + 1)))
>>>>>>> FETCH_HEAD
=======
trackList <- mListCopy$getTrackList(sin=FALSE, fi=2, ff=0.01, tAll=seq(0,500,length.out=max(as.numeric(names(maximaList$maxima)) + 1)))
>>>>>>> parent of 0cfcc93... Create comments to switch between Liz and Jay

# Get rid of short tracks (i.e., cells that are hard too hard to follow for long periods of time like ones that go on and off screen)
trackList$filterTracks(fun = trackLengthFilter, min=50, max=1000000)

# We could get rid of other tracks based on when they start or stop with this filter but for now we'll include all tracks no matter when they start or stop by skipping this step.
# trackList$filterTracks(fun = trackFrameFilter, startMin=0, startMax=1000000, endMin=maximaList$length()-1, endMax=1000000)

# Plot the tracks to get a sense of things
trackList$plotTrackList()

# Fit all the data points with a single curve to determine the phaseShift of the cells
bestFit <- getBulkPhaseShift(trackList)

# Smooth the velocities (i.e., average over multiple frames) to get a more accurate measure, especially for slow moving cells
# We need the bestFit information to estimate the speed of cells over time
# 'dist' refers to the distance in pixels that we want to see the cell travel before quantifying its velocity
# 'maxWidth' refers to the maximum number of frames to average together to estimate velocity (i.e., it's great to try and use more frames but we don't want to use too many)
trackList$smoothVelocities(fit=bestFit, dist=10, maxWidth=25)

# Given the flow switches direction and attached cells can 'wobble', we only want to gage whether as cell is adhered
# after it is done wobbling. To do this we set which frames are 'valid' for determining whether the cell is adhered or not.
# We need the 'bestFit' to estimate where the chages in flow direction are.
# We then want to wait before quantifying because the cell could be wobbling so 'validStart' dictates when
# after the change in flow direction we can see if the cell is adhered or not. 'validStart' is expressed as a 
# fraction of the interval between changed in flow direction (i.e., validStart=0.1 says to wait 10% of the time between direction changes to mark time frames as valid)
# 'validEnd' is also a fraction of the same interval. Frames between 'validStart' and 'validEnd' within each interval bounded by changes in flow direction
# are recorded as being valid.
trackList$setValidFrames(fit=bestFit, validStart=0.15, validEnd=0.9)

# Get and plot the percent of cells adhered over time
results = trackList$getPercentAdhered(velocityThreshold=3)
plot(results$time, results$percentAdhered, xlab='Time [s]', ylab='Percent Adhered [%]')


