arg <- commandArgs(TRUE)[1]

library(methods)
library(curl)

sourceGitHubFile <- function(user, repo, branch, file)
{
	require(curl)
	destfile <- tempfile()
	fileToGet <- paste0("https://raw.githubusercontent.com/", user, "/", repo, "/", branch, "/", file)
	curl_download(url=fileToGet, destfile)
	source(destfile)
}

sourceGitHubFile('jaywarrick','R-General','master','.Rprofile')
sourceGitHubFile('jaywarrick','R-Adhesion','master','CellTracking/R/LAPJV.R')
sourceGitHubFile('jaywarrick','R-Adhesion','master','CellTracking/R/MaximaList.R')
sourceGitHubFile('jaywarrick','R-Adhesion','master','CellTracking/R/Maxima.R')
sourceGitHubFile('jaywarrick','R-Adhesion','master','CellTracking/R/PackageFunctions.R')
sourceGitHubFile('jaywarrick','R-Adhesion','master','CellTracking/R/TrackFilters.R')
sourceGitHubFile('jaywarrick','R-Adhesion','master','CellTracking/R/TrackFitting.R')
sourceGitHubFile('jaywarrick','R-Adhesion','master','CellTracking/R/Tracking.R')
sourceGitHubFile('jaywarrick','R-Adhesion','master','CellTracking/R/Track.R')
sourceGitHubFile('jaywarrick','R-Adhesion','master','CellTracking/R/TrackList.R')

infile <- "maxima.jxd"
outfile1 <- "maximalist.Rdata"
outfile2 <- "maximaROI.jxd"
outfile3 <- "maximatracklist.Rdata"
outfile4 <- "maximatracklistFiltered.Rdata"
outfile5 <- "maximaRTrackListOS.Rdata"
outfile6 <- "maximaBulkFitResults.csv"


#Step 1 filter track oscillatory
maximaList <- new('MaximaList')
maximaList$initializeWithJEXROIFile(path=infile, timeDimName='Time')

startFrame <- max(as.numeric(names(maximaList$maxima)))

maximaList$trackBack(startFrame=startFrame, endFrame=1, maxDist=150, direction=c(1,0,0), directionality=10, uniformityDistThresh=2, digits=1)

maximaList$saveROI(file=outfile2)
trackList <- maximaList$getTrackList(t0_Frame=1, timePerFrame=0.035)
#save(list=c('trackList'), file=outfile3)

#save(list=c('maximaList'), file=outfile1)

#Step 2 Filter Tracklist
trackList$filterTracks(fun = trackLengthFilter, min=3, max=10000000)
trackList$filterTracks(fun = trackFrameFilter, startMin=0, startMax=1000000, endMin=0, endMax=1000000)
#save(list=c('trackList'), file=outfile4)

#Step 3 Apply Oscillatory Track Metadata
t0_Frame <- trackList$meta$t0_Frame
timePerFrame <- trackList$meta$timePerFrame
trackList$setOscillatoryMeta(sin=FALSE, fi=1.0, ff=0.01, sweepDuration=300,t0_Frame=t0_Frame, timePerFrame=timePerFrame)
bestFit <- getBulkPhaseShiftGS(trackList, cores=4)
trackList$meta$bestFit <- bestFit
write.csv(as.data.frame(as.list(bestFit$par)), file=outfile6, row.names=FALSE)
trackList$smoothVelocities(fit=bestFit, dist=10, maxWidth=15)
trackList$calculateValidFrames(fit=bestFit, validStart=0.2, validEnd=0.8)
save(list=c('trackList'), file=outfile5)