# Adhesion Analysis with Liz
rm(list=ls())
source('~/.Rprofile')
source('/Users/jaywarrick/Public/DropBox/GitHub/R-Adhesion/AdhesionAnalysis_HelperFunctions.R')

##### Analysis #####

# Remember cell 328
path <- '/Volumes/BeebeBig/Jay/JEX Databases/Adhesion FACS/Test/Cell_x0_y0/Roi-Tracks/x0_y0.jxd'
path2 <- '/Volumes/BeebeBig/Jay/JEX Databases/Adhesion FACS/RPMI P-Sel 5Hz-100mHz/Cell_x0_y0/Roi-Tracks Roi/x0_y0.jxd'
path3 <- '/Users/jaywarrick/Documents/JEX/Raw Data/LNCaP.arff'

# Determine samplingFreq if necessary
if(is.na(samplingFreq))
{
     temp <- seq(ti, tf, length.out=length.out)
     samplingFreq <- 1 / (temp[2]-temp[1])
}

trackListImport <- trackListFromFile(file=path3, sin=TRUE, fi=1, ff=0.1, t=seq(0, 515, 1))
trackList <- filterTrackList(trackListImport, trackLengthFilter, min=10)
trackList <- filterTrackList(trackList, trackFrameFilter, endMin=length(trackList@tAll)-1)
trackList <- sortTrackList(trackList)
plot(trackList, slot='vx')
trackList <- setTrackListPhaseShift(object=trackList)
aTrack <- getTrack(trackList, 2319)
plot(aTrack, slotY='vx', relY=FALSE)
plotFit(aTrack, xlab='Time [s]', ylab='Velocity [pixel/s]')
trackListF <- filterTrackList(trackList, trackLengthFilter, min=4)

validTrack <- setValidTimeIndices(aTrack, applyWobbleFilter=T)
plotFit(validTrack)

setwd('/Users/jaywarrick/Desktop/AAASandbox/plots')
for(track in trackList@tracks)
{
    print(track@id)
    pdf(paste(getwd(),'/',as.character(track@id),'.pdf',sep=''), width=4, height=3)
    fitFixedPhase <- getTrackSinFitFixedPhase(object=aTrack, slot='vx', phaseShift=trackList@phaseShift, guess=getTrackParamGuess(aTrack))
    plotFit(track, fitFixedPhase)
    dev.off()
}

aTrack <- getTrack(trackList, 328)
fitAll <- getTrackSinFitAll(object=aTrack, guess=getTrackParamGuess(aTrack))
fitFixedPhase <- getTrackSinFitFixedPhase(object=aTrack, phaseShift=fitAll$par['phaseShift'], guess=getTrackParamGuess(aTrack))
plotFit(track=aTrack, fit=fitAll)
plotFit(track=aTrack, fit=fitFixedPhase)


aTrack <- getTrack(trackList, 103)
dxdt <- getDerivative(aTrack@x, 2)
plot(aTrack@t, normalize(dxdt, abs=T), type='l', col='red')
lines(aTrack@t, normalize(aTrack@x), col='black')
#filteredTrack <- filterWandering(aTrack)


##### Test Stuff #####


# library(foreign)
# tracksFile <- read.arff(path)
# tracksFile2 <- reorganizeTable(tracksFile, nameCol='Metadata')
# row=6076
# duh <- new('track',tracksFile2[row,'Track'], tracksFile2[row,'polygonPts'], tracksFile2[row,'patternPts'], fi=0.2, ff=0.01, ti=0, tf=500, length.out=2361)


