# Adhesion Analysis with Liz
rm(list=ls())
source('~/.Rprofile')
source('/Users/jaywarrick/Documents/MMB_Lab_2/Projects/Adhesion/TestAnalysis/AdhesionAnalysis_HelperFunctions.R')

##### Analysis #####

# Remember cell 328
path <- '/Volumes/BeebeBig/Jay/JEX Databases/Adhesion FACS/Test/Cell_x0_y0/Roi-Tracks/x0_y0.jxd'
path2 <- '/Volumes/BeebeBig/Jay/JEX Databases/Adhesion FACS/RPMI P-Sel 5Hz-100mHz/Cell_x0_y0/Roi-Tracks Roi/x0_y0.jxd'

trackListImport <- new('trackList', tracksFilePath=path, fi=5, ff=5, ti=0, tf=500, samplingFreq=0.5)
trackList <- filterTrackList(trackListImport, trackLengthFilter, min=4)
plot(trackList, slot='vx')
trackList <- setTrackListSinPhaseShift(object=trackList)
aTrack <- trackList@tracks[[5001]]
plot(aTrack, slotY='vx')

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



duh <- list(a='hi', b='there')
duh2 <- 'c'
duh <- append(duh, list(duh2='you'))
duh

length(trackList)




