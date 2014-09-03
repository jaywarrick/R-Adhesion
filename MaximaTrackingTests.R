source('~/.Rprofile')
library(foreign)
path1 <- '/Users/jaywarrick/Documents/JEX/LocalTest/temp/JEXData0000000001.arff'
path2 <- '/Users/jaywarrick/Documents/JEX/LocalTest/PC3 vs LNCaP/Cell_x0_y0/Roi-Maxima Upper/x0_y0.jxd'
fileTable <- read.arff(path1)
# read.arff(fileTable$Value[1])
#maxima <- read.arff('/Users/jaywarrick/Documents/JEX/LocalTest/PC3 vs LNCaP/Cell_x0_y0/Roi-Maxima Upper/x0_y0.jxd')

maximaList <- new('maximaList', maximaFilePath=path2, timeDimName='Time')
maxima <- getMaxima(maximaList, frame=0)
plot(maxima)
cbind(maxima@points, data.frame(z=0))

startMaxima <- getMaxima(maximaList, 501)
newList <- new('maximaList')
newList <- setMaxima(newList, startMaxima)
newList@maxID <- max(startMaxima@points$id)

done <- FALSE
count = 1;
# for(i in 1:18)
# {
#      results <- getNextLinkedMaxima(newList, maximaList)
while(!is.null(results <- getNextLinkedMaxima(newList, maximaList)))
{
     print(' ')
     print(paste('Frame = ', results$newMaxima@frame, sep=''))
     print(' ')
     newList@maxima[[as.character(results$newMaxima@frame)]] <- results$newMaxima
     newList@maxID <- results$maxID
     pdf(file=paste('/Users/jaywarrick/Documents/MMB/Projects/Adhesion/R/Analysis/TrackPlots/Frame', results$newMaxima@frame, '.pdf', sep=''), width=9, height=4)
     plot(results$newMaxima, IDs=TRUE, xlim=c(0,1300), ylim=c(400,100), main=as.character(results$newMaxima@frame))
     dev.off()
     count <- count + 1
}

tAll=seq(0,515,1)
trackList <- new('trackList', sin=FALSE, fi=2, ff=0.1, tAll=tAll)
# trackList@tracks <- vector('list', 10001)
# names(trackList@tracks) <- as.character(0:10000)
for(maxima in maximaList@maxima[1:3])
{
     for(i in 1:nrow(maxima@points))
     {
          data <- maxima@points[i,]
          trackList <- addTrackPoint(trackList, id=data$id, x=data$x, y=data$y, t=tAll[maxima@frame+1], frame=maxima@frame)
     }
}
lapply()

print('YAY!')
