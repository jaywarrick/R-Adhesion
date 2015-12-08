#rm(list=ls())
path <- '/Users/jaywarrick/Documents/MMB/Projects/Adhesion/MM AMD/20150724/'
#path <- '/Users/jaywarrick/Documents/MMB/Projects/Adhesion/MM FACS/20150806/'
#path <- '/Users/jaywarrick/Documents/MMB/Projects/Adhesion/MM FACS/20150812/'

#wellNames <- c('x0_y0', 'x0_y1', 'x0_y2', 'x0_y3', 'x0_y4', 'x0_y5', 'x0_y6', 'x0_y7')
#wellNames <- c('x0_y0', 'x0_y1', 'x0_y2')
wellNames <- c('x0_y8')

wellNames2 <- paste0(wellNames, '_c')

dataNames <- wellNames
col <- c(rgb(1,0,0), rgb(1,0,0), rgb(1,0,0), gray(0), gray(0), gray(0), rgb(1,0,0), rgb(1,0,0))
lty <- c(1,1,2,1,2,1,1,1)
plotData <- function(dataNames, xlim=c(0,300), ylim=c(0,70), xlab='Shear Stress [Pa]', ylab='Percent Adhered [%]', col=gray(seq(0.2,0.8,0.6/length(dataNames))), lty=rep(1, length(dataNames)), PEG=TRUE, ...)
{
     if(PEG)
     {
          mu <- 0.0052
     }
     else
     {
          mu <- 0.00078
     }
     tANDf <- getFrequencies(t=xlim)
     tau <- getShearStress(f=tANDf$f, pixelAmplitude=500, mu=mu)
     plot(c(), c(), xlim=tau, ylim=ylim, xlab=xlab, ylab=ylab, log='x', ...)
     for(dataName in dataNames)
     {
          amplitude <- max(as.numeric(get(dataName)$getProp(fun=function(x){r <- x$range('x', rel=TRUE); r <- (r[2]-r[1])/2; return(r)})))
          results <- get(dataName)$getPercentAdhered()
          tANDf <- getFrequencies(t=results$time)
          tau <- getShearStress(f=tANDf$f, pixelAmplitude=amplitude, mu=mu)
          temp <- smooth.spline(tau, results$percentAdhered, spar=0.5)

          points(tau, results$percentAdhered, pch=20, cex=0.75)
          lines(temp$x, temp$y, type='l', col=col[which(dataNames %in% dataName)], lty=lty[which(dataNames %in% dataName)], lwd=3)
     }
}

plotData(dataNames, col=col, lty=lty, PEG=FALSE)
x0_y8 <- trackList
results <- list()
for(wellName in wellNames)
{
     wellName2 <- wellName#paste0(wellName, '_c')
     load(paste0(path, wellName2, '.RData'))
     assign(wellName2, get(wellName)$copy())
     results[[wellName2]] <- get(wellName2)$getPercentAdhered()
     pdf(file=paste0(path,wellName2,'_PercentAdhered.pdf'), width=6, height=5)
     plot(results[[wellName2]]$time, results[[wellName2]]$percentAdhered, xlab='Time [s]', ylab='Percent Adhered [%]', main=wellName2, pch=20, cex=0.75, ylim=c(0,70))
     dev.off()

     load(paste0(path, wellName, '.RData'))
     assign(wellName, get(wellName)$copy())
     results[[wellName]] <- get(wellName)$getPercentAdhered()
     pdf(file=paste0(path,wellName,'_PercentAdhered.pdf'), width=6, height=5)
     plot(results[[wellName]]$time, results[[wellName]]$percentAdhered, xlab='Time [s]', ylab='Percent Adhered [%]', main=wellName, pch=20, cex=0.75, ylim=c(0,70))
     dev.off()
}

# path <- '/Users/jaywarrick/Documents/MMB/Projects/Adhesion/MM FACS/'
# load(paste0(path, 'x0_y0', '.RData'))
# trackList <- x0_y0
# bestFit <- getBulkPhaseShift2(trackList, tiGuess=0)
