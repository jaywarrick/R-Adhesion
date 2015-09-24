#path <- '/Users/jaywarrick/Documents/MMB/Projects/Adhesion/MM AMD/'
#path <- '/Users/jaywarrick/Documents/MMB/Projects/Adhesion/MM FACS/20150806/'
path <- '/Users/jaywarrick/Documents/MMB/Projects/Adhesion/MM FACS/20150812/'

#wellNames <- c('x0_y0', 'x0_y1', 'x0_y2', 'x0_y3', 'x0_y4', 'x0_y5', 'x0_y6', 'x0_y7')
wellNames <- c('x0_y0', 'x0_y1', 'x0_y2')

wellNames2 <- paste0(wellNames, '_c')

results <- list()
for(wellName in wellNames[3])
{
     rm(list=setdiff(ls(), c("wellName","wellNames","path")))
     well <- wellName
     source('~/Public/DropBox/GitHub/R-Adhesion/CellTracking/RunTestAnalysis2.R', echo=TRUE)
     #      results[[wellName]] <- get(wellName)$getPercentAdhered()
     #      pdf(file=paste0(path,wellName,'_PercentAdhered.pdf'), width=6, height=3.5)
     #      plot(results[[wellName]]$time, results[[wellName]]$percentAdhered, xlab='Time [s]', ylab='Percent Adhered [%]', main=wellName, pch=20, cex=0.75, ylim=c(0,70))
     #      dev.off()
}