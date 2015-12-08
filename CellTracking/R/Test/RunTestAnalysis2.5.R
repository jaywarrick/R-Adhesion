

x <- 0
y <- 0:8
dbDir <- '/Users/jaywarrick/Documents/MMB/Projects/Adhesion/MM AMD/'

results <- list()
for(xi in x)
{
     for(yi in y)
     {
          entry <- paste0('Cell_x',xi,'y',yi)
          path <- dBDir #file.path(dbDir, entry)
          rm(list=setdiff(ls(), c("path")))
          well <- paste0('x', xi, '_y', yi)
          source('~/Public/DropBox/GitHub/R-Adhesion/CellTracking/RunTestAnalysis2.R', echo=TRUE)
          #      results[[wellName]] <- get(wellName)$getPercentAdhered()
          #      pdf(file=paste0(path,wellName,'_PercentAdhered.pdf'), width=6, height=3.5)
          #      plot(results[[wellName]]$time, results[[wellName]]$percentAdhered, xlab='Time [s]', ylab='Percent Adhered [%]', main=wellName, pch=20, cex=0.75, ylim=c(0,70))
          #      dev.off()
     }
}