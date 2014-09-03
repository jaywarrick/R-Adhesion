results <- comparison()
library(foreign)
write.arff(results, file='~/Documents/MMB/Projects/Adhesion/R/Testing/TrackingMethodsResults.txt')
results <- read.arff(file='~/Documents/MMB/Projects/Adhesion/R/Testing/TrackingMethodsResults.txt')

plotMethod(results, 'error1Mean', normalize=TRUE)
plotMethod(results, 'error2Mean', normalize=TRUE)
plotMethod(results, 'error3Mean', normalize=TRUE)
plotMethod(results, 'error4Mean', normalize=TRUE)
plotMethod(results, 'error5Mean', normalize=TRUE)
plotMethod(results, 'error6Mean', normalize=TRUE)


directionality <- 10
points <- getPoints(n=300, movement=c(0.1,0,0), direction=c(1,0,0), parallelNoise=0.1/4, perpendicularNoise=(0.1/4)/directionality)
px <- DirectionalLinearAssignment(points, maxDist <- 0.5, direction=c(1,0,0), directionality=directionality, uniformityDistThres=0.03, digits=4)

duh <- testLAPAccuracy_outOfFrame2(n=300)

duh <- testLAPAccuracy_outOfFrame(n=100)

n<-10
mags<-seq(0,0.1,length.out=3)
m <- mags[3]
getCostMatrix(points=newPoints, pointByPoint=FALSE, linkCostFunction=getDirectionalCost2, digits=4, maxDist=m+2*0.1, direction=c(1,0,0), directionality=10)

duh <- testLAPAccuracy_directionality()

duh <- createMovements(n=100, movement=c(0,0,0), direction=c(1,1,0), parallelNoise=0.1/4, perpendicularNoise=0.01)
plot(duh$x,duh$y)

directionality <- 10
points <- getPoints(n=300, movement=c(0,0,0), direction=c(1,0,0), parallelNoise=0.1/4, perpendicularNoise=(0.1/4)/directionality)
duh <- plotAssignments(points, direction=c(1,0,0), directionality=10)
duh2 <- plotAssignments(points, direction=c(1,0,0), directionality=1)
r <- nrow(points$t0)
duh$errorRate
duh2$errorRate
duh$px[1:nrow(t0)] - 1:nrow(t0)
which(diag(duh$costMat[,duh$px])==getBlockingCost())
duh$costMat

points <- getPoints(n=5)
getDirectionalCost2(points)
getCostMatrix(t0=points$t0, t1=points$t1, linkCostFunction=getDirectionalCost, digits=4, maxDist=15, direction=c(1,0,0), directionality=10)

# # Ran into issue with integer solver in going over the .Machine$integer.max within the C code that cause erroneous results
# # This code allowed me to trouble shoot how to avoid this.
# # Essentially, the "blocking distance" is set to 10 times the max distance to consider
# duh <- matrix(getBlockingCost(),5,5)
# for(i in 1:5)
# {
#      duh[i,6-i] <- 6-i
# }
# for(i in 1:5)
# {
#      duh[i,i] <- i
# }
# duh
# LinearAssignment(duh)

## TEST B 

# movement=c(0.1,0,0), direction=c(1,0,0), parallelNoise=0.1/4, perpendicularNoise=(0.1/4)/directionality, directionality=10, maxDist=0.5
## TEST C
# movement=c(0,0,0), direction=c(1,0,0), parallelNoise=0.1/4, perpendicularNoise=(0.1/4)/directionality, directionality=10, maxDist=0.4
