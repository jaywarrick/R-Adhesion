library(foreign)
library(plyr)

#' Take an arff file and reorganize it into a more standard 'table' format. Specifically this is used to
#' import an arff file from JEX as JEX uses a column called 'Measurement' to define the type of measurment
#' or property being stored and 'Value', the value of that property.
#'
#' @param data An object that is the result of using foreign::read.arff(file) on an arff file
#' @param baseName An optional basename to add to whatever label is in the \code{nameCol} portion of each row entry
#' @param convertToNumeric An option to convert the columns of information within \code{data} leading up to
#' \code{nameCol} and \code{valueCol} to numeric or to leave as text. Default is to convert to numeric (i.e., TRUE)
#' @param nameCol The name of the column that describes the nature of the value in the \code{valueCol}
#' @param valueCol The name of the column with the values of the properties listed in the \code{nameCol}
reorganizeTable <- function(data, baseName=NA, convertToNumeric=TRUE, nameCol='Measurement', valueCol='Value')
{
     require(plyr)
     idCols <- names(data)
     idCols <- idCols[-which(idCols %in% c(nameCol,valueCol))]
     newData <- data.frame(stringsAsFactors=FALSE)
     measurements <- unique(data[,nameCol])

     for(m in measurements)
     {
          if(is.na(baseName))
          {
               newColName <- m
               newColName <- gsub(' ','.',newColName, fixed=TRUE) # Get rid of extraneous spaces
          }else
          {
               newColName <- paste(baseName,'.',m, sep='')
               newColName <- gsub(' ','.',newColName, fixed=TRUE) # Get rid of extraneous spaces
          }

          temp <- data[data[,nameCol]==m,]
          temp2 <- temp[,c(idCols,valueCol)]
          if (length(idCols) == 0) {
               temp2 <- data.frame(ReallyRandomNameYo = temp2)
               names(temp2) <- newColName
          }
          else {
               names(temp2)[names(temp2) == valueCol] <- newColName
          }
          if(nrow(newData) == 0)
          {
               newData <- temp2
          }else
          {
               newData <- merge(newData, temp2, by=idCols)
          }
     }

     if(convertToNumeric)
     {
          for(n in idCols)
          {
               newData[,n] <- as.numeric(as.character(newData[,n]))
          }
     }

     return(newData)
}

#' Return information and data for a log frequency sweep with the given parameters
#' @param amplitude numeric
#' @param phaseShift numeric - default=0
#' @param offset numeric - defulat=0
#' @param sin boolean T=sinsoid and F=triangular - default=FALSE
#' @param fi numeric initial frequency of the sweep
#' @param ff numeric final frequency of the sweep
#' @param sweepDurtion numeric
#' @param tAll numeric vector All possible times for this
getSweep <- function(amplitude=1, phaseShift=0, offset=0, sin=FALSE, fi=2, ff=0.1, sweepDuration, t, guess=NULL)
{
     ti <- 0
     tf <- sweepDuration
     N <- log(ff/fi)/log(2)
     R <- N / (tf-ti)

     if(is.null(guess))
     {
          A <- amplitude
          phi <- phaseShift
          b <- offset
     }
     else
     {
          A <- guess['amplitude']
          phi <- guess['phaseShift']
          b <- guess['offset']
     }

     offsetInflections <- (  (-1*(phi/(pi/2))) %/% 1  )
     nf <- -((4*fi-4*ff)*pi*tf+2*log(ff/fi)*phi)/(log(ff/fi)*pi)
     nf <- nf + 1 # for good measure.
     ni <- seq(0,nf,1) + offsetInflections
     suppressWarnings(inflections <- (log((log(ff/fi)*phi)/(2*fi*pi*tf)+(log(ff/fi)*ni)/(4*fi*tf)+1)*tf)/log(ff/fi))
     inflections <- inflections[!is.nan(inflections)]
     inflectionNums <- (  seq(0,length(inflections)-1,1) + (  offsetInflections %% 4  )  ) %% 4
     startTimeI <- which(inflections >= t[1])[1]
     endTimeI <- which(inflections <= t[length(t)])
     endTimeI <- endTimeI[length(endTimeI)]
     if(is.na(endTimeI))
     {
          stop(paste(t[length(t)], "is greater than last time of last inflection point. Therefore, there must be an error in fi/ff or tAll as there are inflections guaranteed in during times ti to tf given fi and ff."))
     }
     inflections <- inflections[startTimeI:endTimeI]
     inflectionNums <- inflectionNums[startTimeI:endTimeI]

     if(sin)
     {
          predicted <- A*sin(2*pi*((fi*(-1+2^(R*t)))/(R*log(2))) - phi) + b
          v=getDerivative(x=predicted, t=t)
          return(list(A=A, t=t, x=predicted, line=line, v=v, inflections=inflections, inflectionNums=inflectionNums))
     }
     else
     {
          predicted <- (2*A/pi)*asin(sin(2*pi*((fi*(-1+2^(R*t)))/(R*log(2))) - phi)) + b
          v=getDerivative(x=predicted, t=t)
          return(list(A=A, t=t, x=predicted, line=line, v=v, inflections=inflections, inflectionNums=inflectionNums))
     }
}
