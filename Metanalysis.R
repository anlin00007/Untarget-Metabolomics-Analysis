#-----------------------------------------------------------------------------------------------#
# Name - metanalysis.r													#
# Desc - Untargeted metabolomics profiling.									#
# Author - Difei Wang (dw670@georgetown.edu); Lin An (la512@georgetown.edu)				#
#-----------------------------------------------------------------------------------------------#

#-----------------------------------------------------------------------------------------------#
#						Name&Path Setting								#
#-----------------------------------------------------------------------------------------------#

## Setting Directories
myDir = "C:/Users/la512/Desktop/xcmstest/DMDneg"
setwd(myDir)

## CONTROL GROUP
myControl = "NORMAL"

## CASE GROUP
myCase = "DMD"

## Folder Name for Image Storage
boxname= "DMDvsNORMALnew"

## Technical Replicate Information
techrep="techrep.txt"


#-----------------------------------------------------------------------------------------------#
#						Parameter Setting								#
#-----------------------------------------------------------------------------------------------#

## Experiment Polarity
polar = "negative"

## Feature Detection Filter
ppmF <- 10
prefpeak <- 3
prefint <- 100

## Statistics Threshold
pthr <- 0.01

## Annotation Parameter
ppmA <- 5


#-----------------------------------------------------------------------------------------------#
#						Raw Data Processing							#
#-----------------------------------------------------------------------------------------------#

## Calling Libraries
require(xcms)
require(CAMERA)

## Load Raw Data
files <- list.files(myDir, pattern="*.CDF", recursive=TRUE, full.names=TRUE)
paste(files)
NumofSample <- length(files)
pd <- xcms:::phenoDataFromPaths(files)

## Feature detection (CentWave method)
xset <- xcmsSet(files,method="centWave", prefilter=c(prefpeak,prefint),
 ppm=ppmF, peakwidth=c(5,10), snthr=8, mzdiff=0.01, noise=0)
xset <- group(xset)

## Retention time correction (obiwarp method)
xsetR <- retcor(xset, method="obiwarp",profStep=1)

## Alignment (density method)
xsetA <- group(xsetR,bw=6, mzwid=0.025, minfrac=0.5, minsamp=1)

## FillPeaks
xsetF <- fillPeaks(xsetA)

## Diffreport
Diffre <- diffreport(xsetF, sortpval=FALSE, myControl, myCase, filebase=boxname, eicmax = 10, eicwidth = 200
, value = "into", classeic = c(myControl,myCase), test="t")
write.csv(Diffre,"Diffreport1.csv",row.names=FALSE)

## Annotated Diffreport
Annodiffre <- annotateDiffreport(xsetF, sample=c(1:NumofSample), polarity=polar, 
sigma=6, perfwhm=0.6, maxcharge=3, maxiso=4, ppm=ppmA, mzabs=0.015, pval_th = pthr, intval="into", quick=TRUE)
write.csv(Annodiffre,"Annotationdiff.csv",row.names=FALSE)

#-----------------------------------------------------------------------------------------------#
#						Technical Replicates Correction					#
#-----------------------------------------------------------------------------------------------#

## Remove Technical Replicates
avgduplicates <- function(reporttab, annot.file, NumofSample, data.start.column  = 14) {
	browser()
	data.end.column = data.start.column + NumofSample - 1
	ends = data.end.column + 1
	endc = data.end.column + 3
	dataMatrix=as.matrix(reporttab[, data.start.column:data.end.column])
	dataMatrix_names=as.matrix(reporttab[, 1])
	dataMatrix_mzrt=as.matrix(reporttab[, 5:13])
	dataMatrix_iap=as.matrix(reporttab[, ends:endc])
	storage.mode(dataMatrix) = 'double'
	col.names = colnames(dataMatrix)
	annot = as.matrix(read.delim(annot.file, header = TRUE))
	idx = match(col.names, annot[,1])
	col.names[!is.na(idx)] = annot[idx[!is.na(idx)], 2]
	colnames(dataMatrix)=col.names
	unique.col.names = unique(col.names)
	avg.dataMatrix = matrix(0, nrow(dataMatrix), length(unique.col.names))
	colnames(avg.dataMatrix) = unique.col.names
	for(i in 1:length(unique.col.names)) {
		cur.subMatrix = dataMatrix[, col.names == unique.col.names[i]]
		avg.dataMatrix[,i] =  rowMeans(cur.subMatrix, na.rm = TRUE)
	}
	avgcol.names = colnames(avg.dataMatrix)
	avgidx = match(avgcol.names, annot[,2])
	c1 = annot[avgidx, 2][annot[avgidx,3] == '0']
	c2 = annot[avgidx, 2][annot[avgidx,3] == '1']
      ## Check against missing Values
      if (any(is.na(avg.dataMatrix[,c(c1,c2)]))) {
          stop("NA values in xcmsSet. Use fillPeaks()")
      }
	mean1 <- rowMeans(avg.dataMatrix[,c1,drop=FALSE], na.rm = TRUE)
	mean2 <- rowMeans(avg.dataMatrix[,c2,drop=FALSE], na.rm = TRUE)
	## Calculate fold change.
      ## For foldchange <1 set fold to 1/fold
      ## See tstat to check which was higher
      fold <- mean2 / mean1
      fold[!is.na(fold) & fold < 1] <- 1/fold[!is.na(fold) & fold < 1]
      ## Calculate tstat
	testval <- avg.dataMatrix[,c(c1,c2)]
      testclab <- c(rep(0,length(c1)),rep(1,length(c2)))

      if (min(length(c1), length(c2)) >= 2) {
          tstat <- mt.teststat(testval, testclab)
          pvalue <- pval(testval, testclab, tstat)
      } else {
          message("Too few samples per class, skipping t-test.")
          tstat <- pvalue <- rep(NA,nrow(testval))
      }
      stat <- data.frame(fold = fold, tstat = tstat, pvalue = pvalue)
	save(avg.dataMatrix, file = 'avg.dataMatrix.Rdata')
	report <- do.call(cbind,list(dataMatrix_names, stat, dataMatrix_mzrt, avg.dataMatrix, dataMatrix_iap))
	write.csv(report, file = 'avg.report.csv', row.names=FALSE)
}

## P-Value Calculation
pval <- function(X, classlabel, teststat) {

    n1 <- rowSums(!is.na(X[,classlabel == 0]))
    n2 <- rowSums(!is.na(X[,classlabel == 1]))
    A <- apply(X[,classlabel == 0], 1, sd, na.rm=TRUE)^2/n1 ## sd(t(X[,classlabel == 0]), na.rm = TRUE)^2/n1
    B <- apply(X[,classlabel == 1], 1, sd, na.rm=TRUE)^2/n2 ## sd(t(X[,classlabel == 1]), na.rm = TRUE)^2/n2
    df <- (A+B)^2/(A^2/(n1-1)+B^2/(n2-1))

    pvalue <- 2 * (1 - pt(abs(teststat), df))
    invisible(pvalue)
}

Postrep <- avgduplicates(Annodiffre, techrep, NumofSample)
