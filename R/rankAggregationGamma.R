#!/usr/bin/Rscript

#===========================================================================
# Command-line execution statement:
# ---------------------------------
#
# Examples:
# ---------
# To run:
# Rscript rankAggregationGamma.R --inputfile="synthetic-samples.csv" --resultsfile="results/synth-results.csv" --outputfolder=results --debug --plot
#
#===========================================================================

rm(list=ls())

suppressMessages(library("optparse"))

srcpath <- "./R"
# Model-based Rank aggregation using mixtures of Gamma Distributions
source(file.path(srcpath,"fnRankAggregationGamma.R"))

option_list <- list(
	make_option("--inputfile", action="store", default="",
		help="Input file that contains 'label' in the first column, 'diff' in second column, and algo scores in following columns"),
	make_option("--resultsfile", action="store", default="",
		help="(optional) Output file which will store the inferred probability of being anomalous"),
	make_option("--outputfolder", action="store", default="",
		help="(optional) Output folder where the plots or parameter files are written to"),
	make_option("--saveparams", action="store_true", default=FALSE,
		help="(optional) specifies whether the parameters should be written to file."),
	make_option("--plot", action="store_true", default=FALSE,
		help="(optional) whether the distributions should be plotted"),
	make_option("--test", action="store_true", default=FALSE,
		help="(optional) whether to just use a test dataset for debug"),
	make_option("--debug", action="store_true", default=FALSE,
		help="(optional) whether to print debug statements")
)

if (F) {
  source("./R/rankAggregationGamma.R")
}

if (F) {
  args <- parse_args(OptionParser(option_list = option_list), args = commandArgs(trailingOnly = TRUE))
} else {
  # use these options for debugging only
  args <- parse_args(OptionParser(option_list = option_list), 
                     args = c(
                       sprintf("--inputfile=%s", file.path(srcpath,"synth1.csv")),
                       "--resultsfile=./results/synth-results.csv",
                       "--outputfolder=./results",
                       #"--test", 
                       "--saveparams", "--plot", "--debug")
  )
}
#print(args)

if (args$inputfile == "" && !args$test) {
  stop("No datafile specified")
}

# Data that either should be generated from test or should be read from file
anomind <- c()
if (args$test) {
  data <- read.csv(file.path(srcpath,"synthetic-samples.csv"),header=T)
  anomind <- rep(0,nrow(data))
  samples <- as.matrix(data)
} else {
  # read data from file
  data <- read.csv(args$inputfile, header=T)
  if (ncol(data) == 1) {
    anomind <- rep(0,nrow(data))
    samples <- as.matrix(data)
  } else {
    anomind <- ifelse(data[,1]=="nominal",0,1)
    data <- data[,2:ncol(data)]
    samples <- as.matrix(data)
  }
}

n <- nrow(samples)
D <- ncol(samples)

# make sure that all values in each detector are positive
min_samples <- apply(samples, 2, min)
for (d in 1:D) {
  # NOTE: Make sure that the samples are never ZERO, 
  # else, assign a very small value
  if (min_samples[d] < 1e-4) {
    samples[,d] = samples[,d]-min_samples[d]+1e-4
  }
}

# make sure results are reproducible
set.seed(32567)
params <- fnRankAggregationGamma(samples,F)
Pi <- params[["pi"]]
Z <- params[["z"]]
alpha_kd <- params[["alphaKd"]]
beta_kd <- params[["betaKd"]]

converged <- params[["converged"]]
epochs <- params[["epochs"]]
diff_ll <- params[["diffLoglik"]]
if (args$debug) cat(ifelse(converged,"converged in","did not converge in"),epochs,"; diff Log lik",diff_ll,"\n")

if (args$resultsfile != "") {
  if (args$test) {
    write.table(Z[,1], file=args$resultsfile, 
              col.names=c("gammamix"),row.names=F,sep=",",quote=F)
  } else {
    # Append the label,diff columns
    write.table(Z[,1], file=args$resultsfile, 
              col.names=c("gammamix"),row.names=F,sep=",",quote=F)
  }
}

#===========================================================================
# All the code below is optional and to be used only for debug.
# This is used only if one (or more) of the following switches is turned on:
# --saveparams --plot
#===========================================================================

if (args$outputfolder=="" && (args$plot || args$saveparams)) {
  stop("Output path needs to be specified")
}

if (args$saveparams || args$plot) {
  outpdfname <- paste("plot-fit.pdf",sep="")
  outpdfpath <- file.path(args$outputfolder, outpdfname)
}

# Parameters alpha_kd, beta_kd, Pi will be stored at location args$outputfolder
# Outputs:
#   Z.csv - Z variables
#   alpha_kd.csv - alpha parameters of Gamma distributions
#   beta_kd.csv - beta parameters of Gamma distributions
#   Pi.csv - Proportion of anomalies
if (args$saveparams) {
  write.table(Z, file=file.path(args$outputfolder,"Z.csv"), 
              col.names=c("P_Anomaly","P_Normal"),row.names=F,sep=",",quote=F)
  write.table(alpha_kd, file=file.path(args$outputfolder,"alpha_kd.csv"), 
              col.names=paste("Detector",as.character(1:ncol(alpha_kd)),sep=""),
              row.names=F,sep=",",quote=F)
  write.table(beta_kd, file=file.path(args$outputfolder,"beta_kd.csv"), 
              col.names=paste("Detector",as.character(1:ncol(beta_kd)),sep=""),
              row.names=F,sep=",",quote=F)
  write.table(matrix(c(Pi,1-Pi),ncol=2), file=file.path(args$outputfolder,"Pi.csv"), 
              col.names=c("P_Anomaly","P_Normal"),row.names=F,sep=",",quote=F)
}

#---------------
# Below commands are only for visualization for debug purposes.
# These should be disabled by default.
# These superimpose the inferred distributions over the score histogram 
if (args$plot) {
  pdffile <- file.path(args$outputfolder, outpdfname)
  algonames <- colnames(samples)
  anomidxs <- which(anomind==1)
  nanoms <- length(anomidxs)
  cat("No. anoms",nanoms,"\n")
  pdf(pdffile)
  par(mfrow=c(1,1))
  for (d in 1:D) {
    smpl <- samples[,d]
    #par(mfrow=c(1,1))
    mx <- max(smpl)
    #mx <- min(20,mx)
    x <- seq(from=1e-4,to=mx,length.out=100)
    y_anom <- dgamma(x,shape=alpha_kd[1,d], rate=beta_kd[1,d])
    y_noml <- dgamma(x,shape=alpha_kd[2,d], rate=beta_kd[2,d])
    hist(smpl, prob=T,
         breaks=seq(from=min(smpl),to=max(smpl),length.out=500),
         xlim=c(0,mx), ylim=c(0,max(c(y_anom, y_noml))), 
         main=paste(algonames[d],"(Pi=",as.character(round(Pi,3)),")",sep=""), 
         xlab="scores", ylab="density",cex.main=1.0)
    # plot anomaly distribution
    lines(x,y_anom*Pi,col="red",lwd=1)
    #lines(x,y,col="red",lwd=1)
    # plot regular distribution
    lines(x,y_noml*(1-Pi),col="green",lwd=1)
    #lines(x,y,col="green",lwd=1)
    if (nanoms > 0) {
      points(smpl[anomidxs],rep(0.25,nanoms),pch=".",col="red")
    }
  }
  dev.off()
}
