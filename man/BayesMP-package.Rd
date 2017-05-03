\name{BayesMP-package}
\alias{BayesMP-package}
\alias{BayesMP}
\docType{package}
\title{
BayesMP
}
\description{
Implementation for Zhiguang Huo, Chi Song and George Tseng. (2016) Bayesian latent hierarchical model for transcriptomic meta-analysis to detect biomarkers with
clustered meta-patterns of differential expression signals. }
\details{
\tabular{ll}{
Package: \tab BayesMP\cr
Type: \tab Package\cr
Version: \tab 0.0.2\cr
Date: \tab 2015-08-16\cr
License: \tab GPL-3\cr
}
~~ An overview of how to use the package, including the most important functions ~~
}
\author{
Zhiguang Huo, Chi Song and George Tseng
Maintainer: Zhiguang Huo <zhh18@pitt.edu>
}
\references{
Zhiguang Huo, Chi Song and George Tseng. (2016) Bayesian latent hierarchical model for transcriptomic meta-analysis to detect biomarkers with
clustered meta-patterns of differential expression signals.
}
\keyword{ package }
\examples{

######################################
## generate data
n1 <- 500
n2 <- 50
S <- 3
truth <- c(rep('nonDE',n1),rep('DE',n2))

set.seed(15213)
aZ <- c(rnorm(n1),rnorm(n2,3)*sample(c(1,-1),n2,replace=TRUE))
set.seed(15214)
bZ <- c(rnorm(n1),rnorm(n2,3)*sample(c(1,-1),n2,replace=TRUE))
set.seed(15215)
cZ <- c(rnorm(n1),rnorm(n2,3)*sample(c(1,-1),n2,replace=TRUE))
Z <- cbind(aZ, bZ, cZ)
G <- nrow(Z)

w <- locfdr(Z, plot=0)
estGamma <- 1 - w$fp0["mlest", "p0"]

niter <- 100
burnin <- 50
nsample <- niter - burnin

mcmc(Z, estGamma, niter=niter, burnin=burnin, trunc=0.1)

HSallRes <- read.table('BayesMP_mcmcHSall.txt')


## Bayesian inference.
## pos=1: HSb. pos=S: HSa. pos=r: HSr.
pos <- 3
HSb_belief <- HSallRes[,3]/nsample
HSb_qvalue <- BayesianFDR(HSb_belief)
sum(HSb_qvalue<0.05)
table(truth, HSb_qvalue<0.05)


## MetaPattern
fileNameFull <- 'BayesMP_mcmcFull.txt'
con  <- file(fileNameFull, open = "r")

pi_g.index <- 1:G
delta_g.index <- 1:G + G
Yindex <- 2*G + 1:(S*G)

resYplus <- matrix(0,G,S)
resYminus <- matrix(0,G,S)


i = 1
while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0) {
  print(i)
  if(i>burnin){
	  seven = strsplit(oneLine, "\t")[[1]]
	  thisY <- matrix(as.numeric(seven[Yindex]),G,S)
  	
	  ## for individual studies
	  resYplus[thisY>0] <- resYplus[thisY>0] + 1
	  resYminus[thisY<0] <- resYminus[thisY<0] + 1
	  
	  ## for meta-analysis
  
	  DEevidence = apply(thisY,1,function(x) sum(x!=0))
	  for(j in 1:S){
		  curLogic = DEevidence == j
	  }
	  
  }    
  i = i + 1
} 

close(con)


## tight clustering
dissimilarity <- distance(resYplus, resYminus, nsample)
tightClustResult <- tightClustPam(dissimilarity, target=2, k.min=10)


}