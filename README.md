# BayesMP
R package for Bayes Meta Pattern


## Required Package
* qvalue
* cluster
* truncnorm
* invgamma

## Install This Package from github

* From GitHub
```R
library(devtools)
install_github("Caleb-Huo/BayesMP")
```

* From Bioconductor
```R
"to be filled in"
```

## Full tutorial

https://bayesmp.github.io


## Short tutorial
```R

## Simulate data
set.seed(15213)
G <- 2000
S <- 4
alpha <- 200
X0 <- matrix(rnorm(G * S), G, S)
Xplus <- matrix(rnorm(G * S, 2), G, S)
Xminus <- matrix(rnorm(G * S, -2), G, S)
piall <- rbeta(G, alpha/G, 1)
delta <- rbeta(G, 1/2, 1/2)
p0 <- 1 - piall
p1 <- piall * delta
p2 <- piall * (1 - delta)
Y <- replicate(S, apply(cbind(p0, p1, p2),1,function(x) sample(c(0,1,-1),1,prob = x)))
Z <- X0 * (Y == 0) + Xplus * (Y == 1) + Xminus * (Y == -1)

## Perform MCMC with Emperical Bayes (EB) method.
niter=200
burnin=50
system.time(BayesMP_EB(Z,niter=niter, burnin=burnin, writeY=T, writeHSall=T))

HSallRes <- read.table('BayesMP_EB_HSall.txt')


## Bayesian inference.
## pos=1: HSb. pos=S: HSa. pos=r (1<r<S): HSr.
HSb_belief <- HSallRes[,1]/(niter - burnin)
HSb_qvalue <- BayesianFDR(HSb_belief)
sum(HSb_qvalue<0.05)


## MetaPattern
fileNameFull <- 'BayesMP_EB_Y.txt'
con  <- file(fileNameFull, open = "r")

resYplus <- matrix(0,G,S)
resYminus <- matrix(0,G,S)


i = 1
while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0) {
  if(i>burnin){
	  print(i)
	  seven = strsplit(oneLine, "\t")[[1]]
	  thisY <- matrix(as.numeric(seven),G,S)
  	
	  ## for individual studies
	  resYplus[thisY>0] <- resYplus[thisY>0] + 1
	  resYminus[thisY<0] <- resYminus[thisY<0] + 1
  }    
  i = i + 1
} 

close(con)

resYplus_DE <- resYplus[HSb_qvalue<0.05,]
resYminus_DE <- resYminus[HSb_qvalue<0.05,]

## tight clustering
dissimilarity <- distance(resYplus_DE, resYminus_DE, niter - burnin)
tightClustResult <- tightClustPam(dissimilarity, target=2, k.min=10)
```

## Citation
* Zhiguang Huo, Chi Song, George C. Tseng. (2018) Bayesian latent hierarchical model for transcriptomic meta-analysis to detect biomarkers with clustered meta-patterns of differential expression signals. Annals of Applied Statistics (Accepted).
