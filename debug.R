setwd("~/Desktop/BayesMP")
devtools::document()
devtools::install()

devtools::use_vignette("BayesMP")


setwd("Desktop/")
library(BayesMP)

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

niter=200
burnin=50
system.time(BayesMP(Z,niter=niter, burnin=burnin, writeY=T, writeHSall=T))

HSallRes <- read.table('BayesMP_HSall.txt')


## Bayesian inference.
## pos=1: HSb. pos=S: HSa. pos=r (1<r<S): HSr.
HSb_belief <- HSallRes[,1]/(niter - burnin)
HSb_qvalue <- BayesianFDR(HSb_belief)
sum(HSb_qvalue<0.05)


## MetaPattern
fileNameFull <- 'BayesMP_Y.txt'
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




set.seed(15213)
G <- 5000
K <- 10
alpha <- 200
X0 <- matrix(rnorm(G * K), G, K)
Xplus <- matrix(rnorm(G * K, 2), G, K)
Xminus <- matrix(rnorm(G * K, -2), G, K)
piall <- rbeta(G, alpha/G, 1)
delta <- rbeta(G, 1/2, 1/2)
p0 <- 1 - piall
p1 <- piall * delta
p2 <- piall * (1 - delta)
Y <- replicate(K, apply(cbind(p0, p1, p2),1,function(x) sample(c(0,1,-1),1,prob = x)))
Z <- X0 * (Y == 0) + Xplus * (Y == 1) + Xminus * (Y == -1)
niter <- 2000

system.time(BayesMP_DP(Z, writeY = T, writePi = T, writeDelta = T, writeGamma = T, writeHSall = T))
system.time(obj_DP <- BayesMP_DP(Z, writeHSall = F))
system.time(obj_EB <- BayesMP(Z, writeHSall = F))
system.time(BayesMP(Z))


system.time(BayesMP(Z,writeY = 1, writePi = 1, writeGamma = 1, writeMu = 1, writeS2 = 1, writeHSall = 1, silence = T))
system.time(BayesMP(Z,writeY = 1, writePi = 1, writeGamma = 1, writeMu = 1, writeS2 = 1, writeHSall = 1))
system.time(BayesMP(Z))




## debug EB

gamma=NULL
updateGamma = TRUE
beta=1/2
empMu = rep(0,ncol(Z))
empSD = rep(1,ncol(Z))
niter=100
burnin=50
silence = F
logDotsPerLine = 50
fileName='BayesMP_'
writeY = 1
writePi = 1
writeGamma = 1
writeMu = 1
writeS2 = 1
writeHSall = 1

