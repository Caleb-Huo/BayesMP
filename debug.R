setwd("~/Desktop/BayesMP")
devtools::document()
devtools::install()




setwd("Desktop/")
library(BayesMP)

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
system.time(obj_EB <- BayesMP_EB(Z, writeHSall = F))
system.time(BayesMP_EB(Z))


system.time(BayesMP_EB(Z,writeY = 1, writePi = 1, writeGamma = 1, writeMu = 1, writeS2 = 1, writeHSall = 1, silence = T))
system.time(BayesMP_EB(Z,writeY = 1, writePi = 1, writeGamma = 1, writeMu = 1, writeS2 = 1, writeHSall = 1))
system.time(BayesMP_EB(Z))




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
fileName='BayesMP_EB_'
writeY = 1
writePi = 1
writeGamma = 1
writeMu = 1
writeS2 = 1
writeHSall = 1

