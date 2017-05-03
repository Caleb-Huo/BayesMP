rm(list=ls())

WD <- "/Users/caleb/Box Sync/chuck/bayesianmixturemodelinmetaanalysis/MCMC20141015/gitHub/BayesMP/debug"
WD1 <- "/Users/caleb/Box Sync/chuck/bayesianmixturemodelinmetaanalysis/MCMC20141015/gitHub/BayesMP/src"

compile <- "g++ -I/Library/Frameworks/R.framework/Resources/include  -DNDEBUG   -std=c++0x  -fpic  -g -O2 -fstack-protector --param=ssp-buffer-size=4 -Wformat -Wformat-security -Werror=format-security -D_FORTIFY_SOURCE=2 -g  -c BayesMP.cpp -o BayesMP.o"
build <- "R CMD SHLIB BayesMP.cpp"

setwd(WD1)
system(compile)
system(build)

dyn.load('BayesMP.so')

## generate data
setwd(WD)
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

niter <- 1000
burnin <- 50
nsample <- niter - burnin
if(F){
  mcmc(Z,niter = 1000, burnin=50)
}

beta=1/2
alpha=1
mu0=0
sigma0=10
sigma=1
trunc=0.1
Pi=NULL
delta=NULL
Y=NULL
niter=1000
burnin=50
fileName='BayesMP_mcmc'
fullRes=1
HSall=1

G <- nrow(Z) ## number of genes
S <- ncol(Z) ## number of studies

if(is.null(Pi)){
  Pi <- rbeta(G, gamma, 1-gamma)		
}
if(is.null(delta)){	
  delta <- rbeta(G, beta, beta)
}
if(is.null(Y)){
  Y0 <- matrix(0,G,S)		
  for(s in 1:S){
    Y_s <- numeric(G)
    Y_s[Z[,s]>2] = 1
    Y_s[Z[,s]< -2] = -1
    Y0[,s] <- Y_s
  }
  Y <- Y0		
}

obj <- .C('mcmc_R2',G=as.integer(G),S=as.integer(S),Z=as.double(Z),gamma=as.double(gamma),beta=as.double(beta),alpha=as.double(alpha),mu0=as.double(mu0),
          sigma0=as.double(sigma0),sigma=as.double(sigma),trunc=as.double(trunc),pi=as.double(Pi),delta=as.double(delta),Y=as.integer(Y0),
          niter=as.integer(niter),burnin=as.integer(burnin),fileName=as.character(fileName),fullRes=as.integer(fullRes),HSall=as.integer(HSall))

if(F){
  BayesMP::mcmc(Z, estGamma, niter=niter, burnin=burnin)
}

