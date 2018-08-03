##' BayesMP_EB
##'
##' implementation for BayesMP, MCMC part. This is a full Bayesian model, and the alternative distribution is modeled via dirichlet process.
##' @title MCMC for BayesMP_EB
##' @param Z Z statistics. Z should be a p*n matrix.
##' @param gamma initial gamma Estimated null proportation, by default, gamma will be estimated by the emperical null method.
##' @param updateGamma If TRUE, will update gamma by MH method. If FALSE, will keep the gamma as constant.
##' @param beta Non-informative prior: given a gene is DE, the prior probablity this gene is up-regulated, default=1/2.
##' @param empMu a vector of mean parameter for the null component, default is 0. Alternatively, this vector can by estimated by the emperical null method.
##' @param empSD a vector of sd parameter for the null component, default is 1. Alternatively, this vector can by estimated by the emperical null method.
##' @param niter Number of iterations. Default 100, suggest to be 10,000
##' @param burnin Number of burnin period. Default 50, suggest to be 500
##' @param silence If FALSE (default), will print the MCMC progress in the console.
##' @param logDotsPerLine Number of dots printed perline in the console when silence is FALSE.
##' @param fileName Base fileName for saving fulll mcmc results.
##' @param writeY If TRUE, will save all (niter) posterior samples of Y from MCMC. 
##' @param writePi If TRUE, will save all (niter) posterior samples of Pi from MCMC. 
##' @param writeDelta If TRUE, will save all (niter) posterior samples of Delta from MCMC. 
##' @param writeGamma If TRUE, will save all (niter) posterior samples of Gamma from MCMC. 
##' @param writeHSall If TRUE, will save the HSall matrix (niter - burnin). Each row of HSall represent a input gene (feature). If the number in the ith column and jth row equals m, 
##' it represents there are m posterior samples of Y for feature j that are DE in at least i studies. This matrix will be the input matrix to calculate the Bayesian FDR.
##' For example, for each column i, after normalized (divided) by total number effective samples (niter - burnin), it is the Bayesian belief that the genes are significiant in at least i studies.
##' The Bayesian FDR can be calculated by the BayesianFDR function.
##' @return The MCMC object for the last iteration.
##' @useDynLib BayesMP
##' @author Zhiguang Huo <zhuo@ufl.edu>
##' @export
##' @examples
##' set.seed(15213)
##' G <- 2000
##' K <- 10
##' alpha <- 200
##' X0 <- matrix(rnorm(G * K), G, K)
##' Xplus <- matrix(rnorm(G * K, 2), G, K)
##' Xminus <- matrix(rnorm(G * K, -2), G, K)
##' piall <- rbeta(G, alpha/G, 1)
##' delta <- rbeta(G, 1/2, 1/2)
##' p0 <- 1 - piall
##' p1 <- piall * delta
##' p2 <- piall * (1 - delta)
##' Y <- replicate(K, apply(cbind(p0, p1, p2),1,function(x) sample(c(0,1,-1),1,prob = x)))
##' Z <- X0 * (Y == 0) + Xplus * (Y == 1) + Xminus * (Y == -1)
##' system.time(BayesMP_EB(Z,writeHSall=F))


BayesMP_EB <- function(Z, gamma=NULL, updateGamma = TRUE, beta=1/2, empMu = rep(0,ncol(Z)), empSD = rep(1,ncol(Z)), 
				niter=100, burnin=50, silence = F, logDotsPerLine = 50, fileName='BayesMP_EB_', 
				writeY = F, writePi = F, writeGamma = F, writeMu = F, writeS2 = F, writeHSall = T) {
  
  ## initialization	
  G <- nrow(Z) ## number of genes
  S <- ncol(Z) ## number of studies

  if(is.null(gamma)){
      pp <- pmin(pnorm(-abs(Z)) * 2, 1)
	  gamma <- max(1 - pi0est(pp)$pi0, 0.01)		
  }	
  
  log = FALSE
  f0 <- function(x) dnorm(x, 0, 1, log = log)
  f1 <- function(x) dnorm(x, 1, 1, log = log)
  f2 <- function(x) dnorm(x, -1, 1, log = log)
  attr(f0, "mu") <- 0
  attr(f0, "var") <- 1
  attr(f1, "mu") <- 1
  attr(f1, "var") <- 1
  attr(f2, "mu") <- -1
  attr(f2, "var") <- 1
  f <- list(replicate(S, f0), replicate(S, f1), replicate(S, f2))
  for(s in 1:S){
	  attr(f[[1]][[s]], "mu") <- empMu[s]
	  attr(f[[1]][[s]], "var") <- empSD[s]
  }
      
  obj <- list(Y = NULL, pi = rep(1/2, G), delta = rep(1/2, G),f = f, gamma = gamma, beta = beta, log = FALSE, update.f = TRUE, update.gamma=updateGamma )
  if(writeHSall){
  	 obj$YHSall <- matrix(0,nrow=nrow(Z),ncol=ncol(Z))
  }
  
  if(updateGamma){
  	obj$MHsd <- 0.1
	obj$MHcount <- 0;
  }
  
  
  
  if(writeY){
	  outfile_Y <- file(paste0(fileName, "Y.txt"))
	  open(outfile_Y, "w")
  }
  if(writePi){
	  outfile_Pi <- file(paste0(fileName, "Pi.txt"))
	  open(outfile_Pi, "w")
  }
  if(writeGamma){
	  outfile_Gamma <- file(paste0(fileName, "Gamma.txt"))
	  open(outfile_Gamma, "w")
  }
  if(writeMu){
	  outfile_Mu <- file(paste0(fileName, "Mu.txt"))
	  open(outfile_Mu, "w")
  }
  if(writeS2){
	  outfile_S2 <- file(paste0(fileName, "S2.txt"))
	  open(outfile_S2, "w")
  }

  
  if (!silence) 
          cat("performing MCMC, i = 1,2,..., niter (= ", niter, ")  [one \".\" per sample]:\n", sep = "")
    
  for(i in 1:niter) {
    if (!silence) 
		cat(".", if (i%%logDotsPerLine == 0) 
        	paste(i, "\n"))
    obj <- iterate.one(obj, Z)
    M <- length(obj$f)
    K <- length(obj$f[[1]])
    mu <- sapply(obj$f, function(x) sapply(x, function(y) attr(y, "mu")))
    s2 <- sapply(obj$f, function(x) sapply(x, function(y) attr(y, "var")))
    if(obj$update.f & writeHSall & i>burnin) obj$YHSall <- update.HSall(obj)
	if(writeY){
		thisY <- obj$Y - 1
		thisY[thisY==2] <- -1
    	writeLines(paste(as.vector(thisY), collapse = "\t"), con = outfile_Y)			
	}
	if(writePi){
    	writeLines(paste(as.vector(obj$pi), collapse = "\t"), con = outfile_Pi)			
	}
	if(writeGamma){
    	writeLines(paste(as.vector(obj$gamma), collapse = "\t"), con = outfile_Gamma)			
	}
	if(writeMu){
    	writeLines(paste(as.vector(mu), collapse = "\t"), con = outfile_Mu)			
	}
	if(writeS2){
    	writeLines(paste(as.vector(s2), collapse = "\t"), con = outfile_S2)						
	}				
  }
  
  if(writeY){
      close(outfile_Y)  	
  }
  if(writePi){
      close(outfile_Pi)  	
  }
  if(writeGamma){
      close(outfile_Gamma)  	
  }
  if(writeMu){
      close(outfile_Mu)  	
  }
  if(writeS2){
      close(outfile_S2)  	
  }
  
  if(writeHSall){
	 outfile_HSall <- paste0(fileName, "HSall.txt")
	 write.table(obj$YHSall, outfile_HSall, quote = FALSE, sep="\t", col.names = FALSE, row.names = FALSE)
  }
    
  return(obj)
}


