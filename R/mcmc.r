##' BayesMP
##'
##' implementation for BayesMP, MCMC part
##' @title MCMC for BayesMP
##' @param Z Z statistics. Z should be a p*n matrix.
##' @param initial gamma Estimated null proportation.
##' @param beta Non-informative prior: given a gene is DE, the prior probablity this gene is up-regulated, default=1/2
##' @param alpha Concentration parameter for DPs, default=1
##' @param mu0 Mean parameter for base function, default=0
##' @param sigma0 sqrt root of variance parameter for base function, default=10
##' @param sigma sqrt root of variance parameter for DP mixture component, default=1
##' @param trunc truncation parameter for base function (For both positive component and negative component), default=0
##' @param Pi initial pi (DE probablity vector), default rbeta(G, estGamma/(G-estGamma), 1)
##' @param delta initial delta (DE direction probablity vector), default rbeta(G, beta, beta)
##' @param Y initial component indicator. Y should be a p*n matrix. Y = ..., -2, -1, 0, 1, 2, ...
##' @param niter Number of iterations. Default 100, suggest to be 10,000
##' @param burnin Number of burnin period. Default 50, suggest to be 500
##' @param fileName Base fileName for saving fulll mcmc results, hypothesis HS results.
##' @param fullRes binary: 0: do not save full mcmc results; 1: save full mcmc results.
##' @param HSall binary: 0: do not save hypothesis HS results; 1: save hypothesis HS results.
##' @return computing time
##' @author Zhiguang Huo
##' @export
##' @examples
##' n1 <- 1000
##' n2 <- 100
##' S <- 3
##' set.seed(15213)
##' aZ <- c(rnorm(n1),rnorm(n2,3)*sample(c(1,-1),n2,replace=TRUE))
##' set.seed(15214)
##' bZ <- c(rnorm(n1),rnorm(n2,3)*sample(c(1,-1),n2,replace=TRUE))
##' set.seed(15215)
##' cZ <- c(rnorm(n1),rnorm(n2,3)*sample(c(1,-1),n2,replace=TRUE))
##' Z <- cbind(aZ, bZ, cZ)
##' G <- nrow(Z)
##' mcmc(Z)


mcmc <- function(Z, gamma=NULL, randomGamma = TRUE, beta=1/2, alpha=1, mu0=0, sigma0=10, sigma=1, trunc=0, empMu = rep(0,ncol(Z)), empSD = rep(0,ncol(Z)), Pi=NULL, delta=NULL, Y=NULL, niter=100, burnin=50, fileName='BayesMP_mcmc', fullRes=1, HSall=1){
	G <- nrow(Z) ## number of genes
	S <- ncol(Z) ## number of studies

	if(is.null(gamma)){
		pp <- pmin(pnorm(-abs(Z)) * 2, 1)
		gamma <- max(1 - pi0est(pp)$pi0, 0.01)		
	}	
	cat("Initial gamma: ",gamma, "\n")

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
		
	keepTime <- system.time(
		obj <- .C('mcmc_R3',G=as.integer(G),S=as.integer(S),Z=as.double(Z),gamma=as.double(gamma),randomGamma=as.integer(randomGamma),
				empMu=as.double(empMu), empSD=as.double(empSD),
				beta=as.double(beta),alpha=as.double(alpha),mu0=as.double(mu0),
				sigma0=as.double(sigma0),sigma=as.double(sigma),trunc=as.double(trunc),pi=as.double(Pi),delta=as.double(delta),Y=as.integer(Y0),
				niter=as.integer(niter),burnin=as.integer(burnin),fileName=as.character(fileName),fullRes=as.integer(fullRes),HSall=as.integer(HSall))
		)
	return(keepTime)
}

