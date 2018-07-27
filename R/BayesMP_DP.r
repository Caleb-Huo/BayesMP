##' BayesMP_DP
##'
##' implementation for BayesMP, MCMC part. This is a full Bayesian model, and the alternative distribution is modeled via dirichlet process.
##' @title MCMC for BayesMP_DP
##' @param Z Z statistics. Z should be a p*n matrix.
##' @param gamma initial gamma Estimated null proportation.
##' @param randomGamma update gamma if true.
##' @param beta Non-informative prior: given a gene is DE, the prior probablity this gene is up-regulated, default=1/2
##' @param alpha Concentration parameter for DPs, default=1
##' @param mu0 Mean parameter for base function, default=0
##' @param sigma0 sqrt root of variance parameter for base function, default=10
##' @param sigma sqrt root of variance parameter for DP mixture component, default=1
##' @param trunc truncation parameter for base function (For both positive component and negative component), default=0
##' @param empMu a vector of mean parameter for the null component, default is 0.
##' @param empSD a vector of sd parameter for the null component, default is 1.
##' @param Pi initial pi (DE probablity vector), default rbeta(G, estGamma/(G-estGamma), 1)
##' @param delta initial delta (DE direction probablity vector), default rbeta(G, beta, beta)
##' @param Y initial component indicator. Y should be a p*n matrix. Y = ..., -2, -1, 0, 1, 2, ...
##' @param niter Number of iterations. Default 100, suggest to be 10,000
##' @param burnin Number of burnin period. Default 50, suggest to be 500
##' @param fileName Base fileName for saving fulll mcmc results, hypothesis HS results.
##' @param fullRes binary: 0: do not save full mcmc results; 1: save full mcmc results.
##' @param HSall binary: 0: do not save hypothesis HS results; 1: save hypothesis HS results.
##' @return computing time
##' @useDynLib BayesMP
##' @author Zhiguang Huo
##' @export
##' @examples
##' G <- 5000
##' K <- 10
##' alpha <- 200
##' X0 <- matrix(rnorm(G * K), G, K)
##' X1 <- matrix(rnorm(G * K, 2), G, K)
##' X2 <- matrix(rnorm(G * K, -2), G, K)
##' piall <- rbeta(G, alpha/G, 1)
##' delta <- rbeta(G, 1/2, 1/2)
##' p0 <- 1 - piall
##' p1 <- piall * delta
##' p2 <- piall * (1 - delta)
##' Y <- replicate(K, rMultinom(cbind(p0, p1, p2)))
##' Z <- X0 * (Y == 1) + X1 * (Y == 2) + X2 * (Y == 3)
##' niter <- 2000
##' system.time(BayesMP_DP(Z,fullRes=1))


BayesMP_DP <- function(Z, gamma=NULL, updateGamma = TRUE, beta=1/2, alpha=1, mu0=0, sigma0=10, sigma=1, trunc=0, empMu = rep(0,ncol(Z)), empSD = rep(1,ncol(Z)), 
		niter=100, burnin=50, silence = F, logDotsPerLine = 50, fileName='BayesMP_DP_', 
		writeY = F, writePi = F, writeDelta = F, writeGamma = F, writeHSall = T){
			
	G <- nrow(Z) ## number of genes
	S <- ncol(Z) ## number of studies

	if(is.null(gamma)){
		pp <- pmin(pnorm(-abs(Z)) * 2, 1)
		gamma <- max(1 - pi0est(pp)$pi0, 0.01)		
	}	
	#cat("Initial gamma: ",gamma, "\n")

	Pi <- rbeta(G, gamma, 1-gamma)
	delta <- rbeta(G, beta, beta)
	Y0 <- matrix(0,G,S)		
	for(s in 1:S){
		Y_s <- numeric(G)
		Y_s[Z[,s]>2] = 1
		Y_s[Z[,s]< -2] = -1
		Y0[,s] <- Y_s
	}
	Y <- Y0		
		
	fileName_Y <- paste0(fileName,"Y.txt")
	fileName_Pi <- paste0(fileName,"Pi.txt")
	fileName_Delta <- paste0(fileName,"Delta.txt")
	fileName_Gamma <- paste0(fileName,"Gamma.txt")
	fileName_HSall <- paste0(fileName,"HSall.txt")
		
    if (!silence) 
            cat("performing MCMC, i = 1,2,..., niter (= ", niter, ")  [one \".\" per sample]:\n", sep = "")

	obj <- .C('mcmc_R3',G=as.integer(G),S=as.integer(S),Z=as.double(Z),gamma=as.double(gamma),randomGamma=as.integer(updateGamma),
			empMu=as.double(empMu), empSD=as.double(empSD),
			beta=as.double(beta),alpha=as.double(alpha),mu0=as.double(mu0),
			sigma0=as.double(sigma0),sigma=as.double(sigma),trunc=as.double(trunc),pi=as.double(Pi),delta=as.double(delta),Y=as.integer(Y0),
			niter=as.integer(niter),burnin=as.integer(burnin),silence=as.integer(silence),logDotsPerLine=as.integer(logDotsPerLine),			
			fileName_Y=as.character(fileName_Y),fileName_Pi=as.character(fileName_Pi),fileName_Delta=as.character(fileName_Delta),fileName_Gamma=as.character(fileName_Gamma),fileName_HSall=as.character(fileName_HSall),
			writeY = as.integer(writeY), writePi = as.integer(writePi), writeDelta = as.integer(writeDelta), writeGamma = as.integer(writeGamma), writeHSall = as.integer(writeHSall)
			)
	return(obj)
}
