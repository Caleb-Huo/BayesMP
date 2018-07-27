##' Baysian FDR
##'
##' This function is used to calculate Baysian FDR.
##' The result is similar to q-value.
##' @title Bayesian FDR
##' @param belief A posterior belief vector.
##' @return Bayesian q-value, which is analog to frequentists' q-value
##' @author Zhiguang Huo <zhuo@ufl.edu>
##' @export
##' @examples
##' set.seed(15213)
##' aBeliefVec <- runif(1000)
##' bayesianQValue <- BayesianFDR(aBeliefVec)
BayesianFDR <- function(belief){
	pvalue <- 1 - belief
	pvalue_order <- order(pvalue)
	sortedP <- pvalue[pvalue_order]
	sortedQ <- cumsum(sortedP)/(1:length(sortedP))
	qvalue <- sortedQ[match(1:length(pvalue_order),pvalue_order)]
	return(qvalue)
}
