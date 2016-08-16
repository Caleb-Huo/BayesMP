##' The proportion of null component is estimated following the procedure by Storey (2002).
##'
##' Estimate null proportion for multiple studies.
##' @title Estimate null proportion for multiple studies.
##' @param X Matrix of Z statistics (transformed p value) p*S, p is number of genes and S is number of studies.
##' @param delta a neighborhood (-delta, delta) around 0 to be used to estiamte null proportion.
##' @return Null proportion estimation of multiple studies.
##' @author Zhiguang Huo
##' @export
##' @examples
##' n1 <- 1000
##' n2 <- 100
##' S <- 3
##' set.seed(15213)
##' aZ <- c(rnorm(n1),rnorm(n2,3)*sample(n2))
##' set.seed(15214)
##' bZ <- c(rnorm(n1),rnorm(n2,3)*sample(n2))
##' set.seed(15215)
##' cZ <- c(rnorm(n1),rnorm(n2,3)*sample(n2))
##' Z <- cbind(aZ, bZ, cZ)
##' delta <- 0.4
##' estimateNull_Efron_delta(Z, delta)

estimatePrior <- function(X,delta){
	priorsss <- numeric(ncol(X))
	for(i in 1:ncol(X)){
		priorsss[i] <- 	estimateNull_Efron_delta(X[,i],d)
	}
	priorsss
}
