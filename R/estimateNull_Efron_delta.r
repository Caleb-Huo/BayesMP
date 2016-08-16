##' The proportion of null component is estimated following the procedure by Storey (2002).
##'
##' Estimate null proportion in single study.
##' @title Estimate null proportion in single study.
##' @param x Vector of Z statistics (transformed p value) in single study.
##' @param delta a neighborhood (-delta, delta) around 0 to be used to estiamte null proportion.
##' @return Null proportion in single study.
##' @author Zhiguang Huo
##' @export
##' @examples
##' set.seed(15213)
##' n1 <- 1000
##' n2 <- 100
##' aZ <- c(rnorm(n1),runif(n2,2,3)*sample(n2))
##' delta <- 0.4
##' estimateNull_Efron_delta(aZ, delta)

estimateNull_Efron_delta <- function(x,delta){
	res = numeric(length(delta))
	n = length(x)
	for(i in 1:length(delta)){
		adelta <- delta[i]
		p = 2*pnorm(adelta) - 1
		res[i] <- sum(abs(x)<adelta)/n/p
	}
	return(res)
}
