##' Dissimilarity matrix calculation.
##'
##' Based on DE posterior probablity of gene*study matrix,
##' calculate gene-gene distance
##' @title Distance
##' @param plus p*n matrix. Each element represents the posterior probablity of a gene to be positive DE.
##' @param minus p*n matrix. Each element represents the posterior probablity of a gene to be negative DE.
##' @param nsample Total number of posterior samples excluding burning period.
##' @param alpha Forge parameter
##' @param method Choose one from those: cos, KL, L1, L2, L2_2D, Hellinger.
##' @return Dissimilarity matrix. p*p
##' @author Zhiguang Huo
##' @export
##' @examples
##' G<-10
##' S<-3
##' nsample <- 10000
##' plus <- matrix(c(sample(nsample,G*S,replace = TRUE),rep(0,G*S)),2*G,S)
##' minus <- matrix(c(rep(0,G*S),sample(nsample,G*S,replace = TRUE)),2*G,S)
##' dissimilarity <- distance(plus, minus, nsample)

distance <- function(plus,minus,nsample,alpha=0,method='cos'){
	## alpha is for laplace smoothing
	## each row represent a gene

	if(method!='cos' & method!='KL' & method!='L1' & method!='L2' & method!='L2_2D' & method!='Hellinger'){
		stop('please specify correct distance cos, KL, L1, L2, L2_2D, Hellinger')
	}

	cosineDist <- function(p,q){
		1 - sum(p*q)/sqrt(sum(p^2)*sum(q^2))
	}

	KLDist <- function(p,q){
		sum(p*log(p/q)+q*log(q/p)) /2
	}

	L1Dist <- function(p,q){
		sum(abs(p-q))
	}

	L2Dist <- function(p,q){
		sum((p-q)^2)
	}

	L2Dist2D <- function(p,q){
		sum((p[-3]-q[-3])^2)
	}

	HellingerDist <- function(p,q){
		sqrt(sum((sqrt(p)-sqrt(q))^2))/sqrt(2)
	}

	simpleDist <- function(p,q){
		a <- sum( abs(p[1:2] - q[1:2]))
		if(a==0) return	(0)
		res <- sum( abs(p[1:2] - q[1:2]) / (p[1:2] + q[1:2]) )
	}


	selfDistance <- function(a,b,fun){
	  outer(a, b, function(x,y) vapply(seq_along(x), function(i) fun(x[[i]], y[[i]]), numeric(1)))
	}

	K <- ncol(plus)
	netro = nsample - plus - minus
	Pplus = (plus + alpha)/(nsample + alpha*K)
	Pminus = (minus + alpha)/(nsample + alpha*K)
	Pnetro = (netro + alpha)/(nsample + alpha*K)

	resdist = 0
	for(k in 1:K){
		tmpData = cbind(Pplus[,k],Pminus[,k],Pnetro[,k])
		tmpList <- lapply(seq_len(nrow(tmpData)), function(i) tmpData[i,]/sum(tmpData[i,]))
		if(method=='cos'){
			adist = selfDistance(tmpList,tmpList,cosineDist)
			}
		else if(method=='KL'){
			if(alpha==0)
				warnings('for KL distance, it is better to specify non zero alpha!')
			adist = selfDistance(tmpList,tmpList,KLDist)
			}
		else if(method=='L1'){
			adist = selfDistance(tmpList,tmpList,L1Dist)
			}
		else if(method=='L2'){
			adist = selfDistance(tmpList,tmpList,L2Dist)
			}
		else if(method=='L2_2D'){
			adist = selfDistance(tmpList,tmpList,L2Dist2D)
			}
		else if(method=='Hellinger'){
			adist = selfDistance(tmpList,tmpList,HellingerDist)
			}
		else {stop('there is a bug for function distance 1.')}
		resdist = resdist + adist/K
	}
	resdist
}
