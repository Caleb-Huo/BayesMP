#### Step 0. Utilities
rMultinom <- function (p) 
{
    if (any(p < 0)) 
        stop("non-positive probability")
    if (!isTRUE(all.equal(rowSums(p), rep(1, nrow(p))))) 
        stop("The sums of probabilities have to be 1.")
    k <- ncol(p)
	if (k!=3)
		stop("ncol of the input has to be three")
    n <- nrow(p)
    csm <- ifelse(outer(1:k, 1:(k - 1), FUN = `<=`), 1, 0)
    cp <- p %*% csm
    U <- runif(n)
    rowSums(sweep(cp, 1, U, `<=`)) + 1
}

#### Step 1. Update Y

update.y <- function(X, f, pi, delta, log = FALSE) {
  K <- ncol(X)
  G <- nrow(X)
  M <- length(f)
  
  Pi <- cbind(1 - pi, pi * delta, pi * (1 - delta))
  lik <- lapply(1:M, function(m) sapply(1:K, function(k) f[[m]][[k]](X[, k])))
  if(log) {
    post1 <- lapply(1:M, function(m) lik[[m]] + log(Pi[, m]))
    const1 <- do.call(pmax, post1)
    post2 <- lapply(post1, function(x) exp(x - const1))
    const2 <- Reduce("+", post2)
    post <- lapply(post2, function(x) x/const2)
  } else {
    post1 <- lapply(1:M, function(m) lik[[m]] * Pi[, m])
    const1 <- Reduce("+", post1)
    post <- lapply(post1, function(x) x/const1)
  }
  postv <- sapply(post, as.vector)
  Y <- matrix(rMultinom(postv), nrow = G, ncol = K)  ## 1 = 0, 2 = +1, 3 = -1
  return(Y)
}

#### Step 2. Update pi

update.pi <- function(Y, gamma) {
  G <- nrow(Y)
  K <- ncol(Y)
  Y2 <- rowSums(Y == 2)
  Y3 <- rowSums(Y == 3)
  rbeta(G, gamma + Y2 + Y3, 1 + K - Y2 -Y3)  
}

update.delta <- function(Y, beta = 1/2) {
  G <- nrow(Y)
  K <- ncol(Y)
  Y2 <- rowSums(Y == 2)
  Y3 <- rowSums(Y == 3)
  rbeta(G, beta + Y2, beta + Y3)
}

#### Step 3. Update f (optional)

update.fnorm.one <- function(f, xx, yy, m, log = FALSE, a = -Inf, b = Inf) {
  mu.old <- attr(f, "mu")
  var.old <- attr(f, "var")
  ny <- sum(yy == m)
  Xbar <- mean(xx[yy == m])
  mu <- rtruncnorm(1, a = a, b = b, mean = Xbar, sd = sqrt(var.old/ny))
  SS <- sum((xx[yy == m] - mu)^2)
  var <- rinvgamma(1, ny/2, SS/2)
  fout <- function(x) dnorm(x, mu, sqrt(var), log = log)
  attr(fout, "mu") <- mu
  attr(fout, "var") <- var
  return(fout)
}

update.f <- function(f, X, Y, log = FALSE) {
  M <- length(f)
  K <- ncol(X)
  a <- c(0, -Inf)
  b <- c(Inf, 0)
  for(m in 2:M) {
    for(k in 1:K) {
      if (M == 3) {
        f[[m]][[k]] <- update.fnorm.one(f[[m]][[k]], X[, k], Y[, k], m, log, a[m - 1], b[m - 1])
      } else {
        f[[m]][[k]] <- update.fnorm.one(f[[m]][[k]], X[, k], Y[, k], m, log)
      }
    }
  }
  return(f)
}

#### Step 4. Update gamma (optional)
update.gamma <- function(piDE, preGamma, MHsd){
	aprop = rnorm(1,preGamma,MHsd);
	if(aprop <= 0 | aprop >= 1){
      gamma <- preGamma;
  	  MHsd <- MHsd / 1.01;	
  	  MHcount <- 0;		
	} else if(runif(1) < exp(loglikelihood(aprop, piDE) - loglikelihood(preGamma, piDE))){
	  gamma <- aprop;
	  MHcount <- 1;
	  MHsd <- MHsd * 1.01;
	} else {
        gamma <- preGamma;
    	MHsd <- MHsd / 1.01;	
   	 	MHcount <- 0;		
	}
	return(list(gamma=gamma, MHsd=MHsd, MHcount=MHcount))	
}

loglikelihood <- function(gamma, piDE){
	sum(dbeta(piDE,gamma,1-gamma, log=T))
}


update.HSall <- function(obj){
	rowSumY <- rowSums(obj$Y>1)
	YHSall <- obj$YHSall
	for(i in 1:ncol(YHSall)){
		thisSig <- rowSumY>=i
		YHSall[thisSig,i] <- YHSall[thisSig,i] + 1
	}
	YHSall
}

#### One iteration
iterate.one <- function(obj, X) {
  M <- length(obj$f)
  obj$Y <- update.y(X, obj$f, obj$pi, obj$delta, log = obj$log)
  obj$pi <- update.pi(obj$Y, obj$gamma)
  obj$delta <- update.pi(obj$Y, obj$beta)
  if(obj$update.f) obj$f <- update.f(obj$f, X, obj$Y, obj$log)
  if(obj$update.gamma){
	gammaUpdateList <- update.gamma(piDE=obj$pi, preGamma=obj$gamma, MHsd = obj$MHsd)	  
  	obj$gamma <- gammaUpdateList$gamma
  	obj$MHsd <- gammaUpdateList$MHsd
  	obj$MHcount <- obj$MHcount + gammaUpdateList$MHcount 
  } 
  return(obj)
}

