rm(list=ls())
    
    set.seed(2)
    
    n <- 12
    n.u <- 3
    n.beta <- 3
    
    X <- cbind(1, rep(c(1,1,0,0),3), rep(c(0,0,1,1),3)) 
    Z <- matrix(0, n, n.u)
    
    Z[1:4,1] <- 1
    Z[5:8,2] <- 1
    Z[9:12,3] <- 1
    
    beta <- rep(42, n.beta)
    u <- rnorm(n.u, sd=2)
    y <- X%*%beta + Z%*%u + rnorm(n, sd=2)
    
    loglik_fun <- function(theta,y,X,Z) { 
      n <- length(y) # length of the data vector 
      u <- theta[1:ncol(X)] # parameter vector b
      sigmas <- theta[-(1:ncol(X))] # sig_u and sig
      V <- Z%*%t(Z)*exp(sigmas[1]*2) + diag(n)*exp(sigmas[2]*2) 
      
      L <- chol(V) # cholesky decomposition
      quad <- forwardsolve(t(L),y-X%*%u) 
      l <- -n*log(2*pi)/2 - sum(log(diag(L))) - sum(quad*quad)/2 # eq(2.8)
      -l } 
    initia <- rep(41, n.beta)
    theta <- c(initia,0,0)
    fit <- optim(theta,loglik_fun,method="BFGS",y=y,X=X,Z=Z) 
    fit$par[4:5] <- exp(fit$par[4:5])
    fit$par
  
  