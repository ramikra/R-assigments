## set.seed for reproducibility 
set.seed(1)
n <- 100
## generate random values from the uniform distribution
x <- sort(runif(n))
## the formula to generate data
## y = a1*x1 + a2*x1^2 + e 
y <- 1+ 0.1*x+3.2*x^2+rnorm(n)*0.1

## create the model matrix 
## poly() generates orthogonal polynomials
## recall > sqrt((1/sqrt(n))^2 *n) = 1
## orthonormal, if <e,e>=1
X <- cbind(n^-.5,poly(x,3))
## finding the eigenvalues
d <- svd(X)$d 
## computing kappa
max(d)/min(d) 
## comparing to the machine precision
1/.Machine$double.eps
## fitting using lm
b = lm(y~X-1)
plot(x,y,cex.lab=1,cex.axis=1) 
lines(x,fitted(b),lwd=2, col = "darkred") 
## summary of the model
summary(b)


###############
###############
###############
rm(list=ls())

set.seed(1)

n <- 12
n.u <- 3
n.beta <- 3

X <- cbind(1, rep(c(1,1,0,0),3), rep(c(0,0,1,1),3)) 
Z <- matrix(0, n, n.u)

Z[1:4,1] <- 1
Z[5:8,2] <- 1
Z[9:12,3] <- 1

beta <- rep(42, n.beta)
u <- rnorm(n.u)
y <- X%*%beta + Z%*%u + rnorm(n)

loglik_fun <- function(theta,y,X,Z) { 
  n <- length(y) # length of the data vector 
  u <- theta[1:ncol(X)] # parameter vector b
  sigmas <- theta[-(1:ncol(X))] # sig_u and sig
  V <- Z%*%t(Z)*exp(sigmas[1]*2) + diag(n)*exp(sigmas[2]*2) 
  
  L <- chol(V) # cholesky decomposition
  quad <- forwardsolve(t(L),y-X%*%u) 
  l <- -n*log(2*pi)/2 - sum(log(diag(L))) - sum(quad*quad)/2 # eq(2.8)
  -l } 
theta <- c(beta,0,0)
fit <- optim(theta,loglik_fun,method="BFGS",y=y,X=X,Z=Z) 
fit$par


#####
#####

response <- y
analyst <- rep(c(1,1,0,0),3)
part <- c(rep(1,4),rep(2,4),rep(3,4))
df <- data.frame(response,analyst,part)


####

###############
rm(list=ls())

#set.seed(1)


set.seed(10)
n <- 100
n.b <- 10
n.beta <- 5

X <- cbind(1, matrix(runif(n*n.beta-n), n, n.beta-1)) 

Z <- matrix(runif(n*n.b), n, n.b)

beta <- rep(2, n.beta)
b <- rnorm(n.b)

y <- X%*%beta + Z%*%b + rnorm(n)



loglik_fun <- function(theta,y,X,Z) { 
  n <- length(y) # length of the data vector 
  u <- theta[1:ncol(X)] # parameter vector b
  sigmas <- theta[-(1:ncol(X))] # sig_u and sig
  V <- Z%*%t(Z)*exp(sigmas[1]*2) + diag(n)*exp(sigmas[2]*2) 
  
  L <- chol(V) # cholesky decomposition
  quad <- forwardsolve(t(L),y-X%*%u) 
  l <- -n*log(2*pi)/2 - sum(log(diag(L))) - sum(quad*quad)/2 # eq(2.8)
  -l } 
theta <- c(beta,0,0)
fit <- optim(theta,loglik_fun,method="BFGS",y=y,X=X,Z=Z) 
fit$par[1:3]
exp(fit$par)[4:5]


#####
#####

