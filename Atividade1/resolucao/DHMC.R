library(mvtnorm)
U <- function(theta,r=12,ni){
  N <- theta[1]
  pi <- theta[2:length(theta)]
  Ux <- log(gamma(N+1))-log(gamma(N-r+1)) + sum(ni*log(pi)+(N-ni)*log(1-pi)) - .01*(N-50)^2
  return(-Ux)
}
K <- function(M,rho){
  .5*t(rho)%*%M%*%rho
}
dUdt <- function(theta,r=12,ni){
  N <- theta[1]
  pi <- theta[2:length(theta)]
  dUdp <- (N*p-ni)/(pi*(1-pi))
  dUdN <- -(digamma(N+1)/gamma(N+1)) + (digamma(N-r+1)/gamma(N-r+1)) - sum(log(1-pi)) + .02*(N-50)
  grad <- c(dUdN,dUdp)
  return(grad)
}
dKdrho <- function(M,rho){
  
}
iter=500
sample <- matrix(0,nrow = iter,ncol=15)
colnames(sample) <- c("N",paste0("p",1:14))
sample[1,] <- c(rpois(1,50),runif(14))
M <- cor(gbs%>%select(N,p1:p14))
L=50;delta=.01
for (i in 2:iter) {
  rho0 <- theta0 <- matrix(0,nrow = L,ncol = 15)
  rho0[1,] <- rmvnorm(1, mean = rep(0,15), sigma = M)
  theta0[1,] <- sample[i-1,]
  #leapfrog
  for (j in 1:L) {
    
  }
  
}