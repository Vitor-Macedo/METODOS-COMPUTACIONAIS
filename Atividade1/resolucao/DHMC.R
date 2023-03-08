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
  dUdp <- (N*pi-ni)/(pi*(1-pi))
  dUdN <- -(digamma(N+1)/gamma(N+1)) + (digamma(N-r+1)/gamma(N-r+1)) - sum(log(1-pi)) + .02*(N-50)
  grad <- c(dUdN,dUdp)
  return(grad)
}

iter=2000
sample <- matrix(0,nrow = iter,ncol=15)
colnames(sample) <- c("N",paste0("p",1:14))
sample[1,] <- c(15,runif(14))
M <- cor(gbs%>%select(N,p1:p14))
L=50;delta=.001
for (i in 2:iter) {
  rho0 <- theta0 <- matrix(0,nrow = L,ncol = 15)
  rho0[1,] <- rmvnorm(1, mean = rep(0,15), sigma = M)
  theta0[1,] <- sample[i-1,]
  #leapfrog
  for (j in 2:L) {
    rhotemp <- rho0[j-1,] - (delta/2)*dUdt(theta = theta0[j-1,],ni=nj)
    theta0[j,] <- theta0[j-1,] + delta*rhotemp
    rho0[j,] <- rhotemp - (delta/2)*dUdt(theta = theta0[j,],ni=nj)
  }
  alpha=min(1,exp( -U(theta = theta0[L,],ni=nj) + U(theta = theta0[1,],ni=nj)- K(M,rho0[L,]) + K(M,rho0[1,])))
  if(runif(1)<alpha){
    sample[i,] <- theta0[L,]
  }else{
    sample[i,] <- sample[i-1,]
  }
}
sample%>%
  as.data.frame()%>%
  mutate(iteration=1:iter)%>%
  select(p1:p14,N,iteration)%>%
  pivot_longer(cols = -iteration,names_to = 'param',values_to = 'valor')%>%
  mutate(param=factor(param,levels = c('N',paste0('p',1:14))))%>%
  group_by(param)%>%
  summarise(media=mean(valor),
            mediana=median(valor),
            lower=quantile(valor,.025)%>%round(3),
            upper=quantile(valor,.975)%>%round(3))%>%
  mutate(ic=paste0('[',lower,';',upper,']'))%>%
  select(param,media,mediana,ic)%>%
  xtable::xtable(digits = 3,rownames=F)
sample%>%
  as.data.frame()%>%
  mutate(iteration=1:iter)%>%
  ggplot(aes(x=iteration,y=p8))+
  geom_line()
for (i in 1:14) {
  sample%>%
    as.data.frame()%>%
    mutate(iteration=1:iter)%>%
    select(-N)%>%
    pivot_longer(cols = c(-iteration),names_to = 'param',values_to = 'pi')%>%
    filter(param==paste0('p',i))%>%
    ggplot(aes(x=iteration,y=pi))+
    geom_line(alpha=.5)+
    labs(y="p")+
    theme_bw()
  ggsave(filename = paste0('./imgs/HMC_R',i,'.pdf'),width = 4,height = 2.5)
}
# Gráficos do stan
load('mcmc2.Rdata')
cadeias <- rstan::extract(MCMC)
stan <- cbind(cadeias$N,cadeias$p)
colnames(stan) <- c("N",paste0("p",1:14))
stan%>%
  as.data.frame()%>%
  mutate(iteration=rep(1:2000,4))%>%
  select(p1:p14,N,iteration)%>%
  pivot_longer(cols = -iteration,names_to = 'param',values_to = 'valor')%>%
  mutate(param=factor(param,levels = c('N',paste0('p',1:14))))%>%
  group_by(param)%>%
  summarise(media=mean(valor),
            mediana=median(valor),
            lower=quantile(valor,.025)%>%round(3),
            upper=quantile(valor,.975)%>%round(3))%>%
  mutate(ic=paste0('[',lower,';',upper,']'))%>%
  select(param,media,mediana,ic)%>%
  xtable::xtable(digits = 3)
rstan::traceplot(MCMC,pars=c('N','p'))



