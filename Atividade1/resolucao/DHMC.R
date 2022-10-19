library(tidyverse)
setwd("D:/UNICAMP/Disciplinas/COMPUTACIONAL/Atividade1/resolucao")
theme_set(theme_bw())
lambda=50
df <- data.frame(grid=25:75)%>%
  mutate(poi=dpois(x = grid,lambda = lambda),
         gaus=dnorm(x = grid,mean = lambda,sd = sqrt(lambda)))
ggplot(df,aes(x=grid,y=poi))+
  geom_bar(stat='identity',fill='#ff6666')+
  geom_point(size=.1)+
  labs(y='f(x)',x='x')+
  geom_function(fun=dnorm,args=list(mean=lambda,sd=sqrt(lambda)),size=.1)
ggsave(filename = './imgs/g1.pdf',width = 4,height = 2.5)
df <- data.frame(q = ppoints(20))%>%
  mutate(pois=qpois(q,lambda),
         gauss=qnorm(q,lambda,sqrt(lambda)))
ggplot(df,aes(y=pois,x=gauss))+
  geom_point()+
  labs(x='Gaussian quantiles',y='Poisson quantiles')+
  geom_abline(slope = 1,intercept = 0,col='blue')
ggsave(filename = './imgs/g2.pdf',width = 4,height = 2.5)

mod <- rstan::stan_model('cap_recap.stan')
nj <- c(2,4,2,4,4,2,3,4,2,3,5,10,2,5)
mj <- c(0,0,2,4,4,2,3,2,2,3,5,6,2,5)
MCMC <- rstan::sampling(object=mod,data=list(j=14,ni=nj,mi=mj),chains=4,iter=2000,control=list(max_treedepth=15))
chains <- rstan::extract(MCMC)