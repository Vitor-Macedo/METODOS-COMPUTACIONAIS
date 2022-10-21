library(readr)
library(tidyverse)
setwd("D:/UNICAMP/Disciplinas/COMPUTACIONAL/Atividade1/resolucao")
dados <- read_table("./dados.txt", 
                    col_names = FALSE)
colnames(dados) <- c('onca',paste0('armadilha',1:14))
nj <- dados%>%
  pivot_longer(-onca,values_to = "capturada",
               names_to = "armadilha",names_prefix = 'armadilha',
               names_transform = list(armadilha = as.numeric))%>%
  group_by(armadilha)%>%
  summarise(nj=sum(capturada))
mj <- dados%>%
  pivot_longer(-onca,values_to = "capturada",names_to = "armadilha",
               names_prefix = 'armadilha',
               names_transform = list(armadilha = as.numeric))%>%
  group_by(onca)%>%
  mutate(mj=case_when(cumsum(capturada)*capturada>1 ~ 1 , TRUE ~ 0))%>%
  group_by(armadilha)%>%
  summarise(mj=sum(mj))
nj <- nj$nj 
Ngibbs=10e4
iter0 <- c(rpois(1,50),runif(14))
chain <- matrix(0,ncol = Ngibbs,nrow=length(iter0))
chain[,1] <- iter0
rownames(chain) <- c("N",paste0("p",1:14))
lambda=50
r <- sum(nj$nj)-sum(mj$mj)
nj <- nj$nj 
for (i in 2:Ngibbs) {
  lambda_pos <- lambda * prod(1-chain[-1,i-1])
  N_i_plus_1 <- rpois(n=1,lambda = lambda_pos) + r
  chain[1,i] <- N_i_plus_1 
  for (j in 1:length(nj)) {
    aux <- rbeta(n = 1,shape1 = nj[j]+1,shape2 = N_i_plus_1-nj[j]+1)
    chain[j+1,i] <- aux
  }
}

gibbs_popSize <- function(iter=500,chains=4,r=12,nj,lambda=50){
  sample <- NULL
  aux <- matrix(0,nrow = iter,ncol=15)
  colnames(aux) <- c("N",paste0("p",1:14))
  for (k in 1:chains) {
    iter0 <- c(rpois(1,50),runif(14))
    aux[1,] <- iter0
    for (i in 2:iter) {
      lambda_pos <- lambda * prod(1-aux[i-1,-1])
      N_i_plus_1 <- rpois(n=1,lambda = lambda_pos) + r
      aux[i,1] <- N_i_plus_1 
      for (j in 1:length(nj)) {
        aux2 <- rbeta(n = 1,shape1 = nj[j]+1,shape2 = N_i_plus_1-nj[j]+1)
        aux[i,j+1] <- aux2
      }
    }
    sample <- bind_rows(sample,data.frame(aux,chain=as.character(k),iteration=1:iter))
  }
  return(sample)
}

gbs <- gibbs_popSize(iter = 1000,chains = 4,r = 12,nj = nj,lambda = 50)
gbs%>%
  ggplot(aes(x=iteration,y=N,color=chain))+
  geom_line(alpha=.5)+
  theme_bw()
ggsave(filename = paste0('./imgs/N','.pdf'),width = 4,height = 2.5)

for (i in 1:14) {
  gbs%>%
    select(-N)%>%
    pivot_longer(cols = c(-iter,-chain),names_to = 'param',values_to = 'pi')%>%
    filter(param==paste0('p',i))%>%
    ggplot(aes(x=iter,y=pi,color=chain))+
    geom_line(alpha=.5)+
    labs(y="p")+
    theme_bw()
ggsave(filename = paste0('./imgs/p',i,'.pdf'),width = 4,height = 2.5)
}

# ACF N
Nacf <- acf(gbs%>%filter(chain==1)%$%N)
Nacfdf <- with(Nacf, data.frame(lag, acf))
ggplot(data = Nacfdf, mapping = aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) +
  geom_segment(mapping = aes(xend = lag, yend = 0))+
  geom_hline(yintercept = c(-.025,.025),linetype=2,col='blue')+
  theme_bw()
ggsave(filename = paste0('./imgs/Nacf','.pdf'),width = 4,height = 2.5)
gbs%>%
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
