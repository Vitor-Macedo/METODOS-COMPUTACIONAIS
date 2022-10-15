library(readr)
library(tidyverse)
dados <- read_table("D:/UNICAMP/Disciplinas/COMPUTACIONAL/Atividade1/resolucao/dados.txt", 
                    col_names = FALSE)
colnames(dados) <- c('Individuo',paste0('armadilha',1:14))
nj <- dados%>%
  select(-Individuo)%>%
  colSums()
Ngibbs=50e4
x0 <- c(rpois(1,50),runif(14))
chain <- matrix(0,ncol = Ngibbs,nrow=length(x0))
chain[,1] <- x0
rownames(chain) <- c("N",paste0("p",1:14))
lambda=50
r <- 12
for (i in 2:Ngibbs) {
  lambda_pos <- lambda * prod(1-chain[-1,i-1])
  N_i_plus_1 <- rpois(n=1,lambda = lambda_pos) + r
  chain[1,i] <- N_i_plus_1 
  for (j in 1:length(nj)) {
    aux <- rbeta(n = 1,shape1 = nj[j]+1,shape2 = N_i_plus_1-nj[j]+1)
    chain[j+1,i] <- aux
  }
}
final_sample <- chain[,seq(Ngibbs/2,Ngibbs,by=10)]
plot(final_sample[1,],type='l')
rowMeans(final_sample)
