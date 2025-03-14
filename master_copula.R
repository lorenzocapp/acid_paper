rm(list=ls())
source("utilities.R")

#arg<-commandArgs(TRUE)
#idx <- as.integer(arg[1])
#idx <- 1
for (idx in 1:100){
n <- 200
y <- scan(paste("data/data_n",n,"_",idx,".csv",sep=""),sep = ",")
true <- scan(paste("data_true_d/data_true_d_n",n,"_",idx,".csv",sep=""),sep = ",")

#Load output from python
library(reticulate)

np <- import("numpy")
name <- paste("output_py/out_python_n",n,"_",idx,".npy",sep="")
mat <- np$load(name)
mat <- np$matrix(exp(mat))/sd(y)

grid <- seq(-4,3.95,0.05)*sd(y)+mean(y)

q2.5 <- apply(data.frame(mat), 2, quantile, probs=0.025)
q97.5 <- apply(data.frame(mat), 2, quantile, probs=0.975)
q50 <- apply(data.frame(mat), 2, quantile, probs=0.5)
mean <- apply(data.frame(mat), 2, mean)
posterior <- list(grid=grid,mean=mean,q2.5=q2.5,q97.5=q97.5,q50=q50)



#evaluation_metric

df <- data.frame(n=n,data=idx,method="copula",env=envelope(posterior,true),
                  bias=bias(posterior,true),rwd=relwid(posterior,true),
                  dev=dev(posterior,true))

name <- paste("out/out_copula_n",n,"_",idx,".RData",sep="")
save(df,file=name)
}
