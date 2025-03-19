rm(list=ls())
source("../functions_kernel.R")
source("../utilities.R")
library(diagL1)
arg<-commandArgs(TRUE)
idx <- as.integer(arg[1])
#idx <- 1
#sample data
for (n in c(50,200)){

#Load toy dataset
y <- scan(paste("data/data_n",n,"_",idx,".csv",sep=""),sep = ",")
true <- scan(paste("data_true_d/data_true_d_n",n,"_",idx,".csv",sep=""),sep = ",")

#check
grid <- seq(-4,3.95,0.05)*sd(y)+mean(y)


#Run
M <- 1000
B <- 400

df <- data.frame(n=integer(),type=character(),data=integer(),method=character(),
                 env=numeric(),abs_bias=numeric(),rel_bias=numeric(),abs_rwd=numeric(),rel_rwd=numeric(),
                 abs_dev=numeric(),rel_dev=numeric(),MSE_med=numeric(),MSE_mean=numeric())
for (k in c("gaussian","laplace","uniform")){
  for (al in c("standard","recursive")){
    for (b in c("SJ","ucv")){
      for (d in c(0.5,0.7,0.9,1,1.1)){
        h <- d*initialization_kernel_resampling(y,M,bw=b)
        posterior <- ps_posterior(y=y,grid=grid,h=h,M=M,B=B,kernel=k,algo=al,bw=b)
        method_name <- paste("kernel_",substr(k,1,1),"_",substr(al,1,1),"_",d,b,sep="")
        df <- rbind(df,data.frame(n=n,type="v1",data=idx,method=method_name,env=envelope(posterior,true),
                                  abs_bias=bias(posterior,true)$abs,rel_bias=bias(posterior,true)$rel,
                                  abs_rwd=relwid(posterior,true)$abs,rel_rwd=relwid(posterior,true)$rel,
                                  abs_dev=dev(posterior,true)$abs,rel_dev=dev(posterior,true)$rel,
                                  MSE_med=NA,MSE_mean=NA))
      }
    }
  }
}

name <- paste("out/out_kernel_n",n,"_",idx,".RData",sep="")
save(df,file=name)
}
