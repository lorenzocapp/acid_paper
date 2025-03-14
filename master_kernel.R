rm(list=ls())
source("functions_kernel.R")
source("utilities.R")
library(diagL1)
arg<-commandArgs(TRUE)
idx <- as.integer(arg[1])
#idx <- 1
#sample data
n <- 200

#Load toy dataset
y <- scan(paste("data/data_n",n,"_",idx,".csv",sep=""),sep = ",")
true <- scan(paste("data_true_d/data_true_d_n",n,"_",idx,".csv",sep=""),sep = ",")

#check
grid <- seq(-4,3.95,0.05)*sd(y)+mean(y)


#Run
M <- 1000
B <- 400
df <- data.frame(n=integer(),data=integer(),method=character(),
                 env=numeric(),bias=numeric(),rwd=numeric(),
                 dev=numeric())
for (k in c("gaussian","laplace","uniform")){
  for (al in c("standard","recursive")){
    for (b in c("SJ","nrd","nrd0","ucv")){
      h <- initialization_kernel_resampling(y,M,bw=b)
      posterior <- ps_posterior(y=y,grid=grid,h=h,M=M,B=B,kernel=k,algo=al,bw=b)
      method_name <- paste("kernel_",substr(k,1,1),"_",substr(al,1,1),"_",b,sep="")
      df <- rbind(df,data.frame(n=n,data=idx,method=method_name,env=envelope(posterior,true),
                       bias=bias(posterior,true),rwd=relwid(posterior,true),
                       dev=dev(posterior,true)))
    }
  }
}

name <- paste("out/out_kernel_n",n,"_",idx,".RData",sep="")
save(df,file=name)

