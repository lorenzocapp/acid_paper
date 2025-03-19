rm(list=ls())
source("../utilities.R")

#arg<-commandArgs(TRUE)
#idx <- as.integer(arg[1])
#idx <- 1
df <- data.frame(n=integer(),type=character(),data=integer(),method=character(),
                 env=numeric(),abs_bias=numeric(),rel_bias=numeric(),abs_rwd=numeric(),rel_rwd=numeric(),
                 abs_dev=numeric(),rel_dev=numeric(),MSE_med=numeric(),MSE_mean=numeric())
for (idx in 1:100){
  for (n in c(50,200)){
    for (rho in c(0.4, 0.6, 0.8)){
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
      method_name <- paste("copula_",rho,sep="")
      df <- rbind(df,data.frame(n=n,type="v1",data=idx,method=method_name,env=envelope(posterior,true),
                                abs_bias=bias(posterior,true)$abs,rel_bias=bias(posterior,true)$rel,
                                abs_rwd=relwid(posterior,true)$abs,rel_rwd=relwid(posterior,true)$rel,
                                abs_dev=dev(posterior,true)$abs,rel_dev=dev(posterior,true)$rel,
                                MSE_med=NA,MSE_mean=NA))
      

    }
    name <- paste("out/out_copula_n",n,"_",idx,".RData",sep="")
    save(df,file=name)
 }
}