rm(list=ls())
source("../functions_kernel.R")
source("../functions_kernel_regr.R")
source("../utilities.R")

library(diagL1)
library(plgp)
library(ks)


arg<-commandArgs(TRUE)
idx <- as.integer(arg[1])
#data
#idx <- 1
for (n in c(100,200)){
for (t in c("g","t","k")){

#Load data set 
train <- read.csv(paste("data_train/data_train_",t,"_n_",n,"_",idx,".csv",sep=""),sep = ",",header = FALSE)
test <- read.csv(paste("data_test/data_test_",t,"_n_",n,"_",idx,".csv",sep=""),sep = ",",header = FALSE)
#Put variables in the required format
#Data
ytrain <- train[,1]
Xtrain <- train[,2]
data <- cbind(ytrain,Xtrain)
#test
ytest <- test[,1]
Xtest <- matrix(test[,2],ncol=1)
mXtest <- test[,3]
d <- dim(data)[2]


df <- data.frame(n=integer(),type=character(),data=integer(),method=character(),
                 env=numeric(),abs_bias=numeric(),rel_bias=numeric(),abs_rwd=numeric(),rel_rwd=numeric(),
                 abs_dev=numeric(),rel_dev=numeric(),MSE_med=numeric(),MSE_mean=numeric())
#Run
M <- 1000
B <- 400
for (k in c("gaussian","laplace","uniform")){
  for (d in c(0.5,0.7,0.9,1,1.1)){
    for (b in c("pi","scv")){
        h <- d*initialization_kernel_regr_resampling(data,M,bw=b)
        posterior <- ps_posterior_regr(data=data,grid=Xtest,h=h,M=M,B=B,kernel=k)
        method_name <- paste("kernel_",substr(k,1,1),"_",d,b,sep="")
        df <- rbind(df,data.frame(n=n,type=t,data=idx,method=method_name,env=envelope(posterior,mXtest),
                                  abs_bias=bias(posterior,mXtest)$abs,rel_bias=bias(posterior,mXtest)$rel,
                                  abs_rwd=relwid(posterior,mXtest)$abs,rel_rwd=relwid(posterior,mXtest)$rel,
                                  abs_dev=dev(posterior,mXtest)$abs,rel_dev=dev(posterior,mXtest)$rel,
                                  MSE_med=mean((posterior$q50-ytest)^2),MSE_mean=mean((posterior$mean-ytest)^2)))
}}}

#plot_ps_posterior_regr(grid=Xtest,posterior,mXtest,data)

name <- paste("out/out_kernel_",t,"_n_",n,"_",idx,".RData",sep="")
save(df,file=name)
}}



