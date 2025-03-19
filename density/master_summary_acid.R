

#define matrix
d<-data.frame(n=integer(),type=character(),data=integer(),method=character(),
              env=numeric(),abs_bias=numeric(),rel_bias=numeric(),abs_rwd=numeric(),rel_rwd=numeric(),
              abs_dev=numeric(),rel_dev=numeric(),MSE_med=numeric(),MSE_mean=numeric())
r<-data.frame(n=integer(),type=character(),data=integer(),method=character(),
              env=numeric(),abs_bias=numeric(),rel_bias=numeric(),abs_rwd=numeric(),rel_rwd=numeric(),
              abs_dev=numeric(),rel_dev=numeric(),MSE_med=numeric(),MSE_mean=numeric())

for (n in c(50,200)){
for (i in 1:100){
    file=paste("out/out_dpmm_n",n,"_",i,".RData",sep="")
    if (file.exists(file)){
      load(file)
      dp <-df
      d <- rbind(d,df)
    }
    file=paste("out/out_kernel_n",n,"_",i,".RData",sep="")
    if (file.exists(file)){
      load(file)
      dk <-df
      d <- rbind(d,df)
      r <- rbind(r,cbind(dk[,1:3],sweep(as.matrix(dk[,4:7]), 1, as.numeric(dp[4:7]), "/")))
    }
    file=paste("out/out_copula_n",n,"_",i,".RData",sep="")
    if (file.exists(file)){
      load(file)
      dc <-df
      d <- rbind(d,df)
      r <- rbind(r,c(dc[1:3],dc[4:7]/dp[4:7]))
    }
    
    
}} 

a <- list(d=d,r=r)
file=paste("summary_acid.rda",sep="")
save(a,file=file)


