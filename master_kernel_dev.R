rm(list=ls())
source("functions_kernel.R")
#sample data
n <- 200

u <- runif(n)
y <- rnorm(n,-2,1)
y1 <- rnorm(n,2,1)
y[u>0.8] <- y1[u>0.8]
#y <- (y-mean(y))/sd(y)


#grid
xmin=-8
xmax=8
ngrid=200
grid <- seq(xmin,xmax,(xmax-xmin)/ngrid)
true <- 0.8*dnorm(grid,-2,1)+0.2*dnorm(grid,2,1)



#Run
M <- 500
h <- initialization_kernel_resampling(y,M,bw="SJ")
posterior <- ps_posterior(y,grid,h,M,B=400,kernel="gaussian")
plot_ps_posterior(grid,posterior,true)





#R functions to select bandwidth with Gaussian kernel
library(MASS)
library(ks)
density(y,kernel="gaussian",bw="nrd0")$bw 
density(y,kernel="gaussian",bw="ucv")$bw 

bw.nrd0(y) #silverman rule
bw.nrd(y) #Scott's variation of Silverman rule
bw.SJ(x)    # Sheather-Jones plug-in method
bw.ucv(y) #unbiased CV


plot(density(y,kernel="gaussian",bw=bw.ucv(y)))
lines(grid,true)


