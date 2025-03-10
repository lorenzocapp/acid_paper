rm(list=ls())
source("functions_kernel.R")

library(diagL1)

idx <- 1
#sample data
n <- 150

u <- runif(n)
y <- rnorm(n,-2,1)
y1 <- rnorm(n,2,1)
y[u>0.8] <- y1[u>0.8]
#y <- (y-mean(y))/sd(y)
name <- paste("data/data_",idx,".csv",sep="")
write.csv(y,"data/data.csv",row.names = FALSE)
M <- 500
B <- 400
characteristics <- data.frame(M=M,B=B)


#grid
#xmin=-8
#xmax=8
#ngrid=200
grid <- seq(-4,3.95,0.05)*sd(y)+mean(y) #it is an automatic grid, so the Gaussian copula is identical
true <- 0.8*dnorm(grid,-2,1)+0.2*dnorm(grid,2,1)



#Run
M <- 500
h <- initialization_kernel_resampling(y,M,bw="SJ")
posterior <- ps_posterior(y=y,grid=grid,h=h,M=M,B=B,kernel="laplace")
plot_ps_posterior(grid,posterior,true)

#h0=h[1]
#hist(y,breaks=50,freq=FALSE)
#lines(grid,laplace_kde(y,grid,h0))


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



###Python practice
library(reticulate)
#py_require("pr_copula")
#py_require("numpy")
#py_install("pr_copula")
#py_require("jax")#this runs and check which packages are necessary
#source_python("python_codes/1_univariate_copula_modified.py")
#source_python("python_codes/bivariate_copula.py")
#source_python("python_codes/copula_density_functions.py")
#config <- reticulate::py_config()
#system2(config$python, c("-m", "pip", "install", "--quiet", shQuote("./")))

source_python("master_python.py")

y <- as.matrix(y)
check <- run_copula_algo(y,100,100,n)

np <- import("numpy")
mat <- np$load("output/out_python_1.npy")
mat <- np$matrix(exp(mat))/sd(y)


grid_fong <- seq(-4,3.95,0.05)*sd(y)+mean(y)#taken from experiments_plots notebook
plot_ps_fong(mat,grid,grid_fong,true,ylim=c(0,0.4),xlim=c(-6,6))
