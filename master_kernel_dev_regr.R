rm(list=ls())
source("functions_kernel.R")
source("functions_kernel_regr.R")

library(diagL1)
library(plgp)
library(ks)


idx <- 1
#sample data
n <- 100
X0 <- matrix(seq(0, 5, length=n), ncol=1)
D <- distance(X0)
eps <- sqrt(.Machine$double.eps) 
Sigma <- exp(-D) + diag(eps, n)
mX  <- rmvnorm(1, sigma=Sigma)
y <- matrix(mX + rnorm(n,0,sd=0.5),ncol=1)
plot(X0, mX, type="l")
points(X0,y)



data <- cbind(y,X0)
d <- dim(data)[2]

#create grid 
grid <- matrix(0,nrow=length(seq(min(X0),max(X0),0.05)),ncol=(d-1))
for (i in 2:d){grid[,i-1] <- seq(min(X0),max(X0),0.05)} 





Hpi.diag(data)
Hscv.diag(data)


#Run
M <- 500
B <- 200
h <- initialization_kernel_regr_resampling(data,M,bw="scv")
posterior <- ps_posterior_regr(data=data,grid=grid,h=h,M=M,B=B,kernel="gaussian")
plot_ps_posterior_regr(grid,posterior,mX,data)


#h0=h[1]
#hist(y,breaks=50,freq=FALSE)
#lines(grid,laplace_kde(y,grid,h0))


#R functions to select bandwidth with Gaussian kernel
# library(MASS)
# library(ks)
# density(y,kernel="gaussian",bw="nrd0")$bw 
# density(y,kernel="gaussian",bw="ucv")$bw 
# 
# bw.nrd0(y) #silverman rule
# bw.nrd(y) #Scott's variation of Silverman rule
# bw.SJ(x)    # Sheather-Jones plug-in method
# bw.ucv(y) #unbiased CV


# plot(density(y,kernel="gaussian",bw=bw.ucv(y)))
# lines(grid,true)



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
name <- paste("output_py/out_python_n",n,"_",idx,".npy",sep="")
mat <- np$load(name)
mat <- np$matrix(exp(mat))/sd(y)


grid_fong <- seq(-4,3.95,0.05)*sd(y)+mean(y)#taken from experiments_plots notebook
plot_ps_fong(mat,grid,grid_fong,true)

             