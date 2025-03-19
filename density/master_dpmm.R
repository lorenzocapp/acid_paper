library(dirichletprocess)
source("../utilities.R")
#The following code has been taken from https://github.com/edfong/MP 
# In particular the folder run_expts
#Small modifications has been added to conform to the current simulation study. 
arg<-commandArgs(TRUE)     
idx <- as.integer(arg[1])
#idx <- 1
for (n in c(50,200)){

#Load toy dataset
y <- scan(paste("data/data_n",n,"_",idx,".csv",sep=""),sep = ",")
true <- scan(paste("data_true_d/data_true_d_n",n,"_",idx,".csv",sep=""),sep = ",")
#Elicit prior and carry out MCMC
n_samp = 4000
k = 1
system.time(dp_gmm <- DirichletProcessGaussian(y,g0Priors = c(0,k,0.1,0.1))) #g0Priors = c(0,1,0.1,0.1)
system.time(dp_gmm <- Fit(dp_gmm, n_samp))

#Initialize y_plot
grid <- seq(-4,3.95,0.05)*sd(y)+mean(y)
dy = grid[2] - grid[1]

#Compute posterior samples of pdf
ind = 1:n_samp
system.time(pdf_samp_gmm <- t(sapply(ind, function(i) PosteriorFunction(dp_gmm,i)(grid))))
#discard burn-in
posterior <- pdf_samp_gmm[2001:n_samp,]

#Compute posterior mean and quantiles of pdf samples
q2.5 <- apply(data.frame(posterior), 2, quantile, probs=0.025)
q97.5 <- apply(data.frame(posterior), 2, quantile, probs=0.975)
q50 <- apply(data.frame(posterior), 2, quantile, probs=0.5)
mean <- apply(data.frame(posterior), 2, mean)
posterior <- list(grid=grid,mean=mean,q2.5=q2.5,q97.5=q97.5,q50=q50)


#plot_ps_posterior(grid,posterior,true,ylim=c(),xlim=c())


#evaluation_metric

df <- data.frame(n=n,type="v1",data=idx,method="dpmm",env=envelope(posterior,true),
           abs_bias=bias(posterior,true)$abs,rel_bias=bias(posterior,true)$rel,
           abs_rwd=relwid(posterior,true)$abs,rel_rwd=relwid(posterior,true)$rel,
           abs_dev=dev(posterior,true)$abs,rel_dev=dev(posterior,true)$rel,
           MSE_med=NA,MSE_mean=NA)

name <- paste("out/out_dpmm_n",n,"_",idx,".RData",sep="")
save(df,file=name)

}