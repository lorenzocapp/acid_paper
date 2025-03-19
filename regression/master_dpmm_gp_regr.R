
rm(list=ls())
library(dirichletprocess)
library(laGP)
source("../utilities.R")


arg<-commandArgs(TRUE)     
idx <- as.integer(arg[1])
#data
#idx <- 1
for (n in c(100,200)){
for (t in c("g","k","t")){

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


########################
######## DP ############
########################
z <- cbind(Xtrain,ytrain)


# Sample with MCMC
n_samples = 4000 #3000 for timing
g0Priors <- list(mu0 = rep_len(0, length.out = ncol(z)), 
                 Lambda = diag(ncol(z)), kappa0 = ncol(z), nu = ncol(z))

dp <- DirichletProcessMvnormal(z,g0Priors = g0Priors)
system.time(dp <- Fit(dp, n_samples))




n_plot = 100
y_grid = seq(min(c(ytrain,ytest)),max(c(ytrain,ytest)),length.out = n_plot)
dy = y_grid[2] - y_grid[1]

# Preallocate a matrix for efficiency
reg_means <- matrix(NA, nrow = 2000, ncol = length(Xtest)) 

for (iter in 2001:4000) {
  dp_sample <- PosteriorFunction(dp, iter)  # Get function for specific posterior sample
  
  # Create matrix of test points (Xtest, y_grid)
  xy_pairs <- cbind(rep(Xtest, each = n_plot), rep(y_grid, times = length(Xtest)))
  
  # Evaluate f(x, y) for all test pairs
  f_xy <- matrix(dp_sample(as.matrix(xy_pairs)), nrow = length(Xtest), byrow = TRUE)
  
  # Compute marginal f(x) by integrating over y
  f_x <- rowSums(f_xy) * dy
  
  # Compute conditional density f(y | x)
  f_y_given_x <- sweep(f_xy, 1, f_x, "/")
  
  # Compute conditional expectation E[Y | X=x]
  E_Y_given_X <- rowSums(f_y_given_x * matrix(y_grid, nrow = length(Xtest), ncol = length(y_grid), byrow = TRUE)) * dy
  
  # Store results efficiently
  reg_means[iter - 2000, ] <- E_Y_given_X
}


q2.5 <- apply(data.frame(reg_means), 2, quantile, probs=0.025)
q97.5 <- apply(data.frame(reg_means), 2, quantile, probs=0.975)
q50 <- apply(data.frame(reg_means), 2, quantile, probs=0.5)
mean <- apply(data.frame(reg_means), 2, mean)
posterior <- list(grid=Xtest,mean=mean,q2.5=q2.5,q97.5=q97.5,q50=q50)


# Plot with uncertainty
##plot(x_values, mean_reg, type = "l", col = "blue", ylim = range(lower, upper))
#polygon(c(x_values, rev(x_values)), c(lower, rev(upper)), col = rgb(0, 0, 1, 0.2), border = NA)
df <- data.frame(n=n,type=t,data=idx,method="dpmm",env=envelope(posterior,mXtest),
                          abs_bias=bias(posterior,mXtest)$abs,rel_bias=bias(posterior,mXtest)$rel,
                          abs_rwd=relwid(posterior,mXtest)$abs,rel_rwd=relwid(posterior,mXtest)$rel,
                          abs_dev=dev(posterior,mXtest)$abs,rel_dev=dev(posterior,mXtest)$rel,
                          MSE_med=mean((posterior$q50-ytest)^2),MSE_mean=mean((posterior$mean-ytest)^2))
name <- paste("out/out_dpmm_",t,"_n_",n,"_",idx,".RData",sep="")
save(df,file=name)


#########################
####  GP process @#######
#########################

###standard GP
# Define training data
Xtrain <- matrix(Xtrain,ncol=1)

# Define test points for the regression function
Xtest <- matrix(Xtest, ncol = 1)

eps <- sqrt(.Machine$double.eps) 

gp_model <- newGPsep(Xtrain, ytrain, d = 0.1, g = 1e-6,dK=TRUE)  # d = lengthscale, g = nugget
mle <- mleGPsep(gp_model, param="both", tmin=c(eps, eps), tmax=c(10, var(ytrain)))
# Compute posterior mean and variance (removing noise variance)
mean_reg <- predGPsep(gp_model, Xtest, lite = TRUE)$mean  # Mean function E[Y | X]
sd_reg <- sqrt(predGPsep(gp_model, Xtest, lite = TRUE)$s2)  # Posterior std dev (no noise)


q2.5 <- mean_reg - 1.96 * sd_reg
q97.5 <- mean_reg + 1.96 * sd_reg
mean <- mean_reg 
posterior <- list(grid=Xtest,mean=mean,q2.5=q2.5,q97.5=q97.5,q50=mean)
#plot_ps_posterior_regr(grid=Xtest,posterior,mXtest,data)


df <- data.frame(n=n,type=t,data=idx,method="gp",env=envelope(posterior,mXtest),
                 abs_bias=bias(posterior,mXtest)$abs,rel_bias=bias(posterior,mXtest)$rel,
                 abs_rwd=relwid(posterior,mXtest)$abs,rel_rwd=relwid(posterior,mXtest)$rel,
                 abs_dev=dev(posterior,mXtest)$abs,rel_dev=dev(posterior,mXtest)$rel,
                 MSE_med=mean((posterior$mean-ytest)^2),MSE_mean=mean((posterior$mean-ytest)^2))
name <- paste("out/out_gp_",t,"_n_",n,"_",idx,".RData",sep="")
save(df,file=name)

# Clean up model
deleteGPsep(gp_model)


}}


          