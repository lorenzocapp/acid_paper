###############################################
######## FUNCTIONS KERNEL PredResampl ########
###############################################

#the final density estimation does not need to be done with the recursive kernel estimator


initialization_kernel_resampling <- function(y,M,bw="nrd"){

 'Determines a sequence of bandwidth
 - y: data
 - M: desidered number of synthetic samples'

 'The function goes as follows:
 1: bandwdith from 1 to n are chosen to be equal to the value given by the chosen criteria:
    - ucv: unbiased cross validation
    - nrd0: Silverman rule
    - nrd: Scott s variation of Silverman rule
 2: Sample 10 synthetic datasets of size M (10 is arbitrary)
 3: Find new bandwidth of the 10 new datasets (y+M synthetic obs) and take the average
 4: Interpolated the bandwidth from n to n+M using an exponential function
    - there are many choices at the moment I fit a function of the form a1 e^(-a2 n)
 5: define the sequence of bandwidth:
    - the first n are identical (not decaying equal to the quantity at begininng
    - then the n+1:M are decaying expoentially all the way to the precoomputed quantity'
  
  n <- length(y)
  bw_function <- match.fun(paste0("bw.", bw))
  #1.  
  b1 <- bw_function(y)
  #2 & 3
  b2 <- c()
  for (i in 1:10){
    yB <- y[sample(n,M,replace=T)]
    b2 <- c(b2,bw_function(c(y,yB)))
  }
  b2 <- mean(b2)
  #4
  a1 <- b1
  a2 <- log(b1/b2)/M
  #5
  h <- c(rep(b1,n),a1*exp(-a2*(seq(1,M))))
  return(h)
}


ps_kernel <- function(y,M,h,kernel="gaussian"){
  
  'Implement a predictive resampling with kernels (single run!)
  - y: data
  - M: # synthetic data
  - kernel type'
  
  n <- length(y)
  Y <- y
  if (kernel=="gaussian"){
    for (i in 1:M){
      j <- sample(n+i-1,1) 
      Y <- c(Y,rnorm(1,mean=Y[j],sd=h[j]))
    }
  } else if (kernel=="uniform"){
    for (i in 1:M){
      j <- sample(n+i-1,1) 
      Y <- c(Y,runif(1,min=Y[j]-h[j],max=Y[j]+h[j]))
    }
  } else if (kernel=="laplace"){
    #h <- 1.1 *h #inflate kernel for Laplace
    for (i in 1:M){
      j <- sample(n+i-1,1) 
      Y <- c(Y,rlaplace(1,location=Y[j],scale=h[j]))
    }
  } else {
    print("Kernel not yet implemented!")
  }
    boot <- Y[n+1:M]
  
  return(boot)
}

                          
recursive_kernel <- function(y,grid,h,kernel="gaussian",ps=TRUE,M=1000){
  
  'Implement the kernel recursive estimator, it can generate predictive resampling samples!'
  
  n <- length(y)
  if (ps==TRUE){
    yB <- ps_kernel(y,M,h,kernel)
    y  <- c(y,yB)
  }
  
  
  if (kernel=="gaussian"){
    f <- dnorm((grid-y[1])/h[1])/h[1]
    for (i in 2:length(y)){
      f <- (1-i^(-1))*f + i^(-1)/h[i]* dnorm((grid-y[i])/h[i])
      #Z <- (grid[2]-grid[1])/2*(2*sum(f)-f[1]-f[length(f)])
      #f <- f/Z
    }
  } else if (kernel=="uniform"){
    f <- dunif(grid,min=y[1]-h[1],max=y[1]+h[1])
    for (i in 2:length(y)){
      f <- (1-i^(-1))*f + i^(-1)*dunif(grid,min=y[i]-h[i],max=y[i]+h[i])
      #Z <- (grid[2]-grid[1])/2*(2*sum(f)-f[1]-f[length(f)])
      #f <- f/Z
    }
  } else if (kernel=="laplace"){
    #h <- 1.1*h #inflate the kernel for Laplace
    f <- dlaplace(grid,location=y[1],scale=h[1])
    for (i in 2:length(y)){
      f <- (1-i^(-1))*f + i^(-1)*dlaplace(grid,location=y[i],scale=h[i])
    }
    #we adjust to make it sum to one at the end
    #Z <- (grid[2]-grid[1])/2*(2*sum(f)-f[1]-f[length(f)])
    #f <- f/Z
  } else {
    print("Kernel not yet implemented!")
  }  
  return(f)
}



standard_kernel <- function(y,grid,h,kernel="gaussian",bw="nrd",ps=TRUE,M=1000){
  
  'Implement the kernel standard estimator (unless it is Gaussian, since it already exists): it can generate predictive resampling samples!'
  
  #1.Predictive resampling
  n <- length(y)
  if (ps==TRUE){
    yB <- ps_kernel(y,M,h,kernel)
    y  <- c(y,yB)
  }
  #2. Determines optimal bw at the end
  bw_function <- match.fun(paste0("bw.", bw))
  h0 <- bw_function(y)
  
  #3. Do kde with optimal bw and data+synthetic data
  if (kernel=="gaussian"){
    f <- gaussian_kde(y,grid,h0=h0)
  } else if (kernel=="uniform"){
    f <- uniform_kde(y,grid,h0=h0)
  } else if (kernel=="laplace"){
    f <- laplace_kde(y,grid,h0=h0)
  } else {
    print("Kernel not yet implemented!")
  }  
  return(f)
}



laplace_kde <- function(y,grid,h0) {
  
  "Implements a standard Kde"
  
  # laplace kernel function
  f <- dlaplace(grid,location=y[1],scale=h0)
  for (i in 2:length(y)){
    f <- (1-i^(-1))*f + i^(-1)*dlaplace(grid,location=y[i],scale=h0)
  }
  #adjust to make it sum to one at the end
  Z <- (grid[2]-grid[1])/2*(2*sum(f)-f[1]-f[length(f)])
  f <- f/Z

  return(f)
}

uniform_kde <- function(y,grid,h0) {
  
  "Implements a standard Kde"
  
  # uniform kernel function
  f <- dunif(grid,min=y[1]-h0,max=y[1]+h0)
  for (i in 2:length(y)){
    f <- (1-i^(-1))*f + i^(-1)*dunif(grid,min=y[i]-h0,max=y[i]+h0)
    Z <- (grid[2]-grid[1])/2*(2*sum(f)-f[1]-f[length(f)])
    f <- f/Z
  }
  return(f)
}


gaussian_kde <- function(y,grid,h0) {
  
  "Implements a standard Kde"
  
  # Gaussian kernel function
  f <- dnorm((grid-y[1])/h0)/h0
  for (i in 2:length(y)){
    f <- (1-i^(-1))*f + i^(-1)/h0* dnorm((grid-y[i])/h0)
    Z <- (grid[2]-grid[1])/2*(2*sum(f)-f[1]-f[length(f)])
    f <- f/Z
  }
  return(f)
}

ps_posterior <- function(y,grid,h,M=1000,B=100,kernel="gaussian",algo="recursive",bw="ndr"){
  
  "Implement Predictive resampling posterior
  Inputs:
  - algo: recursive it should use the recursive kernel estimator, standard if uses the standard kernel density estimatr
  Returns:
  - mean
  - median
  - 97.5% and 2.5% bands"
  
  posterior <- matrix(0,ncol=length(grid),nrow=B)
  if (algo=="recursive") {
    for (b in 1:B){
      posterior[b,] <- recursive_kernel(y,grid,h,kernel,ps=TRUE,M)
    }
  } else if (algo=="standard"){
    for (b in 1:B){
      posterior[b,] <- standard_kernel(y,grid,h,kernel,M,bw,ps=TRUE,)
    }
  }
  
  
  q2.5 <- apply(data.frame(posterior), 2, quantile, probs=0.025)
  q97.5 <- apply(data.frame(posterior), 2, quantile, probs=0.975)
  q50 <- apply(data.frame(posterior), 2, quantile, probs=0.5)
  mean <- apply(data.frame(posterior), 2, mean)
  posterior <- list(grid=grid,mean=mean,q2.5=q2.5,q97.5=q97.5,q50=q50)
  return(posterior)
}






plot_ps_posterior <- function(grid,posterior,true,ylim=c(),xlim=c()){

  if (length(ylim)==0){
    ylim = c(0,max(c(posterior$q97.5,true)))
  }
  if (length(xlim)==0){
    xlim <- c(min(grid),max(grid))
  }
  plot(grid,posterior$mean,type="l",ylim=ylim,xlim=xlim,lwd=2)
  polygon(c(grid, rev(grid)), c(posterior$q2.5, rev(posterior$q97.5)),
          col = "#D3D3D3")
  lines(grid,posterior$mean,lwd=2)
  lines(grid,true, lty=2,lwd=2)
  lines(grid,posterior$q2.5,col="black")
  lines(grid,posterior$q97.5,col="black")
}


plot_ps_fong <- function(mat,grid,grid_fong,true,ylim=c(),xlim=c()){
  
  q2.5 <- apply(data.frame(mat), 2, quantile, probs=0.025)
  q97.5 <- apply(data.frame(mat), 2, quantile, probs=0.975)
  q50 <- apply(data.frame(mat), 2, quantile, probs=0.5)
  mean <- apply(data.frame(mat), 2, mean)
  posterior <- list(grid_fong=grid_fong,mean=mean,q2.5=q2.5,q97.5=q97.5,q50=q50)
  if (length(ylim)==0){
    ylim = c(0,max(c(posterior$q97.5,true)))
  }
  if (length(xlim)==0){
    xlim <- c(min(grid_fong),max(grid_fong))
  }
  plot(grid_fong,posterior$mean,type="l",ylim=ylim,xlim=xlim,lwd=2)
  polygon(c(grid_fong, rev(grid_fong)), c(posterior$q2.5, rev(posterior$q97.5)),
          col = "#D3D3D3")
  lines(grid_fong,posterior$mean,lwd=2)
  lines(grid,true, lty=2,lwd=2)
  lines(grid_fong,posterior$q2.5,col="black")
  lines(grid_fong,posterior$q97.5,col="black")
}




##################################
######## UNUSED (FOR NOW) ########
##################################

#Not sure if necessary for now
recursive_kernel_one <- function(f,grid,y,h,i,kernel="gaussian"){
  "it does one-step-ahead kernel update"
  '- f: last density estimate
   - grid: x
   - y,h,i: data, bandwidth corresponding to the i-th update'
  if (kernel=="gaussian"){
    f <- (1-i^(-1))*f + i^(-1)/h* dnorm((grid-y)/h)
    #Renormalize for numerical prediction
    Z <- (grid[2]-grid[1])/2*(2*sum(f)-f[1]-f[length(f)])
    f <- f/Z
  }
  return(f)
}

