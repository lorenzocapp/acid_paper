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
  for (i in 1:M){
    j <- sample(n+i-1,1) #determine from which kernel to sample
    Y <- c(Y,rnorm(1,mean=Y[j],sd=h[j]))
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
    Z <- (grid[2]-grid[1])/2*(2*sum(f)-f[1]-f[length(f)])
    f <- f/Z
  }}
  return(f)
}


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

ps_posterior <- function(y,grid,h,M=1000,B=100,kernel="gaussian"){
  
  "Implement Predictive resampling posterior
  Returns:
  - mean
  - median
  - 97.5% and 2.5% bands"
  
  posterior <- matrix(0,ncol=length(grid),nrow=B)
  for (b in 1:B){
    posterior[b,] <- recursive_kernel(y,grid,h,kernel="gaussian",ps=TRUE,M)
  }
  
  q2.5 <- apply(data.frame(posterior), 2, quantile, probs=0.025)
  q97.5 <- apply(data.frame(posterior), 2, quantile, probs=0.975)
  q50 <- apply(data.frame(posterior), 2, quantile, probs=0.5)
  mean <- apply(data.frame(posterior), 2, mean)
  posterior <- list(grid=grid,mean=mean,q2.5=q2.5,q97.5=q97.5,q50=q50)
  return(posterior)
}


plot_ps_posterior <- function(grid,posterior,true){

  ymax <- max(c(posterior$q97.5,true))
  plot(grid,posterior$mean,type="l",ylim=c(0,ymax),lwd=2)
  polygon(c(grid, rev(grid)), c(posterior$q2.5, rev(posterior$q97.5)),
          col = "#D3D3D3")
  lines(grid,posterior$mean,lwd=2)
  lines(grid,true, lty=2,lwd=2)
  lines(grid,posterior$q2.5,col="black")
  lines(grid,posterior$q97.5,col="black")
}
