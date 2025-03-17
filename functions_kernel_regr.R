##########################################################
######## FUNCTIONS KERNEL  REGRESSION PredResampl ########
#########################################################



initialization_kernel_regr_resampling <- function(data,M,bw="pi"){
  
  
  'Determines a sequence of bandwidth
 - data: the first column is y,then the covariates. It is n x (p+1)
 - M: desidered number of synthetic samples'
  
  'The function goes as follows:
 1: bandwdith from 1 to n are chosen to be equal to the value given by the chosen criteria:
    - scv: SCV selector of Jones Marron Park (1991)
    - pi : plug in method of Wand and Jones (1994)
 2- 5 proceeds as the univariate function'
  
  n <- dim(data)[1]
  bw_function <- match.fun(paste0("H", bw,".diag"))
  #1.  
  b1 <- diag(bw_function(data))
  #2 & 3
  b2 <- c()
  for (i in 1:10){
    dataB <- data[sample(n,M,replace=T),]
    b2 <- cbind(b2,diag(bw_function(rbind(data,dataB))))
  }
  b2 <- rowMeans(b2)
  #4
  a1 <- b1
  a2 <- log(b1/b2)/M
  #5
  h <- bandwidth_matrix(a1,a2,M)
  return(h)
}


bandwidth_matrix <- function(a1,a2,M){
  
  "Generate a bandwidth sequence in the multivariate case"
  
  h0 <- matrix(rep(a1,n),nrow=n,byrow=T)
  h1 <- exp(-outer(seq(1,M),a2))* matrix(a1,nrow=M,ncol=dim(h0)[2],byrow=TRUE)
  return(rbind(h0,h1))
}




ps_posterior_regr <- function(data,grid,h,M=1000,B=100,kernel="gaussian",algo="recursive",bw="ndr"){
  
  "Implement Predictive resampling posterior for a regression function
  Inputs:
  - algo: recursive it should use the recursive kernel estimator, standard if uses the standard kernel density estimatr
  Returns:
  - mean
  - median
  - 97.5% and 2.5% bands"
  
  #it is still 1d
  posterior <- matrix(0,ncol=nrow(grid),nrow=B)
  if (algo=="recursive") {
    for (b in 1:B){
      posterior[b,] <- recursive_kernel_regr(data,grid,h,kernel,ps=TRUE,M)
    }
  } else if (algo=="standard"){
    for (b in 1:B){               
      posterior[b,] <- standard_kernel(y,grid,h,kernel,ps=TRUE,M,bw)
    }
  } else {
    print("which estimator you want to use?")
  }
  
  
  q2.5 <- apply(data.frame(posterior), 2, quantile, probs=0.025)
  q97.5 <- apply(data.frame(posterior), 2, quantile, probs=0.975)
  q50 <- apply(data.frame(posterior), 2, quantile, probs=0.5)
  mean <- apply(data.frame(posterior), 2, mean)
  posterior <- list(grid=grid,mean=mean,q2.5=q2.5,q97.5=q97.5,q50=q50)
  return(posterior)
}



recursive_kernel_regr <- function(data,grid,h,kernel="gaussian",ps=TRUE,M=1000){
  
  'Implement the kernel recursive regression estimator, 
  it can generate predictive resampling samples!'
  
  n <- dim(data)[1]
  d <- dim(data)[2]
  if (ps==TRUE){
    dataB <- ps_kernel_mv(data,M,h,kernel)
    data  <- rbind(data,dataB)
  }
  
  y <- data[,1]
  X <- matrix(data[,2:d],ncol=(d-1)) #Check when it is multivariate!
  
  if (kernel=="gaussian"){
    m <- gaussian_kernel_regression(X, y, h, grid) 
  } else if (kernel=="uniform"){
    f <- dunif(grid,min=y[1]-h[1],max=y[1]+h[1])
    for (i in 2:length(y)){
      #f <- (1-i^(-1))*f + i^(-1)*dunif(grid,min=y[i]-h[i],max=y[i]+h[i])
      f <- f + dunif(grid,min=y[i]-h[i],max=y[i]+h[i])
      #Z <- (grid[2]-grid[1])/2*(2*sum(f)-f[1]-f[length(f)])
      #f <- f/Z
    }
  } else if (kernel=="laplace"){
    #h <- 1.1*h #inflate the kernel for Laplace
    #f <- dlaplace(grid,location=y[1],scale=h[1])
    f <- exp(-abs(grid - y[1]) / h[1])/(2 * h[1])
    for (i in 2:length(y)){
      #f <- (1-i^(-1))*f + i^(-1)*dlaplace(grid,location=y[i],scale=h[i])
      #f <- f + dlaplace(grid,location=y[i],scale=h[i])
      f <- f + exp(-abs(grid - y[i]) / h[i])/(2 * h[i])
    }
    #we adjust to make it sum to one at the end
    #Z <- (grid[2]-grid[1])/2*(2*sum(f)-f[1]-f[length(f)])
    #f <- f/Z
  } else {
    print("Kernel not yet implemented!")
  }  
  return(m)
}


ps_kernel_mv <- function(data,M,h,kernel="gaussian"){
  
  'Implement a predictive resampling with kernels multivariate (single run!)
  - data: data
  - M: # synthetic data
  - kernel type'
  
  n <- dim(data)[1]
  d <- dim(data)[2]
  Y <- data
  if (kernel=="gaussian"){
    for (i in 1:M){
      j <- sample(n+i-1,1) 
      Y <- rbind(Y,rmvnorm(1,mean=Y[j,], sigma= diag(h[j,]))) 
    }
  } else if (kernel=="uniform"){
    for (i in 1:M){
      j <- sample(n+i-1,1) 
      Y <- rbind(Y,runif(d,min=Y[j,]-h[j,],max=Y[j,]+h[j,]))
    }
  } else if (kernel=="laplace"){
    for (i in 1:M){
      j <- sample(n+i-1,1) 
      Y <- rbind(Y,mapply(rlaplace, n = 1, location=Y[j,],scale=h[j,]))
    }
  } else {
    print("Kernel not yet implemented!")
  }
  boot <- Y[n+1:M,]
  
  return(boot)
}


gaussian_kernel_regression <- function(X, y, h, grid) {
  
  # X: matrix of predictors (n x p)
  # y: vector of responses (n)
  # h: matrix of bandwidth
  # x_new: new data points to predict (m x p)
  
  # Gaussian Kernel Function
  kernel <- function(x, x_i, h) {
    exp(-sum((x - x_i)^2 / (2 * h^2)))
  }
  
  # Predict for each grid_i (each x in which I am evaluating the grid)
  m <- sapply(1:nrow(grid), function(i) {
    grid_i <- grid[i, ]
    
    # Compute kernel weights 
    weights <- sapply(1:nrow(X), function(j) kernel(grid_i, X[j, ], h[j,]))
    
    # Weighted sum of responses
    weighted_sum <- sum(weights * y) / sum(weights)
    return(weighted_sum)
  })
  
  return(m)
}


plot_ps_posterior_regr <- function(grid,posterior,mX,data,ylim=c(),xlim=c()){
  
  y <- data[,1]
  X0 <- data[,2:ncol(data)]
  if (length(ylim)==0){
    ylim = c(min(c(posterior$q2.5,mX)),max(c(posterior$q97.5,mX)))
  }
  if (length(xlim)==0){
    xlim <- c(min(grid),max(grid))
  }
  plot(grid,posterior$mean,type="l",ylim=ylim,xlim=xlim,lwd=2)
  polygon(c(grid, rev(grid)), c(posterior$q2.5, rev(posterior$q97.5)),
          col = "#D3D3D3")
  lines(grid,posterior$mean,lwd=2)
  lines(X0,mX, lty=2,lwd=2)
  lines(grid,posterior$q2.5,col="black")
  lines(grid,posterior$q97.5,col="black")
  points(X0,y)
}


