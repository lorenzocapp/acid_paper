
rm(list=ls())
library(plgp)
generate_regression_data <- function(n,type){


    #X0 <- matrix(seq(0, 5, length=n), ncol=1)
    Xtrain <- matrix(sort(5*runif(n)), ncol=1)
    Xtest <-  matrix(seq(0,4.95,0.05),ncol=1)
    X <- matrix(sort(rbind(Xtrain,Xtest)),ncol=1)
    D <- distance(X)
    eps <- sqrt(.Machine$double.eps) 
    if (type=="g"){
      Sigma <- exp(-D) + diag(eps, nrow(X))
      mX  <- rmvnorm(1, sigma=Sigma)
      y <- matrix(mX + rnorm(nrow(X)),ncol=1)
    } else if (type=="t"){
      Sigma <- exp(-D) + diag(eps, nrow(X))
      mX  <- rmvnorm(1, sigma=Sigma)
      y <- matrix(mX + rt(nrow(X),df=5),ncol=1)
    } else if (type=="k"){
      D1 <- sqrt(D)/20
      Sigma <- 0.8*exp(-D)+0.2*exp(-D1) + diag(eps, nrow(X))
      mX  <- rmvnorm(1, sigma=Sigma)
      y <- matrix(mX + rnorm(nrow(X)),ncol=1)
    }
    idtrain <- X %in% Xtrain
    ytrain <- y[idtrain]
    ytest <- y[!idtrain]
    mXtest <- mX[!idtrain]
    #plot(X, mX, type="l")
    #points(X,y)
    train <- data.frame(ytrain,Xtrain)
    test <- data.frame(ytest,Xtest,mXtest)
    return(list(train=train,test=test))
}

# Example usage:
set.seed(123)
n <- 200
type <- "g"
for (i in 1:100){
  result <- generate_regression_data(n,type)
  train <- result$train
  test <- result$test
  name <- paste("data_train/data_train_",type,"_n_",n,"_",i,".csv",sep="")
  write.table(train, name, sep = ",", row.names = FALSE, col.names = FALSE)
  name <- paste("data_test/data_test_",type,"_n_",n,"_",i,".csv",sep="")
  write.table(test, name, sep = ",", row.names = FALSE, col.names = FALSE)
}



#hist(y, breaks = 50, col = "lightblue", main = "Generated Mixture Data", freq = FALSE)
#lines(grid, true, col = "red", lwd = 2)

lines(test$Xtest,test$mXtest,col="red")
