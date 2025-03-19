generate_mixture_data <- function(n) {
   # For reproducibility
  num_components <- sample(2:3, 1)  # Randomly choose 2 or 3 components
  
  repeat {
    means <- runif(num_components, -3, 3)  # Generate means in (-3,3)
    if (min(dist(means)) > 1.5) break  # Ensure separation of at least 1.5
  }
  
  sds <- runif(num_components, 0.5, 2)  # Random variances in a reasonable range
  
  kernel_types <- sample(c("gaussian", "student"), num_components, replace = TRUE)
  
  component_probs <- runif(num_components)
  component_probs <- component_probs / sum(component_probs)  # Normalize to sum to 1
  
  y <- numeric(n)
  
  for (i in seq_len(n)) {
    comp <- sample(seq_len(num_components), 1, prob = component_probs)
    if (kernel_types[comp] == "gaussian") {
      y[i] <- rnorm(1, mean = means[comp], sd = sds[comp])
    } else {  # Student's t-distribution
      df <- 10 # Random degrees of freedom between 3 and 10
      y[i] <- means[comp] + sds[comp] * rt(1, df = df)
    }
  }
  
  grid <- seq(-4, 3.95, 0.05) * sd(y) + mean(y)
  true_density <- numeric(length(grid))
  
  for (j in seq_along(grid)) {
    true_density[j] <- sum(component_probs * ifelse(kernel_types == "gaussian",
                                                    dnorm(grid[j], mean = means, sd = sds),
                                                    dt((grid[j] - means) / sds, df = 5) / sds))
  }
  
  return(list(y = y, grid = grid, true_density = true_density))
}

# Example usage:
set.seed(123)
n <- 200
for (i in 1:100){
  result <- generate_mixture_data(n)
  y <- result$y
  true <- result$true_density
  name <- paste("data/data_n",n,"_",i,".csv",sep="")
  write.table(y, name, sep = ",", row.names = FALSE, col.names = FALSE)
  name <- paste("data_true_d/data_true_d_n",n,"_",i,".csv",sep="")
  write.table(true, name, sep = ",", row.names = FALSE, col.names = FALSE)
}



#hist(y, breaks = 50, col = "lightblue", main = "Generated Mixture Data", freq = FALSE)
#lines(grid, true, col = "red", lwd = 2)
