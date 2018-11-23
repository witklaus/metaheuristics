simulatedAnnealing = function(f, x, alpha, t, delta, maxIter)
{# Simulated Annealin Algorithm
  # f - objective function
  # x - inital solution
  # alpha - annealing schedule parameter
  # t - inital temperature
  # delta - neighbourhood radius
  # maxIt - maximum no. of iterations
  algorithmStart <- Sys.time()
  result = list(x_opt = x, f_opt = f(x), x_hist = x, f_hist = f(x), temperature = t, time = c(), cumTime = c())

  currIter <- 1
  finished <- FALSE
  x_s <- x
  
  while(finished == FALSE){
    # x_c - candidat sol. drawn uniformly fron N(x)
    u = runif(length(x_s))
    x_c = x_s + (-delta + 2 * delta * u)
    
    # A - Metropolis activation function
    A = min(1, exp(- (f(x_c) - f(x_s)) / t))
    
    # transition to candidate solution
    if (runif(1) < A)
    {
      x_s <- x_c
    }
    
    # temperature update
    t = alpha * t

    if(currIter<maxIter){
      if(f(x_s)<f(result$x_opt)){
        result$x_opt <- x_s
        result$f_opt <- f(x_s)
      }
      time <- Sys.time()
      cumTime <- time-algorithmStart
      result$time <- rbind(result$time,time)
      result$cumTime <- rbind(result$cumTime,cumTime)
      result$x_hist       <-rbind(result$x_hist, x_s)
      result$f_hist       <- rbind(result$f_hist, f(x_s))
      result$temperature  <- rbind(result$temperature, t)
      result$transProb    <- rbind(result$transProb, A)
    }else{
      finished            <- TRUE
    }
    currIter <- currIter + 1
  }
  return(result)
}