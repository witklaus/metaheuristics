geneticAlgorithm = function(f, x_min, x_max, cel, popSize, pMut, maxIter)
{
  # geneticAlgorithm
  # INPUT
  #      - f objective function
  #      - x_min vector of the minimum values of coordinates
  #      - x_max vector of the maximum values of coordinates
  #      - cel coordinate encryption length 
  #      - popSize size of the population
  #      - pMut probability of single genome mutation
  #      - maxIter number of generations
  algorithmStart <- Sys.time()
  result <- list(x_opt = c(), f_opt = c(), x_hist= c(), f_hist= c(), f_mean = c(), time=c(), cumTime=c())
  
  # Check the number of dimensions
  Dim <- length(x_min) 
    
  # Initialize Population
  population <- matrix(NA, nrow = popSize, ncol = cel*Dim)
  for(i in 1 : popSize){
    population[i,] <- runif(cel*Dim)>1
  }
  coordinates <- getCoordinates(population, cel, x_min, x_max)
  
  # Calculate fittness of individuals
  objFunction <- rep(NA, popSize)
  for(i in 1 : popSize){
    objFunction[i] <- f(coordinates[i,])
  }
  
  # Assign the first population to output 
  result$x_opt <- coordinates[which.min(objFunction),]
  result$f_opt <- f(coordinates[which.min(objFunction),])
  
  # The generational loop
  finished <- FALSE
  currIter <- 1
  while(finished == FALSE){
    # Assign the output
    if(currIter <= maxIter){
      if(result$f_opt > f(coordinates[which.min(objFunction),])){
        result$x_opt <- coordinates[which.min(objFunction),]
        result$f_opt <- f(coordinates[which.min(objFunction),])
      }
      time <- Sys.time()
      cumTime <- time-algorithmStart
      result$time <- rbind(result$time,time)
      result$cumTime <- rbind(result$cumTime,cumTime)
      result$f_hist <- rbind(result$f_hist, result$f_opt) 
      result$x_hist <- rbind(result$x_hist, coordinates[which.min(objFunction),])
      result$f_mean <- rbind(result$f_mean, mean(objFunction)) 
    }else{
      finished <- TRUE
    }
    
    # Translate binary coding into real values  
    coordinates <- getCoordinates(population, cel, x_min, x_max)
    
    # Calculate fittness of the individuals
    objFunction <- rep(NA, popSize)
      for(i in 1 : popSize){
      objFunction[i] <- f(coordinates[i,])
    }
    rFitt <- min(objFunction)/objFunction # Relative Fittness
    nrFitt <- rFitt / sum(rFitt)          # Relative Normalized (sum up to 1) Fittness
    
    # Selection operator (Roulette wheel)
    selectedPool<- rep(0, popSize)
    for(i in 1 : popSize){
      selectedPool[i] <- which.min(runif(1)>cumsum(nrFitt))
    }
    
    # Crossover operator (for selected pool)
    nextGeneration <- matrix(NA, nrow = popSize, ncol = cel*Dim)
    for(i in 1 : popSize){
      parentId <- round(runif(1,1,popSize))
      cutId <- round(runif(1,1,Dim*cel-1)) # Please, do not exceed the matrix sizes
      # Create offspring
      nextGeneration[i, 1 : cutId] <- population[selectedPool[i], 1 : cutId]
      nextGeneration[i, (cutId + 1) : (Dim*cel)] <- population[selectedPool[parentId], (cutId + 1) : (Dim*cel)]
    }
    
    # Mutation operator
    for(i in 1 : popSize){
      genomeMutId <- which(runif(Dim*cel)>pMut) # Draw the genomes that will mutate
      for(j in 1 : length(genomeMutId)){
        nextGeneration[i, genomeMutId[j]] <- !nextGeneration[i, genomeMutId[j]] 
      }
    }
    
    # Replace the old population
    population <- nextGeneration
    currIter <- currIter + 1
  }
  return(result)
}
  
  
intbin = function(x){
  # Translate the binary coding to real values numbers
  return(sum(2^(which(rev(x==1))-1)))
}
  
getCoordinates = function(population, cel, x_min, x_max, pMut){
  # Transform the binary coding into coordinates
  coordinates <- matrix(NA, nrow = dim(population)[1], ncol = 2)
  for(i in 1 : dim(population)[1]){
    for(j in 1 : 2){
      coordinatesTemp <- intbin(population[i, seq(cel*(j-1)+1, j*cel)])
      coordinates[i,j] <- ((x_max[j]-x_min[j])/(2^cel-1))*coordinatesTemp+x_min[j]
    }
  }
  return(coordinates)
}