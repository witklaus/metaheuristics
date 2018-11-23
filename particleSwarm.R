#Algorytm Particle Swarm Optimization
  particleSwarm <- function(f,v=c(0.5,2.05,2.05),swarmSize=20,maxIter=1000,dim=2){
  #parametry funkcji particleSwarm:
    #f - funkcja celu
    #v[1] - omega (wspolczynnik bezwladnosci - odpowiada za inercje)
    #v[2] - phi(p) - waga swiadomosci
    #v[3] - phi(g) - waga myslenia spolecznego
    #swarmSize - liczba czasteczek w roju
    #maxIter - liczba iteracji
    #dim - wymiar
    #oznaczenie czasu rozpoczecia algorytmu
  algorithmStart <- Sys.time()
    #zwracanie listy wynikowej
  result <- list(x_opt=c(),f_opt=c(),x_hist=c(),f_hist=c(),time=c(),cumTime=c())
    #generowanie poczatkowych pozycji czasteczek
  position <- matrix(NA,nrow=swarmSize,ncol=dim)
  for(i in 1:swarmSize){
    for(j in 1:dim){
      position[i,j]=runif(1,-20,20)
    }
  }
    #generowanie poczatkowych wartosci wektora predkoscci
  velocity <- matrix(NA,nrow=swarmSize,ncol=dim)
  for(i in 1:swarmSize){
    for(j in 1:dim){
      velocity[i,j]=runif(1,0,1)
    }
  }
    #ocena rozwiazan poczatkowych
  p_best <- position
  funApp <- apply(position,1,myFun)
  g_best <- p_best[which.min(funApp),]
  fg_best <- f(g_best)
   #aktualizacja najlepszego rozwiazania globalnego g_best
  result$x_opt <- g_best
    #glowna petla algorytmu
  for(currIter in 1:(maxIter)){
    for(i in 1:swarmSize){
      for(j in 1:dim){
        #r_p - zaburzenie na poziomie lokalnym
        #r_g - zaburzenie na poziomie globalnym
        r_p <- runif(1,0,1)
        r_g <- runif(1,0,1)
        #uaktualnienie wektora prêdkoœci
        velocity[i,j] <- v[1]*velocity[i,j]+v[2]*r_p*(p_best[i,j]-position[i,j])+v[3]*r_g*(g_best[j]-position[i,j])
        #uaktualnienie wektora pozycji czasteczek
        position[i,j] <- position[i,j]+velocity[i,j]
      }
        #ocena nowego rozwiazania na poziomie lokalnym
      if(f(position[i,])<f(p_best[i,])){
        p_best[i,] <- position[i,]
        #ocena nowego rozwiazania na poziomie globalnym
        if(f(p_best[i,])<fg_best){
          g_best <- p_best[i,]
          fg_best <- f(g_best)
        }
      }
    }
        #sprawdzenie warunku stopu
    if(currIter<maxIter){
        #aktualizacja optymalnej wartosci funkcji celu
      if(f(g_best)<f(result$x_opt)){
        result$x_opt <- g_best
        result$f_opt <- fg_best
      }
      #przygotowanie listy wynikowej z wartosciami historycznymi i czasem iteracji (skumulowanym)
      time <- Sys.time()
      cumTime <- time-algorithmStart
      result$time <- rbind(result$time,time)
      result$cumTime <- rbind(result$cumTime,cumTime)
      result$x_hist <- rbind(result$x_hist,g_best)
      result$f_hist <- rbind(result$f_hist,fg_best)
    }
    #inkrementacja iteracji
    currIter <- currIter+1
  }
  #zwracanie listy wynikowej
  return(result)
}