#Nieklasyczne metody optymalizacji - case study: Particle Swarm Optimization

#Czyszczenie srodowiska
rm(list=ls())

#Sourcowanie funkcji algorytmow
source('numGradient.R')
source('steepestDescent.R')
source('simulatedAnnealing.R')
source('geneticAlgorithm.R')
source('particleSwarm.R')

#Funkcja celu - funkcja Schaffera:
myFun = function(x){
  return( 0.6 + ((sin(x[1]^2-x[2]^2))^2-0.5)/((1+0.001*(x[1]^2+x[2]^2))^2) )
}

#Parametry
xSeed <- c(3,4)
n_grid  <- 100
maxIter <- 500
set.seed(10)

#Algorytmy
sd <- steepestDescent(myFun, xSeed, 0.01, 10^-13, maxIter)
sa <- simulatedAnnealing(myFun, xSeed, 0.99, 10000, 0.2, maxIter)
ga <- geneticAlgorithm(myFun, c(-20, -20), c(20, 20), cel=50, popSize = 30, maxIter, pMut = 0.05)
pso <- particleSwarm(myFun,v=c(0.5,2.05,2.05),swarmSize=20,maxIter+1,dim=2)

#Zbieznosc funkcji celu dla algorytmow
plot(sd$f_hist[1:maxIter],col='purple',type='l',lwd=1,xlim=c(0,30),ylim=c(0,1.5),xlab="Numer iteracji",ylab="Wartosc f-cji celu")
lines(sa$f_hist[1:maxIter],col='blue',type='l', lwd=1)
lines(ga$f_hist[1:maxIter],col='green',type ='l',lwd=1)
lines(pso$f_hist[1:maxIter],col='red',type='l',lwd=1)

#Wykres sciezki optymalizacji dla algorytmow
x_seq <- seq(-20, 20, length = n_grid)
matrVal <- matrix(0, nrow = n_grid, ncol = n_grid)
for(iRow in 1 : n_grid){
  for(iCol in 1 : n_grid){
    matrVal[iRow, iCol] <- myFun(c(x_seq[iRow], x_seq[iCol]))    
  }
}
contour(x_seq, x_seq, matrVal)
lines(sd$x_hist,col='purple',type='l',lwd=8)
lines(sa$x_hist,col='blue',type='l',lwd=4)
lines(ga$x_hist,col='green',type='p',lwd=3)
lines(pso$x_hist,col='red',type='p',lwd=3)

#Wybor najlepszej metody optymalizacji
algoNames <- c("Steepest Descent", "Simulated Annealing", "Genetic Algorithm","Particle Swarm Optimization")
cat("The best solution was found by: ", algoNames[which.min(c(sd$f_opt, sa$f_opt, ga$f_opt, pso$f_opt))])

#Wykres skumulowanego czasu iteracji algorytmu
plot(sd$cumTime[1:maxIter],col='purple',type='l',lwd=1,ylim=c(0,10))
lines(sa$cumTime[1:maxIter],col='blue',type='l',lwd=1)
lines(ga$cumTime[1:maxIter],col='green',type='l',lwd=1)
lines(pso$cumTime[1:maxIter],col='red',type='l',lwd=1)

#Przy ustawieniu ziarna randomizacji set.seed(10) najefektywniejszym algorytmem okazuje sie byc Particle Swarm Optimization.
#Zarowno PSO, jak i algorytm genetyczny zbiegaja do minimum globalnego. Wyzarzanie symulacyjne utyka w minimum lokalnym (choc
#wyjatkowo blisko minimum globalnego - w przypadku wielu prob "udaje sie" w zlym kierunku. Metoda najszybszego spadku zatrzymuje
#sie w minimum lokalnym blisko punktu startowego niezaleznie od liczby prob.

#PSO ma przewage nad algorytmem genetycznym szczegolnie pod wzgledem czas komputacji algorytmu (jak mozna zauwazyc na powyzszym
#wykresie, czas komputacji PSO jest znacznie mniejszy - wynosi srednio ok. 1 sekundy, natomiast algorytmu genetycznego srednio ok. 12 sekund). Przy
#zalozeniu maxIter rownemu 500 iteracji roznica miedzy czasem komputacji tych algorytmow wynosi:
diff <- ga$cumTime[1:maxIter][500]-pso$cumTime[1:maxIter][500]
#Roznica w czasie komputacji algorytmu PSO i genetycznego wynosi ok. 11 sekund.
#PSO okazuje sie byc lepsza metaheurystyka dla problemu optymalizacyjnego funkcji Schaffera.