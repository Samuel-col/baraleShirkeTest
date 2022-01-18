library(Rcpp)
library(RcppArmadillo)

setwd("/home/samuel/Documentos/U/Análisis/")
sourceCpp("baraleShirke.cpp")


# baraleShirkeTest(muestra1,muestra2,B)
# Las muestras deben ser matrices


baraleShirkeTest(X1,X2,5000)


# Parámetro de centralidad:
n4 <- 50;n5 <- 50;p <- 2
X4 <- matrix(rnorm(n4*p),nrow=n4)
X5 <- matrix(rnorm(n5*p,1),nrow=n5)

baraleShirkeTest(X4,X5,30000)

# Parámetro de dispersión:
n6 <- 50;n7 <- 50;p <- 2
X6 <- matrix(rnorm(n6*p),nrow=n4)
X7 <- matrix(rnorm(n7*p,0,sqrt(2)),nrow=n5)

baraleShirkeTest(X6,X7,30000)


# Ambos parámetros

n1 <- 50;n2 <- 50;p <- 2

X1 <- matrix(rnorm(n1*p),nrow=n1)
X2 <- matrix(rnorm(n2*p,1/2,1.5),nrow=n2)

baraleShirkeTest(X1,X2,30000)

# Tortugas pintadas

tur <- read.csv("turtles.csv",sep=" ")

F_tur <- tur[,1:3]
M_tur <- tur[,4:6]

baraleShirkeTest(as.matrix(F_tur),as.matrix(M_tur),5000)
baraleShirkeTest(as.matrix(M_tur),as.matrix(F_tur),5000)


t <- Sys.time()
baraleShirkeTest(as.matrix(F_tur),as.matrix(M_tur),10^6)
Sys.time()-t
