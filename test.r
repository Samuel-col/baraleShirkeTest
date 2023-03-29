setwd("/home/samuel/Documentos/U/Seminario/baraleShirkeTest/")
source("baraleShirke.R")


# baraleshirke.test <- function(X1,X2,
#                               depth = "Mahalanobis",
#                               NIter = 1000,
#                               alpha = 0.05,
#                               returnDepths = F,
#                               returnSamples = F)
# Las muestras deben ser matrices
# Profundidades disponibles:
# * betaSkeleton
# * halfspace
# * L2
# * Mahalanobis
# * projection
# * potential
# * qhpeeling
# * simplicial
# * simplicialVolume
# * spatial
# * zonoid



n1 <- 70;n2 <- 70;p <- 2
X1 <- matrix(rnorm(n1*p,1,0.5),nrow=n1)
X2 <- matrix(rnorm(n2*p,0,0.32),nrow=n2)

my.test <- baraleshirke.test(X1,X2,depth = "projection",
                             returnDepths = T,
                             returnSamples = T)

my.test



plot(my.test$DepthVals,col = rep(1:2,c(n1,n2)),
     pch = 20)

plot(c(my.test$DepthVals[1:n1],my.test$DepthVals[(n1+1):(n1+n2)]),
     rep(0,n1+n2),
     col = rep(1:2,c(n1,n2)),
     pch = 20,ylim = c(-1,1),cex = 4)
abline(h = 0)

plot(density(my.test$Samples),
     xlim = range(c(my.test$Samples,my.test$Statistic)))
abline(v = my.test$Statistic, h = 0)






# Parámetro de centralidad:
n4 <- 50;n5 <- 50;p <- 2
X4 <- matrix(rnorm(n4*p),nrow=n4)
X5 <- matrix(rnorm(n5*p,1),nrow=n5)

baraleshirke.test(X4,X5,NIter = 30000,alpha = 0.03,
                  returnDepths = T,returnSamples = T)
baraleshirke.test(X4,X5,NIter = 30000,alpha = 0.03,
                  returnDepths = F,returnSamples = F)




# Parámetro de dispersión:
n6 <- 50;n7 <- 50;p <- 2
X6 <- matrix(rnorm(n6*p),nrow=n4)
X7 <- matrix(rnorm(n7*p,0,sqrt(2)),nrow=n5)

baraleshirke.test(X6,X7,NIter = 3000,depth = "simplicial")


# Ambos parámetros

n1 <- 50;n2 <- 50;p <- 2

X1 <- matrix(rnorm(n1*p),nrow=n1)
X2 <- matrix(rnorm(n2*p,1/2,1.5),nrow=n2)

baraleshirke.test(X1,X2,NIter = 3000,depth = "potential")
baraleshirke.test(X1,X2,NIter = 3000,depth = "qhpeeling")

#------ Tortugas pintadas -----------

tur <- read.csv("turtles.csv",sep=" ")
View(tur)


F_tur <- tur[,1:3]
M_tur <- tur[,4:6]

baraleshirke.test(F_tur,M_tur,NIter = 1e4,depth = "zonoid")
baraleshirke.test(M_tur,F_tur,NIter = 1e2,depth = "zonoid")



t <- Sys.time()
baraleshirke.test(F_tur,M_tur,depth = "L2",NIter = 1e4)
Sys.time()-t


my.test <- baraleshirke.test(F_tur,M_tur,
                             depth = "halfspace",
                             NIter = 3e4,
                             returnDepths = T,
                             returnSamples = T)

my.test


plot(my.test$DepthVals,col = rep(1:2,c(n1,n2)),
     pch = c(20,17))

plot(c(my.test$DepthVals[1:n1],my.test$DepthVals[(n1+1):(n1+n2)]),
     rep(0,n1+n2),
     col = rep(1:2,c(n1,n2)),
     pch = 20,ylim = c(-1,1),cex = 4)
abline(h = 0)

plot(density(my.test$Samples),
     xlim = range(c(my.test$Samples,my.test$Statistic)))
abline(v = my.test$Statistic, h = 0)
