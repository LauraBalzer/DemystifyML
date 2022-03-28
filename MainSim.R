# R code for the commentary "Demystifying Statistical Inference 
# When Using Machine Learning in Causal Research" 
# by Laura B. Balzer & Ted Westling

# Code by Laura B. Balzer

rm(list=ls())

gen.data <- function(n, do.pop=F,  do.misp=F, my.transform=F){

  C1 <- rnorm(n,0,1)
  C2 <- rnorm(n,0,1)
  C3 <- rnorm(n,0,1)
  C4 <- rnorm(n,0,1)
  
  pscore <- plogis(-1 +log(1.75)*(C1+C2+C3+C4) )
  X <- rbinom(n, 1, pscore)
  
  eps <- rnorm(n,0,6)
  get.Y <- function(C1,C2,C3,C4,X,eps){
    120+6*X+3*(C1+C2+C3+C4)+eps
  }
  Y1 <- get.Y(C1,C2,C3,C4,X=1,eps)
  Y0 <- get.Y(C1,C2,C3,C4,X=0,eps)
  Y <- get.Y(C1,C2,C3,C4,X,eps)
  
  if(do.pop){
    # return the pscore & counterfactual outcomes
    yay <- data.frame(cbind(pscore, Y1,Y0))
  } else {
    
    if(do.misp & my.transform){
      # BALZER & WESTLING transformation of the confounders 
      C1 <- exp(C1/2)
      C2 <- C2/(1+exp(C1)) + 10
      C3 <- (C1*C3/25 + 0.6)^3
      C4 <- (C2 + C4 + 20)^2
      
    } else if (do.misp & !my.transform){
      # NAIMI et al. transformation of the confounders
      Z1 <- C1
      Z2 <- C2
      Z3 <- C3
      Z4 <- C4
      
      C1 <- exp(Z1/2)
      C2 <- Z2/(1+exp(Z1)) + 10
      C3 <- (Z1*Z3/25 + 0.6)^3
      C4 <- (Z2 + Z4 + 20)^2
    }
    # return observed data
    yay <- data.frame(cbind(C1,C2,C3,C4,X,Y))
    
  }
  yay
}  



get.summary<- function(x, ATE){
  data.frame(cbind(
    pt=mean(as.numeric(x$estimate)),
    bias= mean(as.numeric(x$estimate)-ATE),
    sd=sd(as.numeric(x$estimate)),
    sd.est=mean(as.numeric(x$std.dev)),
    cover= mean(as.numeric(x$CI1)<= ATE & ATE <= as.numeric(x$CI2) )
  ))
}

library('SuperLearner')
library('ltmle')  

set.seed(1)

# sample size
n <- 5000
# number of iterations
R <- 500
# indicator for "complex" setting
do.misp <- T
# indicator of which transformation in the "complex" setting
my.transform <- T
# use Super Learner for estimation 
do.SL <- T
# bound the pscore
do.gbounds <- F
# numbe of folds in Super Learner
V <- 10

# population-level metrics
pop <- gen.data(n=100000, do.pop=T)
ATE <- mean(pop$Y1)-mean(pop$Y0);ATE
summary(pop)
png('Hist.png')
hist(pop$pscore, main='', xlab='True propensity scores') 
dev.off()

if(do.SL){
  SL.library <- c('SL.glm', 'SL.step.interaction', 'SL.earth', 'SL.mean')
  SL.file <- ''
}else{
  SL.library <- NULL
  SL.file <- 'glm'
}

if(do.gbounds){
  gbounds <- c(0.025, 0.975)
  g.file <- 'Gbounds'
}else{
  gbounds <- c(0.01, 1) #default 
  g.file <- ''
}

file.name <- paste('Complex', do.misp,'MyTransform', my.transform,
                   'V', V, SL.file, g.file,
                   'n',n, 'R', R, '.', Sys.Date(), '.Rdata', sep='')
print(file.name)


these <- c('estimate', 'std.dev', 'CI1', 'CI2')
OUT <- data.frame(matrix(NA, nrow=R, ncol=length(these)))
colnames(OUT) <- these
for(r in 1:R){
  O <- gen.data(n, do.misp=do.misp, my.transform=my.transform)
  est  <- ltmle(data=O, Anodes='X', Ynodes='Y', abar=list(1,0), 
                SL.library=SL.library, estimate.time=F,
                SL.cvControl=list(V=V),
                gbounds=gbounds)
 
  OUT[r,] <-  unlist(summary(est)$effect.measures$ATE)[these]
  print(r)              
}


get.summary(OUT, ATE)

save(ATE,n,R,do.misp, SL.library, OUT, file=file.name)

# 
# ###############################################################
# # # Demonstrating ills of overfitting
# # # See discussion in Supplementary Materials
# # # NOT recommended in practice 
# O <- gen.data(n=200, do.misp=F)
# # estimate E(Y|X,W) with observed outcome
# init <- O$Y
# pscore <- predict(glm(X~C1+C2+C3+C4, data=O, family='binomial'), type='response')
# H <- O$X/pscore
# glm(O$Y ~ -1 +offset(init) + H)
