library("MASS")

flag.para <- flag.var <- betaT <- betaR1 <- betaR2 <- gamma1 <- gamma2<- rep(0,1)
theta1 <- theta2 <- theta3 <- ro2 <- ro1 <- ro3 <- se.betaT <- se.betaR1 <- se.betaR2 <- rep(0,1)
se.gamma1 <- se.gamma2 <- se.theta1 <- se.theta2 <- se.theta3 <- se.ro1 <- se.ro2 <- se.ro3 <- rep(0,1)

z1 <- matrix(c(rep(c(rep(1,10),rep(0,200)),19),rep(1,10)),ncol=20)

betaR1.simu <- as.matrix(c(1.6))
gamma1.simu <- as.matrix(c(-0.2))
betaR2.simu <- as.matrix(c(1.0))
gamma2.simu <- as.matrix(c(-0.8))
betaT.simu  <- as.matrix(c(-1.2))

mu <- c(0, 0, 0)
theta1 <- 1
theta2 <- 1
theta3 <- 1
ro1 <- 0.5
ro2 <- 0.5
ro3 <- 0.5

sigma <- matrix(c(theta1,ro1*sqrt(theta1*theta2),ro2*sqrt(theta1*theta3),ro1*sqrt(theta1*theta2),
                  theta2,ro3*sqrt(theta2*theta3),ro2*sqrt(theta1*theta3),ro3*sqrt(theta2*theta3),theta3),nrow=3)

no <- rep(1,10) 
for(k in 2:20)
{
  no <- c(no,rep(k,10))
}

j <- 1
repeat
{
  x.simu <- runif(200,0,1)
  x.I <- rep(1,200)
  x <- mvrnorm(20, mu, sigma)
  U1 <- as.matrix(x[,1])
  U2 <- as.matrix(x[,2])
  V <- as.matrix(x[,3])
  
  yetaR1.simu <- x.simu %*% betaR1.simu + x.I %*% gamma1.simu + z1 %*% U1
  yetaR2.simu <- x.simu %*% betaR2.simu + x.I %*% gamma2.simu  + z1 %*% U2
  yetaT.simu <- x.simu %*% betaT.simu + z1 %*% V
  
  P0 <- 1 / (1 + exp(yetaR1.simu) + exp(yetaR2.simu))
  P1 <- exp(yetaR1.simu) / (1 + exp(yetaR1.simu) + exp(yetaR2.simu))
  P2 <- exp(yetaR2.simu) / (1 + exp(yetaR1.simu) + exp(yetaR2.simu))
  
  u <- runif(200,0,1)
  R.simu <- rep(0,200)
  
  for(l in 1:200)
  {
    if(u[l] > P0[l]+P1[l])
    {R.simu[l] <- 2}
    if( (u[l] <= P0[l]+P1[l]) & (u[l] > P0[l]) )
    {R.simu[l] <- 1}
  }
  
  T.simu <- rep(0,200)
  C.simu <- runif(200,50,500)
  for(l in 1:200)
  {
    T.simu[l] <- -500 * log(1 - runif(1,0,1)) / exp(yetaT.simu[l])
  }
  delta <- ifelse(T.simu <= C.simu,1,0)
  t.simu <- ifelse(delta == 1,T.simu,C.simu)
  
  data <- cbind(no,t.simu,delta,R.simu,x.simu)
  
  cat("j=",j,'\n')
  result<-dependentmix(data=data,theta1=1,theta2=1,theta3=1,ro1=0.1,ro2=0.1,ro3=0.1,itmax=500,epsilon=0.001)
  
  betaR1[j] <- result$betaR1
  gamma1[j] <- result$gamma1
  betaR2[j] <- result$betaR2
  gamma2[j] <- result$gamma2
  betaT[j] <- result$betaT
  theta1[j] <- result$theta1
  theta2[j] <- result$theta2
  theta3[j] <- result$theta3
  ro1[j] <- result$ro1
  ro2[j] <- result$ro2
  ro3[j] <- result$ro3
  
  se.betaR1[j] <- result$se.betaR1
  se.gamma1[j] <- result$se.gamma1
  se.betaR2[j] <- result$se.betaR2
  se.gamma2[j] <- result$se.gamma2
  se.betaT[j] <- result$se.betaT
  se.theta1[j] <- result$se.theta1
  se.theta2[j] <- result$se.theta2
  se.theta3[j] <- result$se.theta3
  se.ro1[j] <- result$se.ro1
  se.ro2[j] <- result$se.ro2
  se.ro3[j] <- result$se.ro3
  flag.para[j] <- result$flag.para
  flag.var[j] <- result$flag.var
  
  if(result$flag.var == 0)
    { j <- j }
  else{j <- j+1}
  
  if(j == 1000) break
}	

cat("betaR1=",mean(betaR1),"se.betaR1=",mean(se.betaR1),"SD.betaR1=",sqrt(var(betaR1)),'\n')
cat("betaR2=",mean(betaR2),"se.betaR2=",mean(se.betaR2),"SD.betaR2=",sqrt(var(betaR2)),'\n')
cat("gamma1=",mean(gamma1),"se.gamma1=",mean(se.gamma1),"SD.gamma1=",sqrt(var(gamma1)),'\n')
cat("gamma2=",mean(gamma2),"se.gamma2=",mean(se.gamma2),"SD.gamma2=",sqrt(var(gamma2)),'\n')
cat("betaT=",mean(betaT),"se.betaT=",mean(se.betaT),"SD.betaT=",sqrt(var(betaT)),'\n')
cat("theta1=",mean(theta1),"se.theta1=",mean(se.theta1),"SD.theta1=",sqrt(var(theta1)),'\n')
cat("theta2=",mean(theta2),"se.theta2=",mean(se.theta2),"SD.theta2=",sqrt(var(theta2)),'\n')
cat("theta3=",mean(theta3),"se.theta3=",mean(se.theta3),"SD.theta3=",sqrt(var(theta3)),'\n')
cat("ro1=",mean(ro1),"se.ro1=",mean(se.ro1),"SD.ro1=",sqrt(var(ro1)),'\n')
cat("ro2=",mean(ro2),"se.ro2=",mean(se.ro2),"SD.ro2=",sqrt(var(ro2)),'\n')
cat("ro3=",mean(ro3),"se.ro3=",mean(se.ro3),"SD.ro3=",sqrt(var(ro3)),'\n')

output <- rbind(c(mean(betaR1),mean(se.betaR1),sqrt(var(betaR1))),
                c(mean(gamma1),mean(se.gamma1),sqrt(var(gamma1))),
                c(mean(betaR2),mean(se.betaR2),sqrt(var(betaR2))),
                c(mean(gamma2),mean(se.gamma2),sqrt(var(gamma2))),
                c(mean(betaT),mean(se.betaT),sqrt(var(betaT))),
                c(mean(theta1),mean(se.theta1),sqrt(var(theta1))),
                c(mean(theta2),mean(se.theta2),sqrt(var(theta2))),
                c(mean(theta3),mean(se.theta3),sqrt(var(theta3))),
                c(mean(ro1),mean(se.ro1),sqrt(var(ro1))),
                c(mean(ro2),mean(se.ro2),sqrt(var(ro2))),
                c(mean(ro3),mean(se.ro3),sqrt(var(ro3))))
