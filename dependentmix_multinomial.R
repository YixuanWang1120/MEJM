
#---------------------------------------------
# dependent censor model with random effects
#---------------------------------------------
dependentmix<-function(data,theta1=1,theta2=1,theta3=1,ro1=0.1,ro2=0.1,ro3=0.1,itmax=500,epsilon=0.001)
{
	data <- as.matrix(data)
	clinic <- data[,1]
  n <- length(data[,1])
	nb <- ncol(data) - 4
	m <- length(unique(data[,1]))
	
#generate Z for random effects
  z <- rep(1,sum(clinic == 1))
	for(i1 in 2:m)
	{
	 z <- c(z,rep(0,n),rep(1,sum(clinic == i1)))	
  }
  z <- matrix(z,ncol = m)    
  r <- cbind(data, z)
  r <- r[sort.list(r[,2]),]
  x <- as.matrix(r[, 5:(4 + nb)])
  t <- as.vector(r[, 2])
  deltaT <- as.vector(r[, 3])
  z <- as.matrix(r[,(5 + nb):(4 + nb + m)])
  R0 <- as.vector(ifelse(r[, 4] == 0, 1, 0))
  R1 <- as.vector(ifelse(r[, 4] == 1, 1, 0))
  R2 <- as.vector(ifelse(r[, 4] == 2, 1, 0))
  M <- lower.tri(matrix(rep(1, n^2), ncol = n), diag = T)
  M <- ifelse(M == T, 1, 0)
  x.I <- as.matrix(rep(1, n))

  one <- rep(1, n)
  tT.event <- cbind(t, deltaT)

for(k2 in 2:n)
{
	if(all(tT.event[k2,] == tT.event[k2-1,]))
	{
		M[,k2] <- M[,k2-1]
	}
}

#initial values
	
  flag.var <- 0
	flag.para <- 0
	Diverge <- 0
	betaR1 <- rep(0.5, nb)
	gamma1 <- rep(-0.2, 1)
	betaR2 <- rep(0.5, nb)
	gamma2 <- rep(-0.8, 1)
	betaT <- rep(-1.2, nb)
	u1 <- rep(0, m)
	u2 <- rep(0, m)
	v <- rep(0, m)
  yetaR1 <- x %*% betaR1 + x.I %*% gamma1  + z %*% u1
  yetaR2 <- x %*% betaR2 + x.I %*% gamma2  + z %*% u2
  yetaT <- x %*% betaT + z %*% v
################################
	for(iter0 in 1:itmax)
	{
		for(iter in 1:itmax)
		{
		  W1 <- as.vector((exp(yetaR1)+exp(yetaR1+yetaR2)) / ((1+exp(yetaR1)+exp(yetaR2))^2))
		  SecdR11 <- diag(W1)
		  
		  W2 <- as.vector((exp(yetaR2)+exp(yetaR1+yetaR2)) / ((1+exp(yetaR1)+exp(yetaR2))^2))
		  SecdR22 <- diag(W2)
		  
	    W3 <- as.vector( -exp(yetaR1+yetaR2) / ((1 + exp(yetaR1) + exp(yetaR2))^2))
	    SecdR12 <- diag(W3)
	    
	    W4 <- as.vector(exp(yetaT))
		  A <- as.vector(deltaT / (t(M) %*% (exp(yetaT))))
	    B <- as.vector(M %*% (A*one))
	    SecdT <- diag(W4 * B)-diag(W4) %*% M %*% diag(A*A) %*% t(M) %*% diag(W4)
	    
	    x.R <- cbind(x,x.I)
	    
	    Info11 <- t(x.R) %*% SecdR11 %*% x.R
	    Info12 <- t(x.R) %*% SecdR12 %*% x.R
	    Info13 <- matrix(rep(0,(nb+1)*nb), ncol=nb)
	    Info14 <- t(x.R) %*% SecdR11 %*% z
	    Info15 <- t(x.R) %*% SecdR12 %*% z
	    Info16 <- matrix(rep(0, (nb+1)*m), ncol=m)
	    
	    Info21 <- t(Info12)
	    Info22 <- t(x.R) %*% SecdR22 %*% x.R
      Info23 <- matrix(rep(0,(nb+1)*nb), ncol=nb)
      Info24 <- t(x.R) %*% SecdR12 %*% z
      Info25 <- t(x.R) %*% SecdR22 %*% z
      Info26 <- matrix(rep(0, (nb+1)*m), ncol=m)
      
      Info31 <- t(Info13)
      Info32 <- t(Info23)
      Info33 <- t(x) %*% SecdT %*% x
      Info34 <- matrix(rep(0, nb*m), ncol=m)
      Info35 <- matrix(rep(0, nb*m), ncol=m)
      Info36 <- t(x) %*% SecdT %*% z
      
      Info41 <- t(Info14)
      Info42 <- t(Info24)
      Info43 <- t(Info34)
      Info44 <- t(z) %*% SecdR11 %*% z + diag(rep((1-ro3^2)/(theta1*(1-(ro1^2+ro2^2+ro3^2-2*ro1*ro2*ro3))),m))
      Info45 <- t(z) %*% SecdR12 %*% z + diag(rep((ro2*ro3-ro1)/(sqrt(theta1*theta2)*(1-(ro1^2+ro2^2+ro3^2-2*ro1*ro2*ro3))),m))
      Info46 <- diag(rep((ro1*ro3-ro2)/(sqrt(theta1*theta3)*(1-(ro1^2+ro2^2+ro3^2-2*ro1*ro2*ro3))),m))
      
      Info51 <- t(Info15)
      Info52 <- t(Info25)
      Info53 <- t(Info35)
      Info54 <- t(Info45)
      Info55 <- t(z) %*% SecdR22 %*% z + diag(rep((1-ro2^2)/(theta2*(1-(ro1^2+ro2^2+ro3^2-2*ro1*ro2*ro3))),m))
      Info56 <- diag(rep((ro1*ro2-ro3)/(sqrt(theta2*theta3)*(1-(ro1^2+ro2^2+ro3^2-2*ro1*ro2*ro3))),m))
      
      Info61 <- t(Info16)
      Info62 <- t(Info26)
      Info63 <- t(Info36)
      Info64 <- t(Info46)
      Info65 <- t(Info56)
      Info66 <- t(z) %*% SecdT %*% z + diag(rep((1-ro1^2)/(theta3*(1-(ro1^2+ro2^2+ro3^2-2*ro1*ro2*ro3))),m))
      
      Info1 <- cbind(Info11, Info12, Info13, Info14, Info15, Info16)
	    Info2 <- cbind(Info21, Info22, Info23, Info24, Info25, Info26)
	    Info3 <- cbind(Info31, Info32, Info33, Info34, Info35, Info36)
	    Info4 <- cbind(Info41, Info42, Info43, Info44, Info45, Info46)
	    Info5 <- cbind(Info51, Info52, Info53, Info54, Info55, Info56)
	    Info6 <- cbind(Info61, Info62, Info63, Info64, Info65, Info66)
	    Info0 <- rbind(Info1, Info2, Info3, Info4, Info5, Info6)
      Info <- solve(Info0)
          
      First.yetaR1 <- as.matrix(R1 - exp(yetaR1) / (1 + exp(yetaR1) + exp(yetaR2)))
      First.yetaR2 <- as.matrix(R2 - exp(yetaR2) / (1 + exp(yetaR1) + exp(yetaR2)))
      First.yetaT <- as.matrix(deltaT - W4 * (M %*% (A*one)))	 
      
      First.betaR1 <- t(x) %*% First.yetaR1
      First.gamma1 <- t(x.I) %*% First.yetaR1
      First.betaR2 <- t(x) %*% First.yetaR2
      First.gamma2 <- t(x.I) %*% First.yetaR2
      First.betaT <- t(x) %*% First.yetaT           
	    
      First.u1 <- t(z)%*%First.yetaR1 - (( u1 * (1-ro3^2) * theta2 * theta3
	                                       - u2 * (ro1-ro2*ro3) * theta3 * sqrt(theta1 * theta2)
	                                       - v  * (ro2-ro1*ro3) * theta2 * sqrt(theta1 * theta3))
	                                       /(theta1*theta2*theta3*(1-(ro1^2+ro2^2+ro3^2-2*ro1*ro2*ro3))))   
	    
	    First.u2 <- t(z)%*%First.yetaR2 - ((- u1 * (ro1-ro2*ro3) * theta3 * sqrt(theta1 * theta2)
	                                        + u2 * (1-ro2^2) * theta1 * theta3 
	                                        - v  * (ro3-ro1*ro2) * theta1 * sqrt(theta2 * theta3))
	                                        /(theta1*theta2*theta3*(1-(ro1^2+ro2^2+ro3^2-2*ro1*ro2*ro3))))    
	    
	    First.v <- t(z)%*%First.yetaT -((- u1 * (ro2-ro1*ro3) * theta2 * sqrt(theta1 * theta3)
	                                     - u2 * (ro3-ro1*ro2) * theta1 * sqrt(theta2 * theta3)
	                                     + v  * (1-ro1^2) * theta1 * theta2 )
	                                     /(theta1*theta2*theta3*(1-(ro1^2+ro2^2+ro3^2-2*ro1*ro2*ro3))))
	
      para <- as.matrix(c(betaR1, gamma1, betaR2, gamma2, betaT, u1, u2, v))
      
      para1 <- para + Info %*% c(First.betaR1,First.gamma1,First.betaR2,First.gamma2,First.betaT,First.u1,First.u2,First.v)
			
      betaR1 <- as.matrix(para1[1:nb, 1])
      gamma1 <- as.matrix(para1[nb+1:1, 1])
			betaR2 <- as.matrix(para1[nb+1+1:nb, 1])
			gamma2 <- as.matrix(para1[2*nb+1+1:1, 1])
			betaT <- as.matrix(para1[2*nb+2+1:nb, 1])
			u1 <- as.matrix(para1[3*nb+2+1:m, 1])
			u2 <- as.matrix(para1[3*nb+2+m+1:m, 1])
			v <- as.matrix(para1[3*nb+2+m+m+1:m, 1])
			
			yetaR1 <- x %*% betaR1 + x.I %*% gamma1 + z %*% u1
			yetaR2 <- x %*% betaR2 + x.I %*% gamma2 + z %*% u2
			yetaT <- x %*% betaT + z %*% v
			
			
			if(max(abs(para1 - para)) > 30)
			{cat("break");
			  Diverge <- 1;break}
			
			if(max(abs(para1 - para)) < epsilon)
			{flag.para <- 1;break}
			
			para <- para1
			#cat("iter=",iter,"betaT=",betaT,"betaR1=",betaR1,"betaR2=",betaR2,'\n')		
		}
	  
		  if(Diverge)
		    {break} 
	 
	  T44 <-Info[3*nb+2+1:(3*m), 3*nb+2+1:(3*m)]
	  J1 <- diag(c(rep(1,m), rep(0,m), rep(0,m)))
	  J2 <- rbind(cbind(diag(rep(0,m)), diag(rep(1,m)), diag(rep(0,m))),
	              cbind(diag(rep(1,m)), diag(rep(0,m)), diag(rep(0,m))),
	              cbind(diag(rep(0,m)), diag(rep(0,m)), diag(rep(0,m))))
	  J3 <- rbind(cbind(diag(rep(0,m)), diag(rep(0,m)), diag(rep(1,m))),
	              cbind(diag(rep(0,m)), diag(rep(0,m)), diag(rep(0,m))),
	              cbind(diag(rep(1,m)), diag(rep(0,m)), diag(rep(0,m))))
	  J4 <- diag(c(rep(0,m), rep(1,m), rep(0,m)))
	  J5 <- rbind(cbind(diag(rep(0,m)), diag(rep(0,m)), diag(rep(0,m))),
	              cbind(diag(rep(0,m)), diag(rep(0,m)), diag(rep(1,m))),
	              cbind(diag(rep(0,m)), diag(rep(1,m)), diag(rep(0,m))))
	  J6 <- diag(c(rep(0,m), rep(0,m), rep(1,m)))
	  L1 <- sum(diag(J1 %*% (T44 + c(u1,u2,v) %*% t(c(u1,u2,v)))))
	  L2 <- sum(diag(J2 %*% (T44 + c(u1,u2,v) %*% t(c(u1,u2,v)))))/2
	  L3 <- sum(diag(J3 %*% (T44 + c(u1,u2,v) %*% t(c(u1,u2,v)))))/2
	  L4 <- sum(diag(J4 %*% (T44 + c(u1,u2,v) %*% t(c(u1,u2,v)))))
	  L5 <- sum(diag(J5 %*% (T44 + c(u1,u2,v) %*% t(c(u1,u2,v)))))/2
	  L6 <- sum(diag(J6 %*% (T44 + c(u1,u2,v) %*% t(c(u1,u2,v)))))
	  
	  theta10 <- L1/m
	  theta20 <- L4/m
	  theta30 <- L6/m
	  
		ro10 <- L2/sqrt(L1*L4)
		ro20 <- L3/sqrt(L1*L6)
		ro30 <- L5/sqrt(L4*L6)
		  
		if(max(abs(c(theta1,theta2,theta3,ro1,ro2,ro3)-c(theta10,theta20,theta30,ro10,ro20,ro30)))>100)
		  {cat("break");Diverge<-1;break}
		if(max(abs(c(theta1,theta2,theta3,ro1,ro2,ro3)-c(theta10,theta20,theta30,ro10,ro20,ro30)))<epsilon)
		  {flag.var<-1;break}
		  theta1 <- theta10
		  theta2 <- theta20
		  theta3 <- theta30
		  ro1 <- ro10
		  ro2 <- ro20
		  ro3 <- ro30
		    
		  #cat("iter0=",iter0,"theta1=",theta1,"theta2=",theta2,'\n')
		}              
      if(flag.var)

########## standard error ###########
  if(Diverge==1 | flag.var==0){
	
    cat('break',Diverge,flag.para,flag.var,'\n')
  
    list(flag.var=flag.var,Diverge=Diverge,betaT=0,betaR1=0,betaR2=0,gamma1=0,gamma2=0,
       se.betaT=0,se.betaR1=0,se.betaR2=0,se.gamma1=0,se.gamma2=0,theta1=0,
       theta2=0,theta3=0,ro1=0,ro2=0,ro3=0,se.theta1=0,se.theta2=0,se.theta3=0,
       se.ro1=0,se.ro2=0,se.ro3=0)
  }
  else{
    se.betaR1 <- sqrt(diag(Info)[1:nb])
    se.gamma1 <- sqrt(diag(Info)[nb+1:1])
    se.betaR2 <- sqrt(diag(Info)[nb+1+1:nb])
    se.gamma2 <- sqrt(diag(Info)[2*nb+1+1:1])
    se.betaT <- sqrt(diag(Info)[2*nb+2+1:nb])
  
    Info.random <- Info[3*nb+2+1:(3*m),3*nb+2+1:(3*m)]
    
    KK <- 1-(ro1^2+ro2^2+ro3^2-2*ro1*ro2*ro3)
    Omega <- rbind(cbind(diag(rep(theta1,m)),diag(rep(ro1*sqrt(theta1*theta2),m)),diag(rep(ro2*sqrt(theta1*theta3),m))),
                   cbind(diag(rep(ro1*sqrt(theta1*theta2),m)),diag(rep(theta2,m)),diag(rep(ro3*sqrt(theta2*theta3),m))),
                   cbind(diag(rep(ro2*sqrt(theta1*theta3),m)),diag(rep(ro3*sqrt(theta2*theta3),m)),diag(rep(theta3,m))))
    
    InO11 <- diag(rep((1-ro3^2)*theta2*theta3, m))
    InO12 <- diag(rep((-ro1+ro2*ro3)*theta3*sqrt(theta1*theta2), m))
    InO13 <- diag(rep((-ro2+ro1*ro3)*theta2*sqrt(theta1*theta3), m))
    InO21 <- t(InO12)
    InO22 <- diag(rep((1-ro2^2)*theta1*theta3, m))
    InO23 <- diag(rep((-ro3+ro1*ro2)*theta1*sqrt(theta2*theta3), m))
    InO31 <- t(InO13)
    InO32 <- t(InO23)
    InO33 <- diag(rep((1-ro1^2)*theta1*theta2, m))
    InO1 <- cbind(InO11, InO12, InO13)
    InO2 <- cbind(InO21, InO22, InO23)
    InO3 <- cbind(InO31, InO32, InO33)
    InO <- rbind(InO1, InO2, InO3)/(theta1*theta2*theta3*KK)
    
    OTh111 <- diag(rep(1, m))
    OTh112 <- diag(rep( ro1*theta2/(2*sqrt(theta1*theta2)), m))
    OTh113 <- diag(rep( ro2*theta3/(2*sqrt(theta1*theta3)), m))
    OTh121 <- t(OTh112)
    OTh122 <- diag(rep(0, m))
    OTh123 <- diag(rep(0, m))
    OTh131 <- t(OTh113)
    OTh132 <- t(OTh123)
    OTh133 <- diag(rep(0, m))
    OTh11 <- cbind(OTh111, OTh112, OTh113)
    OTh12 <- cbind(OTh121, OTh122, OTh123)
    OTh13 <- cbind(OTh131, OTh132, OTh133)
    OTh1 <- rbind(OTh11, OTh12, OTh13)
    
    OTh211 <- diag(rep(0, m))
    OTh212 <- diag(rep( ro1*theta1/(2*sqrt(theta1*theta2)), m))
    OTh213 <- diag(rep(0, m))
    OTh221 <- t(OTh212)
    OTh222 <- diag(rep(1, m))
    OTh223 <- diag(rep( ro3*theta3/(2*sqrt(theta2*theta3)), m))
    OTh231 <- t(OTh213)
    OTh232 <- t(OTh223)
    OTh233 <- diag(rep(0, m))
    OTh21 <- cbind(OTh211, OTh212, OTh213)
    OTh22 <- cbind(OTh221, OTh222, OTh223)
    OTh23 <- cbind(OTh231, OTh232, OTh233)
    OTh2 <- rbind(OTh21, OTh22, OTh23)
    
    OTh311 <- diag(rep(0, m))
    OTh312 <- diag(rep(0, m))
    OTh313 <- diag(rep( ro2*theta1/(2*sqrt(theta1*theta3)), m))
    OTh321 <- t(OTh312)
    OTh322 <- diag(rep(0, m))
    OTh323 <- diag(rep( ro3*theta2/(2*sqrt(theta2*theta3)), m))
    OTh331 <- t(OTh313)
    OTh332 <- t(OTh323)
    OTh333 <- diag(rep(1, m))
    OTh31 <- cbind(OTh311, OTh312, OTh313)
    OTh32 <- cbind(OTh321, OTh322, OTh323)
    OTh33 <- cbind(OTh331, OTh332, OTh333)
    OTh3 <- rbind(OTh31, OTh32, OTh33)
    
    ORo111 <- diag(rep( 0, m))
    ORo112 <- diag(rep( sqrt(theta1*theta2), m))
    ORo113 <- diag(rep( 0, m))
    ORo121 <- t(ORo112)
    ORo122 <- diag(rep(0, m))
    ORo123 <- diag(rep(0, m))
    ORo131 <- t(ORo113)
    ORo132 <- t(ORo123)
    ORo133 <- diag(rep(0, m))
    ORo11 <- cbind(ORo111, ORo112, ORo113)
    ORo12 <- cbind(ORo121, ORo122, ORo123)
    ORo13 <- cbind(ORo131, ORo132, ORo133)
    ORo1 <- rbind(ORo11, ORo12, ORo13)
    
    ORo211 <- diag(rep( 0, m))
    ORo212 <- diag(rep( 0, m))
    ORo213 <- diag(rep( sqrt(theta1*theta3), m))
    ORo221 <- t(ORo212)
    ORo222 <- diag(rep(0, m))
    ORo223 <- diag(rep(0, m))
    ORo231 <- t(ORo213)
    ORo232 <- t(ORo223)
    ORo233 <- diag(rep(0, m))
    ORo21 <- cbind(ORo211, ORo212, ORo213)
    ORo22 <- cbind(ORo221, ORo222, ORo223)
    ORo23 <- cbind(ORo231, ORo232, ORo233)
    ORo2 <- rbind(ORo21, ORo22, ORo23)
    
    ORo311 <- diag(rep( 0, m))
    ORo312 <- diag(rep( 0, m))
    ORo313 <- diag(rep( 0, m))
    ORo321 <- t(ORo312)
    ORo322 <- diag(rep(0, m))
    ORo323 <- diag(rep(sqrt(theta2*theta3), m))
    ORo331 <- t(ORo313)
    ORo332 <- t(ORo323)
    ORo333 <- diag(rep(0, m))
    ORo31 <- cbind(ORo311, ORo312, ORo313)
    ORo32 <- cbind(ORo321, ORo322, ORo323)
    ORo33 <- cbind(ORo331, ORo332, ORo333)
    ORo3 <- rbind(ORo31, ORo32, ORo33)
    
    InOTh111 <- diag(rep( -2*(1-ro3^2)*theta2*theta3, m))
    InOTh112 <- diag(rep( (ro1-ro2*ro3)*theta3*sqrt(theta1*theta2), m))
    InOTh113 <- diag(rep( (ro2-ro1*ro3)*theta2*sqrt(theta1*theta3), m))
    InOTh121 <- t(InOTh112)
    InOTh122 <- diag(rep(0, m))
    InOTh123 <- diag(rep(0, m))
    InOTh131 <- t(InOTh113)
    InOTh132 <- t(InOTh123)
    InOTh133 <- diag(rep(0, m))
    InOTh11 <- cbind(InOTh111, InOTh112, InOTh113)
    InOTh12 <- cbind(InOTh121, InOTh122, InOTh123)
    InOTh13 <- cbind(InOTh131, InOTh132, InOTh133)
    InOTh1 <- rbind(InOTh11, InOTh12, InOTh13)/(2*theta1^2*theta2*theta3*KK)
    
    InOTh211 <- diag(rep( 0, m))
    InOTh212 <- diag(rep( (ro1-ro2*ro3)*theta3*sqrt(theta1*theta2), m))
    InOTh213 <- diag(rep( 0, m))
    InOTh221 <- t(InOTh212)
    InOTh222 <- diag(rep( -2*(1-ro2^2)*theta1*theta3, m))
    InOTh223 <- diag(rep( (ro3-ro1*ro2)*theta1*sqrt(theta2*theta3), m))
    InOTh231 <- t(InOTh213)
    InOTh232 <- t(InOTh223)
    InOTh233 <- diag(rep(0, m))
    InOTh21 <- cbind(InOTh211, InOTh212, InOTh213)
    InOTh22 <- cbind(InOTh221, InOTh222, InOTh223)
    InOTh23 <- cbind(InOTh231, InOTh232, InOTh233)
    InOTh2 <- rbind(InOTh21, InOTh22, InOTh23)/(2*theta1*theta2^2*theta3*KK)
    
    InOTh311 <- diag(rep( 0, m))
    InOTh312 <- diag(rep( 0, m))
    InOTh313 <- diag(rep( (ro2-ro1*ro3)*theta2*sqrt(theta1*theta3), m))
    InOTh321 <- t(InOTh312)
    InOTh322 <- diag(rep( 0, m))
    InOTh323 <- diag(rep( (ro3-ro1*ro2)*theta1*sqrt(theta2*theta3), m))
    InOTh331 <- t(InOTh313)
    InOTh332 <- t(InOTh323)
    InOTh333 <- diag(rep(-2*(1-ro1^2)*theta1*theta2, m))
    InOTh31 <- cbind(InOTh311, InOTh312, InOTh313)
    InOTh32 <- cbind(InOTh321, InOTh322, InOTh323)
    InOTh33 <- cbind(InOTh331, InOTh332, InOTh333)
    InOTh3 <- rbind(InOTh31, InOTh32, InOTh33)/(2*theta1*theta2*theta3^2*KK)
    
    InORo111 <- diag(rep( 2*(ro1-ro2*ro3)*(1-ro3^2)*theta2*theta3, m))
    InORo112 <- diag(rep( -(KK+2*(ro1-ro2*ro3)^2)*theta3*sqrt(theta1*theta2), m))
    InORo113 <- diag(rep( (ro3*(1+ro1^2+ro2^2-ro3^2)-2*ro1*ro2)*theta2*sqrt(theta1*theta3), m))
    InORo121 <- t(InORo112)
    InORo122 <- diag(rep( 2*(ro1-ro2*ro3)*(1-ro2^2)*theta1*theta3, m))
    InORo123 <- diag(rep( (ro2*(1+ro1^2-ro2^2+ro3^2)-2*ro1*ro3)*theta1*sqrt(theta2*theta3), m))
    InORo131 <- t(InORo113)
    InORo132 <- t(InORo123)
    InORo133 <- diag(rep( 2*(ro1*(ro2^2+ro3^2)-ro2*ro3*(ro1^2+1))*theta1*theta2, m))
    InORo11 <- cbind(InORo111, InORo112, InORo113)
    InORo12 <- cbind(InORo121, InORo122, InORo123)
    InORo13 <- cbind(InORo131, InORo132, InORo133)
    InORo1 <- rbind(InORo11, InORo12, InORo13)/(theta1*theta2*theta3*KK^2)
    
    InORo211 <- diag(rep( 2*(ro2-ro1*ro3)*(1-ro3^2)*theta2*theta3, m))
    InORo212 <- diag(rep( (ro3*(1+ro1^2+ro2^2-ro3^2)-2*ro1*ro2)*theta3*sqrt(theta1*theta2), m))
    InORo213 <- diag(rep( -(KK+2*(ro2-ro1*ro3)^2)*theta2*sqrt(theta1*theta3), m))
    InORo221 <- t(InORo212)
    InORo222 <- diag(rep( 2*(ro2*(ro1^2+ro3^2)-ro1*ro3*(ro2^2+1))*theta1*theta3, m))
    InORo223 <- diag(rep( (ro1*(1-ro1^2+ro2^2+ro3^2)-2*ro2*ro3)*theta1*sqrt(theta2*theta3), m))
    InORo231 <- t(InORo213)
    InORo232 <- t(InORo223)
    InORo233 <- diag(rep( 2*(ro2-ro1*ro3)*(1-ro1^2)*theta1*theta2, m))
    InORo21 <- cbind(InORo211, InORo212, InORo213)
    InORo22 <- cbind(InORo221, InORo222, InORo223)
    InORo23 <- cbind(InORo231, InORo232, InORo233)
    InORo2 <- rbind(InORo21, InORo22, InORo23)/(theta1*theta2*theta3*KK^2)
    
    InORo311 <- diag(rep( 2*(ro3*(ro1^2+ro2^2)-ro1*ro2*(ro3^2+1))*theta2*theta3, m))
    InORo312 <- diag(rep( (ro2*(1+ro1^2-ro2^2+ro3^2)-2*ro1*ro3)*theta3*sqrt(theta1*theta2), m))
    InORo313 <- diag(rep( (ro1*(1-ro1^2+ro2^2+ro3^2)-2*ro2*ro3)*theta2*sqrt(theta1*theta3), m))
    InORo321 <- t(InORo312)
    InORo322 <- diag(rep( 2*(ro3-ro1*ro2)*(1-ro2^2)*theta1*theta3, m))
    InORo323 <- diag(rep( -(KK+2*(ro3-ro1*ro2)^2)*theta1*sqrt(theta2*theta3), m))
    InORo331 <- t(InORo313)
    InORo332 <- t(InORo323)
    InORo333 <- diag(rep( 2*(ro3-ro1*ro2)*(1-ro1^2)*theta1*theta2, m))
    InORo31 <- cbind(InORo311, InORo312, InORo313)
    InORo32 <- cbind(InORo321, InORo322, InORo323)
    InORo33 <- cbind(InORo331, InORo332, InORo333)
    InORo3 <- rbind(InORo31, InORo32, InORo33)/(theta1*theta2*theta3*KK^2)
  
    K1 <- T44 %*% InOTh1
    K2 <- Omega %*% InOTh1
    K3 <- T44 %*% InOTh2
    K4 <- Omega %*% InOTh2
    K5 <- T44 %*% InOTh3
    K6 <- Omega %*% InOTh3
    K7 <- T44 %*% InORo1
    K8 <- Omega %*% InORo1
    K9 <- T44 %*% InORo2
    K10 <- Omega %*% InORo2
    K11 <- T44 %*% InORo3
    K12 <- Omega %*% InORo3
  
    a11 <- sum(diag(K1 %*% K1 + K2 %*% K2 - 2 * K1 %*% K2))
    a12 <- sum(diag(K1 %*% K3 + K2 %*% K4 - 2 * K1 %*% K4))
    a13 <- sum(diag(K1 %*% K5 + K2 %*% K6 - 2 * K1 %*% K6))
    a14 <- sum(diag(K1 %*% K7 + K2 %*% K8 - 2 * K1 %*% K8))
    a15 <- sum(diag(K1 %*% K9 + K2 %*% K10 - 2 * K1 %*% K10))
    a16 <- sum(diag(K1 %*% K11 + K2 %*% K12 - 2 * K1 %*% K12))
    a22 <- sum(diag(K3 %*% K3 + K4 %*% K4 - 2 * K3 %*% K4))
    a23 <- sum(diag(K3 %*% K5 + K4 %*% K6 - 2 * K3 %*% K6))
    a24 <- sum(diag(K3 %*% K7 + K4 %*% K8 - 2 * K3 %*% K8))
    a25 <- sum(diag(K3 %*% K9 + K4 %*% K10 - 2 * K3 %*% K10))
    a26 <- sum(diag(K3 %*% K11 + K4 %*% K12 - 2 * K3 %*% K12))
    a33 <- sum(diag(K5 %*% K5 + K6 %*% K6 - 2 * K5 %*% K6))
    a34 <- sum(diag(K5 %*% K7 + K6 %*% K8 - 2 * K5 %*% K8))
    a35 <- sum(diag(K5 %*% K9 + K6 %*% K10 - 2 * K5 %*% K10))
    a36 <- sum(diag(K5 %*% K11 + K6 %*% K12 - 2 * K5 %*% K12))
    a44 <- sum(diag(K7 %*% K7 + K8 %*% K8 - 2 * K7 %*% K8))
    a45 <- sum(diag(K7 %*% K9 + K8 %*% K10 - 2 * K7 %*% K10))
    a46 <- sum(diag(K7 %*% K11 + K8 %*% K12 - 2 * K7 %*% K12))
    a55 <- sum(diag(K9 %*% K9 + K10 %*% K10 - 2 * K9 %*% K10))
    a56 <- sum(diag(K9 %*% K11 + K10 %*% K12 - 2 * K9 %*% K12))
    a66 <- sum(diag(K11 %*% K11 + K12 %*% K12 - 2 * K11 %*% K12))
    a1 <- cbind(a11, a12, a13, a14, a15, a16)
    a2 <- cbind(t(a12), a22, a23, a24, a25, a26)
    a3 <- cbind(t(a13), t(a23), a33, a34, a35, a36)
    a4 <- cbind(t(a14), t(a24), t(a34), a44, a45, a46)
    a5 <- cbind(t(a15), t(a25), t(a35), t(a45), a55, a56)
    a6 <- cbind(t(a16), t(a26), t(a36), t(a46), t(a56), a66)
    a <- rbind(a1, a2, a3, a4, a5, a6)

    var.v <- as.vector(diag(2 * solve(a)))
    #cat("var=",var.v,'\n')

    se.theta1 <- sqrt(var.v[1])
    se.theta2 <- sqrt(var.v[2])
    se.theta3 <- sqrt(var.v[3])
    se.ro1 <- sqrt(var.v[4])
    se.ro2 <- sqrt(var.v[5])
    se.ro3 <- sqrt(var.v[6])
  
    cat('OVER')
    cat("betaR1=",betaR1,"se=",se.betaR1,"betaR1/se=",betaR1/se.betaR1,"p-value=",2*(1-pnorm(abs(betaR1/se.betaR1))),'\n')
    cat("gamma1=",gamma1,"se=",se.gamma1,"gamma1/se=",gamma1/se.gamma1,"p-value=",2*(1-pnorm(abs(gamma1/se.gamma1))),'\n')
    cat("betaR2=",betaR2,"se=",se.betaR2,"betaR2/se=",betaR2/se.betaR2,"p-value=",2*(1-pnorm(abs(betaR2/se.betaR2))),'\n')
    cat("gamma2=",gamma2,"se=",se.gamma2,"gamma2/se=",gamma2/se.gamma2,"p-value=",2*(1-pnorm(abs(gamma2/se.gamma2))),'\n')
    cat("betaT=",betaT,"se=",se.betaT,"betaT/se=",betaT/se.betaT,"p-value=",2*(1-pnorm(abs(betaT/se.betaT))),'\n')
  
    cat("theta1=",theta1,"se=",se.theta1,"theta1/se=",theta1/se.theta1,"p-value=",2*(1-pnorm(abs(theta1/se.theta1))),'\n')
    cat("theta2=",theta2,"se=",se.theta2,"theta2/se=",theta2/se.theta2,"p-value=",2*(1-pnorm(abs(theta2/se.theta2))),'\n')
    cat("theta3=",theta3,"se=",se.theta3,"theta3/se=",theta3/se.theta3,"p-value=",2*(1-pnorm(abs(theta3/se.theta3))),'\n')
  
    cat("ro1=",ro1,"se=",se.ro1,"ro1/se=",ro1/se.ro1,"p-value=",2*(1-pnorm(abs(ro1/se.ro1))),'\n')
    cat("ro2=",ro2,"se=",se.ro2,"ro2/se=",ro2/se.ro2,"p-value=",2*(1-pnorm(abs(ro2/se.ro2))),'\n')
    cat("ro3=",ro3,"se=",se.ro3,"ro3/se=",ro3/se.ro3,"p-value=",2*(1-pnorm(abs(ro3/se.ro3))),'\n')
  
    list(Diverge=Diverge,flag.para=flag.para,flag.var=flag.var,betaR1=betaR1,gamma1=gamma1,
       betaR2=betaR2,gamma2=gamma2,betaT=betaT,theta1=theta1,theta2=theta2,theta3=theta3,
       ro1=ro1,ro2=ro2,ro3=ro3,se.betaR1=se.betaR1,se.gamma1=se.gamma1,se.betaR2=se.betaR2,
       se.gamma2=se.gamma2,se.betaT=se.betaT,se.theta1=se.theta1,se.theta2=se.theta2,
       se.theta3=se.theta3,se.ro1=se.ro1,se.ro2=se.ro2,se.ro3=se.ro3)
}
}
