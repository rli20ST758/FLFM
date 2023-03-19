library(MASS)

GenerateData_L5 <- function(Nsubj=100, argval=101, J=60, L=5, FICC=1/3, SNR=5){
  ############################################################ 
  #Arguments:
  #Nsubj -- number of subjects, default to 100
  #argval - number of points at which functions are observed, default to 101
  #J -- number of dimensions
  #L -- number of reduced dimensions
  #FICC -- functional intra-class correlation; joint variation / total variation
  #SNR -- the signal to noise ratio, default to 5
  ############################################################ 
  # Vector of subject visits
  subjID <- sort(rep(1:Nsubj, J))
  # Time points at which functional responses are collected
  t <- seq(0, 1, length = argval)
  
  ############################################################ 
  # Define mean fucntions for each variate
  ############################################################ 
  mu <- array(0, dim=c(Nsubj, argval, J))
  mu[,,1:15] <- matrix(rep(3*sin(2*pi*t), Nsubj), byrow=TRUE, nr=Nsubj)
  mu[,,16:30] <- matrix(rep(3*cos(2*pi*t), Nsubj), byrow=TRUE, nr=Nsubj)
  mu[,,31:45] <- matrix(rep(3*(t-1)^2, Nsubj), byrow=TRUE, nr=Nsubj)
  mu[,,46:60] <- matrix(rep(3*(t^2-0.5*t), Nsubj), byrow=TRUE, nr=Nsubj)
  # Put mu in a list form
  mean <- lapply(1:J, function(x){mu[,,x]})
  
  ############################################################ 
  # Common part
  ############################################################
  phi_1 <- cbind(sqrt(2)*sin(2*pi*t), sqrt(2)*cos(4*pi*t), sqrt(2)*sin(4*pi*t))
  phi_2 <- cbind(sqrt(2)*sin(pi*t), sqrt(2)*sin(2*pi*t), sqrt(2)*sin(3*pi*t))
  phi_3 <- cbind(sqrt(2)*cos(pi*t), sqrt(2)*sin(3*pi*t), sqrt(2)*sin(5*pi*t))
  phi_4 <- cbind(sqrt(2)*sin(pi*t), sqrt(2)*cos(pi*t), sqrt(2)*sin(3*pi*t))
  phi_5 <- cbind(sqrt(2)*sin(2*pi*t), sqrt(2)*sin(4*pi*t), sqrt(2)*sin(3*pi*t))
  gamma_1 <- diag( c(6, 6, 6) )
  gamma_2 <- diag( c(6, 5, 5) )
  gamma_3 <- diag( c(5, 5, 4) )
  gamma_4 <- diag( c(4, 4, 4) ) 
  gamma_5 <- diag( c(4, 3, 3) ) 
  gamma <- list(c(6, 6, 6),c(6, 5, 5), c(5, 5, 4), c(4, 4, 4), c(4, 3, 3))
  if(FICC==1/3){
    gamma_1 <- gamma_1/2
    gamma_2 <- gamma_2/2
    gamma_3 <- gamma_3/2
    gamma_4 <- gamma_4/2
    gamma_5 <- gamma_5/2
    gamma <- list(c(6, 6, 6)/2,c(6, 5, 5)/2, c(5, 5, 4)/2, c(4, 4, 4)/2, c(4, 3, 3)/2)
  }
  # Covariance
  C_1 <- phi_1 %*% gamma_1 %*% t(phi_1)
  C_2 <- phi_2 %*% gamma_2 %*% t(phi_2)
  C_3 <- phi_3 %*% gamma_3 %*% t(phi_3)
  C_4 <- phi_4 %*% gamma_4 %*% t(phi_4)
  C_5 <- phi_5 %*% gamma_5 %*% t(phi_5)
  # Generate U_il
  xi <- lapply(1:L, function(x){
    do.call("cbind",lapply(1:length(gamma[[x]]),function(y){rnorm(n=Nsubj, mean=0, sd=sqrt(gamma[[x]][y]))}))  })
  U_1 <- xi[[1]] %*% t(phi_1)
  U_2 <- xi[[2]] %*% t(phi_2)
  U_3 <- xi[[3]] %*% t(phi_3)
  U_4 <- xi[[4]] %*% t(phi_4)
  U_5 <- xi[[5]] %*% t(phi_5)
  #Generate beta_jl
  beta <- matrix(0, nc=L, nr=J)
  for(l in 1:L){
    set.seed(l)
    if (l %in% c(1,2)){
      beta[1:20,l] <- runif(20,min=0.5,max=1) 
      beta[,l] <- beta[,l]/norm(beta[,l],"2") 
    } else if(l==3){
      beta[21:40,l] <- runif(20,min=0.5,max=1) 
      beta[,l] <- beta[,l]/norm(beta[,l],"2")
    } else {
      beta[41:60,l] <- runif(20,min=0.5,max=1) 
      beta[,l] <- beta[,l]/norm(beta[,l],"2") 
    }
  }
  beta[,1:2] <- gramSchmidt(beta[,1:2])$Q
  beta[,4:5] <- gramSchmidt(beta[,4:5])$Q
  
 
  W <- lapply(1:J, function(x){
    beta[x,1]*U_1 + beta[x,2]*U_2 + beta[x,3]*U_3 + beta[x,4]*U_4 + beta[x,5]*U_5
  })
  ############################################################ 
  # Self part: Joint/Individual = FICC / (1-FICC)
  ############################################################ 
  psi <- cbind(rep(1,argval), sqrt(3)*(2*t-1), sqrt(5)*(6*t^2-6*t+1))
  delta <- diag(c(0.5, 0.5^2, 0.5^2))
  G <- psi %*% delta %*% t(psi)
  lambda <- rowSums(do.call(cbind, lapply(1:L, function(x){sum(gamma[[x]])*beta[,x]^2}))) / ( (FICC/(1-FICC)) * sum(delta))
  gamma2 <- lapply(1:J, function(x){lambda[x]*c(0.5, 0.5^2, 0.5^2)})
  eta <- lapply(1:J, function(x){
    do.call("cbind",lapply(1:length(gamma2[[x]]),function(y){rnorm(n=Nsubj, mean=0, sd=sqrt(gamma2[[x]][y]))}))  })
  
  V <- lapply(1:J, function(x){eta[[x]] %*% t(psi)})
  
  ############################################################ 
  # Generate random error terms
  ############################################################ 
  sigma_sq <- (rowSums(do.call(cbind, lapply(1:L, function(x){sum(gamma[[x]])*beta[,x]^2}))) + lambda*sum(delta)) / SNR
  epsilon <- lapply(1:J, function(x) {
    matrix(rnorm(Nsubj*argval, mean=0, sd=sqrt(sigma_sq[x])), nr=Nsubj)} )
  
  ###########################################################################################
  #combine to get data Y_ij(t) = mu_j(t) + V_ij(t) + W_ij(t)  + epsilon_ij(t)
  ###########################################################################################
  Y <- lapply(1:J, function(x) {mu[[x]] + W[[x]] + V[[x]] + epsilon[[x]]})
  
  return(list(Y=Y, mean=mean, W=W, V=V, U=list(U_1,U_2,U_3,U_4,U_5), beta=beta, xi=xi, eta=eta,  
              C=list(C_1,C_2,C_3,C_4,C_5), G=G,  lambda=lambda, sigma2=sigma_sq) )
}



