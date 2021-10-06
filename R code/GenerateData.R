library(MASS)

GenerateData <- function(Nsubj=100, argval=101, J=20, L=3, FICC=1/3, SNR=5){
  ############################################################ 
  #Arguments:
  #Nsubj -- number of subjects, default to 100
  #argval - number of points at which functions are observed; default to 101
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
  mu[,,1:5]<- matrix(rep(3*sin(2*pi*t), Nsubj), byrow=TRUE, nr=Nsubj)
  mu[,,6:10] <- matrix(rep(3*cos(2*pi*t), Nsubj), byrow=TRUE, nr=Nsubj)
  mu[,,11:15] <- matrix(rep(3*(t-1)^2, Nsubj), byrow=TRUE, nr=Nsubj)
  mu[,,16:20] <- matrix(rep(3*(t^2-0.5*t), Nsubj), byrow=TRUE, nr=Nsubj)
  # Put mu in a list form
  mean <- lapply(1:J, function(x){mu[,,x]})
  ############################################################ 
  # Common part
  ############################################################ 
  phi_1 <- cbind(sqrt(2)*sin(2*pi*t), sqrt(2)*cos(4*pi*t), sqrt(2)*sin(4*pi*t))
  phi_2 <- cbind(sqrt(2)*cos(pi*t), sqrt(2)*sin(3*pi*t), sqrt(2)*sin(5*pi*t))
  phi_3 <- cbind(sqrt(2)*sin(pi*t), sqrt(2)*sin(2*pi*t), sqrt(2)*sin(3*pi*t))
  gamma_1 <- diag( c(10, 10, 10) )
  gamma_2 <- diag( c(7, 7, 6) )
  gamma_3 <- diag( c(4, 3, 3) )
  gamma <- list(c(10, 10, 10), c(7, 7, 6), c(4, 3, 3))
  if(FICC==1/3){
    gamma_1 <- gamma_1/2
    gamma_2 <- gamma_2/2
    gamma_3 <- gamma_3/2
    gamma <- list(c(10, 10, 10)/2, c(7, 7, 6)/2, c(4, 3, 3)/2)
  }
  # # Covariance
  C_1 <- phi_1 %*% gamma_1 %*% t(phi_1)
  C_2 <- phi_2 %*% gamma_2 %*% t(phi_2)
  C_3 <- phi_3 %*% gamma_3 %*% t(phi_3)
  
  # Generate U_il
  xi <- lapply(1:L, function(x){
    do.call("cbind",lapply(1:length(gamma[[x]]),function(y){rnorm(n=Nsubj, mean=0, sd=sqrt(gamma[[x]][y]))}))  })
  U_1 <- xi[[1]] %*% t(phi_1)
  U_2 <- xi[[2]] %*% t(phi_2)
  U_3 <- xi[[3]] %*% t(phi_3)

  #Generate beta_jl
  #beta_1 <- c(0.4, 0.3, 0.3, 0.4, rep(0,3), 0.2, 0.4, 0.2, 0.3,rep(0,3), 0.2, -0.3, 0.2, rep(0,3))
  #beta_2 <- c(0.3, -0.2, 0, 0, 0.4, 0.3, 0.3, 0.2, rep(0,3), 0.4, 0.2, 0.2, -0.2, 0.2, 0, 0.3, 0.2, 0.2)
  #beta_3 <- c(0.2, rep(0,2), -0.3, -0.4, 0, 0, 0, 0.3, 0.2, -0.4, 0.2, 0.3, -0.3, rep(0,3), -0.2, 0.4, 0)
  beta_1 <- c(0.4, 0, 0.3, 0.1, -0.5, 0, 0, -0.2, -0.2, 0.1,0, 0, 0.4, 0.3, -0.2, 0.3, 0, 0, 0.1, -0.1)
  beta_2 <- c(0.3, sqrt(0.03), -0.1, 0.3, 0, -0.1, 0.1, 0.4, 0.2, 0.2, -sqrt(0.03), -0.3, 0, 0.3, 0.2, -0.3, sqrt(0.05), -0.3, 0.2, 0)
  beta_3 <- c(0.2, -0.2, rep(0,3), sqrt(0.1), sqrt(0.1), 0.4, 0, -0.2, -0.2, 0.3, 0.2, 0, 0.3, 0.1, 0, 0.4, 0, 0.3)
  beta <- cbind(beta_1, beta_2, beta_3)
  W <- lapply(1:J, function(x){
    beta_1[x]*U_1 + beta_2[x]*U_2 + beta_3[x]*U_3
  })

  ############################################################ 
  # Self part: Joint/Individual = FICC / (1-FICC)
  ############################################################ 
  psi <- cbind(rep(1,argval), sqrt(3)*(2*t-1), sqrt(5)*(6*t^2-6*t+1))
  delta <- diag(c(0.5, 0.5^2, 0.5^2))
  G <- psi %*% delta %*% t(psi)
  lambda <- (sum(gamma_1)*beta_1^2+sum(gamma_2)*beta_2^2+sum(gamma_3)*beta_3^2) / ( (FICC/(1-FICC)) * sum(delta))
  gamma2 <- lapply(1:J, function(x){lambda[x]*c(0.5, 0.5^2, 0.5^2)})
  eta <- lapply(1:J, function(x){
    do.call("cbind",lapply(1:length(gamma2[[x]]),function(y){rnorm(n=Nsubj, mean=0, sd=sqrt(gamma2[[x]][y]))}))  })
  V <- lapply(1:J, function(x){eta[[x]] %*% t(psi)})


  ############################################################ 
  # Generate random error terms
  ############################################################ 
  sigma_sq <- (sum(gamma_1)*beta_1^2+sum(gamma_2)*beta_2^2+sum(gamma_3)*beta_3^2 + lambda*sum(delta))/SNR
  epsilon <- lapply(1:J, function(x) {
    matrix(rnorm(Nsubj*argval, mean=0, sd=sqrt(sigma_sq[x])), nr=Nsubj)} )

  ###########################################################################################
  #combine to get data Y_ij(t) = mu_j(t) + V_ij(t) + W_ij(t)  + epsilon_ij(t)
  ###########################################################################################
  Y <- lapply(1:J, function(x) {mu[[x]] + W[[x]] + V[[x]] + epsilon[[x]]})
  
  return(list(Y=Y, mean=mean, W=W, V=V, U=list(U_1, U_2, U_3), beta=beta, xi=xi, eta=eta,
              C=list(C_1, C_2, C_3), G=G, lambda=lambda, sigma2=sigma_sq))
}



