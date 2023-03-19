require(Rcpp)
require(RSpectra)
require(refund)
require(face)
require(splines)
require(MASS)
require(MFPCA)
require(mgcv)
require(Matrix)
require(fields)
require(rARPACK)
library(glmnet)
require(quadprog)
require(elasticnet)
require(pracma)
require(ggplot2)
require(gridExtra)
require(reshape2)

# Load functions
source("./R code/FLFM.R")
source("./R code/Validation_L.R")
source("./R code/predict.FLFM.R")
source("./R code/predict.MFPCA.R")
source("./R code/fbps.cov.R") # for smooth covariance matrices
source("./R code/Renew_beta.R") # for renew Beta
source("./R code/Renew_C.R") # for renew shared covariance fucntions
source("./R code/Renew_G_sigma.R") # for renew outcome-specific terms
source("./R code/Renew_G_gamma_sigma.R") # for renew outcome-specific terms


#########################################################
# Data: EEG data
#########################################################
load("./data/EEG.RData") 
data <- EEG
Nsubj <- length(data)
variables <- unique(data[[1]][,2])
J <- length(variables)
Y <- list()
for(j in 1:J){
  index <- which(data[[1]][,2]==variables[j])
  temp <- lapply(1:Nsubj, function(x){data[[x]][index,4]})
  Y[[j]] <- do.call("rbind", temp)
}
t <- ncol(Y[[1]])
argvals <- seq(0, 1, length=t) 
label <- c(rep(1,65),rep(0,44),rep(1,12),0)
# reorder electrodes
ord0 <- c(1,3,9,10,5,7,14,16,12,25,19,23,21,36,28,32,30,34,39,45,41,43,48,52,
          50,58,54,57,55,63,61,27,4,8,11,17,18,26,2,24,20,35,29,15,13,46,38,
          6,40,44,47,53,31,33,59,60,22,56,62,51,49,42,64,37)
alcho0 <- matrix(data[[25]][,4], nc=J, nr=t, byrow=F)
contr0 <- matrix(data[[82]][,4], nc=J, nr=t, byrow=F)
ord <- order(ord0)



#########################################################
# LFM
#########################################################
knots <- 30 #list of two vectors of knots or number of equidistant knots for all dimensions
p <- 3 #degrees of B-splines
m <- 2 #order of differencing penalty
tol <- 5*1e-3 #tolerance for iteration
Nmax <- 200  #number of maximal iteration
pve <- 0.95


# 1. Select the best value L by cross validation
# seqL = c(1,2,3,4,5)
# L_result <- Validation_L(Y=Y, Nsubj=Nsubj, argvals=argvals, seqL=seqL, sparse=FALSE, rho=0,
#                          kfold=5, knots=knots, p=p, m=m, tol=tol, Nmax=100)
# L <- L_result$L_min_se

# 2. Estimate parameters
L <- 2
LFM_fit <- FLFM(Y=Y, Nsubj=Nsubj, argvals=argvals, L=L, sparse=FALSE, rho=0.1, pve1=0.95, 
                pve2=0.95, knots=knots, p=p, m=m, tol=tol, Nmax=Nmax, fit=T)

# 3. Select the best value rho1 by cross validation
# initial <- list(C=LFM_fit$C, beta=LFM_fit$beta, G=LFM_fit$G, sigma2=LFM_fit$sigma2)
# seqrho <- seq(0.1, 0.2, 0.05)
# rho_result <- Validation_rho(Y=Y, Nsubj=Nsubj, argvals=argvals, L=L, seqrho=seqrho, kfold=5, 
#                              initial=initial, knots=knots, p=p, m=m, tol=tol, Nmax=100)
# rho <- rho_result$rho_min
# 4. Estimate parameters for sparse=TRUE
# LFM_sparse <- FLFM(Y=Y, Nsubj=Nsubj, argvals=argvals, L=L, sparse=TRUE, rho=0.3, 
#                    pve1=0.95, pve2=0.90, initial=initial, knots=knots, p=p, m=m, tol=tol, Nmax=Nmax, fit=T)


LFM_score <- t(rbind(LFM_fit$xi, LFM_fit$zeta))
# LFM_score <- t(rbind(LFM_sparse$xi, LFM_sparse$zeta))
pred_LFM <- c()
for(i in 1:Nsubj){
  Y_train <- LFM_score[-i,]
  label_train <- label[-i]
  Y_test <- LFM_score[i,]
  data_train <- data.frame(group=label_train, data=as.matrix(Y_train))
  data_test <- data.frame(data=t(as.matrix(Y_test)))
  
  # Fit model
  cv.lasso <- cv.glmnet(x=Y_train, y=label_train, alpha=1, family="binomial", nfolds=nrow(Y_train))
  mylogit <- glmnet(x=Y_train, y=label_train, family = "binomial", alpha=1, lambda=cv.lasso$lambda.min)
  fitted.results <- predict(mylogit,newx=t(as.matrix(Y_test)), type='response')
  pred_LFM <- c(pred_LFM, fitted.results)
  print(i)
}
label_pre <- (pred_LFM>0.5)
# Classification accuracy
LFM_acc <- sum((label_pre-label)==0)/Nsubj
LFM_acc

#=============================
# Plot the shared covariance
#=============================
C <- LFM_fit$C
C1 <- C[[1]]
C2 <- C[[2]]

timep = 1:256
x = rep(timep, 256)
y = sort(rep(timep, 256))
df1 <- data.frame(x=x, y=y, Value=as.vector(as.matrix(C1)), label=rep("C1",256*256))
df2 <- data.frame(x=x, y=y, Value=as.vector(as.matrix(C2)), label=rep("C2",256*256))
df <- rbind(df1,df2)

ggplot(df, aes(x=x, y=y, z=Value, fill=Value)) + 
  geom_tile() + coord_equal() + geom_contour(color = "white", alpha = 0.5) + 
  facet_grid(.~label) + theme_bw() +
  theme(strip.text.x=element_blank(),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16)) +
  scale_fill_distiller(palette="Spectral", na.value="white") +
  xlab("Time(s)") + ylab("Time(s)") + scale_x_continuous(breaks = c(1,64,128,192,256), labels=c(0, 0.25, 0.5,0.75, 1)) + 
  scale_y_continuous(breaks = c(1,64,128,192,256), labels=c(0, 0.25, 0.5, 0.75, 1))


#=============================
# plot beta matrix
#=============================
beta <- LFM_fit$beta
beta <- beta[ord,]
colnames(beta) <- 1:2
rownames(beta) <- 1:64
longData <- melt(beta)

my.labs=c(expression(beta[.][1]),expression(beta[.][2]))
ggplot(longData, aes(x=rep(1:64,2), y=rep(1:2,each=64))) + 
  geom_raster(aes(fill=value)) + 
  scale_fill_gradient2(low="black",mid="white",high="red",midpoint=0) +
  labs(x="Electrodes", y=" ", title=" ") + labs(fill="Value") +
  scale_y_discrete(limits=c("beta1","beta2"),
                   labels=my.labs) + scale_x_continuous(breaks = c(1, 20, 40, 64)) +
  theme_bw() + theme(axis.text.x=element_text(size=16, angle=0, vjust=0.3),
                     axis.text.y=element_text(size=16), 
                     axis.title.x = element_text(size = 18),
                     axis.title.y = element_text(size = 22),
                     axis.title=element_text(size=20))


#############################################
# MFPCA
#############################################
sim <- multiFunData(lapply(1:J, function(x){funData(argvals=argvals, X=Y[[x]])}))
# MFPCA based on univariate FPCA
M <- 80
utype <- list(type="uFPCA", pve=pve)
MFPCA_fit <- MFPCA(sim, M=M, uniExpansions=list(utype,utype,utype,utype,utype,utype,utype,utype,utype,utype,
                                                utype,utype,utype,utype,utype,utype,utype,utype,utype,utype,
                                                utype,utype,utype,utype,utype,utype,utype,utype,utype,utype,
                                                utype,utype,utype,utype,utype,utype,utype,utype,utype,utype,
                                                utype,utype,utype,utype,utype,utype,utype,utype,utype,utype,
                                                utype,utype,utype,utype,utype,utype,utype,utype,utype,utype,
                                                utype,utype,utype,utype), fit=TRUE)
MFPCA_score <- MFPCA_fit$scores 

pred_MFPCA <- c()
for(i in 1:Nsubj){
  Y_train <- MFPCA_score[-i,]
  label_train <- as.factor(label[-i])
  Y_test <- MFPCA_score[i,]
  data_train <- data.frame(group=label_train, data=as.matrix(Y_train))
  data_test <- data.frame(data=t(Y_test))
  
  # Fit model
  cv.lasso <- cv.glmnet(x=Y_train, y=label_train, alpha=1, family="binomial", nfolds=nrow(Y_train))
  mylogit <- glmnet(x=Y_train, y=label_train, family = "binomial", alpha=1, lambda=cv.lasso$lambda.min)
  fitted.results <- predict(mylogit, newx=t(as.matrix(Y_test)), type='response')
  pred_MFPCA <- c(pred_MFPCA, fitted.results)
  print(i)
}
label_pre <- (pred_MFPCA>0.5)
# Classification accuracy
MFPCA_acc <- sum((label_pre-label)==0)/Nsubj
MFPCA_acc




