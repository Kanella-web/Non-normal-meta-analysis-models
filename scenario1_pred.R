###################### Number of studies k=14 ####################
####################### Scenario 1: Normal data, k=14, mu=0, tau2= 0.12 #######################
# install.packages("doParallel")
# install.packages("SimDesign")
# install.packages("metafor")
# install.packages("Overlapping")
# install.packages("meta")
library(doParallel)
library(SimDesign)
library(metafor)
library(overlapping)
library(meta)
library(dplyr)

N.sim = 1000
df <- list()
dich_data <- function(k, tau2, mu){
  mati <- matrix(1:k,nrow = k, ncol=2) 
  
  for (i in 1:N.sim){
    ## random effects distribution 
    thetai <- rnorm(k, mu, sqrt(tau2))
    
    ni <-sample(50:500,k,replace=T)
    nti =nci =ni
    ###Obtain the total number of events ci for the control group from a Binomial(nic, pic) distribution
    ##with pic ~ Uniform (0.05, 0.65)
    pci <- runif(k, 0.05, 0.65)
    ci <- rbinom(k, nci, pci)
    ##Obtain the total number of events ti for the treatment group from a Binomial(nit, pit) distribution
    
    pti <- pci*exp(thetai)/(1-pci+pci*exp(thetai))
    ti <- rbinom(k,nti,pti)
    
    dati <- cbind(ci, ti, nci, nti,thetai)
    mati <- cbind(mati, dati) 
  }
  df[[i]] <- mati
  return(df[[i]])
}
######### RUN THE CODE ###########
set.seed(210)
m <- dich_data(14,0.12,0)  
m <- m[,-c(1,2)]

q <- list()
l <- seq(1, ncol(m), by=5)
for(i in 1:length(l)){
  for(j in l[i]){
    q[[i]] <- cbind(m[,j], m[ ,j+1], m[ , j+2], m[ ,j+3], m[ , j+4])
  }
}
p <- list()
for(i in 1:N.sim){
  p[[i]] = as.data.frame(q[[i]])
  colnames(p[[i]]) = c("ci", "ti","nci","nti", "true_eff")
  
}
####### CHECK IF THERE ARE STUDIES WITH 0 IN TWO ARMS ####
N.studies <- nrow(p[[i]])
d <- vector("list", N.sim)
for(i in 1:N.sim){
  d[[i]] <- vector("numeric", N.studies)
  for(j in 1:N.studies){
    if(p[[i]]$ti[j] == 0 & p[[i]]$ci[j] == 0){  
      d[[i]][j] <- 1
    }
    else{
      d[[i]][j] <- 0
    }
  }
}
d

count <- c()
for(i in 1:N.sim){
  count[i] <- sum(d[[i]])
}
count
studies_zero_2arms <- which(count == 1)
#studies_zero_2arms

mu = 0
tau2 = 0.12
tau = sqrt(tau2)
###############  BAYESIAN MODELS ##############
############### NORMAL MODEL #################
############# BINOMIAL-NORMAL(Unif) #############

N.sim = 1000

writeLines("
  model{
  for(i in 1:ns){ 
  m[i] ~ dnorm(0,.0001) 
  
  r1[i] ~ dbin(p1[i], n1[i])
  logit(p1[i]) <- m[i] 
  
  r2[i] ~ dbin(p2[i], n2[i])
  logit(p2[i]) <- m[i]  + delta12[i]
  
  delta12[i] ~ dnorm(mu, tau1)

  }
 
 mu ~ dnorm(0,.0001) # vague prior for treatment effect
 
 tau1 = 1/tau_sqr
 tau_sqr = tau*tau
 tau ~ dunif(0,10)
 
 delta_new ~ dnorm(mu, tau1)
 
} 
           ", con="MA.model.binU.txt")

modfile = 'MA.model.binU.txt'

datU <- list()
for(i in 1:N.sim){
  datU[[i]] <-list(ns = nrow(p[[i]]),
                   r1 = p[[i]]$ci,
                   r2 = p[[i]]$ti,
                   n1 = p[[i]]$nci,
                   n2 = p[[i]]$nti
  ) }

detectCores()
cl = makeCluster(19)
registerDoParallel(cl)
resultsNU <- list()

resultsNU = foreach(i = 1:N.sim,.packages = c("rjags", "R2jags","tidyverse")) %dopar% {
  set.seed(2101)
  run.model = jags(
    data = datU[[i]],
    inits = NULL,
    parameters.to.save = c(
      "mu",
      "tau",
      "delta12",
      "delta_new"
    ),
    n.chains = 2,
    n.iter = 50000,
    
    n.burnin = 10000,
    DIC = T,
    model.file = modfile
  )
}
stopCluster(cl)

resultsU <- list()
posteriorU <- list()
postU <- list()
for(i in 1:N.sim){
  resultsU[[i]] <- as.data.frame(resultsNU[[i]]$BUGSoutput$summary) 
  posteriorU[[i]] = as.data.frame(resultsNU[[i]]$BUGSoutput$sims.matrix)
  postU[[i]] <- posteriorU[[i]]$mu
}

mu_Normaljags.binU <- c()
LBmu.jags.binU <- c()
UBmu.jags.binU <- c()
pred.binU <- c()
LBpred.binU <- c()
UBpred.binU <- c()
tau_Normaljags.binU <- c()
LBtau.jags.binU <- c()
UBtau.jags.binU <- c()
tau2_Normaljags.binU <- c()
LBtau2.jags.binU <- c()
UBtau2.jags.binU <- c()
rel_effU <- list()
LB_rel_effU <- list()
UB_rel_effU <- list()
sd_rel_effU <- list()
precN_muU <- c()
precN_predU <- c()
precN_tau2U <- c()
precN_tauU <- c()
listaU <- list()
Rhat_muNU <- c()
Rhat_predNU <- c()
Rhat_tauNU <- c()
Rhat_deltaNU <- list()
f1U <- list()
m1U <- list()
f3U <- list()
m3U <- list()
f2U <- list()
m2U <- list()
fdU <- list()
mdU <- list()
for(i in 1:N.sim){
  f1U[[i]] <- grepl("mu", row.names(resultsU[[i]]))
  m1U[[i]] <- resultsU[[i]][f1U[[i]],]
  mu_Normaljags.binU[i] <- m1U[[i]]$`50%` 
  LBmu.jags.binU[i] <- m1U[[i]]$`2.5%`
  UBmu.jags.binU[i] <- m1U[[i]]$`97.5%`
  Rhat_muNU[i] <- m1U[[i]]$Rhat
  precN_muU[i] <-  UBmu.jags.binU[i] - LBmu.jags.binU[i]
  
  f3U[[i]] <- grepl("delta_new", row.names(resultsU[[i]]))
  m3U[[i]] <- resultsU[[i]][f3U[[i]],]
  pred.binU[i] <- m3U[[i]]$`50%` 
  LBpred.binU[i] <- m3U[[i]]$`2.5%`
  UBpred.binU[i] <- m3U[[i]]$`97.5%`
  Rhat_predNU[i] <- m3U[[i]]$Rhat
  precN_predU[i] <-  UBpred.binU[i] - LBpred.binU[i]
  
  
  f2U[[i]] <- grepl("tau", row.names(resultsU[[i]]))
  m2U[[i]] <- resultsU[[i]][f2U[[i]],]
  tau_Normaljags.binU[i] <- m2U[[i]]$`50%`
  LBtau.jags.binU[i] <- m2U[[i]]$`2.5%`
  UBtau.jags.binU[i] <- m2U[[i]]$`97.5%`
  tau2_Normaljags.binU[i] <- (m2U[[i]]$`50%`)^2
  LBtau2.jags.binU[i] <- (m2U[[i]]$`2.5%`)^2
  UBtau2.jags.binU[i] <- (m2U[[i]]$`97.5%`)^2
  Rhat_tauNU[i] <- m2U[[i]]$Rhat
  precN_tau2U[i] <-  UBtau2.jags.binU[i] - LBtau2.jags.binU[i]
  precN_tauU[i] <-  UBtau.jags.binU[i] - LBtau.jags.binU[i]
  
  fdU[[i]] <- grepl("delta12", row.names(resultsU[[i]]))
  mdU[[i]] <- resultsU[[i]][fdU[[i]],]
  rel_effU[[i]] <- mdU[[i]]$`50%`
  LB_rel_effU[[i]] <- mdU[[i]]$`2.5%`
  UB_rel_effU[[i]] <- mdU[[i]]$`97.5%`
  sd_rel_effU[[i]] <- mdU[[i]]$sd
  Rhat_deltaNU[[i]] <- mdU[[i]]$Rhat
  listaU[[i]] <- cbind.data.frame(rel_effU[[i]],sd_rel_effU[[i]],LB_rel_effU[[i]], UB_rel_effU[[i]], Rhat_deltaNU[[i]])
}


RhatsNU <- cbind.data.frame(mu_Normaljags.binU, LBmu.jags.binU, UBmu.jags.binU, tau2_Normaljags.binU,
                            LBtau2.jags.binU, UBtau2.jags.binU,Rhat_muNU, 
                            tau_Normaljags.binU, LBtau.jags.binU, UBtau.jags.binU,
                            pred.binU, LBpred.binU, UBpred.binU,
                            Rhat_tauNU,Rhat_predNU,
                            precN_muU, precN_tau2U,precN_tauU, precN_predU)

############REMOVE Rhats > 1.05 ##############

condition1U <- which(RhatsNU$Rhat_muNU > 1.05)
condition2U <- which(RhatsNU$Rhat_tauNU > 1.05)

########## check if condition1 == condition2 ##########
dist.condNU = c(condition1U,condition2U)
dist.condNU = unique(dist.condNU)

RhatN_outU = round((length(dist.condNU)/N.sim), 4)


if (length(dist.condNU)== 0) {
  RhatsNU <- RhatsNU
  N.sim <- nrow(RhatsNU)
  listaNormU <- listaU
  pNU <- p
} else {
  RhatsNU <- RhatsNU[-dist.condNU, ]
  listaNormU <- listaU[-dist.condNU]
  N.sim <- nrow(RhatsNU)
  pNU <- p[-dist.condNU]
}


count_NU = list()
for(i in length(listaNormU)){
  count_NU[[i]] = which(listaNormU[[i]]$`Rhat_deltaNU[[i]]` > 1.05 )
  tell_meNU = which(count_NU[[i]] != 0)
}

tell_meNU

if(length(tell_meNU) == 0){
  RhatsNU <- RhatsNU
  N.sim <- nrow(RhatsNU)
  listaNormU <- listaNormU
  pNU <- pNU
  RhatN_outU <- RhatN_outU
} else {
  RhatsNU <- RhatsNU[-tell_meNU, ]
  listaNormU <- listaNormU[-tell_meNU]
  N.sim = nrow(RhatsNU)
  pNU <- pNU[-tell_meNU]
  RhatN_outU <- RhatN_outU + ((length(tell_meNU))/N.sim)
}

RhatN_outU

########## AVERAGE RELATIVE BIAS OF mu ##########
rel_biasNormaljags.binU <- c()
for (i in 1:N.sim){
  rel_biasNormaljags.binU[i] <- ( RhatsNU$mu_Normaljags.binU[i] - mu) 
}
avg_rel_biasNormaljags.binU <- round(mean(rel_biasNormaljags.binU), 2)

########## AVERAGE ABSOLUTE BIAS OF mu ##########
abs_biasNormaljags.binU <- c()
for (i in 1:N.sim){
  abs_biasNormaljags.binU[i] <- abs( RhatsNU$mu_Normaljags.binU[i] - mu) 
}
avg_abs_biasNormaljags.binU <- round(mean(abs_biasNormaljags.binU), 2)

############### MSE OF mu  ##################
MSE_Normaljags.binU <- c()
for (i in 1:N.sim){
  MSE_Normaljags.binU[i] <- RMSE(estimate =RhatsNU$mu_Normaljags.binU[i],
                                 parameter = mu,
                                 type = "RMSE",
                                 MSE = TRUE,           
                                 percent = FALSE,
                                 unname = FALSE)
}
avg_MSE_Normaljags.binU <-  round(mean(MSE_Normaljags.binU),4)

########## COVERAGE OF mu ##########
interval_contains_true_mean <- function(i) { 
  ( mu >= RhatsNU$LBmu.jags.binU[i]) && ( mu <= RhatsNU$UBmu.jags.binU[i])
}
interval_contains_true_mean(i)

coverageNormaljags.binU <- c()
for (i in 1:N.sim){0
  if (interval_contains_true_mean(i) == TRUE){
    coverageNormaljags.binU[i] = 1
  }
  else {
    coverageNormaljags.binU[i]=0
  }
  print(coverageNormaljags.binU[i])
}
coverageNormaljags.binU
coverNormaljags.binU <- round(mean(coverageNormaljags.binU)*100, 2)

################
#### AVERAGE BIAS OF tau2 ########
bias_tau2Normaljags.binU <- c()
for (i in 1:N.sim){
  bias_tau2Normaljags.binU[i] <- abs(RhatsNU$tau2_Normaljags.binU[i] - tau2)
}
avg_bias_tau2Normaljags.binU <- round(mean(bias_tau2Normaljags.binU),2)

pbias_tau2Normaljags.binU <- round(mean(bias_tau2Normaljags.binU / tau2),2)
########### MSE OF tau2 ################
MSE_Normaljags.tau2.binU <- c()
for (i in 1:N.sim){
  MSE_Normaljags.tau2.binU[i] <- RMSE(estimate =RhatsNU$tau2_Normaljags.binU[i],
                                      parameter = tau2,
                                      type = "RMSE",
                                      MSE = TRUE           
  )
}
avg_MSE_Normaljags.tau2.binU <- round(mean(MSE_Normaljags.tau2.binU),4)

nMSE_Normaljags.tau2.binU <- round(mean(MSE_Normaljags.tau2.binU / (tau2)^2 ),4)

############# COVERAGE OF tau2 ###################
interval_contains_true_mean <- function(i) { 
  ( tau2 >= RhatsNU$LBtau2.jags.binU[i]) && ( tau2 <=  RhatsNU$UBtau2.jags.binU[i])
}
interval_contains_true_mean(i)

coverage_tau2Normaljags.binU <- c()
for (i in 1:N.sim){
  if (interval_contains_true_mean(i) == TRUE){
    coverage_tau2Normaljags.binU[i] = 1
  }
  else  {
    coverage_tau2Normaljags.binU[i] = 0
  }
  print(coverage_tau2Normaljags.binU[i])
}
coverage_tau2Normaljags.binU
cover_tau2Normaljags.binU <- round(mean(coverage_tau2Normaljags.binU)*100, 2)


################
#### AVERAGE relative BIAS OF tau ########
rel_bias_tauNormaljags.binU <- c()
for (i in 1:N.sim){
  rel_bias_tauNormaljags.binU[i] <- (RhatsNU$tau_Normaljags.binU[i] - tau)
}
avg_rel_bias_tauNormaljags.binU <- round(mean(rel_bias_tauNormaljags.binU),2)

prel_bias_tauNormaljags.binU <- round(mean(rel_bias_tauNormaljags.binU / tau)*100,2)

#### AVERAGE absolute BIAS OF tau ########
abs_bias_tauNormaljags.binU <- c()
for (i in 1:N.sim){
  abs_bias_tauNormaljags.binU[i] <- abs(RhatsNU$tau_Normaljags.binU[i] - tau)
}
avg_abs_bias_tauNormaljags.binU <- round(mean(abs_bias_tauNormaljags.binU),2)

pabs_bias_tauNormaljags.binU <- round(mean(abs_bias_tauNormaljags.binU / tau)*100,2)

########### MSE OF tau ################
MSE_Normaljags.tau.binU <- c()
for (i in 1:N.sim){
  MSE_Normaljags.tau.binU[i] <- RMSE(estimate =RhatsNU$tau_Normaljags.binU[i],
                                     parameter = tau,
                                     type = "RMSE",
                                     MSE = TRUE           
  )
}
avg_MSE_Normaljags.tau.binU <- round(mean(MSE_Normaljags.tau.binU),4)

nMSE_Normaljags.tau.binU <- round(mean(MSE_Normaljags.tau.binU / tau2 ),4)

############# COVERAGE OF tau ###################
interval_contains_true_mean <- function(i) { 
  ( tau >= RhatsNU$LBtau.jags.binU[i]) && ( tau <=  RhatsNU$UBtau.jags.binU[i])
}
interval_contains_true_mean(i)

coverage_tauNormaljags.binU <- c()
for (i in 1:N.sim){
  if (interval_contains_true_mean(i) == TRUE){
    coverage_tauNormaljags.binU[i] = 1
  }
  else  {
    coverage_tauNormaljags.binU[i] = 0
  }
  print(coverage_tauNormaljags.binU[i])
}
coverage_tauNormaljags.binU
cover_tauNormaljags.binU <- round(mean(coverage_tauNormaljags.binU)*100, 2)


######################## ABSOLUTE BIAS OF RELATIVE EFFECTS ##################
abs_bias_relNU <- list()
mean_abs_bias_relNU <- list()
for(i in 1:N.sim){
  for(j in 1:nrow(pNU[[i]])){
    abs_bias_relNU[[i]] <- round(abs( listaNormU[[i]]$`rel_effU[[i]]` - pNU[[i]]$true_eff  ),2)
    mean_abs_bias_relNU[[i]] <- mean(abs_bias_relNU[[i]])
  }
}

################## mean bias of the mean biases of each data set for relative effects  ###################
mean_abs_bias_of_releffU <- unlist(mean_abs_bias_relNU)
mean_abs_bias_of_mean_abs_bias_relNU  <- round(mean(mean_abs_bias_of_releffU),2)

######################## RELATIVE BIAS OF RELATIVE EFFECTS ##################
rel_bias_relNU <- list()
mean_rel_bias_relNU <- list()
for(i in 1:N.sim){
  for(j in 1:nrow(pNU[[i]])){
    rel_bias_relNU[[i]] <- round(( listaNormU[[i]]$`rel_effU[[i]]` - pNU[[i]]$true_eff  ),2)
    mean_rel_bias_relNU[[i]] <- mean(rel_bias_relNU[[i]])
  }
}

################## mean bias of the mean biases of each data set for relative effects  ###################
mean_rel_bias_of_releffU <- unlist(mean_rel_bias_relNU)
mean_rel_bias_of_mean_rel_bias_relNU  <- round(mean(mean_rel_bias_of_releffU),2)

########################### FOR THE OVERLAPPING ######################
dataNU <- list()
overlapN1U <- list()
overlapNU <- list()
for(i in 1:N.sim){
  dataNU[[i]] <- list(NORMAL =listaNormU[[i]]$`rel_effU[[i]]` , True =pNU[[i]]$true_eff )
  overlapN1U[[i]] <- overlap(dataNU[[i]], type = "2")
  overlapNU[[i]] <- overlapN1U[[i]]$OV
}
overlapNormU <- unlist(overlapNU)
avg_overlapNormU <- round(mean(overlapNormU),2)

######################### SAVE THE OUTPUTS ####################
lista_rel_effU = cbind( mean_rel_bias_relNU, mean_rel_bias_of_mean_rel_bias_relNU,
                        mean_abs_bias_relNU, mean_abs_bias_of_mean_abs_bias_relNU) 


len_mu_NU = mean(RhatsNU$precN_muU)
len_PI_NU = mean(RhatsNU$precN_predU)

Norm_bin_jagsU <- cbind(RhatsNU , MSE_Normaljags.binU, MSE_Normaljags.tau2.binU,coverageNormaljags.binU,
                        coverage_tau2Normaljags.binU, rel_biasNormaljags.binU, abs_biasNormaljags.binU,
                        bias_tau2Normaljags.binU, pbias_tau2Normaljags.binU,
                        avg_rel_biasNormaljags.binU, avg_abs_biasNormaljags.binU, 
                        avg_bias_tau2Normaljags.binU,
                        avg_MSE_Normaljags.binU, avg_MSE_Normaljags.tau2.binU, nMSE_Normaljags.tau2.binU,
                        coverNormaljags.binU,cover_tau2Normaljags.binU,
                        
                        abs_bias_tauNormaljags.binU, pabs_bias_tauNormaljags.binU,
                        avg_abs_bias_tauNormaljags.binU,
                        
                        rel_bias_tauNormaljags.binU, avg_rel_bias_tauNormaljags.binU,prel_bias_tauNormaljags.binU,
                        avg_MSE_Normaljags.tau.binU, nMSE_Normaljags.tau.binU,
                        coverNormaljags.binU,cover_tauNormaljags.binU,
                        len_mu_NU, 
                        len_PI_NU ,
                        RhatN_outU)

Overalp_NormalU <- cbind(overlapNormU, avg_overlapNormU )
# 

write.csv(Norm_bin_jagsU, "C:\\Users\\Tianqi YU\\Desktop\\SIMULATION SCENARIOS FINAL\\NORMAL DATA\\NORMAL MODEL\\rev_UNIF_ResultsNormal_jags_Normal_bin_data_k=14_mu=0_tau2=0.12.csv",row.names=FALSE )  
write.csv(listaNormU, "C:\\Users\\Tianqi YU\\Desktop\\SIMULATION SCENARIOS FINAL\\NORMAL DATA\\NORMAL MODEL\\rev_UNIF_Relative_eff_ResultsNormal_jags_Normal_bin_data_k=14_mu=0_tau2=0.12.csv",row.names=FALSE )  
write.csv(lista_rel_effU, "C:\\Users\\Tianqi YU\\Desktop\\SIMULATION SCENARIOS FINAL\\NORMAL DATA\\NORMAL MODEL\\rev_UNIF_Bias_Relative_eff_ResultsNormal_jags_Normal_bin_data_k=14_mu=0_tau2=0.12.csv",row.names=FALSE )  
write.csv(Overalp_NormalU, "C:\\Users\\Tianqi YU\\Desktop\\SIMULATION SCENARIOS FINAL\\NORMAL DATA\\NORMAL MODEL\\rev_UNIF_Overalp_Normal_Normal_bin_data_k=14_mu=0_tau2=0.12.csv",row.names=FALSE )  



###############  BAYESIAN MODELS ##############
############### NORMAL MODEL #################
############# BINOMIAL-NORMAL(HN) #############

N.sim = 1000
mu = 0
tau2 = 0.12
tau = sqrt(tau2)

writeLines("
  model{
  for(i in 1:ns){ 
  m[i] ~ dnorm(0,.0001) 
  
  r1[i] ~ dbin(p1[i], n1[i])
  logit(p1[i]) <- m[i] 
  
  r2[i] ~ dbin(p2[i], n2[i])
  logit(p2[i]) <- m[i]  + delta12[i]
  
  delta12[i] ~ dnorm(mu, tau1)

  }
 
 mu ~ dnorm(0,.0001) # vague prior for treatment effect
 
 tau1 = 1/tau_sqr
 tau_sqr = tau*tau
 tau~dnorm(0,1)T(0,)
 
 delta_new ~ dnorm(mu, tau1)


} 
           ", con="MA.model.binHN.txt")

modfile = 'MA.model.binHN.txt'

datHN <- list()
for(i in 1:N.sim){
  datHN[[i]] <-list(ns = nrow(p[[i]]),
                    r1 = p[[i]]$ci,
                    r2 = p[[i]]$ti,
                    n1 = p[[i]]$nci,
                    n2 = p[[i]]$nti
  ) }

detectCores()
cl = makeCluster(19)
registerDoParallel(cl)
resultsNHN <- list()

resultsNHN = foreach(i = 1:N.sim,.packages = c("rjags", "R2jags","tidyverse")) %dopar% {
  set.seed(2102)
  run.model = jags(
    data = datHN[[i]],
    inits = NULL,
    parameters.to.save = c(
      "mu",
      "tau",
      "delta12",
      "delta_new"
    ),
    n.chains = 2,
    n.iter = 50000,
    
    n.burnin = 10000,
    DIC = T,
    model.file = modfile
  )
}
stopCluster(cl)

resultsHN <- list()
posteriorHN <- list()
postHN <- list()
for(i in 1:N.sim){
  resultsHN[[i]] <- as.data.frame(resultsNHN[[i]]$BUGSoutput$summary) 
  posteriorHN[[i]] = as.data.frame(resultsNHN[[i]]$BUGSoutput$sims.matrix)
  postHN[[i]] <- posteriorHN[[i]]$mu
}

mu_Normaljags.binHN <- c()
LBmu.jags.binHN <- c()
UBmu.jags.binHN <- c()
pred.binHN <- c()
LBpred.binHN <- c()
UBpred.binHN <- c()
tau_Normaljags.binHN <- c()
LBtau.jags.binHN <- c()
UBtau.jags.binHN <- c()
tau2_Normaljags.binHN <- c()
LBtau2.jags.binHN <- c()
UBtau2.jags.binHN <- c()
rel_effHN <- list()
LB_rel_effHN <- list()
UB_rel_effHN <- list()
sd_rel_effHN <- list()
precN_muHN <- c()
precN_tau2HN <- c()
precN_predHN <- c()
precN_tauHN <- c()
listaHN <- list()
Rhat_muNHN <- c()
Rhat_tauNHN <- c()
Rhat_deltaNHN <- list()
Rhat_predNHN <- c()
f1HN <- list()
m1HN <- list()
f3HN <- list()
m3HN <- list()
f2HN <- list()
m2HN <- list()
fdHN <- list()
mdHN <- list()
for(i in 1:N.sim){
  f1HN[[i]] <- grepl("mu", row.names(resultsHN[[i]]))
  m1HN[[i]] <- resultsHN[[i]][f1HN[[i]],]
  mu_Normaljags.binHN[i] <- m1HN[[i]]$`50%` 
  LBmu.jags.binHN[i] <- m1HN[[i]]$`2.5%`
  UBmu.jags.binHN[i] <- m1HN[[i]]$`97.5%`
  Rhat_muNHN[i] <- m1HN[[i]]$Rhat
  precN_muHN[i] <-  UBmu.jags.binHN[i] - LBmu.jags.binHN[i]
  
  f3HN[[i]] <- grepl("delta_new", row.names(resultsHN[[i]]))
  m3HN[[i]] <- resultsHN[[i]][f3HN[[i]],]
  pred.binHN[i] <- m3HN[[i]]$`50%` 
  LBpred.binHN[i] <- m3HN[[i]]$`2.5%`
  UBpred.binHN[i] <- m3HN[[i]]$`97.5%`
  Rhat_predNHN[i] <- m3HN[[i]]$Rhat
  precN_predHN[i] <-  UBpred.binHN[i] - LBpred.binHN[i]
  
  f2HN[[i]] <- grepl("tau", row.names(resultsHN[[i]]))
  m2HN[[i]] <- resultsHN[[i]][f2HN[[i]],]
  tau_Normaljags.binHN[i] <- m2HN[[i]]$`50%`
  LBtau.jags.binHN[i] <- m2HN[[i]]$`2.5%`
  UBtau.jags.binHN[i] <- m2HN[[i]]$`97.5%`
  tau2_Normaljags.binHN[i] <- (m2HN[[i]]$`50%`)^2
  LBtau2.jags.binHN[i] <- (m2HN[[i]]$`2.5%`)^2
  UBtau2.jags.binHN[i] <- (m2HN[[i]]$`97.5%`)^2
  Rhat_tauNHN[i] <- m2HN[[i]]$Rhat
  precN_tau2HN[i] <-  UBtau2.jags.binHN[i] - LBtau2.jags.binHN[i]
  precN_tauHN[i] <-  UBtau.jags.binHN[i] - LBtau.jags.binHN[i]
  
  fdHN[[i]] <- grepl("delta12", row.names(resultsHN[[i]]))
  mdHN[[i]] <- resultsHN[[i]][fdHN[[i]],]
  rel_effHN[[i]] <- mdHN[[i]]$`50%`
  LB_rel_effHN[[i]] <- mdHN[[i]]$`2.5%`
  UB_rel_effHN[[i]] <- mdHN[[i]]$`97.5%`
  sd_rel_effHN[[i]] <- mdHN[[i]]$sd
  Rhat_deltaNHN[[i]] <- mdHN[[i]]$Rhat
  listaHN[[i]] <- cbind.data.frame(rel_effHN[[i]],sd_rel_effHN[[i]],LB_rel_effHN[[i]], UB_rel_effHN[[i]], Rhat_deltaNHN[[i]])
}


RhatsNHN <- cbind.data.frame(mu_Normaljags.binHN, LBmu.jags.binHN, UBmu.jags.binHN, tau2_Normaljags.binHN,
                             LBtau2.jags.binHN, UBtau2.jags.binHN,
                             tau_Normaljags.binHN, LBtau.jags.binHN, UBtau.jags.binHN,
                             pred.binHN, LBpred.binHN, UBpred.binHN,
                             Rhat_muNHN, Rhat_tauNHN,Rhat_predNHN,
                             precN_muHN, precN_tau2HN,precN_tauHN, precN_predHN)

############REMOVE Rhats > 1.05 ##############

condition1HN <- which(RhatsNHN$Rhat_muNHN > 1.05)
condition2HN <- which(RhatsNHN$Rhat_tauNHN > 1.05)

########## check if condition1 == condition2 ##########
dist.condNHN = c(condition1HN,condition2HN)
dist.condNHN = unique(dist.condNHN)

RhatN_outHN = round((length(dist.condNHN)/N.sim), 4)


if (length(dist.condNHN)== 0) {
  RhatsNHN <- RhatsNHN
  N.sim <- nrow(RhatsNHN)
  listaNormHN <- listaHN
  pNHN <- p
} else {
  RhatsNHN <- RhatsNHN[-dist.condNHN, ]
  listaNormHN <- listaHN[-dist.condNHN]
  N.sim <- nrow(RhatsNHN)
  pNHN <- p[-dist.condNHN]
}


count_NHN = list()
for(i in length(listaNormHN)){
  count_NHN[[i]] = which(listaNormHN[[i]]$`Rhat_deltaNHN[[i]]` > 1.05 )
  tell_meNHN = which(count_NHN[[i]] != 0)
}

tell_meNHN

if(length(tell_meNHN) == 0){
  RhatsNHN <- RhatsNHN
  N.sim <- nrow(RhatsNHN)
  listaNormHN <- listaNormHN
  pNHN <- pNHN
  RhatN_outHN <- RhatN_outHN
} else {
  RhatsNHN <- RhatsNHN[-tell_meNHN, ]
  listaNormHN <- listaNormHN[-tell_meNHN]
  N.sim = nrow(RhatsNHN)
  pNHN <- pNHN[-tell_meNHN]
  RhatN_outHN <- RhatN_outHN + ((length(tell_meNHN))/N.sim)
}

RhatN_outHN

########## AVERAGE RELATIVE BIAS OF mu ##########
rel_biasNormaljags.binHN <- c()
for (i in 1:N.sim){
  rel_biasNormaljags.binHN[i] <- ( RhatsNHN$mu_Normaljags.binHN[i] - mu) 
}
avg_rel_biasNormaljags.binHN <- round(mean(rel_biasNormaljags.binHN), 2)

########## AVERAGE ABSOLUTE BIAS OF mu ##########
abs_biasNormaljags.binHN <- c()
for (i in 1:N.sim){
  abs_biasNormaljags.binHN[i] <- abs( RhatsNHN$mu_Normaljags.binHN[i] - mu) 
}
avg_abs_biasNormaljags.binHN <- round(mean(abs_biasNormaljags.binHN), 2)

############### MSE OF mu  ##################
MSE_Normaljags.binHN <- c()
for (i in 1:N.sim){
  MSE_Normaljags.binHN[i] <- RMSE(estimate =RhatsNHN$mu_Normaljags.binHN[i],
                                  parameter = mu,
                                  type = "RMSE",
                                  MSE = TRUE,           
                                  percent = FALSE,
                                  unname = FALSE)
}
avg_MSE_Normaljags.binHN <-  round(mean(MSE_Normaljags.binHN),4)

########## COVERAGE OF mu ##########
interval_contains_true_mean <- function(i) { 
  ( mu >= RhatsNHN$LBmu.jags.binHN[i]) && ( mu <= RhatsNHN$UBmu.jags.binHN[i])
}
interval_contains_true_mean(i)

coverageNormaljags.binHN <- c()
for (i in 1:N.sim){0
  if (interval_contains_true_mean(i) == TRUE){
    coverageNormaljags.binHN[i] = 1
  }
  else {
    coverageNormaljags.binHN[i]=0
  }
  print(coverageNormaljags.binHN[i])
}
coverageNormaljags.binHN
coverNormaljags.binHN <- round(mean(coverageNormaljags.binHN)*100, 2)

################
#### AVERAGE BIAS OF tau2 ########
bias_tau2Normaljags.binHN <- c()
for (i in 1:N.sim){
  bias_tau2Normaljags.binHN[i] <- abs(RhatsNHN$tau2_Normaljags.binHN[i] - tau2)
}
avg_bias_tau2Normaljags.binHN <- round(mean(bias_tau2Normaljags.binHN),2)

pbias_tau2Normaljags.binHN <- round(mean(bias_tau2Normaljags.binHN / tau2),2)

########### MSE OF tau2 ################
MSE_Normaljags.tau2.binHN <- c()
for (i in 1:N.sim){
  MSE_Normaljags.tau2.binHN[i] <- RMSE(estimate =RhatsNHN$tau2_Normaljags.binHN[i],
                                       parameter = tau2,
                                       type = "RMSE",
                                       MSE = TRUE           
  )
}
avg_MSE_Normaljags.tau2.binHN <- round(mean(MSE_Normaljags.tau2.binHN),4)

nMSE_Normaljags.tau2.binHN <- round(mean(MSE_Normaljags.tau2.binHN / (tau2)^2 ),4)

############# COVERAGE OF tau2 ###################
interval_contains_true_mean <- function(i) { 
  ( tau2 >= RhatsNHN$LBtau2.jags.binHN[i]) && ( tau2 <=  RhatsNHN$UBtau2.jags.binHN[i])
}
interval_contains_true_mean(i)

coverage_tau2Normaljags.binHN <- c()
for (i in 1:N.sim){
  if (interval_contains_true_mean(i) == TRUE){
    coverage_tau2Normaljags.binHN[i] = 1
  }
  else  {
    coverage_tau2Normaljags.binHN[i] = 0
  }
  print(coverage_tau2Normaljags.binHN[i])
}
coverage_tau2Normaljags.binHN
cover_tau2Normaljags.binHN <- round(mean(coverage_tau2Normaljags.binHN)*100, 2)


################
#### AVERAGE RELATIVE BIAS OF tau ########
rel_bias_tauNormaljags.binHN <- c()
for (i in 1:N.sim){
  rel_bias_tauNormaljags.binHN[i] <- (RhatsNHN$tau_Normaljags.binHN[i] - tau)
}
avg_rel_bias_tauNormaljags.binHN <- round(mean(rel_bias_tauNormaljags.binHN),2)

prel_bias_tauNormaljags.binHN <- round(mean(rel_bias_tauNormaljags.binHN / tau)*100,2)

#### AVERAGE absolute BIAS OF tau ########
abs_bias_tauNormaljags.binHN <- c()
for (i in 1:N.sim){
  abs_bias_tauNormaljags.binHN[i] <- abs(RhatsNHN$tau_Normaljags.binHN[i] - tau)
}
avg_abs_bias_tauNormaljags.binHN <- round(mean(abs_bias_tauNormaljags.binHN),2)

pabs_bias_tauNormaljags.binHN <- round(mean(abs_bias_tauNormaljags.binHN / tau)*100,2)


########### MSE OF tau ################
MSE_Normaljags.tau.binHN <- c()
for (i in 1:N.sim){
  MSE_Normaljags.tau.binHN[i] <- RMSE(estimate =RhatsNHN$tau_Normaljags.binHN[i],
                                      parameter = tau,
                                      type = "RMSE",
                                      MSE = TRUE           
  )
}
avg_MSE_Normaljags.tau.binHN <- round(mean(MSE_Normaljags.tau.binHN),4)

nMSE_Normaljags.tau.binHN <- round(mean(MSE_Normaljags.tau.binHN / tau2 ),4)

############# COVERAGE OF tau ###################
interval_contains_true_mean <- function(i) { 
  ( tau >= RhatsNHN$LBtau.jags.binHN[i]) && ( tau <=  RhatsNHN$UBtau.jags.binHN[i])
}
interval_contains_true_mean(i)

coverage_tauNormaljags.binHN <- c()
for (i in 1:N.sim){
  if (interval_contains_true_mean(i) == TRUE){
    coverage_tauNormaljags.binHN[i] = 1
  }
  else  {
    coverage_tauNormaljags.binHN[i] = 0
  }
  print(coverage_tauNormaljags.binHN[i])
}
coverage_tauNormaljags.binHN
cover_tauNormaljags.binHN <- round(mean(coverage_tauNormaljags.binHN)*100, 2)



######################## ABSOLUTE BIAS OF RELATIVE EFFECTS ##################
abs_bias_relNHN <- list()
mean_abs_bias_relNHN <- list()
for(i in 1:N.sim){
  for(j in 1:nrow(pNHN[[i]])){
    abs_bias_relNHN[[i]] <- round(abs( listaNormHN[[i]]$`rel_effHN[[i]]` - pNHN[[i]]$true_eff  ),2)
    mean_abs_bias_relNHN[[i]] <- mean(abs_bias_relNHN[[i]])
  }
}

################## mean bias of the mean biases of each data set for relative effects  ###################
mean_abs_bias_of_releffHN <- unlist(mean_abs_bias_relNHN)
mean_abs_bias_of_mean_abs_bias_relNHN  <- round(mean(mean_abs_bias_of_releffHN),2)

######################## RELATIVE  BIAS OF RELATIVE EFFECTS ##################
rel_bias_relNHN <- list()
mean_rel_bias_relNHN <- list()
for(i in 1:N.sim){
  for(j in 1:nrow(pNHN[[i]])){
    rel_bias_relNHN[[i]] <- round(( listaNormHN[[i]]$`rel_effHN[[i]]` - pNHN[[i]]$true_eff  ),2)
    mean_rel_bias_relNHN[[i]] <- mean(rel_bias_relNHN[[i]])
  }
}

################## mean bias of the mean biases of each data set for relative effects  ###################
mean_rel_bias_of_releffHN <- unlist(mean_rel_bias_relNHN)
mean_rel_bias_of_mean_rel_bias_relNHN  <- round(mean(mean_rel_bias_of_releffHN),2)

########################### FOR THE OVERLAPPING ######################
dataNHN <- list()
overlapN1HN <- list()
overlapNHN <- list()
for(i in 1:N.sim){
  dataNHN[[i]] <- list(NORMAL =listaNormHN[[i]]$`rel_effHN[[i]]` , True =pNHN[[i]]$true_eff )
  overlapN1HN[[i]] <- overlap(dataNHN[[i]], type = "2")
  overlapNHN[[i]] <- overlapN1HN[[i]]$OV
}
overlapNormHN <- unlist(overlapNHN)
avg_overlapNormHN <- round(mean(overlapNormHN),2)

######################### SAVE THE OUTPUTS ####################
lista_rel_effHN = cbind( mean_rel_bias_relNHN, mean_rel_bias_of_mean_rel_bias_relNHN,
                         mean_abs_bias_relNHN, mean_abs_bias_of_mean_abs_bias_relNHN) 

len_mu_NHN = mean(RhatsNHN$precN_muHN)
len_PI_NHN = mean(RhatsNHN$precN_predHN)

Norm_bin_jagsHN <- cbind(RhatsNHN , MSE_Normaljags.binHN, MSE_Normaljags.tau2.binHN,coverageNormaljags.binHN,
                         coverage_tau2Normaljags.binHN, rel_biasNormaljags.binHN, abs_biasNormaljags.binHN,
                         bias_tau2Normaljags.binHN,
                         avg_rel_biasNormaljags.binHN,avg_abs_biasNormaljags.binHN,
                         avg_bias_tau2Normaljags.binHN,pbias_tau2Normaljags.binHN,
                         avg_MSE_Normaljags.binHN, avg_MSE_Normaljags.tau2.binHN,nMSE_Normaljags.tau2.binHN,
                         coverNormaljags.binHN,
                         cover_tau2Normaljags.binHN, 
                         
                         abs_bias_tauNormaljags.binHN, pabs_bias_tauNormaljags.binHN,
                         avg_abs_bias_tauNormaljags.binHN,
                         
                         rel_bias_tauNormaljags.binHN, avg_rel_bias_tauNormaljags.binHN,prel_bias_tauNormaljags.binHN,
                         avg_MSE_Normaljags.tau.binHN, nMSE_Normaljags.tau.binHN,
                         coverNormaljags.binHN,cover_tauNormaljags.binHN,
                         len_mu_NHN ,
                         len_PI_NHN ,
                         RhatN_outHN)

Overalp_NormalHN <- cbind(overlapNormHN, avg_overlapNormHN )



write.csv(Norm_bin_jagsHN, "C:\\Users\\Tianqi YU\\Desktop\\SIMULATION SCENARIOS FINAL\\NORMAL DATA\\NORMAL MODEL\\rev_HN_ResultsNormal_jags_Normal_bin_data_k=14_mu=0_tau2=0.12.csv",row.names=FALSE )  
write.csv(listaNormHN, "C:\\Users\\Tianqi YU\\Desktop\\SIMULATION SCENARIOS FINAL\\NORMAL DATA\\NORMAL MODEL\\rev_HN_Relative_eff_ResultsNormal_jags_Normal_bin_data_k=14_mu=0_tau2=0.12.csv",row.names=FALSE )  
write.csv(lista_rel_effHN, "C:\\Users\\Tianqi YU\\Desktop\\SIMULATION SCENARIOS FINAL\\NORMAL DATA\\NORMAL MODEL\\rev_HN_Bias_Relative_eff_ResultsNormal_jags_Normal_bin_data_k=14_mu=0_tau2=0.12.csv",row.names=FALSE )  
write.csv(Overalp_NormalHN, "C:\\Users\\Tianqi YU\\Desktop\\SIMULATION SCENARIOS FINAL\\NORMAL DATA\\NORMAL MODEL\\rev_HN_Overalp_Normal_Normal_bin_data_k=14_mu=0_tau2=0.12.csv",row.names=FALSE )  



####################### DIRICHLET PROCESS MODEL ######################
############# BINOMIAL-DP-26(HN/UNIF) #############
cat("model{

  for( i in 1: ns ) {

    m[i] ~ dnorm(0,.0001)
    
    logit(p1[i]) <- m[i]
    r1[i] ~ dbin(p1[i],n1[i])
    
    logit(p2[i]) <- m[i] + delta12[i]
    r2[i] ~ dbin(p2[i],n2[i])
    
    delta12[i] <- theta[Z[i]]
    
    Z[i] ~ dcat(p[]) #Z is an integer variable
  }
  # Constructive DP
  #stick-breaking prior
  p[1] <- r[1]
  for (j in 2:(N-1)){
   p[j] <- r[j]*(1-r[j-1])*p[j-1]/r[j-1]
  }
  for (k in 1:(N-1)){
   r[k] ~ dbeta(1,alpha)T(0,0.99)
   }
  #assumption to ensure sum p[] is 1 Ishwaran truncation
  ps <- sum(p[1:(N-1)])
  for(k in N:N){
   p[k]<-1-ps
   }
  # Baseline distribution
  for(k in 1:N){
   theta[k] ~ dnorm(basemu ,basetau1)
  }
  #### priors 

  basemu~dnorm(0,0.0001)
  basetau1 <- 1/basetau_sqr
  basetau_sqr <- basetau*basetau
  basetau~dnorm(0,1)T(0,)

  # DPP parameter prior
  alpha~dunif(0.3,5)
  
  # Random effects distribution mean#
 
  for(i in 1:N){
   meancl[i]<-p[i]*theta[i]
   }
  poptrue<-sum(meancl[])  ### the E[X]
  
  # Random effects distribution variance #
  for(i in 1:N){
   mom2[i]<-p[i]*theta[i]*theta[i]  ####E[X2]
   }
  mom2.true<-sum(mom2[])
  var.true<-mom2.true-(poptrue*poptrue) ###E[X2] - E[X]2
  
  # Programming for calculating summary statistics #
  for(i in 1:ns){
   for (j in 1:N){
    SC[i,j] <- equals(j,Z[i])
   }
  }
  # total clusters K#
  for (j in 1:N){
   cl[j] <-step(sum(SC[,j])-1)
   }
  K<-sum(cl[])
  
}",file="DPmodel2.bin_HN_U_26.txt")
modfile = 'DPmodel2.bin_HN_U_26.txt'
N.sim =1000
mu = 0
tau2 = 0.12
tau = sqrt(tau2)

cl = makeCluster(19)
registerDoParallel(cl)
DPresults1_HN_U_26 <- list()

DPresults1_HN_U_26 = foreach(i = 1:N.sim,.packages = c("rjags", "R2jags","tidyverse")) %dopar% {
  set.seed(2103)
  run.model = jags(
    data =list(ns = nrow(p[[i]]),
               r1 = p[[i]]$ci,
               r2 = p[[i]]$ti,
               n1 = p[[i]]$nci,
               n2 = p[[i]]$nti,
               N= 26
               
    ) ,  
    inits = NULL,
    parameters.to.save = c(
      "basemu",    ## tau of the base Normal distribution
      "basetau",   ## mu of the base Normal distribution
      "poptrue",   ##the overall mu
      "var.true", ## the heterogeneity 
      "delta12", #### the random effects
      "K",       ### the total number of clusters 
      "p",     ## the weights of the process 
      "alpha"  ## the concentration parameter
    ), 
    
    n.chains = 2,
    n.iter = 50000,
    n.burnin = 10000,
    DIC = T,
    model.file = modfile
    
  )
  
}
stopCluster(cl)


DPresults_HN_U_26 <- list()
posteriorDP_HN_U_26 <- list()
postDP_HN_U_26 <- list()
for(i in 1:N.sim){
  DPresults_HN_U_26[[i]] <- as.data.frame(DPresults1_HN_U_26[[i]]$BUGSoutput$summary) 
  posteriorDP_HN_U_26[[i]] = as.data.frame(DPresults1_HN_U_26[[i]]$BUGSoutput$sims.matrix)
  postDP_HN_U_26[[i]] <- posteriorDP_HN_U_26[[i]]$poptrue
}

DP_results_HN_U_26 <- list()
for(i in 1:N.sim){
  DP_results_HN_U_26[[i]] <- DPresults_HN_U_26[[i]] 
}

############# extraction of the parameters of interest ###########
alpha_HN_U_26 <- c()
LB_alpha_HN_U_26 <- c()
UB_alpha_HN_U_26 <- c()
mu_DP_HN_U_26 <- c()
LB_mu_DP_HN_U_26 <- c()
UB_mu_DP_HN_U_26 <- c()
tau2_DP_HN_U_26 <- c()
LB_tau2_DP_HN_U_26 <- c()
UB_tau2_DP_HN_U_26 <- c()
precDP_mu_HN_U_26 <- c()
precDP_tau2_HN_U_26 <- c()
prec_basemu_HN_U_26 <- c()
prec_basetau2_HN_U_26 <- c()
prec_alpha_HN_U_26 <- c()
base_mu_HN_U_26 <- c()
base_tau_HN_U_26 <- c()
base_tau2_HN_U_26 <- c()
LB_base_mu_HN_U_26 <- c()
UB_base_mu_HN_U_26 <- c()
LB_base_tau2_HN_U_26 <- c()
UB_base_tau2_HN_U_26 <- c()
Rhat_muDP_HN_U_26 <- c()
Rhat_tau2DP_HN_U_26 <- c()
Rhat_basemu_HN_U_26 <- c()
Rhat_basetau_HN_U_26 <- c()
Rhat_alpha_HN_U_26 <- c()
Rhat_deltaDP_HN_U_26 <- list()
median_K_HN_U_26 <- list()
LB_K_HN_U_26 <- list()
UB_K_HN_U_26 <- list()
rel_effDP_HN_U_26 <- list()
LB_rel_effDP_HN_U_26 <- list()
UB_rel_effDP_HN_U_26 <- list()
sd_rel_effDP_HN_U_26 <- list()
pi_HN_U_26 <- list()
f11_HN_U_26 <- list()
m11_HN_U_26 <- list()
f22_HN_U_26 <- list()
m22_HN_U_26 <- list()
fdd_HN_U_26 <- list()
mdd_HN_U_26 <- list()
f33_HN_U_26 <- list()
m33_HN_U_26 <- list()
f44_HN_U_26 <- list()
m44_HN_U_26 <- list()
f55_HN_U_26 <- list()
m55_HN_U_26 <- list()
fcl_HN_U_26 <- list()
mcl_HN_U_26 <- list()
fp_HN_U_26 <- list()
mp_HN_U_26 <- list()
listaDP1_HN_U_26 <- list()
listaDP_weights_HN_U_26 <- list()
for(i in 1:N.sim){
  f11_HN_U_26[[i]] <- grepl("poptrue", row.names(DP_results_HN_U_26[[i]]))
  m11_HN_U_26[[i]] <- DP_results_HN_U_26[[i]][f11_HN_U_26[[i]],]
  mu_DP_HN_U_26[i] <- m11_HN_U_26[[i]]$`50%`
  LB_mu_DP_HN_U_26[i] <- m11_HN_U_26[[i]]$`2.5%`
  UB_mu_DP_HN_U_26[i] <- m11_HN_U_26[[i]]$`97.5%`
  Rhat_muDP_HN_U_26[i] <- m11_HN_U_26[[i]]$Rhat
  precDP_mu_HN_U_26[i] <- UB_mu_DP_HN_U_26[i] - LB_mu_DP_HN_U_26[i]
  
  f22_HN_U_26[[i]] <- grepl("var.true", row.names(DP_results_HN_U_26[[i]]))
  m22_HN_U_26[[i]] <- DP_results_HN_U_26[[i]][f22_HN_U_26[[i]],]
  tau2_DP_HN_U_26[i] <- m22_HN_U_26[[i]]$`50%`
  LB_tau2_DP_HN_U_26[i] <- m22_HN_U_26[[i]]$`2.5%`
  UB_tau2_DP_HN_U_26[i] <- m22_HN_U_26[[i]]$`97.5%`
  Rhat_tau2DP_HN_U_26[i] <- m22_HN_U_26[[i]]$Rhat
  precDP_tau2_HN_U_26[i] <- UB_tau2_DP_HN_U_26[i] -  LB_tau2_DP_HN_U_26[i]
  
  fdd_HN_U_26[[i]] <- grepl("delta12", row.names(DP_results_HN_U_26[[i]]))
  mdd_HN_U_26[[i]] <- DP_results_HN_U_26[[i]][fdd_HN_U_26[[i]],]
  rel_effDP_HN_U_26[[i]] <- mdd_HN_U_26[[i]]$`50%`
  LB_rel_effDP_HN_U_26[[i]] <- mdd_HN_U_26[[i]]$`2.5%`
  UB_rel_effDP_HN_U_26[[i]] <- mdd_HN_U_26[[i]]$`97.5%`
  sd_rel_effDP_HN_U_26[[i]] <- mdd_HN_U_26[[i]]$sd
  Rhat_deltaDP_HN_U_26[[i]] <- mdd_HN_U_26[[i]]$Rhat
  
  f33_HN_U_26[[i]] <- grepl("basemu", row.names(DP_results_HN_U_26[[i]]))
  m33_HN_U_26[[i]] <- DP_results_HN_U_26[[i]][f33_HN_U_26[[i]],]
  base_mu_HN_U_26[i] <- m33_HN_U_26[[i]]$`50%`
  LB_base_mu_HN_U_26[i] <- m33_HN_U_26[[i]]$`2.5%`
  UB_base_mu_HN_U_26[i] <- m33_HN_U_26[[i]]$`97.5%`
  Rhat_basemu_HN_U_26[i] <- m33_HN_U_26[[i]]$Rhat
  prec_basemu_HN_U_26[i] <- UB_base_mu_HN_U_26[i] - LB_base_mu_HN_U_26[i]
  
  f44_HN_U_26[[i]] <- grepl("basetau", row.names(DP_results_HN_U_26[[i]]))
  m44_HN_U_26[[i]] <- DP_results_HN_U_26[[i]][f44_HN_U_26[[i]],]
  base_tau_HN_U_26[i] <- m44_HN_U_26[[i]]$`50%`
  base_tau2_HN_U_26[i] <- (m44_HN_U_26[[i]]$`50%`)^2
  LB_base_tau2_HN_U_26[i] <- (m44_HN_U_26[[i]]$`2.5%`)^2
  UB_base_tau2_HN_U_26[i] <- (m44_HN_U_26[[i]]$`97.5%`)^2
  Rhat_basetau_HN_U_26[i] <- (m44_HN_U_26[[i]]$Rhat)^2
  prec_basetau2_HN_U_26[i] <- UB_base_tau2_HN_U_26[i] - LB_base_tau2_HN_U_26[i]
  
  f55_HN_U_26[[i]] <- grepl("alpha", row.names(DP_results_HN_U_26[[i]]))
  m55_HN_U_26[[i]] <- DP_results_HN_U_26[[i]][f55_HN_U_26[[i]],]
  alpha_HN_U_26[i] <- m55_HN_U_26[[i]]$`50%`
  LB_alpha_HN_U_26[i] <- m55_HN_U_26[[i]]$`2.5%`
  UB_alpha_HN_U_26[i] <- m55_HN_U_26[[i]]$`97.5%`
  Rhat_alpha_HN_U_26[i] <- m55_HN_U_26[[i]]$Rhat
  prec_alpha_HN_U_26[i] <- UB_alpha_HN_U_26[i] - LB_alpha_HN_U_26[i]
  
  fcl_HN_U_26[[i]] <- grepl("K", row.names(DP_results_HN_U_26[[i]]))  
  mcl_HN_U_26[[i]] <- DP_results_HN_U_26[[i]][fcl_HN_U_26[[i]],]
  median_K_HN_U_26[[i]] <- mcl_HN_U_26[[i]]$`50%` 
  LB_K_HN_U_26[[i]] <- mcl_HN_U_26[[i]]$`2.5%`
  UB_K_HN_U_26[[i]] <- mcl_HN_U_26[[i]]$`97.5%`
  
  fp_HN_U_26[[i]] <- grepl("p", row.names(DP_results_HN_U_26[[i]]))
  mp_HN_U_26[[i]] <- DP_results_HN_U_26[[i]][fp_HN_U_26[[i]],]
  mp_HN_U_26[[i]] <- mp_HN_U_26[[i]][!grepl("pop", row.names(mp_HN_U_26[[i]])),]
  mp_HN_U_26[[i]] <- mp_HN_U_26[[i]][!grepl("alpha", row.names(mp_HN_U_26[[i]])),]
  pi_HN_U_26[[i]] <- mp_HN_U_26[[i]]$mean
  
  listaDP1_HN_U_26[[i]] <- cbind.data.frame(rel_effDP_HN_U_26[[i]],sd_rel_effDP_HN_U_26[[i]], LB_rel_effDP_HN_U_26[[i]],UB_rel_effDP_HN_U_26[[i]],Rhat_deltaDP_HN_U_26[[i]])
  listaDP_weights_HN_U_26[[i]] <- cbind.data.frame(pi_HN_U_26[[i]])
  
}

numclus_HN_U_26 <- unlist(median_K_HN_U_26)
LB_K_HN_U_26 <- unlist(LB_K_HN_U_26)
UB_K_HN_U_26 <- unlist(UB_K_HN_U_26)

mean_alpha_HN_U_26 = mean(alpha_HN_U_26)
mean_alpha_LB_HN_U_26 = mean(LB_alpha_HN_U_26)
mean_alpha_UB_HN_U_26 = mean(UB_alpha_HN_U_26)

tau_DP_HN_U_26 = sqrt(tau2_DP_HN_U_26)
LB_tau_DP_HN_U_26 = sqrt(LB_tau2_DP_HN_U_26)
UB_tau_DP_HN_U_26 = sqrt(UB_tau2_DP_HN_U_26)
precDP_tau_HN_U_26 = UB_tau_DP_HN_U_26 -  LB_tau_DP_HN_U_26


RhatsDP_HN_U_26 <- cbind.data.frame(base_mu_HN_U_26,LB_base_mu_HN_U_26,UB_base_mu_HN_U_26,base_tau_HN_U_26,base_tau2_HN_U_26, LB_base_tau2_HN_U_26,UB_base_tau2_HN_U_26,
                                    mu_DP_HN_U_26, LB_mu_DP_HN_U_26, UB_mu_DP_HN_U_26, tau2_DP_HN_U_26, LB_tau2_DP_HN_U_26, UB_tau2_DP_HN_U_26,
                                    tau_DP_HN_U_26, LB_tau_DP_HN_U_26, UB_tau_DP_HN_U_26,
                                    precDP_mu_HN_U_26, precDP_tau2_HN_U_26, prec_basemu_HN_U_26, prec_basetau2_HN_U_26,prec_alpha_HN_U_26,  alpha_HN_U_26 , LB_alpha_HN_U_26, UB_alpha_HN_U_26,
                                    mean_alpha_HN_U_26, mean_alpha_LB_HN_U_26, mean_alpha_UB_HN_U_26,
                                    Rhat_muDP_HN_U_26, Rhat_tau2DP_HN_U_26, Rhat_basemu_HN_U_26, Rhat_basetau_HN_U_26, Rhat_alpha_HN_U_26, numclus_HN_U_26, LB_K_HN_U_26, UB_K_HN_U_26)

##########REMOVE Rhats > 1.05 ##############

condition1_HN_U_26 <- which(RhatsDP_HN_U_26$Rhat_muDP_HN_U_26 > 1.05)
condition2_HN_U_26 <- which(RhatsDP_HN_U_26$Rhat_tau2DP_HN_U_26 > 1.05)
condition3_HN_U_26 <- which(RhatsDP_HN_U_26$Rhat_basemu_HN_U_26 > 1.05)
condition4_HN_U_26 <- which(RhatsDP_HN_U_26$Rhat_basetau_HN_U_26 > 1.05)
condition5_HN_U_26 <- which(RhatsDP_HN_U_26$Rhat_alpha_HN_U_26 > 1.05)


dist.condDP_HN_U_26 = c(condition1_HN_U_26,condition2_HN_U_26,condition3_HN_U_26,condition4_HN_U_26,condition5_HN_U_26)
dist.condDP_HN_U_26 = unique(dist.condDP_HN_U_26)

RhatDP_out_HN_U_26 = round((length(dist.condDP_HN_U_26)/N.sim), 4)

############### Extract and remove the datasets with Rhat > 1.05 #########
if (length(dist.condDP_HN_U_26)== 0) {
  RhatsDP_HN_U_26 <- RhatsDP_HN_U_26
  listaDP_HN_U_26 <- listaDP1_HN_U_26
  listaDP_weights_HN_U_26 <- listaDP_weights_HN_U_26
  N.sim <- nrow(RhatsDP_HN_U_26)
  pDP_HN_U_26 <- p
  
} else {
  RhatsDP_HN_U_26 <- RhatsDP_HN_U_26[-dist.condDP_HN_U_26, ]
  listaDP_HN_U_26 <- listaDP1_HN_U_26[-dist.condDP_HN_U_26]
  listaDP_weights_HN_U_26 <- listaDP_weights_HN_U_26[-dist.condDP_HN_U_26]
  N.sim <- nrow(RhatsDP_HN_U_26)
  pDP_HN_U_26 <- p[-dist.condDP_HN_U_26]
}

count_DP_HN_U_26 = list()
for(i in length(listaDP_HN_U_26)){
  count_DP_HN_U_26[[i]] = which(listaDP_HN_U_26[[i]]$`Rhat_deltaDP_HN_U_26[[i]]` > 1.05 )
  tell_meDP_HN_U_26 = which(count_DP_HN_U_26[[i]] != 0)
}

tell_meDP_HN_U_26

if(length(tell_meDP_HN_U_26) == 0){
  RhatsDP_HN_U_26 <- RhatsDP_HN_U_26
  N.sim <- nrow(RhatsDP_HN_U_26)
  listaDP_HN_U_26 <- listaDP_HN_U_26
  listaDP_weights_HN_U_26 <- listaDP_weights_HN_U_26
  pDP_HN_U_26 <- pDP_HN_U_26
  RhatDP_out_HN_U_26 <- RhatDP_out_HN_U_26
} else {
  RhatsDP_HN_U_26 <- RhatsDP_HN_U_26[-tell_meDP_HN_U_26, ]
  listaDP_HN_U_26 <- listaDP_HN_U_26[-tell_meDP_HN_U_26]
  listaDP_weights_HN_U_26 <- listaDP_weights_HN_U_26[-tell_meDP_HN_U_26]
  N.sim = nrow(RhatsDP_HN_U_26)
  pDP_HN_U_26 <- pDP_HN_U_26[-tell_meDP_HN_U_26]
  RhatDP_out_HN_U_26 <- RhatDP_out_HN_U_26 + ((length(tell_meDP_HN_U_26))/N.sim)
}

RhatDP_out_HN_U_26

#### AVERAGE ABSOLUTE BIAS OF mu ######
abs_biasDPjags.bin_HN_U_26 <- c()
for (i in 1:N.sim){
  abs_biasDPjags.bin_HN_U_26[i] <- abs( RhatsDP_HN_U_26$mu_DP_HN_U_26[i] - mu) 
}
avg_abs_biasDPjags.bin_HN_U_26 <- round(mean(abs_biasDPjags.bin_HN_U_26),2)

#### AVERAGE RELATIVE BIAS OF mu ######
rel_biasDPjags.bin_HN_U_26 <- c()
for (i in 1:N.sim){
  rel_biasDPjags.bin_HN_U_26[i] <- ( RhatsDP_HN_U_26$mu_DP_HN_U_26[i] - mu) 
}
avg_rel_biasDPjags.bin_HN_U_26 <- round(mean(rel_biasDPjags.bin_HN_U_26),2)

######  MSE OF mu ##########
MSE_DPjags.bin_HN_U_26 <- c()
for (i in 1:N.sim){
  MSE_DPjags.bin_HN_U_26[i] <- RMSE(estimate = RhatsDP_HN_U_26$mu_DP_HN_U_26[i],
                                    parameter = mu,
                                    type = "RMSE",
                                    MSE = TRUE )        
}
avg_MSE_DPjags.bin_HN_U_26 <- round(mean(MSE_DPjags.bin_HN_U_26), 4)

############ COVERAGE OF mu ########
interval_contains_true_mean <- function(i) { 
  ( mu >= RhatsDP_HN_U_26$LB_mu_DP_HN_U_26[i]) && ( mu <= RhatsDP_HN_U_26$UB_mu_DP_HN_U_26[i])
}
interval_contains_true_mean(i)

coverageDPjags.bin_HN_U_26 <- c()
for (i in 1:N.sim){
  if (interval_contains_true_mean(i) == TRUE){
    coverageDPjags.bin_HN_U_26[i] = 1
  }
  else {
    coverageDPjags.bin_HN_U_26[i]=0
  }
  print(coverageDPjags.bin_HN_U_26[i])
}
coverageDPjags.bin_HN_U_26
coverDPjags.bin_HN_U_26 <- round(mean(coverageDPjags.bin_HN_U_26)*100, 2)

## BIAS OF base_mu ######
biasDPjags_basemu.bin_HN_U_26 <- c()
for (i in 1:N.sim){
  biasDPjags_basemu.bin_HN_U_26[i] <- abs( RhatsDP_HN_U_26$base_mu_HN_U_26[i] - mu) 
}
avg_biasDPjags_basemu.bin_HN_U_26 <- round(mean(biasDPjags_basemu.bin_HN_U_26),2)

###### MSE OF base_mu ##########
MSE_DPjags_basemu.bin_HN_U_26 <- c()
for (i in 1:N.sim){
  MSE_DPjags_basemu.bin_HN_U_26[i] <- RMSE(estimate =RhatsDP_HN_U_26$base_mu_HN_U_26[i],
                                           parameter = mu,
                                           type = "RMSE",
                                           MSE = TRUE      
  )
}
avg_MSE_DPjags_basemu.bin_HN_U_26 <- round(mean(MSE_DPjags_basemu.bin_HN_U_26), 4)

############ COVERAGE OF base_mu ########
interval_contains_true_mean <- function(i) { 
  ( mu >= RhatsDP_HN_U_26$LB_base_mu_HN_U_26[i]) && ( mu <= RhatsDP_HN_U_26$UB_base_mu_HN_U_26[i])
}
interval_contains_true_mean(i)

coverageDPjags_basemu.bin_HN_U_26 <- c()
for (i in 1:N.sim){
  if (interval_contains_true_mean(i) == TRUE){
    coverageDPjags_basemu.bin_HN_U_26[i] = 1
  }
  else {
    coverageDPjags_basemu.bin_HN_U_26[i]=0
  }
  print(coverageDPjags_basemu.bin_HN_U_26[i])
}
coverageDPjags_basemu.bin_HN_U_26
coverDPjags_basemu.bin_HN_U_26 <- round(mean(coverageDPjags_basemu.bin_HN_U_26)*100, 2)

######## AVERAGE BIAS OF tau2 ###########
bias_tau2DPjags.bin_HN_U_26 <- c()
for (i in 1:N.sim){
  bias_tau2DPjags.bin_HN_U_26[i] <- abs(RhatsDP_HN_U_26$tau2_DP_HN_U_26[i]- tau2)
}
avg_bias_tau2DPjags.bin_HN_U_26 <- round(mean(bias_tau2DPjags.bin_HN_U_26),2)

pbias_tau2DPjags.bin_HN_U_26 <- round(mean(bias_tau2DPjags.bin_HN_U_26 / tau2),2)

######## MSE OF tau2 #############
MSE_DPjags.bin.tau2_HN_U_26 <- c()
for (i in 1:N.sim){
  MSE_DPjags.bin.tau2_HN_U_26[i] <- RMSE(estimate =RhatsDP_HN_U_26$tau2_DP_HN_U_26[i],
                                         parameter = tau2,
                                         type = "RMSE",
                                         MSE = TRUE         
  )
  
}
avg_MSE_DPjags.bin.tau2_HN_U_26 <- round(mean(MSE_DPjags.bin.tau2_HN_U_26),4)

nMSE_DPjags.bin.tau2_HN_U_26 <- round(mean(MSE_DPjags.bin.tau2_HN_U_26/ (tau2)^2),4)

############ COVERAGE OF tau2 ##################
interval_contains_true_mean <- function(i) { 
  ( tau2 >= RhatsDP_HN_U_26$LB_tau2_DP_HN_U_26[i]) && ( tau2 <= RhatsDP_HN_U_26$UB_tau2_DP_HN_U_26[i])
}
interval_contains_true_mean(i)

coverage_tau2DPjags.bin_HN_U_26 <- c()
for (i in 1:N.sim){
  if (interval_contains_true_mean(i) == TRUE){
    coverage_tau2DPjags.bin_HN_U_26[i] = 1
  }
  else  {
    coverage_tau2DPjags.bin_HN_U_26[i] = 0
  }
  print(coverage_tau2DPjags.bin_HN_U_26[i])
}
coverage_tau2DPjags.bin_HN_U_26
cover_tau2DPjags.bin_HN_U_26 <- round(mean(coverage_tau2DPjags.bin_HN_U_26)*100, 2)


######## AVERAGE RELATIVE  BIAS OF tau ###########
rel_bias_tauDPjags.bin_HN_U_26 <- c()
for (i in 1:N.sim){
  rel_bias_tauDPjags.bin_HN_U_26[i] <- (RhatsDP_HN_U_26$tau_DP_HN_U_26[i]- tau)
}
avg_rel_bias_tauDPjags.bin_HN_U_26 <- round(mean(rel_bias_tauDPjags.bin_HN_U_26),2)

prel_bias_tauDPjags.bin_HN_U_26 <- round(mean(rel_bias_tauDPjags.bin_HN_U_26 / tau)*100,2)

######## AVERAGE ABSOLUTE  BIAS OF tau ###########
abs_bias_tauDPjags.bin_HN_U_26 <- c()
for (i in 1:N.sim){
  abs_bias_tauDPjags.bin_HN_U_26[i] <- abs(RhatsDP_HN_U_26$tau_DP_HN_U_26[i]- tau)
}
avg_abs_bias_tauDPjags.bin_HN_U_26 <- round(mean(abs_bias_tauDPjags.bin_HN_U_26),2)

pabs_bias_tauDPjags.bin_HN_U_26 <- round(mean(abs_bias_tauDPjags.bin_HN_U_26 / tau)*100,2)

######## MSE OF tau #############
MSE_DPjags.bin.tau_HN_U_26 <- c()
for (i in 1:N.sim){
  MSE_DPjags.bin.tau_HN_U_26[i] <- RMSE(estimate =RhatsDP_HN_U_26$tau_DP_HN_U_26[i],
                                        parameter = tau,
                                        type = "RMSE",
                                        MSE = TRUE         
  )
  
}
avg_MSE_DPjags.bin.tau_HN_U_26 <- round(mean(MSE_DPjags.bin.tau_HN_U_26),4)

nMSE_DPjags.bin.tau_HN_U_26 <- round(mean(MSE_DPjags.bin.tau_HN_U_26/ tau2),4)

############ COVERAGE OF tau ##################
interval_contains_true_mean <- function(i) { 
  ( tau >= RhatsDP_HN_U_26$LB_tau_DP_HN_U_26[i]) && ( tau <= RhatsDP_HN_U_26$UB_tau_DP_HN_U_26[i])
}
interval_contains_true_mean(i)

coverage_tauDPjags.bin_HN_U_26 <- c()
for (i in 1:N.sim){
  if (interval_contains_true_mean(i) == TRUE){
    coverage_tauDPjags.bin_HN_U_26[i] = 1
  }
  else  {
    coverage_tauDPjags.bin_HN_U_26[i] = 0
  }
  print(coverage_tauDPjags.bin_HN_U_26[i])
}
coverage_tauDPjags.bin_HN_U_26
cover_tauDPjags.bin_HN_U_26 <- round(mean(coverage_tauDPjags.bin_HN_U_26)*100, 2)


######## AVERAGE BIAS OF base_tau2 ###########
bias_basetau2DPjags.bin_HN_U_26 <- c()
for (i in 1:N.sim){
  bias_basetau2DPjags.bin_HN_U_26[i] <- abs(RhatsDP_HN_U_26$base_tau2_HN_U_26[i]- tau2)
}
avg_bias_basetau2DPjags.bin_HN_U_26 <- round(mean(bias_basetau2DPjags.bin_HN_U_26),2)

######## MSE OF base_tau2 #############
MSE_DPjags.bin.basetau2_HN_U_26 <- c()
for (i in 1:N.sim){
  MSE_DPjags.bin.basetau2_HN_U_26[i] <- RMSE(estimate =RhatsDP_HN_U_26$base_tau2_HN_U_26[i],
                                             parameter = tau2,
                                             type = "RMSE",
                                             MSE = TRUE         
  )
  
}
avg_MSE_DPjags.bin.basetau2_HN_U_26 <- round(mean(MSE_DPjags.bin.basetau2_HN_U_26),4)

############ COVERAGE OF base_tau2 ##################
interval_contains_true_mean <- function(i) {
  ( tau2 >= RhatsDP_HN_U_26$LB_base_tau2_HN_U_26[i]) && ( tau2 <= RhatsDP_HN_U_26$UB_base_tau2_HN_U_26[i])
}
interval_contains_true_mean(i)

coverage_basetau2DPjags.bin_HN_U_26 <- c()
for (i in 1:N.sim){
  if (interval_contains_true_mean(i) == TRUE){
    coverage_basetau2DPjags.bin_HN_U_26[i] = 1
  }
  else  {
    coverage_basetau2DPjags.bin_HN_U_26[i] = 0
  }
  print(coverage_basetau2DPjags.bin_HN_U_26[i])
}
coverage_basetau2DPjags.bin_HN_U_26
cover_basetau2DPjags.bin_HN_U_26 <- round(mean(coverage_basetau2DPjags.bin_HN_U_26)*100, 2)

######################## ABSOLUTE BIAS OF RELATIVE EFFECTS ##################
abs_bias_relDP_HN_U_26 <- list()
mean_abs_bias_relDP_HN_U_26 <- list()
for(i in 1:N.sim){
  for(j in 1:nrow(pDP_HN_U_26[[i]])){
    abs_bias_relDP_HN_U_26[[i]] <- round(abs(listaDP_HN_U_26[[i]]$`rel_effDP_HN_U_26[[i]]` - pDP_HN_U_26[[i]]$true_eff ),2)
    mean_abs_bias_relDP_HN_U_26[[i]] <- mean(abs_bias_relDP_HN_U_26[[i]])
  }
}
#mean_bias_relDP

################## mean bias of the mean biases of each data set for relative effects  ###################
mean_abs_bias_of_mean_abs_bias_relDP_HN_U_26 <- unlist(mean_abs_bias_relDP_HN_U_26)
mean_abs_bias_of_mean_abs_bias_relDP_HN_U_26  <- round(mean(mean_abs_bias_of_mean_abs_bias_relDP_HN_U_26),2)

######################## RELATIVE BIAS OF RELATIVE EFFECTS ##################
rel_bias_relDP_HN_U_26 <- list()
mean_rel_bias_relDP_HN_U_26 <- list()
for(i in 1:N.sim){
  for(j in 1:nrow(pDP_HN_U_26[[i]])){
    rel_bias_relDP_HN_U_26[[i]] <- round((listaDP_HN_U_26[[i]]$`rel_effDP_HN_U_26[[i]]` - pDP_HN_U_26[[i]]$true_eff ),2)
    mean_rel_bias_relDP_HN_U_26[[i]] <- mean(rel_bias_relDP_HN_U_26[[i]])
  }
}
#mean_bias_relDP

################## mean bias of the mean biases of each data set for relative effects  ###################
mean_rel_bias_of_mean_rel_bias_relDP_HN_U_26 <- unlist(mean_rel_bias_relDP_HN_U_26)
mean_rel_bias_of_mean_rel_bias_relDP_HN_U_26  <- round(mean(mean_rel_bias_of_mean_rel_bias_relDP_HN_U_26),2)

#####
########################### FOR THE OVERLAPPING ######################
dataDP_HN_U_26 <- list()
overlapDP1_HN_U_26 <- list()
overlapDP_HN_U_26 <- list()
for(i in 1:N.sim){
  dataDP_HN_U_26[[i]] <- list(DP = listaDP_HN_U_26[[i]]$`rel_effDP_HN_U_26[[i]]` , True = pDP_HN_U_26[[i]]$true_eff  )
  overlapDP1_HN_U_26[[i]] <- overlap(dataDP_HN_U_26[[i]], type = "2")
  overlapDP_HN_U_26[[i]] <- overlapDP1_HN_U_26[[i]]$OV
}
overlapDP_HN_U_26 <- unlist(overlapDP_HN_U_26)
avg_overlapDP_HN_U_26 <- round(mean(overlapDP_HN_U_26),2)


listaDP_rel_eff_HN_U_26 <- cbind( mean_abs_bias_relDP_HN_U_26, mean_abs_bias_of_mean_abs_bias_relDP_HN_U_26,
                                  mean_rel_bias_relDP_HN_U_26, mean_rel_bias_of_mean_rel_bias_relDP_HN_U_26) 

######### SAVE THE OUTPUTS ###########
DP_bin_jags_HN_U_26 <- cbind(RhatsDP_HN_U_26,
                             abs_biasDPjags.bin_HN_U_26, rel_biasDPjags.bin_HN_U_26,
                             MSE_DPjags.bin_HN_U_26, biasDPjags_basemu.bin_HN_U_26, MSE_DPjags_basemu.bin_HN_U_26,
                             bias_tau2DPjags.bin_HN_U_26,rel_bias_tauDPjags.bin_HN_U_26, abs_bias_tauDPjags.bin_HN_U_26, 
                             MSE_DPjags.bin.tau2_HN_U_26, MSE_DPjags.bin.tau_HN_U_26,
                             bias_basetau2DPjags.bin_HN_U_26, MSE_DPjags.bin.basetau2_HN_U_26, 
                             avg_abs_biasDPjags.bin_HN_U_26, avg_rel_biasDPjags.bin_HN_U_26, 
                             avg_biasDPjags_basemu.bin_HN_U_26, 
                             avg_bias_tau2DPjags.bin_HN_U_26,pbias_tau2DPjags.bin_HN_U_26,avg_rel_bias_tauDPjags.bin_HN_U_26,avg_abs_bias_tauDPjags.bin_HN_U_26,
                             prel_bias_tauDPjags.bin_HN_U_26,pabs_bias_tauDPjags.bin_HN_U_26,
                             avg_bias_basetau2DPjags.bin_HN_U_26,
                             rel_bias_tauDPjags.bin_HN_U_26,abs_bias_tauDPjags.bin_HN_U_26,
                             avg_MSE_DPjags.bin_HN_U_26,avg_MSE_DPjags_basemu.bin_HN_U_26,avg_MSE_DPjags.bin.tau2_HN_U_26, nMSE_DPjags.bin.tau2_HN_U_26,
                             avg_MSE_DPjags.bin.tau_HN_U_26, nMSE_DPjags.bin.tau_HN_U_26,
                             avg_MSE_DPjags.bin.basetau2_HN_U_26,
                             coverageDPjags.bin_HN_U_26 , coverageDPjags_basemu.bin_HN_U_26, coverage_tau2DPjags.bin_HN_U_26, coverage_basetau2DPjags.bin_HN_U_26, 
                             coverage_tauDPjags.bin_HN_U_26, 
                             coverDPjags.bin_HN_U_26, coverDPjags_basemu.bin_HN_U_26, cover_tau2DPjags.bin_HN_U_26, cover_tauDPjags.bin_HN_U_26,
                             cover_basetau2DPjags.bin_HN_U_26, RhatDP_out_HN_U_26)

Overlap_DPmodel_HN_U_26 <- cbind(overlapDP_HN_U_26, avg_overlapDP_HN_U_26)                    



write.csv(DP_bin_jags_HN_U_26, "C:\\Users\\Tianqi YU\\Desktop\\SIMULATION SCENARIOS FINAL\\NORMAL DATA\\DP MODEL\\rev_HN_UNIF_26_ResultsDP_jags_Normal_bin_data_k=14_mu=0_tau2=0.12.csv",row.names=FALSE )  
write.csv(listaDP_HN_U_26, "C:\\Users\\Tianqi YU\\Desktop\\SIMULATION SCENARIOS FINAL\\NORMAL DATA\\DP MODEL\\rev_HN_UNIF_26_Relative_eff_ResultsDP_jags_Normal_bin_data_k=14_mu=0_tau2=0.12.csv",row.names=FALSE )  
write.csv(listaDP_rel_eff_HN_U_26, "C:\\Users\\Tianqi YU\\Desktop\\SIMULATION SCENARIOS FINAL\\NORMAL DATA\\DP MODEL\\rev_HN_UNIF_26_Bias_Relative_eff_ResultsDP_jags_Normal_bin_data_k=14_mu=0_tau2=0.12.csv",row.names=FALSE )  
write.csv(listaDP_weights_HN_U_26, "C:\\Users\\Tianqi YU\\Desktop\\SIMULATION SCENARIOS FINAL\\NORMAL DATA\\DP MODEL\\rev_HN_UNIF_26_Weights_in_clusters_jags_DP_jags_Normal_bin_data_k=14_mu=0_tau2=0.12.csv",row.names=TRUE )  
write.csv(Overlap_DPmodel_HN_U_26, "C:\\Users\\Tianqi YU\\Desktop\\SIMULATION SCENARIOS FINAL\\NORMAL DATA\\DP MODEL\\rev_Overlap_HN_UNIF_26_DPmodel_Normal_bin_data_k=14_mu=0_tau2=0.12.csv",row.names=TRUE )  
#
###################### DIRICHLET PROCESS MODEL ######################
############# BINOMIAL-DP-51(HN/Unif) #############
cat("model{

  for( i in 1: ns ) {

    m[i] ~ dnorm(0,.0001)
    
    logit(p1[i]) <- m[i]
    r1[i] ~ dbin(p1[i],n1[i])
    
    logit(p2[i]) <- m[i] + delta12[i]
    r2[i] ~ dbin(p2[i],n2[i])
    
    delta12[i] <- theta[Z[i]]
    
    Z[i] ~ dcat(p[]) #Z is an integer variable
  }
  # Constructive DP
  #stick-breaking prior
  p[1] <- r[1]
  for (j in 2:(N-1)){
   p[j] <- r[j]*(1-r[j-1])*p[j-1]/r[j-1]
  }
  for (k in 1:(N-1)){
   r[k] ~ dbeta(1,alpha)T(0,0.99)
   }
  #assumption to ensure sum p[] is 1 Ishwaran truncation
  ps <- sum(p[1:(N-1)])
  for(k in N:N){
   p[k]<-1-ps
   }
  # Baseline distribution
  for(k in 1:N){
   theta[k] ~ dnorm(basemu ,basetau1)
  }
  #### priors 

  basemu~dnorm(0,0.0001)
  basetau1 <- 1/basetau_sqr
  basetau_sqr <- basetau*basetau
  basetau~dnorm(0,1)T(0,)

  # DPP parameter prior
  alpha~dunif(0.3,10)
  
  # Random effects distribution mean#
 
  for(i in 1:N){
   meancl[i]<-p[i]*theta[i]
   }
  poptrue<-sum(meancl[])  ### the E[X]
  
  # Random effects distribution variance #
  for(i in 1:N){
   mom2[i]<-p[i]*theta[i]*theta[i]  ####E[X2]
   }
  mom2.true<-sum(mom2[])
  var.true<-mom2.true-(poptrue*poptrue) ###E[X2] - E[X]2
  
  # Programming for calculating summary statistics #
  for(i in 1:ns){
   for (j in 1:N){
    SC[i,j] <- equals(j,Z[i])
   }
  }
  # total clusters K#
  for (j in 1:N){
   cl[j] <-step(sum(SC[,j])-1)
   }
  K<-sum(cl[])
  
}",file="DPmodel2.bin_HN_U_51.txt")
modfile = 'DPmodel2.bin_HN_U_51.txt'
N.sim =1000
mu = 0
tau2 = 0.12
tau = sqrt(tau2)

cl = makeCluster(19)
registerDoParallel(cl)
DPresults1_HN_U_51 <- list()

DPresults1_HN_U_51 = foreach(i = 1:N.sim,.packages = c("rjags", "R2jags","tidyverse")) %dopar% {
  set.seed(2104)
  run.model = jags(
    data =list(ns = nrow(p[[i]]),
               r1 = p[[i]]$ci,
               r2 = p[[i]]$ti,
               n1 = p[[i]]$nci,
               n2 = p[[i]]$nti,
               N= 51
               
    ) ,  
    inits = NULL,
    parameters.to.save = c(
      "basemu",    ## tau of the base Normal distribution
      "basetau",   ## mu of the base Normal distribution
      "poptrue",   ##the overall mu
      "var.true", ## the heterogeneity 
      "delta12", #### the random effects
      "K",       ### the total number of clusters 
      "p",     ## the weights of the process 
      "alpha"  ## the concentration parameter
    ), 
    
    n.chains = 2,
    n.iter = 50000,
    n.burnin = 10000,
    DIC = T,
    model.file = modfile
    
  )
  
}
stopCluster(cl)


DPresults_HN_U_51 <- list()
posteriorDP_HN_U_51 <- list()
postDP_HN_U_51 <- list()
for(i in 1:N.sim){
  DPresults_HN_U_51[[i]] <- as.data.frame(DPresults1_HN_U_51[[i]]$BUGSoutput$summary) 
  posteriorDP_HN_U_51[[i]] = as.data.frame(DPresults1_HN_U_51[[i]]$BUGSoutput$sims.matrix)
  postDP_HN_U_51[[i]] <- posteriorDP_HN_U_51[[i]]$poptrue
}

DP_results_HN_U_51 <- list()
for(i in 1:N.sim){
  DP_results_HN_U_51[[i]] <- DPresults_HN_U_51[[i]] 
}

############# extraction of the parameters of interest ###########
alpha_HN_U_51 <- c()
LB_alpha_HN_U_51 <- c()
UB_alpha_HN_U_51 <- c()
mu_DP_HN_U_51 <- c()
LB_mu_DP_HN_U_51 <- c()
UB_mu_DP_HN_U_51 <- c()
tau2_DP_HN_U_51 <- c()
LB_tau2_DP_HN_U_51 <- c()
UB_tau2_DP_HN_U_51 <- c()
precDP_mu_HN_U_51 <- c()
precDP_tau2_HN_U_51 <- c()
prec_basemu_HN_U_51 <- c()
prec_basetau2_HN_U_51 <- c()
prec_alpha_HN_U_51 <- c()
base_mu_HN_U_51 <- c()
base_tau_HN_U_51 <- c()
base_tau2_HN_U_51 <- c()
LB_base_mu_HN_U_51 <- c()
UB_base_mu_HN_U_51 <- c()
LB_base_tau2_HN_U_51 <- c()
UB_base_tau2_HN_U_51 <- c()
Rhat_muDP_HN_U_51 <- c()
Rhat_tau2DP_HN_U_51 <- c()
Rhat_basemu_HN_U_51 <- c()
Rhat_basetau_HN_U_51 <- c()
Rhat_alpha_HN_U_51 <- c()
Rhat_deltaDP_HN_U_51 <- list()
median_K_HN_U_51 <- list()
LB_K_HN_U_51 <- list()
UB_K_HN_U_51 <- list()
rel_effDP_HN_U_51 <- list()
LB_rel_effDP_HN_U_51 <- list()
UB_rel_effDP_HN_U_51 <- list()
sd_rel_effDP_HN_U_51 <- list()
pi_HN_U_51 <- list()
f11_HN_U_51 <- list()
m11_HN_U_51 <- list()
f22_HN_U_51 <- list()
m22_HN_U_51 <- list()
fdd_HN_U_51 <- list()
mdd_HN_U_51 <- list()
f33_HN_U_51 <- list()
m33_HN_U_51 <- list()
f44_HN_U_51 <- list()
m44_HN_U_51 <- list()
f55_HN_U_51 <- list()
m55_HN_U_51 <- list()
fcl_HN_U_51 <- list()
mcl_HN_U_51 <- list()
fp_HN_U_51 <- list()
mp_HN_U_51 <- list()
listaDP1_HN_U_51 <- list()
listaDP_weights_HN_U_51 <- list()
for(i in 1:N.sim){
  f11_HN_U_51[[i]] <- grepl("poptrue", row.names(DP_results_HN_U_51[[i]]))
  m11_HN_U_51[[i]] <- DP_results_HN_U_51[[i]][f11_HN_U_51[[i]],]
  mu_DP_HN_U_51[i] <- m11_HN_U_51[[i]]$`50%`
  LB_mu_DP_HN_U_51[i] <- m11_HN_U_51[[i]]$`2.5%`
  UB_mu_DP_HN_U_51[i] <- m11_HN_U_51[[i]]$`97.5%`
  Rhat_muDP_HN_U_51[i] <- m11_HN_U_51[[i]]$Rhat
  precDP_mu_HN_U_51[i] <- UB_mu_DP_HN_U_51[i] - LB_mu_DP_HN_U_51[i]
  
  f22_HN_U_51[[i]] <- grepl("var.true", row.names(DP_results_HN_U_51[[i]]))
  m22_HN_U_51[[i]] <- DP_results_HN_U_51[[i]][f22_HN_U_51[[i]],]
  tau2_DP_HN_U_51[i] <- m22_HN_U_51[[i]]$`50%`
  LB_tau2_DP_HN_U_51[i] <- m22_HN_U_51[[i]]$`2.5%`
  UB_tau2_DP_HN_U_51[i] <- m22_HN_U_51[[i]]$`97.5%`
  Rhat_tau2DP_HN_U_51[i] <- m22_HN_U_51[[i]]$Rhat
  precDP_tau2_HN_U_51[i] <- UB_tau2_DP_HN_U_51[i] -  LB_tau2_DP_HN_U_51[i]
  
  fdd_HN_U_51[[i]] <- grepl("delta12", row.names(DP_results_HN_U_51[[i]]))
  mdd_HN_U_51[[i]] <- DP_results_HN_U_51[[i]][fdd_HN_U_51[[i]],]
  rel_effDP_HN_U_51[[i]] <- mdd_HN_U_51[[i]]$`50%`
  LB_rel_effDP_HN_U_51[[i]] <- mdd_HN_U_51[[i]]$`2.5%`
  UB_rel_effDP_HN_U_51[[i]] <- mdd_HN_U_51[[i]]$`97.5%`
  sd_rel_effDP_HN_U_51[[i]] <- mdd_HN_U_51[[i]]$sd
  Rhat_deltaDP_HN_U_51[[i]] <- mdd_HN_U_51[[i]]$Rhat
  
  f33_HN_U_51[[i]] <- grepl("basemu", row.names(DP_results_HN_U_51[[i]]))
  m33_HN_U_51[[i]] <- DP_results_HN_U_51[[i]][f33_HN_U_51[[i]],]
  base_mu_HN_U_51[i] <- m33_HN_U_51[[i]]$`50%`
  LB_base_mu_HN_U_51[i] <- m33_HN_U_51[[i]]$`2.5%`
  UB_base_mu_HN_U_51[i] <- m33_HN_U_51[[i]]$`97.5%`
  Rhat_basemu_HN_U_51[i] <- m33_HN_U_51[[i]]$Rhat
  prec_basemu_HN_U_51[i] <- UB_base_mu_HN_U_51[i] - LB_base_mu_HN_U_51[i]
  
  f44_HN_U_51[[i]] <- grepl("basetau", row.names(DP_results_HN_U_51[[i]]))
  m44_HN_U_51[[i]] <- DP_results_HN_U_51[[i]][f44_HN_U_51[[i]],]
  base_tau_HN_U_51[i] <- m44_HN_U_51[[i]]$`50%`
  base_tau2_HN_U_51[i] <- (m44_HN_U_51[[i]]$`50%`)^2
  LB_base_tau2_HN_U_51[i] <- (m44_HN_U_51[[i]]$`2.5%`)^2
  UB_base_tau2_HN_U_51[i] <- (m44_HN_U_51[[i]]$`97.5%`)^2
  Rhat_basetau_HN_U_51[i] <- (m44_HN_U_51[[i]]$Rhat)^2
  prec_basetau2_HN_U_51[i] <- UB_base_tau2_HN_U_51[i] - LB_base_tau2_HN_U_51[i]
  
  f55_HN_U_51[[i]] <- grepl("alpha", row.names(DP_results_HN_U_51[[i]]))
  m55_HN_U_51[[i]] <- DP_results_HN_U_51[[i]][f55_HN_U_51[[i]],]
  alpha_HN_U_51[i] <- m55_HN_U_51[[i]]$`50%`
  LB_alpha_HN_U_51[i] <- m55_HN_U_51[[i]]$`2.5%`
  UB_alpha_HN_U_51[i] <- m55_HN_U_51[[i]]$`97.5%`
  Rhat_alpha_HN_U_51[i] <- m55_HN_U_51[[i]]$Rhat
  prec_alpha_HN_U_51[i] <- UB_alpha_HN_U_51[i] - LB_alpha_HN_U_51[i]
  
  fcl_HN_U_51[[i]] <- grepl("K", row.names(DP_results_HN_U_51[[i]]))  
  mcl_HN_U_51[[i]] <- DP_results_HN_U_51[[i]][fcl_HN_U_51[[i]],]
  median_K_HN_U_51[[i]] <- mcl_HN_U_51[[i]]$`50%` 
  LB_K_HN_U_51[[i]] <- mcl_HN_U_51[[i]]$`2.5%`
  UB_K_HN_U_51[[i]] <- mcl_HN_U_51[[i]]$`97.5%`
  
  fp_HN_U_51[[i]] <- grepl("p", row.names(DP_results_HN_U_51[[i]]))
  mp_HN_U_51[[i]] <- DP_results_HN_U_51[[i]][fp_HN_U_51[[i]],]
  mp_HN_U_51[[i]] <- mp_HN_U_51[[i]][!grepl("pop", row.names(mp_HN_U_51[[i]])),]
  mp_HN_U_51[[i]] <- mp_HN_U_51[[i]][!grepl("alpha", row.names(mp_HN_U_51[[i]])),]
  pi_HN_U_51[[i]] <- mp_HN_U_51[[i]]$mean
  
  listaDP1_HN_U_51[[i]] <- cbind.data.frame(rel_effDP_HN_U_51[[i]],sd_rel_effDP_HN_U_51[[i]], LB_rel_effDP_HN_U_51[[i]],UB_rel_effDP_HN_U_51[[i]],Rhat_deltaDP_HN_U_51[[i]])
  listaDP_weights_HN_U_51[[i]] <- cbind.data.frame(pi_HN_U_51[[i]])
  
}

numclus_HN_U_51 <- unlist(median_K_HN_U_51)
LB_K_HN_U_51 <- unlist(LB_K_HN_U_51)
UB_K_HN_U_51 <- unlist(UB_K_HN_U_51)

mean_alpha_HN_U_51 = mean(alpha_HN_U_51)
mean_alpha_LB_HN_U_51 = mean(LB_alpha_HN_U_51)
mean_alpha_UB_HN_U_51 = mean(UB_alpha_HN_U_51)

tau_DP_HN_U_51 = sqrt(tau2_DP_HN_U_51)
LB_tau_DP_HN_U_51 = sqrt(LB_tau2_DP_HN_U_51)
UB_tau_DP_HN_U_51 = sqrt(UB_tau2_DP_HN_U_51)
precDP_tau_HN_U_51 = UB_tau_DP_HN_U_51 -  LB_tau_DP_HN_U_51


RhatsDP_HN_U_51 <- cbind.data.frame(base_mu_HN_U_51,LB_base_mu_HN_U_51,UB_base_mu_HN_U_51,base_tau_HN_U_51,base_tau2_HN_U_51, LB_base_tau2_HN_U_51,UB_base_tau2_HN_U_51,
                                    mu_DP_HN_U_51, LB_mu_DP_HN_U_51, UB_mu_DP_HN_U_51, tau2_DP_HN_U_51, LB_tau2_DP_HN_U_51, UB_tau2_DP_HN_U_51,
                                    tau_DP_HN_U_51, LB_tau_DP_HN_U_51, UB_tau_DP_HN_U_51,
                                    precDP_mu_HN_U_51, precDP_tau2_HN_U_51, prec_basemu_HN_U_51, prec_basetau2_HN_U_51,prec_alpha_HN_U_51,  alpha_HN_U_51 , LB_alpha_HN_U_51, UB_alpha_HN_U_51,
                                    mean_alpha_HN_U_51, mean_alpha_LB_HN_U_51, mean_alpha_UB_HN_U_51,
                                    Rhat_muDP_HN_U_51, Rhat_tau2DP_HN_U_51, Rhat_basemu_HN_U_51, Rhat_basetau_HN_U_51, Rhat_alpha_HN_U_51, numclus_HN_U_51, LB_K_HN_U_51, UB_K_HN_U_51)

##########REMOVE Rhats > 1.05 ##############

condition1_HN_U_51 <- which(RhatsDP_HN_U_51$Rhat_muDP_HN_U_51 > 1.05)
condition2_HN_U_51 <- which(RhatsDP_HN_U_51$Rhat_tau2DP_HN_U_51 > 1.05)
condition3_HN_U_51 <- which(RhatsDP_HN_U_51$Rhat_basemu_HN_U_51 > 1.05)
condition4_HN_U_51 <- which(RhatsDP_HN_U_51$Rhat_basetau_HN_U_51 > 1.05)
condition5_HN_U_51 <- which(RhatsDP_HN_U_51$Rhat_alpha_HN_U_51 > 1.05)


dist.condDP_HN_U_51 = c(condition1_HN_U_51,condition2_HN_U_51,condition3_HN_U_51,condition4_HN_U_51,condition5_HN_U_51)
dist.condDP_HN_U_51 = unique(dist.condDP_HN_U_51)

RhatDP_out_HN_U_51 = round((length(dist.condDP_HN_U_51)/N.sim), 4)

############### Extract and remove the datasets with Rhat > 1.05 #########
if (length(dist.condDP_HN_U_51)== 0) {
  RhatsDP_HN_U_51 <- RhatsDP_HN_U_51
  listaDP_HN_U_51 <- listaDP1_HN_U_51
  listaDP_weights_HN_U_51 <- listaDP_weights_HN_U_51
  N.sim <- nrow(RhatsDP_HN_U_51)
  pDP_HN_U_51 <- p
  
} else {
  RhatsDP_HN_U_51 <- RhatsDP_HN_U_51[-dist.condDP_HN_U_51, ]
  listaDP_HN_U_51 <- listaDP1_HN_U_51[-dist.condDP_HN_U_51]
  listaDP_weights_HN_U_51 <- listaDP_weights_HN_U_51[-dist.condDP_HN_U_51]
  N.sim <- nrow(RhatsDP_HN_U_51)
  pDP_HN_U_51 <- p[-dist.condDP_HN_U_51]
}

count_DP_HN_U_51 = list()
for(i in length(listaDP_HN_U_51)){
  count_DP_HN_U_51[[i]] = which(listaDP_HN_U_51[[i]]$`Rhat_deltaDP_HN_U_51[[i]]` > 1.05 )
  tell_meDP_HN_U_51 = which(count_DP_HN_U_51[[i]] != 0)
}

tell_meDP_HN_U_51

if(length(tell_meDP_HN_U_51) == 0){
  RhatsDP_HN_U_51 <- RhatsDP_HN_U_51
  N.sim <- nrow(RhatsDP_HN_U_51)
  listaDP_HN_U_51 <- listaDP_HN_U_51
  listaDP_weights_HN_U_51 <- listaDP_weights_HN_U_51
  pDP_HN_U_51 <- pDP_HN_U_51
  RhatDP_out_HN_U_51 <- RhatDP_out_HN_U_51
} else {
  RhatsDP_HN_U_51 <- RhatsDP_HN_U_51[-tell_meDP_HN_U_51, ]
  listaDP_HN_U_51 <- listaDP_HN_U_51[-tell_meDP_HN_U_51]
  listaDP_weights_HN_U_51 <- listaDP_weights_HN_U_51[-tell_meDP_HN_U_51]
  N.sim = nrow(RhatsDP_HN_U_51)
  pDP_HN_U_51 <- pDP_HN_U_51[-tell_meDP_HN_U_51]
  RhatDP_out_HN_U_51 <- RhatDP_out_HN_U_51 + ((length(tell_meDP_HN_U_51))/N.sim)
}

RhatDP_out_HN_U_51

#### AVERAGE ABSOLUTE BIAS OF mu ######
abs_biasDPjags.bin_HN_U_51 <- c()
for (i in 1:N.sim){
  abs_biasDPjags.bin_HN_U_51[i] <- abs( RhatsDP_HN_U_51$mu_DP_HN_U_51[i] - mu) 
}
avg_abs_biasDPjags.bin_HN_U_51 <- round(mean(abs_biasDPjags.bin_HN_U_51),2)

#### AVERAGE RELATIVE BIAS OF mu ######
rel_biasDPjags.bin_HN_U_51 <- c()
for (i in 1:N.sim){
  rel_biasDPjags.bin_HN_U_51[i] <- ( RhatsDP_HN_U_51$mu_DP_HN_U_51[i] - mu) 
}
avg_rel_biasDPjags.bin_HN_U_51 <- round(mean(rel_biasDPjags.bin_HN_U_51),2)


######  MSE OF mu ##########
MSE_DPjags.bin_HN_U_51 <- c()
for (i in 1:N.sim){
  MSE_DPjags.bin_HN_U_51[i] <- RMSE(estimate = RhatsDP_HN_U_51$mu_DP_HN_U_51[i],
                                    parameter = mu,
                                    type = "RMSE",
                                    MSE = TRUE   )      
}
avg_MSE_DPjags.bin_HN_U_51 <- round(mean(MSE_DPjags.bin_HN_U_51), 4)

############ COVERAGE OF mu ########
interval_contains_true_mean <- function(i) { 
  ( mu >= RhatsDP_HN_U_51$LB_mu_DP_HN_U_51[i]) && ( mu <= RhatsDP_HN_U_51$UB_mu_DP_HN_U_51[i])
}
interval_contains_true_mean(i)

coverageDPjags.bin_HN_U_51 <- c()
for (i in 1:N.sim){
  if (interval_contains_true_mean(i) == TRUE){
    coverageDPjags.bin_HN_U_51[i] = 1
  }
  else {
    coverageDPjags.bin_HN_U_51[i]=0
  }
  print(coverageDPjags.bin_HN_U_51[i])
}
coverageDPjags.bin_HN_U_51
coverDPjags.bin_HN_U_51 <- round(mean(coverageDPjags.bin_HN_U_51)*100, 2)

## BIAS OF base_mu ######
biasDPjags_basemu.bin_HN_U_51 <- c()
for (i in 1:N.sim){
  biasDPjags_basemu.bin_HN_U_51[i] <- abs( RhatsDP_HN_U_51$base_mu_HN_U_51[i] - mu) 
}
avg_biasDPjags_basemu.bin_HN_U_51 <- round(mean(biasDPjags_basemu.bin_HN_U_51),2)

###### MSE OF base_mu ##########
MSE_DPjags_basemu.bin_HN_U_51 <- c()
for (i in 1:N.sim){
  MSE_DPjags_basemu.bin_HN_U_51[i] <- RMSE(estimate =RhatsDP_HN_U_51$base_mu_HN_U_51[i],
                                           parameter = mu,
                                           type = "RMSE",
                                           MSE = TRUE      
  )
}
avg_MSE_DPjags_basemu.bin_HN_U_51 <- round(mean(MSE_DPjags_basemu.bin_HN_U_51), 4)

############ COVERAGE OF base_mu ########
interval_contains_true_mean <- function(i) { 
  ( mu >= RhatsDP_HN_U_51$LB_base_mu_HN_U_51[i]) && ( mu <= RhatsDP_HN_U_51$UB_base_mu_HN_U_51[i])
}
interval_contains_true_mean(i)

coverageDPjags_basemu.bin_HN_U_51 <- c()
for (i in 1:N.sim){
  if (interval_contains_true_mean(i) == TRUE){
    coverageDPjags_basemu.bin_HN_U_51[i] = 1
  }
  else {
    coverageDPjags_basemu.bin_HN_U_51[i]=0
  }
  print(coverageDPjags_basemu.bin_HN_U_51[i])
}
coverageDPjags_basemu.bin_HN_U_51
coverDPjags_basemu.bin_HN_U_51 <- round(mean(coverageDPjags_basemu.bin_HN_U_51)*100, 2)

######## AVERAGE BIAS OF tau2 ###########
bias_tau2DPjags.bin_HN_U_51 <- c()
for (i in 1:N.sim){
  bias_tau2DPjags.bin_HN_U_51[i] <- abs(RhatsDP_HN_U_51$tau2_DP_HN_U_51[i]- tau2)
}
avg_bias_tau2DPjags.bin_HN_U_51 <- round(mean(bias_tau2DPjags.bin_HN_U_51),2)

pbias_tau2DPjags.bin_HN_U_51 <- round(mean(bias_tau2DPjags.bin_HN_U_51 / tau2),2)

######## MSE OF tau2 #############
MSE_DPjags.bin.tau2_HN_U_51 <- c()
for (i in 1:N.sim){
  MSE_DPjags.bin.tau2_HN_U_51[i] <- RMSE(estimate =RhatsDP_HN_U_51$tau2_DP_HN_U_51[i],
                                         parameter = tau2,
                                         type = "RMSE",
                                         MSE = TRUE         
  )
  
}
avg_MSE_DPjags.bin.tau2_HN_U_51 <- round(mean(MSE_DPjags.bin.tau2_HN_U_51),4)

nMSE_DPjags.bin.tau2_HN_U_51 <- round(mean(MSE_DPjags.bin.tau2_HN_U_51/ (tau2)^2),4)

############ COVERAGE OF tau2 ##################
interval_contains_true_mean <- function(i) { 
  ( tau2 >= RhatsDP_HN_U_51$LB_tau2_DP_HN_U_51[i]) && ( tau2 <= RhatsDP_HN_U_51$UB_tau2_DP_HN_U_51[i])
}
interval_contains_true_mean(i)

coverage_tau2DPjags.bin_HN_U_51 <- c()
for (i in 1:N.sim){
  if (interval_contains_true_mean(i) == TRUE){
    coverage_tau2DPjags.bin_HN_U_51[i] = 1
  }
  else  {
    coverage_tau2DPjags.bin_HN_U_51[i] = 0
  }
  print(coverage_tau2DPjags.bin_HN_U_51[i])
}
coverage_tau2DPjags.bin_HN_U_51
cover_tau2DPjags.bin_HN_U_51 <- round(mean(coverage_tau2DPjags.bin_HN_U_51)*100, 2)


######## AVERAGE RELATIVE  BIAS OF tau ###########
rel_bias_tauDPjags.bin_HN_U_51 <- c()
for (i in 1:N.sim){
  rel_bias_tauDPjags.bin_HN_U_51[i] <- (RhatsDP_HN_U_51$tau_DP_HN_U_51[i]- tau)
}
avg_rel_bias_tauDPjags.bin_HN_U_51 <- round(mean(rel_bias_tauDPjags.bin_HN_U_51),2)

prel_bias_tauDPjags.bin_HN_U_51 <- round(mean(rel_bias_tauDPjags.bin_HN_U_51 / tau)*100,2)

######## AVERAGE ABSOLUTE  BIAS OF tau ###########
abs_bias_tauDPjags.bin_HN_U_51 <- c()
for (i in 1:N.sim){
  abs_bias_tauDPjags.bin_HN_U_51[i] <- abs(RhatsDP_HN_U_51$tau_DP_HN_U_51[i]- tau)
}
avg_abs_bias_tauDPjags.bin_HN_U_51 <- round(mean(abs_bias_tauDPjags.bin_HN_U_51),2)

pabs_bias_tauDPjags.bin_HN_U_51 <- round(mean(abs_bias_tauDPjags.bin_HN_U_51 / tau)*100,2)

######## MSE OF tau #############
MSE_DPjags.bin.tau_HN_U_51 <- c()
for (i in 1:N.sim){
  MSE_DPjags.bin.tau_HN_U_51[i] <- RMSE(estimate =RhatsDP_HN_U_51$tau_DP_HN_U_51[i],
                                        parameter = tau,
                                        type = "RMSE",
                                        MSE = TRUE         
  )
  
}
avg_MSE_DPjags.bin.tau_HN_U_51 <- round(mean(MSE_DPjags.bin.tau_HN_U_51),4)

nMSE_DPjags.bin.tau_HN_U_51 <- round(mean(MSE_DPjags.bin.tau_HN_U_51/ tau2),4)

############ COVERAGE OF tau ##################
interval_contains_true_mean <- function(i) { 
  ( tau >= RhatsDP_HN_U_51$LB_tau_DP_HN_U_51[i]) && ( tau <= RhatsDP_HN_U_51$UB_tau_DP_HN_U_51[i])
}
interval_contains_true_mean(i)

coverage_tauDPjags.bin_HN_U_51 <- c()
for (i in 1:N.sim){
  if (interval_contains_true_mean(i) == TRUE){
    coverage_tauDPjags.bin_HN_U_51[i] = 1
  }
  else  {
    coverage_tauDPjags.bin_HN_U_51[i] = 0
  }
  print(coverage_tauDPjags.bin_HN_U_51[i])
}
coverage_tauDPjags.bin_HN_U_51
cover_tauDPjags.bin_HN_U_51 <- round(mean(coverage_tauDPjags.bin_HN_U_51)*100, 2)


######## AVERAGE BIAS OF base_tau2 ###########
bias_basetau2DPjags.bin_HN_U_51 <- c()
for (i in 1:N.sim){
  bias_basetau2DPjags.bin_HN_U_51[i] <- abs(RhatsDP_HN_U_51$base_tau2_HN_U_51[i]- tau2)
}
avg_bias_basetau2DPjags.bin_HN_U_51 <- round(mean(bias_basetau2DPjags.bin_HN_U_51),2)

######## MSE OF base_tau2 #############
MSE_DPjags.bin.basetau2_HN_U_51 <- c()
for (i in 1:N.sim){
  MSE_DPjags.bin.basetau2_HN_U_51[i] <- RMSE(estimate =RhatsDP_HN_U_51$base_tau2_HN_U_51[i],
                                             parameter = tau2,
                                             type = "RMSE",
                                             MSE = TRUE         
  )
  
}
avg_MSE_DPjags.bin.basetau2_HN_U_51 <- round(mean(MSE_DPjags.bin.basetau2_HN_U_51),4)

############ COVERAGE OF base_tau2 ##################
interval_contains_true_mean <- function(i) {
  ( tau2 >= RhatsDP_HN_U_51$LB_base_tau2_HN_U_51[i]) && ( tau2 <= RhatsDP_HN_U_51$UB_base_tau2_HN_U_51[i])
}
interval_contains_true_mean(i)

coverage_basetau2DPjags.bin_HN_U_51 <- c()
for (i in 1:N.sim){
  if (interval_contains_true_mean(i) == TRUE){
    coverage_basetau2DPjags.bin_HN_U_51[i] = 1
  }
  else  {
    coverage_basetau2DPjags.bin_HN_U_51[i] = 0
  }
  print(coverage_basetau2DPjags.bin_HN_U_51[i])
}
coverage_basetau2DPjags.bin_HN_U_51
cover_basetau2DPjags.bin_HN_U_51 <- round(mean(coverage_basetau2DPjags.bin_HN_U_51)*100, 2)


######################## ABSOLUTE BIAS OF RELATIVE EFFECTS ##################
abs_bias_relDP_HN_U_51 <- list()
mean_abs_bias_relDP_HN_U_51 <- list()
for(i in 1:N.sim){
  for(j in 1:nrow(pDP_HN_U_51[[i]])){
    abs_bias_relDP_HN_U_51[[i]] <- round(abs(listaDP_HN_U_51[[i]]$`rel_effDP_HN_U_51[[i]]` - pDP_HN_U_51[[i]]$true_eff ),2)
    mean_abs_bias_relDP_HN_U_51[[i]] <- mean(abs_bias_relDP_HN_U_51[[i]])
  }
}
#mean_bias_relDP

################## mean bias of the mean biases of each data set for relative effects  ###################
mean_abs_bias_of_mean_abs_bias_relDP_HN_U_51 <- unlist(mean_abs_bias_relDP_HN_U_51)
mean_abs_bias_of_mean_abs_bias_relDP_HN_U_51  <- round(mean(mean_abs_bias_of_mean_abs_bias_relDP_HN_U_51),2)

######################## RELATIVE BIAS OF RELATIVE EFFECTS ##################
rel_bias_relDP_HN_U_51 <- list()
mean_rel_bias_relDP_HN_U_51 <- list()
for(i in 1:N.sim){
  for(j in 1:nrow(pDP_HN_U_51[[i]])){
    rel_bias_relDP_HN_U_51[[i]] <- round((listaDP_HN_U_51[[i]]$`rel_effDP_HN_U_51[[i]]` - pDP_HN_U_51[[i]]$true_eff ),2)
    mean_rel_bias_relDP_HN_U_51[[i]] <- mean(rel_bias_relDP_HN_U_51[[i]])
  }
}
#mean_bias_relDP

################## mean bias of the mean biases of each data set for relative effects  ###################
mean_rel_bias_of_mean_rel_bias_relDP_HN_U_51 <- unlist(mean_rel_bias_relDP_HN_U_51)
mean_rel_bias_of_mean_rel_bias_relDP_HN_U_51  <- round(mean(mean_rel_bias_of_mean_rel_bias_relDP_HN_U_51),2)

#####
########################### FOR THE OVERLAPPING ######################
dataDP_HN_U_51 <- list()
overlapDP1_HN_U_51 <- list()
overlapDP_HN_U_51 <- list()
for(i in 1:N.sim){
  dataDP_HN_U_51[[i]] <- list(DP = listaDP_HN_U_51[[i]]$`rel_effDP_HN_U_51[[i]]` , True = pDP_HN_U_51[[i]]$true_eff  )
  overlapDP1_HN_U_51[[i]] <- overlap(dataDP_HN_U_51[[i]], type = "2")
  overlapDP_HN_U_51[[i]] <- overlapDP1_HN_U_51[[i]]$OV
}
overlapDP_HN_U_51 <- unlist(overlapDP_HN_U_51)
avg_overlapDP_HN_U_51 <- round(mean(overlapDP_HN_U_51),2)


listaDP_rel_eff_HN_U_51 <- cbind( mean_abs_bias_relDP_HN_U_51, mean_abs_bias_of_mean_abs_bias_relDP_HN_U_51,
                                  mean_rel_bias_relDP_HN_U_51, mean_rel_bias_of_mean_rel_bias_relDP_HN_U_51) 

######### SAVE THE OUTPUTS ###########
DP_bin_jags_HN_U_51 <- cbind(RhatsDP_HN_U_51,
                             abs_biasDPjags.bin_HN_U_51, rel_biasDPjags.bin_HN_U_51,
                             MSE_DPjags.bin_HN_U_51, biasDPjags_basemu.bin_HN_U_51, MSE_DPjags_basemu.bin_HN_U_51,
                             bias_tau2DPjags.bin_HN_U_51,rel_bias_tauDPjags.bin_HN_U_51, abs_bias_tauDPjags.bin_HN_U_51, 
                             MSE_DPjags.bin.tau2_HN_U_51, MSE_DPjags.bin.tau_HN_U_51,
                             bias_basetau2DPjags.bin_HN_U_51, MSE_DPjags.bin.basetau2_HN_U_51, 
                             avg_abs_biasDPjags.bin_HN_U_51, avg_rel_biasDPjags.bin_HN_U_51, 
                             avg_biasDPjags_basemu.bin_HN_U_51, 
                             avg_bias_tau2DPjags.bin_HN_U_51,pbias_tau2DPjags.bin_HN_U_51,avg_rel_bias_tauDPjags.bin_HN_U_51,avg_abs_bias_tauDPjags.bin_HN_U_51,
                             prel_bias_tauDPjags.bin_HN_U_51,pabs_bias_tauDPjags.bin_HN_U_51,
                             avg_bias_basetau2DPjags.bin_HN_U_51,
                             rel_bias_tauDPjags.bin_HN_U_51,abs_bias_tauDPjags.bin_HN_U_51,
                             avg_MSE_DPjags.bin_HN_U_51,avg_MSE_DPjags_basemu.bin_HN_U_51,avg_MSE_DPjags.bin.tau2_HN_U_51, nMSE_DPjags.bin.tau2_HN_U_51,
                             avg_MSE_DPjags.bin.tau_HN_U_51, nMSE_DPjags.bin.tau_HN_U_51,
                             avg_MSE_DPjags.bin.basetau2_HN_U_51,
                             coverageDPjags.bin_HN_U_51 , coverageDPjags_basemu.bin_HN_U_51, coverage_tau2DPjags.bin_HN_U_51, coverage_basetau2DPjags.bin_HN_U_51, 
                             coverage_tauDPjags.bin_HN_U_51, 
                             coverDPjags.bin_HN_U_51, coverDPjags_basemu.bin_HN_U_51, cover_tau2DPjags.bin_HN_U_51, cover_tauDPjags.bin_HN_U_51,
                             cover_basetau2DPjags.bin_HN_U_51, RhatDP_out_HN_U_51)

Overlap_DPmodel_HN_U_51 <- cbind(overlapDP_HN_U_51, avg_overlapDP_HN_U_51)                    

write.csv(DP_bin_jags_HN_U_51, "C:\\Users\\Tianqi YU\\Desktop\\SIMULATION SCENARIOS FINAL\\NORMAL DATA\\DP MODEL\\rev_HN-U51_ResultsDP_jags_Normal_bin_data_k=14_mu=0_tau2=0.12.csv",row.names=FALSE )  
write.csv(listaDP_HN_U_51, "C:\\Users\\Tianqi YU\\Desktop\\SIMULATION SCENARIOS FINAL\\NORMAL DATA\\DP MODEL\\rev_HN-U51_Relative_eff_ResultsDP_jags_Normal_bin_data_k=14_mu=0_tau2=0.12.csv",row.names=FALSE )  
write.csv(listaDP_rel_eff_HN_U_51, "C:\\Users\\Tianqi YU\\Desktop\\SIMULATION SCENARIOS FINAL\\NORMAL DATA\\DP MODEL\\rev_HN-U51_Bias_Relative_eff_ResultsDP_jags_Normal_bin_data_k=14_mu=0_tau2=0.12.csv",row.names=FALSE )  
write.csv(listaDP_weights_HN_U_51, "C:\\Users\\Tianqi YU\\Desktop\\SIMULATION SCENARIOS FINAL\\NORMAL DATA\\DP MODEL\\rev_HN-U51_Weights_in_clusters_jags_DP_jags_Normal_bin_data_k=14_mu=0_tau2=0.12.csv",row.names=TRUE )  
write.csv(Overlap_DPmodel_HN_U_51, "C:\\Users\\Tianqi YU\\Desktop\\SIMULATION SCENARIOS FINAL\\NORMAL DATA\\DP MODEL\\rev_HN-U51_Overlap_DPmodel_Normal_bin_data_k=14_mu=0_tau2=0.12.csv",row.names=TRUE )  

####################### DIRICHLET PROCESS MODEL ######################
############# BINOMIAL-DP-26(Unif/Unif) #############
cat("model{

  for( i in 1: ns ) {

    m[i] ~ dnorm(0,.0001)
    
    logit(p1[i]) <- m[i]
    r1[i] ~ dbin(p1[i],n1[i])
    
    logit(p2[i]) <- m[i] + delta12[i]
    r2[i] ~ dbin(p2[i],n2[i])
    
    delta12[i] <- theta[Z[i]]
    
    Z[i] ~ dcat(p[]) #Z is an integer variable
  }
  # Constructive DP
  #stick-breaking prior
  p[1] <- r[1]
  for (j in 2:(N-1)){
   p[j] <- r[j]*(1-r[j-1])*p[j-1]/r[j-1]
  }
  for (k in 1:(N-1)){
   r[k] ~ dbeta(1,alpha)T(0,0.99)
   }
  #assumption to ensure sum p[] is 1 Ishwaran truncation
  ps <- sum(p[1:(N-1)])
  for(k in N:N){
   p[k]<-1-ps
   }
  # Baseline distribution
  for(k in 1:N){
   theta[k] ~ dnorm(basemu ,basetau1)
  }
  #### priors 

  basemu~dnorm(0,0.0001)
  basetau1 <- 1/basetau_sqr
  basetau_sqr <- basetau*basetau
  basetau ~ dunif(0,10)
  
  # DPP parameter prior
  alpha~dunif(0.3,5)
  
  # Random effects distribution mean#
 
  for(i in 1:N){
   meancl[i]<-p[i]*theta[i]
   }
  poptrue<-sum(meancl[])  ### the E[X]
  
  # Random effects distribution variance #
  for(i in 1:N){
   mom2[i]<-p[i]*theta[i]*theta[i]  ####E[X2]
   }
  mom2.true<-sum(mom2[])
  var.true<-mom2.true-(poptrue*poptrue) ###E[X2] - E[X]2
  
  # Programming for calculating summary statistics #
  for(i in 1:ns){
   for (j in 1:N){
    SC[i,j] <- equals(j,Z[i])
   }
  }
  # total clusters K#
  for (j in 1:N){
   cl[j] <-step(sum(SC[,j])-1)
   }
  K<-sum(cl[])
  
}",file="DPmodel2.bin_U_U_26.txt")
modfile = 'DPmodel2.bin_U_U_26.txt'
N.sim =1000
mu = 0
tau2 = 0.12
tau = sqrt(tau2)

cl = makeCluster(19)
registerDoParallel(cl)
DPresults1_U_U_26 <- list()

DPresults1_U_U_26 = foreach(i = 1:N.sim,.packages = c("rjags", "R2jags","tidyverse")) %dopar% {
  set.seed(2105)
  run.model = jags(
    data =list(ns = nrow(p[[i]]),
               r1 = p[[i]]$ci,
               r2 = p[[i]]$ti,
               n1 = p[[i]]$nci,
               n2 = p[[i]]$nti,
               N= 26
               
    ) ,  
    inits = NULL,
    parameters.to.save = c(
      "basemu",    ## tau of the base Normal distribution
      "basetau",   ## mu of the base Normal distribution
      "poptrue",   ##the overall mu
      "var.true", ## the heterogeneity 
      "delta12", #### the random effects
      "K",       ### the total number of clusters 
      "p",     ## the weights of the process 
      "alpha"  ## the concentration parameter
    ), 
    
    n.chains = 2,
    n.iter = 50000,
    n.burnin = 10000,
    DIC = T,
    model.file = modfile
    
  )
  
}
stopCluster(cl)


DPresults_U_U_26 <- list()
posteriorDP_U_U_26 <- list()
postDP_U_U_26 <- list()
for(i in 1:N.sim){
  DPresults_U_U_26[[i]] <- as.data.frame(DPresults1_U_U_26[[i]]$BUGSoutput$summary) 
  posteriorDP_U_U_26[[i]] = as.data.frame(DPresults1_U_U_26[[i]]$BUGSoutput$sims.matrix)
  postDP_U_U_26[[i]] <- posteriorDP_U_U_26[[i]]$poptrue
}

DP_results_U_U_26 <- list()
for(i in 1:N.sim){
  DP_results_U_U_26[[i]] <- DPresults_U_U_26[[i]] 
}

############# extraction of the parameters of interest ###########
alpha_U_U_26 <- c()
LB_alpha_U_U_26 <- c()
UB_alpha_U_U_26 <- c()
mu_DP_U_U_26 <- c()
LB_mu_DP_U_U_26 <- c()
UB_mu_DP_U_U_26 <- c()
tau2_DP_U_U_26 <- c()
LB_tau2_DP_U_U_26 <- c()
UB_tau2_DP_U_U_26 <- c()
precDP_mu_U_U_26 <- c()
precDP_tau2_U_U_26 <- c()
prec_basemu_U_U_26 <- c()
prec_basetau2_U_U_26 <- c()
prec_alpha_U_U_26 <- c()
base_mu_U_U_26 <- c()
base_tau_U_U_26 <- c()
base_tau2_U_U_26 <- c()
LB_base_mu_U_U_26 <- c()
UB_base_mu_U_U_26 <- c()
LB_base_tau2_U_U_26 <- c()
UB_base_tau2_U_U_26 <- c()
Rhat_muDP_U_U_26 <- c()
Rhat_tau2DP_U_U_26 <- c()
Rhat_basemu_U_U_26 <- c()
Rhat_basetau_U_U_26 <- c()
Rhat_alpha_U_U_26 <- c()
Rhat_deltaDP_U_U_26 <- list()
median_K_U_U_26 <- list()
LB_K_U_U_26 <- list()
UB_K_U_U_26 <- list()
rel_effDP_U_U_26 <- list()
LB_rel_effDP_U_U_26 <- list()
UB_rel_effDP_U_U_26 <- list()
sd_rel_effDP_U_U_26 <- list()
pi_U_U_26 <- list()
f11_U_U_26 <- list()
m11_U_U_26 <- list()
f22_U_U_26 <- list()
m22_U_U_26 <- list()
fdd_U_U_26 <- list()
mdd_U_U_26 <- list()
f33_U_U_26 <- list()
m33_U_U_26 <- list()
f44_U_U_26 <- list()
m44_U_U_26 <- list()
f55_U_U_26 <- list()
m55_U_U_26 <- list()
fcl_U_U_26 <- list()
mcl_U_U_26 <- list()
fp_U_U_26 <- list()
mp_U_U_26 <- list()
listaDP1_U_U_26 <- list()
listaDP_weights_U_U_26 <- list()
for(i in 1:N.sim){
  f11_U_U_26[[i]] <- grepl("poptrue", row.names(DP_results_U_U_26[[i]]))
  m11_U_U_26[[i]] <- DP_results_U_U_26[[i]][f11_U_U_26[[i]],]
  mu_DP_U_U_26[i] <- m11_U_U_26[[i]]$`50%`
  LB_mu_DP_U_U_26[i] <- m11_U_U_26[[i]]$`2.5%`
  UB_mu_DP_U_U_26[i] <- m11_U_U_26[[i]]$`97.5%`
  Rhat_muDP_U_U_26[i] <- m11_U_U_26[[i]]$Rhat
  precDP_mu_U_U_26[i] <- UB_mu_DP_U_U_26[i] - LB_mu_DP_U_U_26[i]
  
  f22_U_U_26[[i]] <- grepl("var.true", row.names(DP_results_U_U_26[[i]]))
  m22_U_U_26[[i]] <- DP_results_U_U_26[[i]][f22_U_U_26[[i]],]
  tau2_DP_U_U_26[i] <- m22_U_U_26[[i]]$`50%`
  LB_tau2_DP_U_U_26[i] <- m22_U_U_26[[i]]$`2.5%`
  UB_tau2_DP_U_U_26[i] <- m22_U_U_26[[i]]$`97.5%`
  Rhat_tau2DP_U_U_26[i] <- m22_U_U_26[[i]]$Rhat
  precDP_tau2_U_U_26[i] <- UB_tau2_DP_U_U_26[i] -  LB_tau2_DP_U_U_26[i]
  
  fdd_U_U_26[[i]] <- grepl("delta12", row.names(DP_results_U_U_26[[i]]))
  mdd_U_U_26[[i]] <- DP_results_U_U_26[[i]][fdd_U_U_26[[i]],]
  rel_effDP_U_U_26[[i]] <- mdd_U_U_26[[i]]$`50%`
  LB_rel_effDP_U_U_26[[i]] <- mdd_U_U_26[[i]]$`2.5%`
  UB_rel_effDP_U_U_26[[i]] <- mdd_U_U_26[[i]]$`97.5%`
  sd_rel_effDP_U_U_26[[i]] <- mdd_U_U_26[[i]]$sd
  Rhat_deltaDP_U_U_26[[i]] <- mdd_U_U_26[[i]]$Rhat
  
  f33_U_U_26[[i]] <- grepl("basemu", row.names(DP_results_U_U_26[[i]]))
  m33_U_U_26[[i]] <- DP_results_U_U_26[[i]][f33_U_U_26[[i]],]
  base_mu_U_U_26[i] <- m33_U_U_26[[i]]$`50%`
  LB_base_mu_U_U_26[i] <- m33_U_U_26[[i]]$`2.5%`
  UB_base_mu_U_U_26[i] <- m33_U_U_26[[i]]$`97.5%`
  Rhat_basemu_U_U_26[i] <- m33_U_U_26[[i]]$Rhat
  prec_basemu_U_U_26[i] <- UB_base_mu_U_U_26[i] - LB_base_mu_U_U_26[i]
  
  f44_U_U_26[[i]] <- grepl("basetau", row.names(DP_results_U_U_26[[i]]))
  m44_U_U_26[[i]] <- DP_results_U_U_26[[i]][f44_U_U_26[[i]],]
  base_tau_U_U_26[i] <- m44_U_U_26[[i]]$`50%`
  base_tau2_U_U_26[i] <- (m44_U_U_26[[i]]$`50%`)^2
  LB_base_tau2_U_U_26[i] <- (m44_U_U_26[[i]]$`2.5%`)^2
  UB_base_tau2_U_U_26[i] <- (m44_U_U_26[[i]]$`97.5%`)^2
  Rhat_basetau_U_U_26[i] <- (m44_U_U_26[[i]]$Rhat)^2
  prec_basetau2_U_U_26[i] <- UB_base_tau2_U_U_26[i] - LB_base_tau2_U_U_26[i]
  
  f55_U_U_26[[i]] <- grepl("alpha", row.names(DP_results_U_U_26[[i]]))
  m55_U_U_26[[i]] <- DP_results_U_U_26[[i]][f55_U_U_26[[i]],]
  alpha_U_U_26[i] <- m55_U_U_26[[i]]$`50%`
  LB_alpha_U_U_26[i] <- m55_U_U_26[[i]]$`2.5%`
  UB_alpha_U_U_26[i] <- m55_U_U_26[[i]]$`97.5%`
  Rhat_alpha_U_U_26[i] <- m55_U_U_26[[i]]$Rhat
  prec_alpha_U_U_26[i] <- UB_alpha_U_U_26[i] - LB_alpha_U_U_26[i]
  
  fcl_U_U_26[[i]] <- grepl("K", row.names(DP_results_U_U_26[[i]]))  
  mcl_U_U_26[[i]] <- DP_results_U_U_26[[i]][fcl_U_U_26[[i]],]
  median_K_U_U_26[[i]] <- mcl_U_U_26[[i]]$`50%` 
  LB_K_U_U_26[[i]] <- mcl_U_U_26[[i]]$`2.5%`
  UB_K_U_U_26[[i]] <- mcl_U_U_26[[i]]$`97.5%`
  
  fp_U_U_26[[i]] <- grepl("p", row.names(DP_results_U_U_26[[i]]))
  mp_U_U_26[[i]] <- DP_results_U_U_26[[i]][fp_U_U_26[[i]],]
  mp_U_U_26[[i]] <- mp_U_U_26[[i]][!grepl("pop", row.names(mp_U_U_26[[i]])),]
  mp_U_U_26[[i]] <- mp_U_U_26[[i]][!grepl("alpha", row.names(mp_U_U_26[[i]])),]
  pi_U_U_26[[i]] <- mp_U_U_26[[i]]$mean
  
  listaDP1_U_U_26[[i]] <- cbind.data.frame(rel_effDP_U_U_26[[i]],sd_rel_effDP_U_U_26[[i]], LB_rel_effDP_U_U_26[[i]],UB_rel_effDP_U_U_26[[i]],Rhat_deltaDP_U_U_26[[i]])
  listaDP_weights_U_U_26[[i]] <- cbind.data.frame(pi_U_U_26[[i]])
  
}

numclus_U_U_26 <- unlist(median_K_U_U_26)
LB_K_U_U_26 <- unlist(LB_K_U_U_26)
UB_K_U_U_26 <- unlist(UB_K_U_U_26)

mean_alpha_U_U_26 = mean(alpha_U_U_26)
mean_alpha_LB_U_U_26 = mean(LB_alpha_U_U_26)
mean_alpha_UB_U_U_26 = mean(UB_alpha_U_U_26)

tau_DP_U_U_26 = sqrt(tau2_DP_U_U_26)
LB_tau_DP_U_U_26 = sqrt(LB_tau2_DP_U_U_26)
UB_tau_DP_U_U_26 = sqrt(UB_tau2_DP_U_U_26)
precDP_tau_U_U_26 = UB_tau_DP_U_U_26 -  LB_tau_DP_U_U_26


RhatsDP_U_U_26 <- cbind.data.frame(base_mu_U_U_26,LB_base_mu_U_U_26,UB_base_mu_U_U_26,base_tau_U_U_26,base_tau2_U_U_26, LB_base_tau2_U_U_26,UB_base_tau2_U_U_26,
                                   mu_DP_U_U_26, LB_mu_DP_U_U_26, UB_mu_DP_U_U_26, tau2_DP_U_U_26, LB_tau2_DP_U_U_26, UB_tau2_DP_U_U_26,
                                   tau_DP_U_U_26, LB_tau_DP_U_U_26, UB_tau_DP_U_U_26,
                                   precDP_mu_U_U_26, precDP_tau2_U_U_26, prec_basemu_U_U_26, prec_basetau2_U_U_26,prec_alpha_U_U_26,  alpha_U_U_26 , LB_alpha_U_U_26, UB_alpha_U_U_26,
                                   mean_alpha_U_U_26, mean_alpha_LB_U_U_26, mean_alpha_UB_U_U_26,
                                   Rhat_muDP_U_U_26, Rhat_tau2DP_U_U_26, Rhat_basemu_U_U_26, Rhat_basetau_U_U_26, Rhat_alpha_U_U_26, numclus_U_U_26, LB_K_U_U_26, UB_K_U_U_26)

##########REMOVE Rhats > 1.05 ##############

condition1_U_U_26 <- which(RhatsDP_U_U_26$Rhat_muDP_U_U_26 > 1.05)
condition2_U_U_26 <- which(RhatsDP_U_U_26$Rhat_tau2DP_U_U_26 > 1.05)
condition3_U_U_26 <- which(RhatsDP_U_U_26$Rhat_basemu_U_U_26 > 1.05)
condition4_U_U_26 <- which(RhatsDP_U_U_26$Rhat_basetau_U_U_26 > 1.05)
condition5_U_U_26 <- which(RhatsDP_U_U_26$Rhat_alpha_U_U_26 > 1.05)


dist.condDP_U_U_26 = c(condition1_U_U_26,condition2_U_U_26,condition3_U_U_26,condition4_U_U_26,condition5_U_U_26)
dist.condDP_U_U_26 = unique(dist.condDP_U_U_26)

RhatDP_out_U_U_26 = round((length(dist.condDP_U_U_26)/N.sim), 4)

############### Extract and remove the datasets with Rhat > 1.05 #########
if (length(dist.condDP_U_U_26)== 0) {
  RhatsDP_U_U_26 <- RhatsDP_U_U_26
  listaDP_U_U_26 <- listaDP1_U_U_26
  listaDP_weights_U_U_26 <- listaDP_weights_U_U_26
  N.sim <- nrow(RhatsDP_U_U_26)
  pDP_U_U_26 <- p
  
} else {
  RhatsDP_U_U_26 <- RhatsDP_U_U_26[-dist.condDP_U_U_26, ]
  listaDP_U_U_26 <- listaDP1_U_U_26[-dist.condDP_U_U_26]
  listaDP_weights_U_U_26 <- listaDP_weights_U_U_26[-dist.condDP_U_U_26]
  N.sim <- nrow(RhatsDP_U_U_26)
  pDP_U_U_26 <- p[-dist.condDP_U_U_26]
}

count_DP_U_U_26 = list()
for(i in length(listaDP_U_U_26)){
  count_DP_U_U_26[[i]] = which(listaDP_U_U_26[[i]]$`Rhat_deltaDP_U_U_26[[i]]` > 1.05 )
  tell_meDP_U_U_26 = which(count_DP_U_U_26[[i]] != 0)
}

tell_meDP_U_U_26

if(length(tell_meDP_U_U_26) == 0){
  RhatsDP_U_U_26 <- RhatsDP_U_U_26
  N.sim <- nrow(RhatsDP_U_U_26)
  listaDP_U_U_26 <- listaDP_U_U_26
  listaDP_weights_U_U_26 <- listaDP_weights_U_U_26
  pDP_U_U_26 <- pDP_U_U_26
  RhatDP_out_U_U_26 <- RhatDP_out_U_U_26
} else {
  RhatsDP_U_U_26 <- RhatsDP_U_U_26[-tell_meDP_U_U_26, ]
  listaDP_U_U_26 <- listaDP_U_U_26[-tell_meDP_U_U_26]
  listaDP_weights_U_U_26 <- listaDP_weights_U_U_26[-tell_meDP_U_U_26]
  N.sim = nrow(RhatsDP_U_U_26)
  pDP_U_U_26 <- pDP_U_U_26[-tell_meDP_U_U_26]
  RhatDP_out_U_U_26 <- RhatDP_out_U_U_26 + ((length(tell_meDP_U_U_26))/N.sim)
}

RhatDP_out_U_U_26

#### AVERAGE ABSOLUTE BIAS OF mu ######
abs_biasDPjags.bin_U_U_26 <- c()
for (i in 1:N.sim){
  abs_biasDPjags.bin_U_U_26[i] <- abs( RhatsDP_U_U_26$mu_DP_U_U_26[i] - mu) 
}
avg_abs_biasDPjags.bin_U_U_26 <- round(mean(abs_biasDPjags.bin_U_U_26),2)

#### AVERAGE RELATIVE BIAS OF mu ######
rel_biasDPjags.bin_U_U_26 <- c()
for (i in 1:N.sim){
  rel_biasDPjags.bin_U_U_26[i] <- ( RhatsDP_U_U_26$mu_DP_U_U_26[i] - mu) 
}
avg_rel_biasDPjags.bin_U_U_26 <- round(mean(rel_biasDPjags.bin_U_U_26),2)



######  MSE OF mu ##########
MSE_DPjags.bin_U_U_26 <- c()
for (i in 1:N.sim){
  MSE_DPjags.bin_U_U_26[i] <- RMSE(estimate = RhatsDP_U_U_26$mu_DP_U_U_26[i],
                                   parameter = mu,
                                   type = "RMSE",
                                   MSE = TRUE )         
}
avg_MSE_DPjags.bin_U_U_26 <- round(mean(MSE_DPjags.bin_U_U_26), 4)

############ COVERAGE OF mu ########
interval_contains_true_mean <- function(i) { 
  ( mu >= RhatsDP_U_U_26$LB_mu_DP_U_U_26[i]) && ( mu <= RhatsDP_U_U_26$UB_mu_DP_U_U_26[i])
}
interval_contains_true_mean(i)

coverageDPjags.bin_U_U_26 <- c()
for (i in 1:N.sim){
  if (interval_contains_true_mean(i) == TRUE){
    coverageDPjags.bin_U_U_26[i] = 1
  }
  else {
    coverageDPjags.bin_U_U_26[i]=0
  }
  print(coverageDPjags.bin_U_U_26[i])
}
coverageDPjags.bin_U_U_26
coverDPjags.bin_U_U_26 <- round(mean(coverageDPjags.bin_U_U_26)*100, 2)

## BIAS OF base_mu ######
biasDPjags_basemu.bin_U_U_26 <- c()
for (i in 1:N.sim){
  biasDPjags_basemu.bin_U_U_26[i] <- abs( RhatsDP_U_U_26$base_mu_U_U_26[i] - mu) 
}
avg_biasDPjags_basemu.bin_U_U_26 <- round(mean(biasDPjags_basemu.bin_U_U_26),2)

###### MSE OF base_mu ##########
MSE_DPjags_basemu.bin_U_U_26 <- c()
for (i in 1:N.sim){
  MSE_DPjags_basemu.bin_U_U_26[i] <- RMSE(estimate =RhatsDP_U_U_26$base_mu_U_U_26[i],
                                          parameter = mu,
                                          type = "RMSE",
                                          MSE = TRUE      
  )
}
avg_MSE_DPjags_basemu.bin_U_U_26 <- round(mean(MSE_DPjags_basemu.bin_U_U_26), 4)

############ COVERAGE OF base_mu ########
interval_contains_true_mean <- function(i) { 
  ( mu >= RhatsDP_U_U_26$LB_base_mu_U_U_26[i]) && ( mu <= RhatsDP_U_U_26$UB_base_mu_U_U_26[i])
}
interval_contains_true_mean(i)

coverageDPjags_basemu.bin_U_U_26 <- c()
for (i in 1:N.sim){
  if (interval_contains_true_mean(i) == TRUE){
    coverageDPjags_basemu.bin_U_U_26[i] = 1
  }
  else {
    coverageDPjags_basemu.bin_U_U_26[i]=0
  }
  print(coverageDPjags_basemu.bin_U_U_26[i])
}
coverageDPjags_basemu.bin_U_U_26
coverDPjags_basemu.bin_U_U_26 <- round(mean(coverageDPjags_basemu.bin_U_U_26)*100, 2)

######## AVERAGE BIAS OF tau2 ###########
bias_tau2DPjags.bin_U_U_26 <- c()
for (i in 1:N.sim){
  bias_tau2DPjags.bin_U_U_26[i] <- abs(RhatsDP_U_U_26$tau2_DP_U_U_26[i]- tau2)
}
avg_bias_tau2DPjags.bin_U_U_26 <- round(mean(bias_tau2DPjags.bin_U_U_26),2)

pbias_tau2DPjags.bin_U_U_26 <- round(mean(bias_tau2DPjags.bin_U_U_26 / tau2),2)

######## MSE OF tau2 #############
MSE_DPjags.bin.tau2_U_U_26 <- c()
for (i in 1:N.sim){
  MSE_DPjags.bin.tau2_U_U_26[i] <- RMSE(estimate =RhatsDP_U_U_26$tau2_DP_U_U_26[i],
                                        parameter = tau2,
                                        type = "RMSE",
                                        MSE = TRUE         
  )
  
}
avg_MSE_DPjags.bin.tau2_U_U_26 <- round(mean(MSE_DPjags.bin.tau2_U_U_26),4)

nMSE_DPjags.bin.tau2_U_U_26 <- round(mean(MSE_DPjags.bin.tau2_U_U_26/ (tau2)^2),4)

############ COVERAGE OF tau2 ##################
interval_contains_true_mean <- function(i) { 
  ( tau2 >= RhatsDP_U_U_26$LB_tau2_DP_U_U_26[i]) && ( tau2 <= RhatsDP_U_U_26$UB_tau2_DP_U_U_26[i])
}
interval_contains_true_mean(i)

coverage_tau2DPjags.bin_U_U_26 <- c()
for (i in 1:N.sim){
  if (interval_contains_true_mean(i) == TRUE){
    coverage_tau2DPjags.bin_U_U_26[i] = 1
  }
  else  {
    coverage_tau2DPjags.bin_U_U_26[i] = 0
  }
  print(coverage_tau2DPjags.bin_U_U_26[i])
}
coverage_tau2DPjags.bin_U_U_26
cover_tau2DPjags.bin_U_U_26 <- round(mean(coverage_tau2DPjags.bin_U_U_26)*100, 2)


######## AVERAGE RELATIVE  BIAS OF tau ###########
rel_bias_tauDPjags.bin_U_U_26 <- c()
for (i in 1:N.sim){
  rel_bias_tauDPjags.bin_U_U_26[i] <- (RhatsDP_U_U_26$tau_DP_U_U_26[i]- tau)
}
avg_rel_bias_tauDPjags.bin_U_U_26 <- round(mean(rel_bias_tauDPjags.bin_U_U_26),2)

prel_bias_tauDPjags.bin_U_U_26 <- round(mean(rel_bias_tauDPjags.bin_U_U_26 / tau)*100,2)

######## AVERAGE ABSOLUTE  BIAS OF tau ###########
abs_bias_tauDPjags.bin_U_U_26 <- c()
for (i in 1:N.sim){
  abs_bias_tauDPjags.bin_U_U_26[i] <- abs(RhatsDP_U_U_26$tau_DP_U_U_26[i]- tau)
}
avg_abs_bias_tauDPjags.bin_U_U_26 <- round(mean(abs_bias_tauDPjags.bin_U_U_26),2)

pabs_bias_tauDPjags.bin_U_U_26 <- round(mean(abs_bias_tauDPjags.bin_U_U_26 / tau)*100,2)

######## MSE OF tau #############
MSE_DPjags.bin.tau_U_U_26 <- c()
for (i in 1:N.sim){
  MSE_DPjags.bin.tau_U_U_26[i] <- RMSE(estimate =RhatsDP_U_U_26$tau_DP_U_U_26[i],
                                       parameter = tau,
                                       type = "RMSE",
                                       MSE = TRUE         
  )
  
}
avg_MSE_DPjags.bin.tau_U_U_26 <- round(mean(MSE_DPjags.bin.tau_U_U_26),4)

nMSE_DPjags.bin.tau_U_U_26 <- round(mean(MSE_DPjags.bin.tau_U_U_26/ tau2),4)

############ COVERAGE OF tau ##################
interval_contains_true_mean <- function(i) { 
  ( tau >= RhatsDP_U_U_26$LB_tau_DP_U_U_26[i]) && ( tau <= RhatsDP_U_U_26$UB_tau_DP_U_U_26[i])
}
interval_contains_true_mean(i)

coverage_tauDPjags.bin_U_U_26 <- c()
for (i in 1:N.sim){
  if (interval_contains_true_mean(i) == TRUE){
    coverage_tauDPjags.bin_U_U_26[i] = 1
  }
  else  {
    coverage_tauDPjags.bin_U_U_26[i] = 0
  }
  print(coverage_tauDPjags.bin_U_U_26[i])
}
coverage_tauDPjags.bin_U_U_26
cover_tauDPjags.bin_U_U_26 <- round(mean(coverage_tauDPjags.bin_U_U_26)*100, 2)


######## AVERAGE BIAS OF base_tau2 ###########
bias_basetau2DPjags.bin_U_U_26 <- c()
for (i in 1:N.sim){
  bias_basetau2DPjags.bin_U_U_26[i] <- abs(RhatsDP_U_U_26$base_tau2_U_U_26[i]- tau2)
}
avg_bias_basetau2DPjags.bin_U_U_26 <- round(mean(bias_basetau2DPjags.bin_U_U_26),2)

######## MSE OF base_tau2 #############
MSE_DPjags.bin.basetau2_U_U_26 <- c()
for (i in 1:N.sim){
  MSE_DPjags.bin.basetau2_U_U_26[i] <- RMSE(estimate =RhatsDP_U_U_26$base_tau2_U_U_26[i],
                                            parameter = tau2,
                                            type = "RMSE",
                                            MSE = TRUE         
  )
  
}
avg_MSE_DPjags.bin.basetau2_U_U_26 <- round(mean(MSE_DPjags.bin.basetau2_U_U_26),4)

############ COVERAGE OF base_tau2 ##################
interval_contains_true_mean <- function(i) {
  ( tau2 >= RhatsDP_U_U_26$LB_base_tau2_U_U_26[i]) && ( tau2 <= RhatsDP_U_U_26$UB_base_tau2_U_U_26[i])
}
interval_contains_true_mean(i)

coverage_basetau2DPjags.bin_U_U_26 <- c()
for (i in 1:N.sim){
  if (interval_contains_true_mean(i) == TRUE){
    coverage_basetau2DPjags.bin_U_U_26[i] = 1
  }
  else  {
    coverage_basetau2DPjags.bin_U_U_26[i] = 0
  }
  print(coverage_basetau2DPjags.bin_U_U_26[i])
}
coverage_basetau2DPjags.bin_U_U_26
cover_basetau2DPjags.bin_U_U_26 <- round(mean(coverage_basetau2DPjags.bin_U_U_26)*100, 2)

######################## ABSOLUTE BIAS OF RELATIVE EFFECTS ##################
abs_bias_relDP_U_U_26 <- list()
mean_abs_bias_relDP_U_U_26 <- list()
for(i in 1:N.sim){
  for(j in 1:nrow(pDP_U_U_26[[i]])){
    abs_bias_relDP_U_U_26[[i]] <- round(abs(listaDP_U_U_26[[i]]$`rel_effDP_U_U_26[[i]]` - pDP_U_U_26[[i]]$true_eff ),2)
    mean_abs_bias_relDP_U_U_26[[i]] <- mean(abs_bias_relDP_U_U_26[[i]])
  }
}
#mean_bias_relDP

################## mean bias of the mean biases of each data set for relative effects  ###################
mean_abs_bias_of_mean_abs_bias_relDP_U_U_26 <- unlist(mean_abs_bias_relDP_U_U_26)
mean_abs_bias_of_mean_abs_bias_relDP_U_U_26  <- round(mean(mean_abs_bias_of_mean_abs_bias_relDP_U_U_26),2)

######################## RELATIVE BIAS OF RELATIVE EFFECTS ##################
rel_bias_relDP_U_U_26 <- list()
mean_rel_bias_relDP_U_U_26 <- list()
for(i in 1:N.sim){
  for(j in 1:nrow(pDP_U_U_26[[i]])){
    rel_bias_relDP_U_U_26[[i]] <- round((listaDP_U_U_26[[i]]$`rel_effDP_U_U_26[[i]]` - pDP_U_U_26[[i]]$true_eff ),2)
    mean_rel_bias_relDP_U_U_26[[i]] <- mean(rel_bias_relDP_U_U_26[[i]])
  }
}
#mean_bias_relDP

################## mean bias of the mean biases of each data set for relative effects  ###################
mean_rel_bias_of_mean_rel_bias_relDP_U_U_26 <- unlist(mean_rel_bias_relDP_U_U_26)
mean_rel_bias_of_mean_rel_bias_relDP_U_U_26  <- round(mean(mean_rel_bias_of_mean_rel_bias_relDP_U_U_26),2)

#####
########################### FOR THE OVERLAPPING ######################
dataDP_U_U_26 <- list()
overlapDP1_U_U_26 <- list()
overlapDP_U_U_26 <- list()
for(i in 1:N.sim){
  dataDP_U_U_26[[i]] <- list(DP = listaDP_U_U_26[[i]]$`rel_effDP_U_U_26[[i]]` , True = pDP_U_U_26[[i]]$true_eff  )
  overlapDP1_U_U_26[[i]] <- overlap(dataDP_U_U_26[[i]], type = "2")
  overlapDP_U_U_26[[i]] <- overlapDP1_U_U_26[[i]]$OV
}
overlapDP_U_U_26 <- unlist(overlapDP_U_U_26)
avg_overlapDP_U_U_26 <- round(mean(overlapDP_U_U_26),2)


listaDP_rel_eff_U_U_26 <- cbind( mean_abs_bias_relDP_U_U_26, mean_abs_bias_of_mean_abs_bias_relDP_U_U_26,
                                 mean_rel_bias_relDP_U_U_26, mean_rel_bias_of_mean_rel_bias_relDP_U_U_26) 

######### SAVE THE OUTPUTS ###########
DP_bin_jags_U_U_26 <- cbind(RhatsDP_U_U_26,
                            abs_biasDPjags.bin_U_U_26, rel_biasDPjags.bin_U_U_26,
                            MSE_DPjags.bin_U_U_26, biasDPjags_basemu.bin_U_U_26, MSE_DPjags_basemu.bin_U_U_26,
                            bias_tau2DPjags.bin_U_U_26,rel_bias_tauDPjags.bin_U_U_26, abs_bias_tauDPjags.bin_U_U_26, 
                            MSE_DPjags.bin.tau2_U_U_26, MSE_DPjags.bin.tau_U_U_26,
                            bias_basetau2DPjags.bin_U_U_26, MSE_DPjags.bin.basetau2_U_U_26, 
                            avg_abs_biasDPjags.bin_U_U_26, avg_rel_biasDPjags.bin_U_U_26, 
                            avg_biasDPjags_basemu.bin_U_U_26, 
                            avg_bias_tau2DPjags.bin_U_U_26,pbias_tau2DPjags.bin_U_U_26,avg_rel_bias_tauDPjags.bin_U_U_26,avg_abs_bias_tauDPjags.bin_U_U_26,
                            prel_bias_tauDPjags.bin_U_U_26,pabs_bias_tauDPjags.bin_U_U_26,
                            avg_bias_basetau2DPjags.bin_U_U_26,
                            rel_bias_tauDPjags.bin_U_U_26,abs_bias_tauDPjags.bin_U_U_26,
                            avg_MSE_DPjags.bin_U_U_26,avg_MSE_DPjags_basemu.bin_U_U_26,avg_MSE_DPjags.bin.tau2_U_U_26, nMSE_DPjags.bin.tau2_U_U_26,
                            avg_MSE_DPjags.bin.tau_U_U_26, nMSE_DPjags.bin.tau_U_U_26,
                            avg_MSE_DPjags.bin.basetau2_U_U_26,
                            coverageDPjags.bin_U_U_26 , coverageDPjags_basemu.bin_U_U_26, coverage_tau2DPjags.bin_U_U_26, coverage_basetau2DPjags.bin_U_U_26, 
                            coverage_tauDPjags.bin_U_U_26, 
                            coverDPjags.bin_U_U_26, coverDPjags_basemu.bin_U_U_26, cover_tau2DPjags.bin_U_U_26, cover_tauDPjags.bin_U_U_26,
                            cover_basetau2DPjags.bin_U_U_26, RhatDP_out_U_U_26)

Overlap_DPmodel_U_U_26 <- cbind(overlapDP_U_U_26, avg_overlapDP_U_U_26)                    



write.csv(DP_bin_jags_U_U_26, "C:\\Users\\Tianqi YU\\Desktop\\SIMULATION SCENARIOS FINAL\\NORMAL DATA\\DP MODEL\\rev_UNIF_UNIF26_ResultsDP_jags_Normal_bin_data_k=14_mu=0_tau2=0.12.csv",row.names=FALSE )  
write.csv(listaDP_U_U_26, "C:\\Users\\Tianqi YU\\Desktop\\SIMULATION SCENARIOS FINAL\\NORMAL DATA\\DP MODEL\\rev_UNIF_UNIF26_Relative_eff_ResultsDP_jags_Normal_bin_data_k=14_mu=0_tau2=0.12.csv",row.names=FALSE )  
write.csv(listaDP_rel_eff_U_U_26, "C:\\Users\\Tianqi YU\\Desktop\\SIMULATION SCENARIOS FINAL\\NORMAL DATA\\DP MODEL\\rev_UNIF_UNIF26_Bias_Relative_eff_ResultsDP_jags_Normal_bin_data_k=14_mu=0_tau2=0.12.csv",row.names=FALSE )  
write.csv(listaDP_weights_U_U_26, "C:\\Users\\Tianqi YU\\Desktop\\SIMULATION SCENARIOS FINAL\\NORMAL DATA\\DP MODEL\\rev_UNIF_UNIF26_Weights_in_clusters_jags_DP_jags_Normal_bin_data_k=14_mu=0_tau2=0.12.csv",row.names=TRUE )  
write.csv(Overlap_DPmodel_U_U_26, "C:\\Users\\Tianqi YU\\Desktop\\SIMULATION SCENARIOS FINAL\\NORMAL DATA\\DP MODEL\\rev_UNIF_UNIF26_Overlap_DPmodel_Normal_bin_data_k=14_mu=0_tau2=0.12.csv",row.names=TRUE )  

####################### DIRICHLET PROCESS MODEL ######################
############# BINOMIAL-DP-n(Unif/Gamma) #############
cat("model{

  for( i in 1: ns ) {

    m[i] ~ dnorm(0,.0001)
    
    logit(p1[i]) <- m[i]
    r1[i] ~ dbin(p1[i],n1[i])
    
    logit(p2[i]) <- m[i] + delta12[i]
    r2[i] ~ dbin(p2[i],n2[i])
    
    delta12[i] <- theta[Z[i]]
    
    Z[i] ~ dcat(p[]) #Z is an integer variable
  }
  # Constructive DP
  #stick-breaking prior
  p[1] <- r[1]
  for (j in 2:(N-1)){
   p[j] <- r[j]*(1-r[j-1])*p[j-1]/r[j-1]
  }
  for (k in 1:(N-1)){
   r[k] ~ dbeta(1,alpha)T(0,0.99)
   }
  #assumption to ensure sum p[] is 1 Ishwaran truncation
  ps <- sum(p[1:(N-1)])
  for(k in N:N){
   p[k]<-1-ps
   }
  # Baseline distribution
  for(k in 1:N){
   theta[k] ~ dnorm(basemu ,basetau1)
  }
  #### priors 

  basemu~dnorm(0,0.0001)
  basetau1 <- 1/basetau_sqr
  basetau_sqr <- basetau*basetau
  basetau ~ dunif(0,10)
  
  # DPP parameter prior
  alpha ~ dgamma(1,1)
  
  # Random effects distribution mean#
 
  for(i in 1:N){
   meancl[i]<-p[i]*theta[i]
   }
  poptrue<-sum(meancl[])  ### the E[X]
  
  # Random effects distribution variance #
  for(i in 1:N){
   mom2[i]<-p[i]*theta[i]*theta[i]  ####E[X2]
   }
  mom2.true<-sum(mom2[])
  var.true<-mom2.true-(poptrue*poptrue) ###E[X2] - E[X]2
  
  # Programming for calculating summary statistics #
  for(i in 1:ns){
   for (j in 1:N){
    SC[i,j] <- equals(j,Z[i])
   }
  }
  # total clusters K#
  for (j in 1:N){
   cl[j] <-step(sum(SC[,j])-1)
   }
  K<-sum(cl[])
  
}",file="DPmodel2.bin_U_G_n.txt")
modfile = 'DPmodel2.bin_U_G_n.txt'
N.sim =1000
mu = 0
tau2 = 0.12
tau = sqrt(tau2)

cl = makeCluster(19)
registerDoParallel(cl)
DPresults1_U_G_n <- list()

DPresults1_U_G_n = foreach(i = 1:N.sim,.packages = c("rjags", "R2jags","tidyverse")) %dopar% {
  set.seed(2106)
  run.model = jags(
    data =list(ns = nrow(p[[i]]),
               r1 = p[[i]]$ci,
               r2 = p[[i]]$ti,
               n1 = p[[i]]$nci,
               n2 = p[[i]]$nti,
               N= 14
               
    ) ,  
    inits = NULL,
    parameters.to.save = c(
      "basemu",    ## tau of the base Normal distribution
      "basetau",   ## mu of the base Normal distribution
      "poptrue",   ##the overall mu
      "var.true", ## the heterogeneity 
      "delta12", #### the random effects
      "K",       ### the total number of clusters 
      "p",     ## the weights of the process 
      "alpha"  ## the concentration parameter
    ), 
    
    n.chains = 2,
    n.iter = 50000,
    n.burnin = 10000,
    DIC = T,
    model.file = modfile
    
  )
  
}
stopCluster(cl)


DPresults_U_G_n <- list()
posteriorDP_U_G_n <- list()
postDP_U_G_n <- list()
for(i in 1:N.sim){
  DPresults_U_G_n[[i]] <- as.data.frame(DPresults1_U_G_n[[i]]$BUGSoutput$summary) 
  posteriorDP_U_G_n[[i]] = as.data.frame(DPresults1_U_G_n[[i]]$BUGSoutput$sims.matrix)
  postDP_U_G_n[[i]] <- posteriorDP_U_G_n[[i]]$poptrue
}

DP_results_U_G_n <- list()
for(i in 1:N.sim){
  DP_results_U_G_n[[i]] <- DPresults_U_G_n[[i]] 
}

############# extraction of the parameters of interest ###########
alpha_U_G_n <- c()
LB_alpha_U_G_n <- c()
UB_alpha_U_G_n <- c()
mu_DP_U_G_n <- c()
LB_mu_DP_U_G_n <- c()
UB_mu_DP_U_G_n <- c()
tau2_DP_U_G_n <- c()
LB_tau2_DP_U_G_n <- c()
UB_tau2_DP_U_G_n <- c()
precDP_mu_U_G_n <- c()
precDP_tau2_U_G_n <- c()
prec_basemu_U_G_n <- c()
prec_basetau2_U_G_n <- c()
prec_alpha_U_G_n <- c()
base_mu_U_G_n <- c()
base_tau_U_G_n <- c()
base_tau2_U_G_n <- c()
LB_base_mu_U_G_n <- c()
UB_base_mu_U_G_n <- c()
LB_base_tau2_U_G_n <- c()
UB_base_tau2_U_G_n <- c()
Rhat_muDP_U_G_n <- c()
Rhat_tau2DP_U_G_n <- c()
Rhat_basemu_U_G_n <- c()
Rhat_basetau_U_G_n <- c()
Rhat_alpha_U_G_n <- c()
Rhat_deltaDP_U_G_n <- list()
median_K_U_G_n <- list()
LB_K_U_G_n <- list()
UB_K_U_G_n <- list()
rel_effDP_U_G_n <- list()
LB_rel_effDP_U_G_n <- list()
UB_rel_effDP_U_G_n <- list()
sd_rel_effDP_U_G_n <- list()
pi_U_G_n <- list()
f11_U_G_n <- list()
m11_U_G_n <- list()
f22_U_G_n <- list()
m22_U_G_n <- list()
fdd_U_G_n <- list()
mdd_U_G_n <- list()
f33_U_G_n <- list()
m33_U_G_n <- list()
f44_U_G_n <- list()
m44_U_G_n <- list()
f55_U_G_n <- list()
m55_U_G_n <- list()
fcl_U_G_n <- list()
mcl_U_G_n <- list()
fp_U_G_n <- list()
mp_U_G_n <- list()
listaDP1_U_G_n <- list()
listaDP_weights_U_G_n <- list()
for(i in 1:N.sim){
  f11_U_G_n[[i]] <- grepl("poptrue", row.names(DP_results_U_G_n[[i]]))
  m11_U_G_n[[i]] <- DP_results_U_G_n[[i]][f11_U_G_n[[i]],]
  mu_DP_U_G_n[i] <- m11_U_G_n[[i]]$`50%`
  LB_mu_DP_U_G_n[i] <- m11_U_G_n[[i]]$`2.5%`
  UB_mu_DP_U_G_n[i] <- m11_U_G_n[[i]]$`97.5%`
  Rhat_muDP_U_G_n[i] <- m11_U_G_n[[i]]$Rhat
  precDP_mu_U_G_n[i] <- UB_mu_DP_U_G_n[i] - LB_mu_DP_U_G_n[i]
  
  f22_U_G_n[[i]] <- grepl("var.true", row.names(DP_results_U_G_n[[i]]))
  m22_U_G_n[[i]] <- DP_results_U_G_n[[i]][f22_U_G_n[[i]],]
  tau2_DP_U_G_n[i] <- m22_U_G_n[[i]]$`50%`
  LB_tau2_DP_U_G_n[i] <- m22_U_G_n[[i]]$`2.5%`
  UB_tau2_DP_U_G_n[i] <- m22_U_G_n[[i]]$`97.5%`
  Rhat_tau2DP_U_G_n[i] <- m22_U_G_n[[i]]$Rhat
  precDP_tau2_U_G_n[i] <- UB_tau2_DP_U_G_n[i] -  LB_tau2_DP_U_G_n[i]
  
  fdd_U_G_n[[i]] <- grepl("delta12", row.names(DP_results_U_G_n[[i]]))
  mdd_U_G_n[[i]] <- DP_results_U_G_n[[i]][fdd_U_G_n[[i]],]
  rel_effDP_U_G_n[[i]] <- mdd_U_G_n[[i]]$`50%`
  LB_rel_effDP_U_G_n[[i]] <- mdd_U_G_n[[i]]$`2.5%`
  UB_rel_effDP_U_G_n[[i]] <- mdd_U_G_n[[i]]$`97.5%`
  sd_rel_effDP_U_G_n[[i]] <- mdd_U_G_n[[i]]$sd
  Rhat_deltaDP_U_G_n[[i]] <- mdd_U_G_n[[i]]$Rhat
  
  f33_U_G_n[[i]] <- grepl("basemu", row.names(DP_results_U_G_n[[i]]))
  m33_U_G_n[[i]] <- DP_results_U_G_n[[i]][f33_U_G_n[[i]],]
  base_mu_U_G_n[i] <- m33_U_G_n[[i]]$`50%`
  LB_base_mu_U_G_n[i] <- m33_U_G_n[[i]]$`2.5%`
  UB_base_mu_U_G_n[i] <- m33_U_G_n[[i]]$`97.5%`
  Rhat_basemu_U_G_n[i] <- m33_U_G_n[[i]]$Rhat
  prec_basemu_U_G_n[i] <- UB_base_mu_U_G_n[i] - LB_base_mu_U_G_n[i]
  
  f44_U_G_n[[i]] <- grepl("basetau", row.names(DP_results_U_G_n[[i]]))
  m44_U_G_n[[i]] <- DP_results_U_G_n[[i]][f44_U_G_n[[i]],]
  base_tau_U_G_n[i] <- m44_U_G_n[[i]]$`50%`
  base_tau2_U_G_n[i] <- (m44_U_G_n[[i]]$`50%`)^2
  LB_base_tau2_U_G_n[i] <- (m44_U_G_n[[i]]$`2.5%`)^2
  UB_base_tau2_U_G_n[i] <- (m44_U_G_n[[i]]$`97.5%`)^2
  Rhat_basetau_U_G_n[i] <- (m44_U_G_n[[i]]$Rhat)^2
  prec_basetau2_U_G_n[i] <- UB_base_tau2_U_G_n[i] - LB_base_tau2_U_G_n[i]
  
  f55_U_G_n[[i]] <- grepl("alpha", row.names(DP_results_U_G_n[[i]]))
  m55_U_G_n[[i]] <- DP_results_U_G_n[[i]][f55_U_G_n[[i]],]
  alpha_U_G_n[i] <- m55_U_G_n[[i]]$`50%`
  LB_alpha_U_G_n[i] <- m55_U_G_n[[i]]$`2.5%`
  UB_alpha_U_G_n[i] <- m55_U_G_n[[i]]$`97.5%`
  Rhat_alpha_U_G_n[i] <- m55_U_G_n[[i]]$Rhat
  prec_alpha_U_G_n[i] <- UB_alpha_U_G_n[i] - LB_alpha_U_G_n[i]
  
  fcl_U_G_n[[i]] <- grepl("K", row.names(DP_results_U_G_n[[i]]))  
  mcl_U_G_n[[i]] <- DP_results_U_G_n[[i]][fcl_U_G_n[[i]],]
  median_K_U_G_n[[i]] <- mcl_U_G_n[[i]]$`50%` 
  LB_K_U_G_n[[i]] <- mcl_U_G_n[[i]]$`2.5%`
  UB_K_U_G_n[[i]] <- mcl_U_G_n[[i]]$`97.5%`
  
  fp_U_G_n[[i]] <- grepl("p", row.names(DP_results_U_G_n[[i]]))
  mp_U_G_n[[i]] <- DP_results_U_G_n[[i]][fp_U_G_n[[i]],]
  mp_U_G_n[[i]] <- mp_U_G_n[[i]][!grepl("pop", row.names(mp_U_G_n[[i]])),]
  mp_U_G_n[[i]] <- mp_U_G_n[[i]][!grepl("alpha", row.names(mp_U_G_n[[i]])),]
  pi_U_G_n[[i]] <- mp_U_G_n[[i]]$mean
  
  listaDP1_U_G_n[[i]] <- cbind.data.frame(rel_effDP_U_G_n[[i]],sd_rel_effDP_U_G_n[[i]], LB_rel_effDP_U_G_n[[i]],UB_rel_effDP_U_G_n[[i]],Rhat_deltaDP_U_G_n[[i]])
  listaDP_weights_U_G_n[[i]] <- cbind.data.frame(pi_U_G_n[[i]])
  
}

numclus_U_G_n <- unlist(median_K_U_G_n)
LB_K_U_G_n <- unlist(LB_K_U_G_n)
UB_K_U_G_n <- unlist(UB_K_U_G_n)

mean_alpha_U_G_n = mean(alpha_U_G_n)
mean_alpha_LB_U_G_n = mean(LB_alpha_U_G_n)
mean_alpha_UB_U_G_n = mean(UB_alpha_U_G_n)

tau_DP_U_G_n = sqrt(tau2_DP_U_G_n)
LB_tau_DP_U_G_n = sqrt(LB_tau2_DP_U_G_n)
UB_tau_DP_U_G_n = sqrt(UB_tau2_DP_U_G_n)
precDP_tau_U_G_n = UB_tau_DP_U_G_n -  LB_tau_DP_U_G_n


RhatsDP_U_G_n <- cbind.data.frame(base_mu_U_G_n,LB_base_mu_U_G_n,UB_base_mu_U_G_n,base_tau_U_G_n,base_tau2_U_G_n, LB_base_tau2_U_G_n,UB_base_tau2_U_G_n,
                                  mu_DP_U_G_n, LB_mu_DP_U_G_n, UB_mu_DP_U_G_n, tau2_DP_U_G_n, LB_tau2_DP_U_G_n, UB_tau2_DP_U_G_n,
                                  tau_DP_U_G_n, LB_tau_DP_U_G_n, UB_tau_DP_U_G_n,
                                  precDP_mu_U_G_n, precDP_tau2_U_G_n, prec_basemu_U_G_n, prec_basetau2_U_G_n,prec_alpha_U_G_n,  alpha_U_G_n , LB_alpha_U_G_n, UB_alpha_U_G_n,
                                  mean_alpha_U_G_n, mean_alpha_LB_U_G_n, mean_alpha_UB_U_G_n,
                                  Rhat_muDP_U_G_n, Rhat_tau2DP_U_G_n, Rhat_basemu_U_G_n, Rhat_basetau_U_G_n, Rhat_alpha_U_G_n, numclus_U_G_n, LB_K_U_G_n, UB_K_U_G_n)

##########REMOVE Rhats > 1.05 ##############

condition1_U_G_n <- which(RhatsDP_U_G_n$Rhat_muDP_U_G_n > 1.05)
condition2_U_G_n <- which(RhatsDP_U_G_n$Rhat_tau2DP_U_G_n > 1.05)
condition3_U_G_n <- which(RhatsDP_U_G_n$Rhat_basemu_U_G_n > 1.05)
condition4_U_G_n <- which(RhatsDP_U_G_n$Rhat_basetau_U_G_n > 1.05)
condition5_U_G_n <- which(RhatsDP_U_G_n$Rhat_alpha_U_G_n > 1.05)


dist.condDP_U_G_n = c(condition1_U_G_n,condition2_U_G_n,condition3_U_G_n,condition4_U_G_n,condition5_U_G_n)
dist.condDP_U_G_n = unique(dist.condDP_U_G_n)

RhatDP_out_U_G_n = round((length(dist.condDP_U_G_n)/N.sim), 4)

############### Extract and remove the datasets with Rhat > 1.05 #########
if (length(dist.condDP_U_G_n)== 0) {
  RhatsDP_U_G_n <- RhatsDP_U_G_n
  listaDP_U_G_n <- listaDP1_U_G_n
  listaDP_weights_U_G_n <- listaDP_weights_U_G_n
  N.sim <- nrow(RhatsDP_U_G_n)
  pDP_U_G_n <- p
  
} else {
  RhatsDP_U_G_n <- RhatsDP_U_G_n[-dist.condDP_U_G_n, ]
  listaDP_U_G_n <- listaDP1_U_G_n[-dist.condDP_U_G_n]
  listaDP_weights_U_G_n <- listaDP_weights_U_G_n[-dist.condDP_U_G_n]
  N.sim <- nrow(RhatsDP_U_G_n)
  pDP_U_G_n <- p[-dist.condDP_U_G_n]
}

count_DP_U_G_n = list()
for(i in length(listaDP_U_G_n)){
  count_DP_U_G_n[[i]] = which(listaDP_U_G_n[[i]]$`Rhat_deltaDP_U_G_n[[i]]` > 1.05 )
  tell_meDP_U_G_n = which(count_DP_U_G_n[[i]] != 0)
}

tell_meDP_U_G_n

if(length(tell_meDP_U_G_n) == 0){
  RhatsDP_U_G_n <- RhatsDP_U_G_n
  N.sim <- nrow(RhatsDP_U_G_n)
  listaDP_U_G_n <- listaDP_U_G_n
  listaDP_weights_U_G_n <- listaDP_weights_U_G_n
  pDP_U_G_n <- pDP_U_G_n
  RhatDP_out_U_G_n <- RhatDP_out_U_G_n
} else {
  RhatsDP_U_G_n <- RhatsDP_U_G_n[-tell_meDP_U_G_n, ]
  listaDP_U_G_n <- listaDP_U_G_n[-tell_meDP_U_G_n]
  listaDP_weights_U_G_n <- listaDP_weights_U_G_n[-tell_meDP_U_G_n]
  N.sim = nrow(RhatsDP_U_G_n)
  pDP_U_G_n <- pDP_U_G_n[-tell_meDP_U_G_n]
  RhatDP_out_U_G_n <- RhatDP_out_U_G_n + ((length(tell_meDP_U_G_n))/N.sim)
}

RhatDP_out_U_G_n

#### AVERAGE ABSOLUTE BIAS OF mu ######
abs_biasDPjags.bin_U_G_n <- c()
for (i in 1:N.sim){
  abs_biasDPjags.bin_U_G_n[i] <- abs( RhatsDP_U_G_n$mu_DP_U_G_n[i] - mu) 
}
avg_abs_biasDPjags.bin_U_G_n <- round(mean(abs_biasDPjags.bin_U_G_n),2)

#### AVERAGE RELATIVE BIAS OF mu ######
rel_biasDPjags.bin_U_G_n <- c()
for (i in 1:N.sim){
  rel_biasDPjags.bin_U_G_n[i] <- ( RhatsDP_U_G_n$mu_DP_U_G_n[i] - mu) 
}
avg_rel_biasDPjags.bin_U_G_n <- round(mean(rel_biasDPjags.bin_U_G_n),2)

######  MSE OF mu ##########
MSE_DPjags.bin_U_G_n <- c()
for (i in 1:N.sim){
  MSE_DPjags.bin_U_G_n[i] <- RMSE(estimate = RhatsDP_U_G_n$mu_DP_U_G_n[i],
                                  parameter = mu,
                                  type = "RMSE",
                                  MSE = TRUE   )       
}
avg_MSE_DPjags.bin_U_G_n <- round(mean(MSE_DPjags.bin_U_G_n), 4)

############ COVERAGE OF mu ########
interval_contains_true_mean <- function(i) { 
  ( mu >= RhatsDP_U_G_n$LB_mu_DP_U_G_n[i]) && ( mu <= RhatsDP_U_G_n$UB_mu_DP_U_G_n[i])
}
interval_contains_true_mean(i)

coverageDPjags.bin_U_G_n <- c()
for (i in 1:N.sim){
  if (interval_contains_true_mean(i) == TRUE){
    coverageDPjags.bin_U_G_n[i] = 1
  }
  else {
    coverageDPjags.bin_U_G_n[i]=0
  }
  print(coverageDPjags.bin_U_G_n[i])
}
coverageDPjags.bin_U_G_n
coverDPjags.bin_U_G_n <- round(mean(coverageDPjags.bin_U_G_n)*100, 2)

## BIAS OF base_mu ######
biasDPjags_basemu.bin_U_G_n <- c()
for (i in 1:N.sim){
  biasDPjags_basemu.bin_U_G_n[i] <- abs( RhatsDP_U_G_n$base_mu_U_G_n[i] - mu) 
}
avg_biasDPjags_basemu.bin_U_G_n <- round(mean(biasDPjags_basemu.bin_U_G_n),2)

###### MSE OF base_mu ##########
MSE_DPjags_basemu.bin_U_G_n <- c()
for (i in 1:N.sim){
  MSE_DPjags_basemu.bin_U_G_n[i] <- RMSE(estimate =RhatsDP_U_G_n$base_mu_U_G_n[i],
                                         parameter = mu,
                                         type = "RMSE",
                                         MSE = TRUE      
  )
}
avg_MSE_DPjags_basemu.bin_U_G_n <- round(mean(MSE_DPjags_basemu.bin_U_G_n), 4)

############ COVERAGE OF base_mu ########
interval_contains_true_mean <- function(i) { 
  ( mu >= RhatsDP_U_G_n$LB_base_mu_U_G_n[i]) && ( mu <= RhatsDP_U_G_n$UB_base_mu_U_G_n[i])
}
interval_contains_true_mean(i)

coverageDPjags_basemu.bin_U_G_n <- c()
for (i in 1:N.sim){
  if (interval_contains_true_mean(i) == TRUE){
    coverageDPjags_basemu.bin_U_G_n[i] = 1
  }
  else {
    coverageDPjags_basemu.bin_U_G_n[i]=0
  }
  print(coverageDPjags_basemu.bin_U_G_n[i])
}
coverageDPjags_basemu.bin_U_G_n
coverDPjags_basemu.bin_U_G_n <- round(mean(coverageDPjags_basemu.bin_U_G_n)*100, 2)

######## AVERAGE BIAS OF tau2 ###########
bias_tau2DPjags.bin_U_G_n <- c()
for (i in 1:N.sim){
  bias_tau2DPjags.bin_U_G_n[i] <- abs(RhatsDP_U_G_n$tau2_DP_U_G_n[i]- tau2)
}
avg_bias_tau2DPjags.bin_U_G_n <- round(mean(bias_tau2DPjags.bin_U_G_n),2)

pbias_tau2DPjags.bin_U_G_n <- round(mean(bias_tau2DPjags.bin_U_G_n / tau2),2)

######## MSE OF tau2 #############
MSE_DPjags.bin.tau2_U_G_n <- c()
for (i in 1:N.sim){
  MSE_DPjags.bin.tau2_U_G_n[i] <- RMSE(estimate =RhatsDP_U_G_n$tau2_DP_U_G_n[i],
                                       parameter = tau2,
                                       type = "RMSE",
                                       MSE = TRUE         
  )
  
}
avg_MSE_DPjags.bin.tau2_U_G_n <- round(mean(MSE_DPjags.bin.tau2_U_G_n),4)

nMSE_DPjags.bin.tau2_U_G_n <- round(mean(MSE_DPjags.bin.tau2_U_G_n/ (tau2)^2),4)

############ COVERAGE OF tau2 ##################
interval_contains_true_mean <- function(i) { 
  ( tau2 >= RhatsDP_U_G_n$LB_tau2_DP_U_G_n[i]) && ( tau2 <= RhatsDP_U_G_n$UB_tau2_DP_U_G_n[i])
}
interval_contains_true_mean(i)

coverage_tau2DPjags.bin_U_G_n <- c()
for (i in 1:N.sim){
  if (interval_contains_true_mean(i) == TRUE){
    coverage_tau2DPjags.bin_U_G_n[i] = 1
  }
  else  {
    coverage_tau2DPjags.bin_U_G_n[i] = 0
  }
  print(coverage_tau2DPjags.bin_U_G_n[i])
}
coverage_tau2DPjags.bin_U_G_n
cover_tau2DPjags.bin_U_G_n <- round(mean(coverage_tau2DPjags.bin_U_G_n)*100, 2)


######## AVERAGE RELATIVE  BIAS OF tau ###########
rel_bias_tauDPjags.bin_U_G_n <- c()
for (i in 1:N.sim){
  rel_bias_tauDPjags.bin_U_G_n[i] <- (RhatsDP_U_G_n$tau_DP_U_G_n[i]- tau)
}
avg_rel_bias_tauDPjags.bin_U_G_n <- round(mean(rel_bias_tauDPjags.bin_U_G_n),2)

prel_bias_tauDPjags.bin_U_G_n <- round(mean(rel_bias_tauDPjags.bin_U_G_n / tau)*100,2)

######## AVERAGE ABSOLUTE  BIAS OF tau ###########
abs_bias_tauDPjags.bin_U_G_n <- c()
for (i in 1:N.sim){
  abs_bias_tauDPjags.bin_U_G_n[i] <- abs(RhatsDP_U_G_n$tau_DP_U_G_n[i]- tau)
}
avg_abs_bias_tauDPjags.bin_U_G_n <- round(mean(abs_bias_tauDPjags.bin_U_G_n),2)

pabs_bias_tauDPjags.bin_U_G_n <- round(mean(abs_bias_tauDPjags.bin_U_G_n / tau)*100,2)

######## MSE OF tau #############
MSE_DPjags.bin.tau_U_G_n <- c()
for (i in 1:N.sim){
  MSE_DPjags.bin.tau_U_G_n[i] <- RMSE(estimate =RhatsDP_U_G_n$tau_DP_U_G_n[i],
                                      parameter = tau,
                                      type = "RMSE",
                                      MSE = TRUE         
  )
  
}
avg_MSE_DPjags.bin.tau_U_G_n <- round(mean(MSE_DPjags.bin.tau_U_G_n),4)

nMSE_DPjags.bin.tau_U_G_n <- round(mean(MSE_DPjags.bin.tau_U_G_n/ tau2),4)

############ COVERAGE OF tau ##################
interval_contains_true_mean <- function(i) { 
  ( tau >= RhatsDP_U_G_n$LB_tau_DP_U_G_n[i]) && ( tau <= RhatsDP_U_G_n$UB_tau_DP_U_G_n[i])
}
interval_contains_true_mean(i)

coverage_tauDPjags.bin_U_G_n <- c()
for (i in 1:N.sim){
  if (interval_contains_true_mean(i) == TRUE){
    coverage_tauDPjags.bin_U_G_n[i] = 1
  }
  else  {
    coverage_tauDPjags.bin_U_G_n[i] = 0
  }
  print(coverage_tauDPjags.bin_U_G_n[i])
}
coverage_tauDPjags.bin_U_G_n
cover_tauDPjags.bin_U_G_n <- round(mean(coverage_tauDPjags.bin_U_G_n)*100, 2)


######## AVERAGE BIAS OF base_tau2 ###########
bias_basetau2DPjags.bin_U_G_n <- c()
for (i in 1:N.sim){
  bias_basetau2DPjags.bin_U_G_n[i] <- abs(RhatsDP_U_G_n$base_tau2_U_G_n[i]- tau2)
}
avg_bias_basetau2DPjags.bin_U_G_n <- round(mean(bias_basetau2DPjags.bin_U_G_n),2)

######## MSE OF base_tau2 #############
MSE_DPjags.bin.basetau2_U_G_n <- c()
for (i in 1:N.sim){
  MSE_DPjags.bin.basetau2_U_G_n[i] <- RMSE(estimate =RhatsDP_U_G_n$base_tau2_U_G_n[i],
                                           parameter = tau2,
                                           type = "RMSE",
                                           MSE = TRUE         
  )
  
}
avg_MSE_DPjags.bin.basetau2_U_G_n <- round(mean(MSE_DPjags.bin.basetau2_U_G_n),4)

############ COVERAGE OF base_tau2 ##################
interval_contains_true_mean <- function(i) {
  ( tau2 >= RhatsDP_U_G_n$LB_base_tau2_U_G_n[i]) && ( tau2 <= RhatsDP_U_G_n$UB_base_tau2_U_G_n[i])
}
interval_contains_true_mean(i)

coverage_basetau2DPjags.bin_U_G_n <- c()
for (i in 1:N.sim){
  if (interval_contains_true_mean(i) == TRUE){
    coverage_basetau2DPjags.bin_U_G_n[i] = 1
  }
  else  {
    coverage_basetau2DPjags.bin_U_G_n[i] = 0
  }
  print(coverage_basetau2DPjags.bin_U_G_n[i])
}
coverage_basetau2DPjags.bin_U_G_n
cover_basetau2DPjags.bin_U_G_n <- round(mean(coverage_basetau2DPjags.bin_U_G_n)*100, 2)

######################## ABSOLUTE BIAS OF RELATIVE EFFECTS ##################
abs_bias_relDP_U_G_n <- list()
mean_abs_bias_relDP_U_G_n <- list()
for(i in 1:N.sim){
  for(j in 1:nrow(pDP_U_G_n[[i]])){
    abs_bias_relDP_U_G_n[[i]] <- round(abs(listaDP_U_G_n[[i]]$`rel_effDP_U_G_n[[i]]` - pDP_U_G_n[[i]]$true_eff ),2)
    mean_abs_bias_relDP_U_G_n[[i]] <- mean(abs_bias_relDP_U_G_n[[i]])
  }
}
#mean_bias_relDP

################## mean bias of the mean biases of each data set for relative effects  ###################
mean_abs_bias_of_mean_abs_bias_relDP_U_G_n <- unlist(mean_abs_bias_relDP_U_G_n)
mean_abs_bias_of_mean_abs_bias_relDP_U_G_n  <- round(mean(mean_abs_bias_of_mean_abs_bias_relDP_U_G_n),2)

######################## RELATIVE BIAS OF RELATIVE EFFECTS ##################
rel_bias_relDP_U_G_n <- list()
mean_rel_bias_relDP_U_G_n <- list()
for(i in 1:N.sim){
  for(j in 1:nrow(pDP_U_G_n[[i]])){
    rel_bias_relDP_U_G_n[[i]] <- round((listaDP_U_G_n[[i]]$`rel_effDP_U_G_n[[i]]` - pDP_U_G_n[[i]]$true_eff ),2)
    mean_rel_bias_relDP_U_G_n[[i]] <- mean(rel_bias_relDP_U_G_n[[i]])
  }
}
#mean_bias_relDP

################## mean bias of the mean biases of each data set for relative effects  ###################
mean_rel_bias_of_mean_rel_bias_relDP_U_G_n <- unlist(mean_rel_bias_relDP_U_G_n)
mean_rel_bias_of_mean_rel_bias_relDP_U_G_n  <- round(mean(mean_rel_bias_of_mean_rel_bias_relDP_U_G_n),2)

#####
########################### FOR THE OVERLAPPING ######################
dataDP_U_G_n <- list()
overlapDP1_U_G_n <- list()
overlapDP_U_G_n <- list()
for(i in 1:N.sim){
  dataDP_U_G_n[[i]] <- list(DP = listaDP_U_G_n[[i]]$`rel_effDP_U_G_n[[i]]` , True = pDP_U_G_n[[i]]$true_eff  )
  overlapDP1_U_G_n[[i]] <- overlap(dataDP_U_G_n[[i]], type = "2")
  overlapDP_U_G_n[[i]] <- overlapDP1_U_G_n[[i]]$OV
}
overlapDP_U_G_n <- unlist(overlapDP_U_G_n)
avg_overlapDP_U_G_n <- round(mean(overlapDP_U_G_n),2)


listaDP_rel_eff_U_G_n <- cbind( mean_abs_bias_relDP_U_G_n, mean_abs_bias_of_mean_abs_bias_relDP_U_G_n,
                                mean_rel_bias_relDP_U_G_n, mean_rel_bias_of_mean_rel_bias_relDP_U_G_n) 

######### SAVE THE OUTPUTS ###########
DP_bin_jags_U_G_n <- cbind(RhatsDP_U_G_n,
                           abs_biasDPjags.bin_U_G_n, rel_biasDPjags.bin_U_G_n,
                           MSE_DPjags.bin_U_G_n, biasDPjags_basemu.bin_U_G_n, MSE_DPjags_basemu.bin_U_G_n,
                           bias_tau2DPjags.bin_U_G_n,rel_bias_tauDPjags.bin_U_G_n, abs_bias_tauDPjags.bin_U_G_n, 
                           MSE_DPjags.bin.tau2_U_G_n, MSE_DPjags.bin.tau_U_G_n,
                           bias_basetau2DPjags.bin_U_G_n, MSE_DPjags.bin.basetau2_U_G_n, 
                           avg_abs_biasDPjags.bin_U_G_n, avg_rel_biasDPjags.bin_U_G_n, 
                           avg_biasDPjags_basemu.bin_U_G_n, 
                           avg_bias_tau2DPjags.bin_U_G_n,pbias_tau2DPjags.bin_U_G_n,avg_rel_bias_tauDPjags.bin_U_G_n,avg_abs_bias_tauDPjags.bin_U_G_n,
                           prel_bias_tauDPjags.bin_U_G_n,pabs_bias_tauDPjags.bin_U_G_n,
                           avg_bias_basetau2DPjags.bin_U_G_n,
                           rel_bias_tauDPjags.bin_U_G_n,abs_bias_tauDPjags.bin_U_G_n,
                           avg_MSE_DPjags.bin_U_G_n,avg_MSE_DPjags_basemu.bin_U_G_n,avg_MSE_DPjags.bin.tau2_U_G_n, nMSE_DPjags.bin.tau2_U_G_n,
                           avg_MSE_DPjags.bin.tau_U_G_n, nMSE_DPjags.bin.tau_U_G_n,
                           avg_MSE_DPjags.bin.basetau2_U_G_n,
                           coverageDPjags.bin_U_G_n , coverageDPjags_basemu.bin_U_G_n, coverage_tau2DPjags.bin_U_G_n, coverage_basetau2DPjags.bin_U_G_n, 
                           coverage_tauDPjags.bin_U_G_n, 
                           coverDPjags.bin_U_G_n, coverDPjags_basemu.bin_U_G_n, cover_tau2DPjags.bin_U_G_n, cover_tauDPjags.bin_U_G_n,
                           cover_basetau2DPjags.bin_U_G_n, RhatDP_out_U_G_n)

Overlap_DPmodel_U_G_n <- cbind(overlapDP_U_G_n, avg_overlapDP_U_G_n)                    


write.csv(DP_bin_jags_U_G_n, "C:\\Users\\Tianqi YU\\Desktop\\SIMULATION SCENARIOS FINAL\\NORMAL DATA\\DP MODEL\\rev_UNIF_gamma_ns_ResultsDP_jags_Normal_bin_data_k=14_mu=0_tau2=0.12.csv",row.names=FALSE )  
write.csv(listaDP_U_G_n, "C:\\Users\\Tianqi YU\\Desktop\\SIMULATION SCENARIOS FINAL\\NORMAL DATA\\DP MODEL\\rev_UNIF_gamma_ns_Relative_eff_ResultsDP_jags_Normal_bin_data_k=14_mu=0_tau2=0.12.csv",row.names=FALSE )  
write.csv(listaDP_rel_eff_U_G_n, "C:\\Users\\Tianqi YU\\Desktop\\SIMULATION SCENARIOS FINAL\\NORMAL DATA\\DP MODEL\\rev_UNIF_gamma_ns_Bias_Relative_eff_ResultsDP_jags_Normal_bin_data_k=14_mu=0_tau2=0.12.csv",row.names=FALSE )  
write.csv(listaDP_weights_U_G_n, "C:\\Users\\Tianqi YU\\Desktop\\SIMULATION SCENARIOS FINAL\\NORMAL DATA\\DP MODEL\\rev_UNIF_gamma_ns_Weights_in_clusters_jags_DP_jags_Normal_bin_data_k=14_mu=0_tau2=0.12.csv",row.names=TRUE )  
write.csv(Overlap_DPmodel_U_G_n, "C:\\Users\\Tianqi YU\\Desktop\\SIMULATION SCENARIOS FINAL\\NORMAL DATA\\DP MODEL\\rev_UNIF_gamma_ns_Overlap_DPmodel_Normal_bin_data_k=14_mu=0_tau2=0.12.csv",row.names=TRUE )  

####################### DIRICHLET PROCESS MODEL ######################
############# BINOMIAL-DP-51(Unif/Unif) #############
cat("model{

  for( i in 1: ns ) {

    m[i] ~ dnorm(0,.0001)
    
    logit(p1[i]) <- m[i]
    r1[i] ~ dbin(p1[i],n1[i])
    
    logit(p2[i]) <- m[i] + delta12[i]
    r2[i] ~ dbin(p2[i],n2[i])
    
    delta12[i] <- theta[Z[i]]
    
    Z[i] ~ dcat(p[]) #Z is an integer variable
  }
  # Constructive DP
  #stick-breaking prior
  p[1] <- r[1]
  for (j in 2:(N-1)){
   p[j] <- r[j]*(1-r[j-1])*p[j-1]/r[j-1]
  }
  for (k in 1:(N-1)){
   r[k] ~ dbeta(1,alpha)T(0,0.99)
   }
  #assumption to ensure sum p[] is 1 Ishwaran truncation
  ps <- sum(p[1:(N-1)])
  for(k in N:N){
   p[k]<-1-ps
   }
  # Baseline distribution
  for(k in 1:N){
   theta[k] ~ dnorm(basemu ,basetau1)
  }
  #### priors 

  basemu~dnorm(0,0.0001)
  basetau1 <- 1/basetau_sqr
  basetau_sqr <- basetau*basetau
  basetau ~ dunif(0,10)
  
  # DPP parameter prior
  alpha~dunif(0.3,10)
  
  # Random effects distribution mean#
 
  for(i in 1:N){
   meancl[i]<-p[i]*theta[i]
   }
  poptrue<-sum(meancl[])  ### the E[X]
  
  # Random effects distribution variance #
  for(i in 1:N){
   mom2[i]<-p[i]*theta[i]*theta[i]  ####E[X2]
   }
  mom2.true<-sum(mom2[])
  var.true<-mom2.true-(poptrue*poptrue) ###E[X2] - E[X]2
  
  # Programming for calculating summary statistics #
  for(i in 1:ns){
   for (j in 1:N){
    SC[i,j] <- equals(j,Z[i])
   }
  }
  # total clusters K#
  for (j in 1:N){
   cl[j] <-step(sum(SC[,j])-1)
   }
  K<-sum(cl[])
  
}",file="DPmodel2.bin_U_U_51.txt")
modfile = 'DPmodel2.bin_U_U_51.txt'
N.sim =1000
mu = 0
tau2 = 0.12
tau = sqrt(tau2)

cl = makeCluster(19)
registerDoParallel(cl)
DPresults1_U_U_51 <- list()

DPresults1_U_U_51 = foreach(i = 1:N.sim,.packages = c("rjags", "R2jags","tidyverse")) %dopar% {
  set.seed(21012)
  run.model = jags(
    data =list(ns = nrow(p[[i]]),
               r1 = p[[i]]$ci,
               r2 = p[[i]]$ti,
               n1 = p[[i]]$nci,
               n2 = p[[i]]$nti,
               N= 51
               
    ) ,  
    inits = NULL,
    parameters.to.save = c(
      "basemu",    ## tau of the base Normal distribution
      "basetau",   ## mu of the base Normal distribution
      "poptrue",   ##the overall mu
      "var.true", ## the heterogeneity 
      "delta12", #### the random effects
      "K",       ### the total number of clusters 
      "p",     ## the weights of the process 
      "alpha"  ## the concentration parameter
    ), 
    
    n.chains = 2,
    n.iter = 50000,
    n.burnin = 10000,
    DIC = T,
    model.file = modfile
    
  )
  
}
stopCluster(cl)


DPresults_U_U_51 <- list()
posteriorDP_U_U_51 <- list()
postDP_U_U_51 <- list()
for(i in 1:N.sim){
  DPresults_U_U_51[[i]] <- as.data.frame(DPresults1_U_U_51[[i]]$BUGSoutput$summary) 
  posteriorDP_U_U_51[[i]] = as.data.frame(DPresults1_U_U_51[[i]]$BUGSoutput$sims.matrix)
  postDP_U_U_51[[i]] <- posteriorDP_U_U_51[[i]]$poptrue
}

DP_results_U_U_51 <- list()
for(i in 1:N.sim){
  DP_results_U_U_51[[i]] <- DPresults_U_U_51[[i]] 
}

############# extraction of the parameters of interest ###########
alpha_U_U_51 <- c()
LB_alpha_U_U_51 <- c()
UB_alpha_U_U_51 <- c()
mu_DP_U_U_51 <- c()
LB_mu_DP_U_U_51 <- c()
UB_mu_DP_U_U_51 <- c()
tau2_DP_U_U_51 <- c()
LB_tau2_DP_U_U_51 <- c()
UB_tau2_DP_U_U_51 <- c()
precDP_mu_U_U_51 <- c()
precDP_tau2_U_U_51 <- c()
prec_basemu_U_U_51 <- c()
prec_basetau2_U_U_51 <- c()
prec_alpha_U_U_51 <- c()
base_mu_U_U_51 <- c()
base_tau_U_U_51 <- c()
base_tau2_U_U_51 <- c()
LB_base_mu_U_U_51 <- c()
UB_base_mu_U_U_51 <- c()
LB_base_tau2_U_U_51 <- c()
UB_base_tau2_U_U_51 <- c()
Rhat_muDP_U_U_51 <- c()
Rhat_tau2DP_U_U_51 <- c()
Rhat_basemu_U_U_51 <- c()
Rhat_basetau_U_U_51 <- c()
Rhat_alpha_U_U_51 <- c()
Rhat_deltaDP_U_U_51 <- list()
median_K_U_U_51 <- list()
LB_K_U_U_51 <- list()
UB_K_U_U_51 <- list()
rel_effDP_U_U_51 <- list()
LB_rel_effDP_U_U_51 <- list()
UB_rel_effDP_U_U_51 <- list()
sd_rel_effDP_U_U_51 <- list()
pi_U_U_51 <- list()
f11_U_U_51 <- list()
m11_U_U_51 <- list()
f22_U_U_51 <- list()
m22_U_U_51 <- list()
fdd_U_U_51 <- list()
mdd_U_U_51 <- list()
f33_U_U_51 <- list()
m33_U_U_51 <- list()
f44_U_U_51 <- list()
m44_U_U_51 <- list()
f55_U_U_51 <- list()
m55_U_U_51 <- list()
fcl_U_U_51 <- list()
mcl_U_U_51 <- list()
fp_U_U_51 <- list()
mp_U_U_51 <- list()
listaDP1_U_U_51 <- list()
listaDP_weights_U_U_51 <- list()
for(i in 1:N.sim){
  f11_U_U_51[[i]] <- grepl("poptrue", row.names(DP_results_U_U_51[[i]]))
  m11_U_U_51[[i]] <- DP_results_U_U_51[[i]][f11_U_U_51[[i]],]
  mu_DP_U_U_51[i] <- m11_U_U_51[[i]]$`50%`
  LB_mu_DP_U_U_51[i] <- m11_U_U_51[[i]]$`2.5%`
  UB_mu_DP_U_U_51[i] <- m11_U_U_51[[i]]$`97.5%`
  Rhat_muDP_U_U_51[i] <- m11_U_U_51[[i]]$Rhat
  precDP_mu_U_U_51[i] <- UB_mu_DP_U_U_51[i] - LB_mu_DP_U_U_51[i]
  
  f22_U_U_51[[i]] <- grepl("var.true", row.names(DP_results_U_U_51[[i]]))
  m22_U_U_51[[i]] <- DP_results_U_U_51[[i]][f22_U_U_51[[i]],]
  tau2_DP_U_U_51[i] <- m22_U_U_51[[i]]$`50%`
  LB_tau2_DP_U_U_51[i] <- m22_U_U_51[[i]]$`2.5%`
  UB_tau2_DP_U_U_51[i] <- m22_U_U_51[[i]]$`97.5%`
  Rhat_tau2DP_U_U_51[i] <- m22_U_U_51[[i]]$Rhat
  precDP_tau2_U_U_51[i] <- UB_tau2_DP_U_U_51[i] -  LB_tau2_DP_U_U_51[i]
  
  fdd_U_U_51[[i]] <- grepl("delta12", row.names(DP_results_U_U_51[[i]]))
  mdd_U_U_51[[i]] <- DP_results_U_U_51[[i]][fdd_U_U_51[[i]],]
  rel_effDP_U_U_51[[i]] <- mdd_U_U_51[[i]]$`50%`
  LB_rel_effDP_U_U_51[[i]] <- mdd_U_U_51[[i]]$`2.5%`
  UB_rel_effDP_U_U_51[[i]] <- mdd_U_U_51[[i]]$`97.5%`
  sd_rel_effDP_U_U_51[[i]] <- mdd_U_U_51[[i]]$sd
  Rhat_deltaDP_U_U_51[[i]] <- mdd_U_U_51[[i]]$Rhat
  
  f33_U_U_51[[i]] <- grepl("basemu", row.names(DP_results_U_U_51[[i]]))
  m33_U_U_51[[i]] <- DP_results_U_U_51[[i]][f33_U_U_51[[i]],]
  base_mu_U_U_51[i] <- m33_U_U_51[[i]]$`50%`
  LB_base_mu_U_U_51[i] <- m33_U_U_51[[i]]$`2.5%`
  UB_base_mu_U_U_51[i] <- m33_U_U_51[[i]]$`97.5%`
  Rhat_basemu_U_U_51[i] <- m33_U_U_51[[i]]$Rhat
  prec_basemu_U_U_51[i] <- UB_base_mu_U_U_51[i] - LB_base_mu_U_U_51[i]
  
  f44_U_U_51[[i]] <- grepl("basetau", row.names(DP_results_U_U_51[[i]]))
  m44_U_U_51[[i]] <- DP_results_U_U_51[[i]][f44_U_U_51[[i]],]
  base_tau_U_U_51[i] <- m44_U_U_51[[i]]$`50%`
  base_tau2_U_U_51[i] <- (m44_U_U_51[[i]]$`50%`)^2
  LB_base_tau2_U_U_51[i] <- (m44_U_U_51[[i]]$`2.5%`)^2
  UB_base_tau2_U_U_51[i] <- (m44_U_U_51[[i]]$`97.5%`)^2
  Rhat_basetau_U_U_51[i] <- (m44_U_U_51[[i]]$Rhat)^2
  prec_basetau2_U_U_51[i] <- UB_base_tau2_U_U_51[i] - LB_base_tau2_U_U_51[i]
  
  f55_U_U_51[[i]] <- grepl("alpha", row.names(DP_results_U_U_51[[i]]))
  m55_U_U_51[[i]] <- DP_results_U_U_51[[i]][f55_U_U_51[[i]],]
  alpha_U_U_51[i] <- m55_U_U_51[[i]]$`50%`
  LB_alpha_U_U_51[i] <- m55_U_U_51[[i]]$`2.5%`
  UB_alpha_U_U_51[i] <- m55_U_U_51[[i]]$`97.5%`
  Rhat_alpha_U_U_51[i] <- m55_U_U_51[[i]]$Rhat
  prec_alpha_U_U_51[i] <- UB_alpha_U_U_51[i] - LB_alpha_U_U_51[i]
  
  fcl_U_U_51[[i]] <- grepl("K", row.names(DP_results_U_U_51[[i]]))  
  mcl_U_U_51[[i]] <- DP_results_U_U_51[[i]][fcl_U_U_51[[i]],]
  median_K_U_U_51[[i]] <- mcl_U_U_51[[i]]$`50%` 
  LB_K_U_U_51[[i]] <- mcl_U_U_51[[i]]$`2.5%`
  UB_K_U_U_51[[i]] <- mcl_U_U_51[[i]]$`97.5%`
  
  fp_U_U_51[[i]] <- grepl("p", row.names(DP_results_U_U_51[[i]]))
  mp_U_U_51[[i]] <- DP_results_U_U_51[[i]][fp_U_U_51[[i]],]
  mp_U_U_51[[i]] <- mp_U_U_51[[i]][!grepl("pop", row.names(mp_U_U_51[[i]])),]
  mp_U_U_51[[i]] <- mp_U_U_51[[i]][!grepl("alpha", row.names(mp_U_U_51[[i]])),]
  pi_U_U_51[[i]] <- mp_U_U_51[[i]]$mean
  
  listaDP1_U_U_51[[i]] <- cbind.data.frame(rel_effDP_U_U_51[[i]],sd_rel_effDP_U_U_51[[i]], LB_rel_effDP_U_U_51[[i]],UB_rel_effDP_U_U_51[[i]],Rhat_deltaDP_U_U_51[[i]])
  listaDP_weights_U_U_51[[i]] <- cbind.data.frame(pi_U_U_51[[i]])
  
}

numclus_U_U_51 <- unlist(median_K_U_U_51)
LB_K_U_U_51 <- unlist(LB_K_U_U_51)
UB_K_U_U_51 <- unlist(UB_K_U_U_51)

mean_alpha_U_U_51 = mean(alpha_U_U_51)
mean_alpha_LB_U_U_51 = mean(LB_alpha_U_U_51)
mean_alpha_UB_U_U_51 = mean(UB_alpha_U_U_51)

tau_DP_U_U_51 = sqrt(tau2_DP_U_U_51)
LB_tau_DP_U_U_51 = sqrt(LB_tau2_DP_U_U_51)
UB_tau_DP_U_U_51 = sqrt(UB_tau2_DP_U_U_51)
precDP_tau_U_U_51 = UB_tau_DP_U_U_51 -  LB_tau_DP_U_U_51


RhatsDP_U_U_51 <- cbind.data.frame(base_mu_U_U_51,LB_base_mu_U_U_51,UB_base_mu_U_U_51,base_tau_U_U_51,base_tau2_U_U_51, LB_base_tau2_U_U_51,UB_base_tau2_U_U_51,
                                   mu_DP_U_U_51, LB_mu_DP_U_U_51, UB_mu_DP_U_U_51, tau2_DP_U_U_51, LB_tau2_DP_U_U_51, UB_tau2_DP_U_U_51,
                                   tau_DP_U_U_51, LB_tau_DP_U_U_51, UB_tau_DP_U_U_51,
                                   precDP_mu_U_U_51, precDP_tau2_U_U_51, prec_basemu_U_U_51, prec_basetau2_U_U_51,prec_alpha_U_U_51,  alpha_U_U_51 , LB_alpha_U_U_51, UB_alpha_U_U_51,
                                   mean_alpha_U_U_51, mean_alpha_LB_U_U_51, mean_alpha_UB_U_U_51,
                                   Rhat_muDP_U_U_51, Rhat_tau2DP_U_U_51, Rhat_basemu_U_U_51, Rhat_basetau_U_U_51, Rhat_alpha_U_U_51, numclus_U_U_51, LB_K_U_U_51, UB_K_U_U_51)

##########REMOVE Rhats > 1.05 ##############

condition1_U_U_51 <- which(RhatsDP_U_U_51$Rhat_muDP_U_U_51 > 1.05)
condition2_U_U_51 <- which(RhatsDP_U_U_51$Rhat_tau2DP_U_U_51 > 1.05)
condition3_U_U_51 <- which(RhatsDP_U_U_51$Rhat_basemu_U_U_51 > 1.05)
condition4_U_U_51 <- which(RhatsDP_U_U_51$Rhat_basetau_U_U_51 > 1.05)
condition5_U_U_51 <- which(RhatsDP_U_U_51$Rhat_alpha_U_U_51 > 1.05)


dist.condDP_U_U_51 = c(condition1_U_U_51,condition2_U_U_51,condition3_U_U_51,condition4_U_U_51,condition5_U_U_51)
dist.condDP_U_U_51 = unique(dist.condDP_U_U_51)

RhatDP_out_U_U_51 = round((length(dist.condDP_U_U_51)/N.sim), 4)

############### Extract and remove the datasets with Rhat > 1.05 #########
if (length(dist.condDP_U_U_51)== 0) {
  RhatsDP_U_U_51 <- RhatsDP_U_U_51
  listaDP_U_U_51 <- listaDP1_U_U_51
  listaDP_weights_U_U_51 <- listaDP_weights_U_U_51
  N.sim <- nrow(RhatsDP_U_U_51)
  pDP_U_U_51 <- p
  
} else {
  RhatsDP_U_U_51 <- RhatsDP_U_U_51[-dist.condDP_U_U_51, ]
  listaDP_U_U_51 <- listaDP1_U_U_51[-dist.condDP_U_U_51]
  listaDP_weights_U_U_51 <- listaDP_weights_U_U_51[-dist.condDP_U_U_51]
  N.sim <- nrow(RhatsDP_U_U_51)
  pDP_U_U_51 <- p[-dist.condDP_U_U_51]
}

count_DP_U_U_51 = list()
for(i in length(listaDP_U_U_51)){
  count_DP_U_U_51[[i]] = which(listaDP_U_U_51[[i]]$`Rhat_deltaDP_U_U_51[[i]]` > 1.05 )
  tell_meDP_U_U_51 = which(count_DP_U_U_51[[i]] != 0)
}

tell_meDP_U_U_51

if(length(tell_meDP_U_U_51) == 0){
  RhatsDP_U_U_51 <- RhatsDP_U_U_51
  N.sim <- nrow(RhatsDP_U_U_51)
  listaDP_U_U_51 <- listaDP_U_U_51
  listaDP_weights_U_U_51 <- listaDP_weights_U_U_51
  pDP_U_U_51 <- pDP_U_U_51
  RhatDP_out_U_U_51 <- RhatDP_out_U_U_51
} else {
  RhatsDP_U_U_51 <- RhatsDP_U_U_51[-tell_meDP_U_U_51, ]
  listaDP_U_U_51 <- listaDP_U_U_51[-tell_meDP_U_U_51]
  listaDP_weights_U_U_51 <- listaDP_weights_U_U_51[-tell_meDP_U_U_51]
  N.sim = nrow(RhatsDP_U_U_51)
  pDP_U_U_51 <- pDP_U_U_51[-tell_meDP_U_U_51]
  RhatDP_out_U_U_51 <- RhatDP_out_U_U_51 + ((length(tell_meDP_U_U_51))/N.sim)
}

RhatDP_out_U_U_51

#### AVERAGE ABSOLUTE BIAS OF mu ######
abs_biasDPjags.bin_U_U_51 <- c()
for (i in 1:N.sim){
  abs_biasDPjags.bin_U_U_51[i] <- abs( RhatsDP_U_U_51$mu_DP_U_U_51[i] - mu) 
}
avg_abs_biasDPjags.bin_U_U_51 <- round(mean(abs_biasDPjags.bin_U_U_51),2)

#### AVERAGE RELATIVE BIAS OF mu ######
rel_biasDPjags.bin_U_U_51 <- c()
for (i in 1:N.sim){
  rel_biasDPjags.bin_U_U_51[i] <- ( RhatsDP_U_U_51$mu_DP_U_U_51[i] - mu) 
}
avg_rel_biasDPjags.bin_U_U_51 <- round(mean(rel_biasDPjags.bin_U_U_51),2)


######  MSE OF mu ##########
MSE_DPjags.bin_U_U_51 <- c()
for (i in 1:N.sim){
  MSE_DPjags.bin_U_U_51[i] <- RMSE(estimate = RhatsDP_U_U_51$mu_DP_U_U_51[i],
                                   parameter = mu,
                                   type = "RMSE",
                                   MSE = TRUE   )       
}
avg_MSE_DPjags.bin_U_U_51 <- round(mean(MSE_DPjags.bin_U_U_51), 4)

############ COVERAGE OF mu ########
interval_contains_true_mean <- function(i) { 
  ( mu >= RhatsDP_U_U_51$LB_mu_DP_U_U_51[i]) && ( mu <= RhatsDP_U_U_51$UB_mu_DP_U_U_51[i])
}
interval_contains_true_mean(i)

coverageDPjags.bin_U_U_51 <- c()
for (i in 1:N.sim){
  if (interval_contains_true_mean(i) == TRUE){
    coverageDPjags.bin_U_U_51[i] = 1
  }
  else {
    coverageDPjags.bin_U_U_51[i]=0
  }
  print(coverageDPjags.bin_U_U_51[i])
}
coverageDPjags.bin_U_U_51
coverDPjags.bin_U_U_51 <- round(mean(coverageDPjags.bin_U_U_51)*100, 2)

## BIAS OF base_mu ######
biasDPjags_basemu.bin_U_U_51 <- c()
for (i in 1:N.sim){
  biasDPjags_basemu.bin_U_U_51[i] <- abs( RhatsDP_U_U_51$base_mu_U_U_51[i] - mu) 
}
avg_biasDPjags_basemu.bin_U_U_51 <- round(mean(biasDPjags_basemu.bin_U_U_51),2)

###### MSE OF base_mu ##########
MSE_DPjags_basemu.bin_U_U_51 <- c()
for (i in 1:N.sim){
  MSE_DPjags_basemu.bin_U_U_51[i] <- RMSE(estimate =RhatsDP_U_U_51$base_mu_U_U_51[i],
                                          parameter = mu,
                                          type = "RMSE",
                                          MSE = TRUE      
  )
}
avg_MSE_DPjags_basemu.bin_U_U_51 <- round(mean(MSE_DPjags_basemu.bin_U_U_51), 4)

############ COVERAGE OF base_mu ########
interval_contains_true_mean <- function(i) { 
  ( mu >= RhatsDP_U_U_51$LB_base_mu_U_U_51[i]) && ( mu <= RhatsDP_U_U_51$UB_base_mu_U_U_51[i])
}
interval_contains_true_mean(i)

coverageDPjags_basemu.bin_U_U_51 <- c()
for (i in 1:N.sim){
  if (interval_contains_true_mean(i) == TRUE){
    coverageDPjags_basemu.bin_U_U_51[i] = 1
  }
  else {
    coverageDPjags_basemu.bin_U_U_51[i]=0
  }
  print(coverageDPjags_basemu.bin_U_U_51[i])
}
coverageDPjags_basemu.bin_U_U_51
coverDPjags_basemu.bin_U_U_51 <- round(mean(coverageDPjags_basemu.bin_U_U_51)*100, 2)

######## AVERAGE BIAS OF tau2 ###########
bias_tau2DPjags.bin_U_U_51 <- c()
for (i in 1:N.sim){
  bias_tau2DPjags.bin_U_U_51[i] <- abs(RhatsDP_U_U_51$tau2_DP_U_U_51[i]- tau2)
}
avg_bias_tau2DPjags.bin_U_U_51 <- round(mean(bias_tau2DPjags.bin_U_U_51),2)

pbias_tau2DPjags.bin_U_U_51 <- round(mean(bias_tau2DPjags.bin_U_U_51 / tau2),2)

######## MSE OF tau2 #############
MSE_DPjags.bin.tau2_U_U_51 <- c()
for (i in 1:N.sim){
  MSE_DPjags.bin.tau2_U_U_51[i] <- RMSE(estimate =RhatsDP_U_U_51$tau2_DP_U_U_51[i],
                                        parameter = tau2,
                                        type = "RMSE",
                                        MSE = TRUE         
  )
  
}
avg_MSE_DPjags.bin.tau2_U_U_51 <- round(mean(MSE_DPjags.bin.tau2_U_U_51),4)

nMSE_DPjags.bin.tau2_U_U_51 <- round(mean(MSE_DPjags.bin.tau2_U_U_51/ (tau2)^2),4)

############ COVERAGE OF tau2 ##################
interval_contains_true_mean <- function(i) { 
  ( tau2 >= RhatsDP_U_U_51$LB_tau2_DP_U_U_51[i]) && ( tau2 <= RhatsDP_U_U_51$UB_tau2_DP_U_U_51[i])
}
interval_contains_true_mean(i)

coverage_tau2DPjags.bin_U_U_51 <- c()
for (i in 1:N.sim){
  if (interval_contains_true_mean(i) == TRUE){
    coverage_tau2DPjags.bin_U_U_51[i] = 1
  }
  else  {
    coverage_tau2DPjags.bin_U_U_51[i] = 0
  }
  print(coverage_tau2DPjags.bin_U_U_51[i])
}
coverage_tau2DPjags.bin_U_U_51
cover_tau2DPjags.bin_U_U_51 <- round(mean(coverage_tau2DPjags.bin_U_U_51)*100, 2)


######## AVERAGE RELATIVE  BIAS OF tau ###########
rel_bias_tauDPjags.bin_U_U_51 <- c()
for (i in 1:N.sim){
  rel_bias_tauDPjags.bin_U_U_51[i] <- (RhatsDP_U_U_51$tau_DP_U_U_51[i]- tau)
}
avg_rel_bias_tauDPjags.bin_U_U_51 <- round(mean(rel_bias_tauDPjags.bin_U_U_51),2)

prel_bias_tauDPjags.bin_U_U_51 <- round(mean(rel_bias_tauDPjags.bin_U_U_51 / tau)*100,2)

######## AVERAGE ABSOLUTE  BIAS OF tau ###########
abs_bias_tauDPjags.bin_U_U_51 <- c()
for (i in 1:N.sim){
  abs_bias_tauDPjags.bin_U_U_51[i] <- abs(RhatsDP_U_U_51$tau_DP_U_U_51[i]- tau)
}
avg_abs_bias_tauDPjags.bin_U_U_51 <- round(mean(abs_bias_tauDPjags.bin_U_U_51),2)

pabs_bias_tauDPjags.bin_U_U_51 <- round(mean(abs_bias_tauDPjags.bin_U_U_51 / tau)*100,2)

######## MSE OF tau #############
MSE_DPjags.bin.tau_U_U_51 <- c()
for (i in 1:N.sim){
  MSE_DPjags.bin.tau_U_U_51[i] <- RMSE(estimate =RhatsDP_U_U_51$tau_DP_U_U_51[i],
                                       parameter = tau,
                                       type = "RMSE",
                                       MSE = TRUE         
  )
  
}
avg_MSE_DPjags.bin.tau_U_U_51 <- round(mean(MSE_DPjags.bin.tau_U_U_51),4)

nMSE_DPjags.bin.tau_U_U_51 <- round(mean(MSE_DPjags.bin.tau_U_U_51/ tau2),4)

############ COVERAGE OF tau ##################
interval_contains_true_mean <- function(i) { 
  ( tau >= RhatsDP_U_U_51$LB_tau_DP_U_U_51[i]) && ( tau <= RhatsDP_U_U_51$UB_tau_DP_U_U_51[i])
}
interval_contains_true_mean(i)

coverage_tauDPjags.bin_U_U_51 <- c()
for (i in 1:N.sim){
  if (interval_contains_true_mean(i) == TRUE){
    coverage_tauDPjags.bin_U_U_51[i] = 1
  }
  else  {
    coverage_tauDPjags.bin_U_U_51[i] = 0
  }
  print(coverage_tauDPjags.bin_U_U_51[i])
}
coverage_tauDPjags.bin_U_U_51
cover_tauDPjags.bin_U_U_51 <- round(mean(coverage_tauDPjags.bin_U_U_51)*100, 2)


######## AVERAGE BIAS OF base_tau2 ###########
bias_basetau2DPjags.bin_U_U_51 <- c()
for (i in 1:N.sim){
  bias_basetau2DPjags.bin_U_U_51[i] <- abs(RhatsDP_U_U_51$base_tau2_U_U_51[i]- tau2)
}
avg_bias_basetau2DPjags.bin_U_U_51 <- round(mean(bias_basetau2DPjags.bin_U_U_51),2)

######## MSE OF base_tau2 #############
MSE_DPjags.bin.basetau2_U_U_51 <- c()
for (i in 1:N.sim){
  MSE_DPjags.bin.basetau2_U_U_51[i] <- RMSE(estimate =RhatsDP_U_U_51$base_tau2_U_U_51[i],
                                            parameter = tau2,
                                            type = "RMSE",
                                            MSE = TRUE         
  )
  
}
avg_MSE_DPjags.bin.basetau2_U_U_51 <- round(mean(MSE_DPjags.bin.basetau2_U_U_51),4)

############ COVERAGE OF base_tau2 ##################
interval_contains_true_mean <- function(i) {
  ( tau2 >= RhatsDP_U_U_51$LB_base_tau2_U_U_51[i]) && ( tau2 <= RhatsDP_U_U_51$UB_base_tau2_U_U_51[i])
}
interval_contains_true_mean(i)

coverage_basetau2DPjags.bin_U_U_51 <- c()
for (i in 1:N.sim){
  if (interval_contains_true_mean(i) == TRUE){
    coverage_basetau2DPjags.bin_U_U_51[i] = 1
  }
  else  {
    coverage_basetau2DPjags.bin_U_U_51[i] = 0
  }
  print(coverage_basetau2DPjags.bin_U_U_51[i])
}
coverage_basetau2DPjags.bin_U_U_51
cover_basetau2DPjags.bin_U_U_51 <- round(mean(coverage_basetau2DPjags.bin_U_U_51)*100, 2)

######################## ABSOLUTE BIAS OF RELATIVE EFFECTS ##################
abs_bias_relDP_U_U_51 <- list()
mean_abs_bias_relDP_U_U_51 <- list()
for(i in 1:N.sim){
  for(j in 1:nrow(pDP_U_U_51[[i]])){
    abs_bias_relDP_U_U_51[[i]] <- round(abs(listaDP_U_U_51[[i]]$`rel_effDP_U_U_51[[i]]` - pDP_U_U_51[[i]]$true_eff ),2)
    mean_abs_bias_relDP_U_U_51[[i]] <- mean(abs_bias_relDP_U_U_51[[i]])
  }
}
#mean_bias_relDP

################## mean bias of the mean biases of each data set for relative effects  ###################
mean_abs_bias_of_mean_abs_bias_relDP_U_U_51 <- unlist(mean_abs_bias_relDP_U_U_51)
mean_abs_bias_of_mean_abs_bias_relDP_U_U_51  <- round(mean(mean_abs_bias_of_mean_abs_bias_relDP_U_U_51),2)

######################## RELATIVE BIAS OF RELATIVE EFFECTS ##################
rel_bias_relDP_U_U_51 <- list()
mean_rel_bias_relDP_U_U_51 <- list()
for(i in 1:N.sim){
  for(j in 1:nrow(pDP_U_U_51[[i]])){
    rel_bias_relDP_U_U_51[[i]] <- round((listaDP_U_U_51[[i]]$`rel_effDP_U_U_51[[i]]` - pDP_U_U_51[[i]]$true_eff ),2)
    mean_rel_bias_relDP_U_U_51[[i]] <- mean(rel_bias_relDP_U_U_51[[i]])
  }
}
#mean_bias_relDP

################## mean bias of the mean biases of each data set for relative effects  ###################
mean_rel_bias_of_mean_rel_bias_relDP_U_U_51 <- unlist(mean_rel_bias_relDP_U_U_51)
mean_rel_bias_of_mean_rel_bias_relDP_U_U_51  <- round(mean(mean_rel_bias_of_mean_rel_bias_relDP_U_U_51),2)

#####
########################### FOR THE OVERLAPPING ######################
dataDP_U_U_51 <- list()
overlapDP1_U_U_51 <- list()
overlapDP_U_U_51 <- list()
for(i in 1:N.sim){
  dataDP_U_U_51[[i]] <- list(DP = listaDP_U_U_51[[i]]$`rel_effDP_U_U_51[[i]]` , True = pDP_U_U_51[[i]]$true_eff  )
  overlapDP1_U_U_51[[i]] <- overlap(dataDP_U_U_51[[i]], type = "2")
  overlapDP_U_U_51[[i]] <- overlapDP1_U_U_51[[i]]$OV
}
overlapDP_U_U_51 <- unlist(overlapDP_U_U_51)
avg_overlapDP_U_U_51 <- round(mean(overlapDP_U_U_51),2)


listaDP_rel_eff_U_U_51 <- cbind( mean_abs_bias_relDP_U_U_51, mean_abs_bias_of_mean_abs_bias_relDP_U_U_51,
                                 mean_rel_bias_relDP_U_U_51, mean_rel_bias_of_mean_rel_bias_relDP_U_U_51) 

######### SAVE THE OUTPUTS ###########
DP_bin_jags_U_U_51 <- cbind(RhatsDP_U_U_51,
                            abs_biasDPjags.bin_U_U_51, rel_biasDPjags.bin_U_U_51,
                            MSE_DPjags.bin_U_U_51, biasDPjags_basemu.bin_U_U_51, MSE_DPjags_basemu.bin_U_U_51,
                            bias_tau2DPjags.bin_U_U_51,rel_bias_tauDPjags.bin_U_U_51, abs_bias_tauDPjags.bin_U_U_51, 
                            MSE_DPjags.bin.tau2_U_U_51, MSE_DPjags.bin.tau_U_U_51,
                            bias_basetau2DPjags.bin_U_U_51, MSE_DPjags.bin.basetau2_U_U_51, 
                            avg_abs_biasDPjags.bin_U_U_51, avg_rel_biasDPjags.bin_U_U_51, 
                            avg_biasDPjags_basemu.bin_U_U_51, 
                            avg_bias_tau2DPjags.bin_U_U_51,pbias_tau2DPjags.bin_U_U_51,avg_rel_bias_tauDPjags.bin_U_U_51,avg_abs_bias_tauDPjags.bin_U_U_51,
                            prel_bias_tauDPjags.bin_U_U_51,pabs_bias_tauDPjags.bin_U_U_51,
                            avg_bias_basetau2DPjags.bin_U_U_51,
                            rel_bias_tauDPjags.bin_U_U_51,abs_bias_tauDPjags.bin_U_U_51,
                            avg_MSE_DPjags.bin_U_U_51,avg_MSE_DPjags_basemu.bin_U_U_51,avg_MSE_DPjags.bin.tau2_U_U_51, nMSE_DPjags.bin.tau2_U_U_51,
                            avg_MSE_DPjags.bin.tau_U_U_51, nMSE_DPjags.bin.tau_U_U_51,
                            avg_MSE_DPjags.bin.basetau2_U_U_51,
                            coverageDPjags.bin_U_U_51 , coverageDPjags_basemu.bin_U_U_51, coverage_tau2DPjags.bin_U_U_51, coverage_basetau2DPjags.bin_U_U_51, 
                            coverage_tauDPjags.bin_U_U_51, 
                            coverDPjags.bin_U_U_51, coverDPjags_basemu.bin_U_U_51, cover_tau2DPjags.bin_U_U_51, cover_tauDPjags.bin_U_U_51,
                            cover_basetau2DPjags.bin_U_U_51, RhatDP_out_U_U_51)

Overlap_DPmodel_U_U_51 <- cbind(overlapDP_U_U_51, avg_overlapDP_U_U_51)                    

write.csv(DP_bin_jags_U_U_51, "C:\\Users\\Tianqi YU\\Desktop\\SIMULATION SCENARIOS FINAL\\NORMAL DATA\\DP MODEL\\rev_U_U_51_ResultsDP_jags_Normal_bin_data_k=14_mu=0_tau2=0.12.csv",row.names=FALSE )  
write.csv(listaDP_U_U_51, "C:\\Users\\Tianqi YU\\Desktop\\SIMULATION SCENARIOS FINAL\\NORMAL DATA\\DP MODEL\\rev_U_U_51_Relative_eff_ResultsDP_jags_Normal_bin_data_k=14_mu=0_tau2=0.12.csv",row.names=FALSE )  
write.csv(listaDP_rel_eff_U_U_51, "C:\\Users\\Tianqi YU\\Desktop\\SIMULATION SCENARIOS FINAL\\NORMAL DATA\\DP MODEL\\rev_U_U_51_Bias_Relative_eff_ResultsDP_jags_Normal_bin_data_k=14_mu=0_tau2=0.12.csv",row.names=FALSE )  
write.csv(listaDP_weights_U_U_51, "C:\\Users\\Tianqi YU\\Desktop\\SIMULATION SCENARIOS FINAL\\NORMAL DATA\\DP MODEL\\rev_U_U_51_Weights_in_clusters_jags_DP_jags_Normal_bin_data_k=14_mu=0_tau2=0.12.csv",row.names=TRUE )  
write.csv(Overlap_DPmodel_U_U_51, "C:\\Users\\Tianqi YU\\Desktop\\SIMULATION SCENARIOS FINAL\\NORMAL DATA\\DP MODEL\\rev_U_U_51_Overlap_DPmodel_Normal_bin_data_k=14_mu=0_tau2=0.12.csv",row.names=TRUE )  

##############################  SKEW NORMAL MODEL #####################
######################## Skew normal(HN) ##############
# install.packages("parallel")
# install.packages("cmdstanr")
# install.packages("bayesplot")

library(parallel)
library(cmdstanr)
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
library(bayesplot)
color_scheme_set("brightblue")
check_cmdstan_toolchain()

no_cores <- detectCores()-1
cl <- makeCluster(no_cores)
N.sim=1000
mu = 0
tau2 = 0.12
tau = sqrt(tau2)

X1=list()
for (i in 1:N.sim){ X1[[i]]=list("data"=p[[i]],"logOR"=mu,"tau_2"=tau2)}

################################################################

SN=function(X)
{
  delta = c()
  median <- c()
  LBmu <- c()
  UBmu <- c()
  pred <- c()
  LB_pred <- c()
  UB_pred <- c()
  tau <- c()
  LBtau <- c()
  UBtau <- c()
  tau2 <- c()
  LBtau2 <- c()
  UBtau2 <- c()
  skew <- c()
  UBskew <- c()
  LBskew <- c()
  Rhat_muSN <- c()
  Rhat_tau2SN <- c()
  Rhat_predSN <- c()
  Rhat_skewSN <- c()
  
  fileSN <- file.path(cmdstan_path(), "examples", "Stanmodels", "Skew_normal_model_bin_HNpred.stan")
  modSN <- cmdstan_model(fileSN)
  
  fitSN <- modSN$sample(
    data =list(ns = nrow(X$data),
               r1 = X$data$ci,
               r2 = X$data$ti,
               n1 = X$data$nci,
               n2 = X$data$nti) , 
    seed = 21014, 
    chains = 2, 
    parallel_chains = 2,
    refresh = 500,
    iter_warmup = 10000,
    iter_sampling = 50000,
    adapt_delta = 0.99
    # print update every 500 iters
  )
  
  N.studies <- nrow(X$data)
  
  sn.res <- as.data.frame(fitSN$summary())
  
  d <- grepl("delta",sn.res$variable)
  
  delta <- sn.res[d,]
  
  
  dm <- grepl("mu",sn.res$variable)
  
  median <- sn.res[dm,]$median
  LBmu <- sn.res[dm,]$q5
  UBmu <- sn.res[dm,]$q95
  Rhat_muSN <-  sn.res[dm,]$rhat
  
  
  tau2 <- c()
  dtau2 <- grepl("tau_sqr",sn.res$variable)
  
  tau2 <- sn.res[dtau2,]$median
  LBtau2 <- sn.res[dtau2,]$q5
  UBtau2 <- sn.res[dtau2,]$q95
  Rhat_tau2SN <- sn.res[dtau2,]$rhat
  
  tau <- sqrt(sn.res[dtau2,]$median)
  LBtau <- sqrt(sn.res[dtau2,]$q5)
  UBtau <- sqrt(sn.res[dtau2,]$q95)
  
  
  dxn <- grepl("pred1",sn.res$variable)
  pred <- sn.res[dxn,]$median
  LBpred <- sn.res[dxn,]$q5
  UBpred <- sn.res[dxn,]$q95
  Rhat_predSN <-  sn.res[dxn,]$rhat
  
  
  skew <- c()
  dskew <- grepl("skew",sn.res$variable)
  
  skew <- sn.res[dskew,]$median
  LBskew <- sn.res[dskew,]$q5
  UBskew <- sn.res[dskew,]$q95
  prec_skewSN <-  LBskew -  UBskew
  Rhat_skewSN <- sn.res[dskew,]$rhat
  
  sn.res1_HN<-data.frame(median=median, 
                         lowerCI=LBmu,
                         upperCI=UBmu,
                         pred=pred, 
                         LBpred=LBpred,
                         UBpred=UBpred,
                         tau2=tau2,
                         l_tau2 = LBtau2,
                         u_tau2 = UBtau2,
                         tau=tau,
                         l_tau = LBtau,
                         u_tau = UBtau,
                         skew = skew,
                         l_skew = LBskew,
                         u_skew = UBskew,
                         Rhat_muSN = Rhat_muSN,
                         Rhat_predSN = Rhat_predSN,
                         Rhat_tau2SN = Rhat_tau2SN,
                         Rhat_skewSN = Rhat_skewSN
                         
  )
  
  
  return(list("res"= sn.res,
              "res1"= sn.res1_HN,
              "delta_SN"=delta,
              "skew" = skew,
              "Rhat_muSN" = Rhat_muSN,
              "Rhat_tau2SN" = Rhat_tau2SN,
              "Rhat_skewSN" = Rhat_skewSN,
              "Rhat_predSN" = Rhat_predSN
              
              
  ))
  
}

clusterExport(cl,"X1")

clusterExport(cl,"N.sim")

clusterExport(cl, "SN")

clusterEvalQ(cl, {library(cmdstanr)}) 

clusterEvalQ(cl, {library(SimDesign)})

l2_HN <- parLapply(cl,1:N.sim, function(x) SN(X1[[x]]))


Rhat_muSN_HN = c()
Rhat_tau2SN_HN = c()
Rhat_skewSN_HN = c()
Rhat_predSN_HN = c()
skewSN_HN = c()
LBskewSN_HN = c()
UBskewSN_HN= c()
mu_SN_HN = c()
LBmu_SN_HN = c()
UBmu_SN_HN = c()
pred_SN_HN = c()
LBpred_SN_HN = c()
UBpred_SN_HN = c()
tau_SN_HN = c()
LBtau_SN_HN = c()
UBtau_SN_HN = c()
tau2_SN_HN = c()
LBtau2_SN_HN = c()
UBtau2_SN_HN = c()
prec_skewSN_HN = c()
prec_muSN_HN = c()
prec_predSN_HN = c()
prec_tauSN_HN = c()
prec_tau2SN_HN = c()
delta_SN_HN = list()

for(i in 1:N.sim){
  Rhat_muSN_HN = c(Rhat_muSN_HN, l2_HN[[i]]$Rhat_muSN)
  Rhat_tau2SN_HN = c(Rhat_tau2SN_HN, l2_HN[[i]]$Rhat_tau2SN)
  Rhat_skewSN_HN = c(Rhat_skewSN_HN, l2_HN[[i]]$Rhat_skewSN)
  Rhat_predSN_HN = c(Rhat_predSN_HN, l2_HN[[i]]$Rhat_predSN)
  
  skewSN_HN=c(skewSN_HN, l2_HN[[i]]$skew)
  LBskewSN_HN =c(LBskewSN_HN, l2_HN[[i]]$res1$l_skew)
  UBskewSN_HN = c(UBskewSN_HN, l2_HN[[i]]$res1$u_skew)
  prec_skewSN_HN = c(prec_skewSN_HN,(l2_HN[[i]]$res1$u_skew - l2_HN[[i]]$res1$l_skew) )
  
  mu_SN_HN = c(mu_SN_HN,l2_HN[[i]]$res1$median)
  LBmu_SN_HN = c(LBmu_SN_HN,l2_HN[[i]]$res1$lowerCI)
  UBmu_SN_HN = c(UBmu_SN_HN,l2_HN[[i]]$res1$upperCI)
  prec_muSN_HN = c(prec_muSN_HN,(l2_HN[[i]]$res1$upperCI - l2_HN[[i]]$res1$lowerCI) )
  
  pred_SN_HN = c(pred_SN_HN,l2_HN[[i]]$res1$pred)
  LBpred_SN_HN = c(LBpred_SN_HN,l2_HN[[i]]$res1$LBpred)
  UBpred_SN_HN = c(UBpred_SN_HN,l2_HN[[i]]$res1$UBpred)
  prec_predSN_HN = c(prec_predSN_HN,(l2_HN[[i]]$res1$UBpred - l2_HN[[i]]$res1$LBpred) )
  
  tau2_SN_HN = c(tau2_SN_HN,l2_HN[[i]]$res1$tau2)
  LBtau2_SN_HN = c(LBtau2_SN_HN,l2_HN[[i]]$res1$l_tau2)
  UBtau2_SN_HN = c(UBtau2_SN_HN,l2_HN[[i]]$res1$u_tau2)
  prec_tau2SN_HN = c(prec_tau2SN_HN,(l2_HN[[i]]$res1$u_tau2 - l2_HN[[i]]$res1$l_tau2) )
  
  tau_SN_HN = c(tau_SN_HN,sqrt(l2_HN[[i]]$res1$tau2))
  LBtau_SN_HN = c(LBtau_SN_HN,sqrt(l2_HN[[i]]$res1$l_tau2))
  UBtau_SN_HN = c(UBtau_SN_HN,sqrt(l2_HN[[i]]$res1$u_tau2))
  prec_tauSN_HN = c(prec_tauSN_HN,(sqrt(l2_HN[[i]]$res1$u_tau2) - sqrt(l2_HN[[i]]$res1$l_tau2)) )
  
  delta_SN_HN[[i]] = (l2_HN[[i]]$delta_SN)
}

data_SN_HN = cbind.data.frame(mu_SN_HN, LBmu_SN_HN, UBmu_SN_HN,
                              pred_SN_HN, LBpred_SN_HN, UBpred_SN_HN,
                              tau_SN_HN, LBtau_SN_HN, UBtau_SN_HN,
                              tau2_SN_HN, LBtau2_SN_HN, UBtau2_SN_HN, skewSN_HN, LBskewSN_HN, UBskewSN_HN, Rhat_muSN_HN,
                              Rhat_predSN_HN,
                              Rhat_tau2SN_HN , Rhat_skewSN_HN ,prec_muSN_HN, prec_tau2SN_HN,prec_predSN_HN, prec_tauSN_HN, prec_skewSN_HN)

############### REMOVE DATA SETS WITH RHATS > 1.05 ############

condition1_HN <- which(data_SN_HN$Rhat_muSN_HN > 1.05)
condition2_HN <- which(data_SN_HN$Rhat_tau2SN_HN > 1.05)
condition3_HN <- which(data_SN_HN$Rhat_skewSN_HN > 1.05)

dist.condSN_HN = c(condition1_HN,condition2_HN,condition3_HN)
dist.condSN_HN = unique(dist.condSN_HN)

RhatSN_out_HN = round((length(dist.condSN_HN)/N.sim),4)

if (length(dist.condSN_HN)== 0) {
  data_SN_HN <- data_SN_HN
  delta_SN_HN <- delta_SN_HN
  N.sim <- nrow(data_SN_HN)
  pSN_HN <- p
} else {
  data_SN_HN = data_SN_HN[-dist.condSN_HN, ]
  delta_SN_HN <- delta_SN_HN[-dist.condSN_HN]
  N.sim = nrow(data_SN_HN)
  pSN_HN <- p[-dist.condSN_HN]
}

count_SN_HN = list()
for(i in length(delta_SN_HN)){
  count_SN_HN[[i]] = which(delta_SN_HN[[i]]$rhat > 1.05 )
  tell_meSN_HN = which(count_SN_HN[[i]] != 0)
}

tell_meSN_HN

if(length(tell_meSN_HN) == 0){
  data_SN_HN <- data_SN_HN
  delta_SN_HN <- delta_SN_HN
  N.sim <- nrow(data_SN_HN)
  pSN_HN <- pSN_HN
  RhatSN_out_HN <- RhatSN_out_HN
} else {
  data_SN_HN = data_SN_HN[-tell_meSN_HN, ]
  delta_SN_HN <- delta_SN_HN[-tell_meSN_HN]
  N.sim = nrow(data_SN_HN)
  pSN_HN <- pSN_HN[-tell_meSN_HN]
  RhatSN_out_HN <- RhatSN_out_HN + ((length(tell_meSN_HN))/N.sim)
}

RhatSN_out_HN
############ EXTRACTION OF RELATIVE EFFECTS ########
deltaSN_HN <- data.frame()
for (i in 1:N.sim) {
  temp_df_HN<- data.frame(simulation = i, deltaSN_HN = delta_SN_HN[[i]])
  deltaSN_HN <- rbind(deltaSN_HN, temp_df_HN)
}

###################### SIMULATION MEASURES ###############
########## AVERAGE ABSOLUTE BIAS OF mu ##########
abs_biasSNjags.bin_HN <- c()
for (i in 1:N.sim){
  abs_biasSNjags.bin_HN[i] <- abs( data_SN_HN$mu_SN_HN[i] - mu) 
}
avg_abs_biasSNjags.bin_HN <- round(mean(abs_biasSNjags.bin_HN), 2)


########## RELATIVE ABSOLUTE BIAS OF mu ##########
rel_biasSNjags.bin_HN <- c()
for (i in 1:N.sim){
  rel_biasSNjags.bin_HN[i] <- ( data_SN_HN$mu_SN_HN[i] - mu) 
}
avg_rel_biasSNjags.bin_HN <- round(mean(rel_biasSNjags.bin_HN), 2)

############### MSE OF mu  ##################
MSE_SNjags.bin_HN <- c()
for (i in 1:N.sim){
  MSE_SNjags.bin_HN[i] <- RMSE(estimate =data_SN_HN$mu_SN_HN[i],
                               parameter = mu,
                               type = "RMSE",
                               MSE = TRUE,           
                               percent = FALSE,
                               unname = FALSE)
}
avg_MSE_SNjags.bin_HN <-  round(mean(MSE_SNjags.bin_HN),4)

########## COVERAGE OF mu ##########
interval_contains_true_mean <- function(i) { 
  ( mu >= data_SN_HN$LBmu_SN_HN[i]) && ( mu <=  data_SN_HN$UBmu_SN_HN[i])
}
interval_contains_true_mean(i)

coverageSNjags.bin_HN <- c()
for (i in 1:N.sim){0
  if (interval_contains_true_mean(i) == TRUE){
    coverageSNjags.bin_HN[i] = 1
  }
  else {
    coverageSNjags.bin_HN[i]=0
  }
  print(coverageSNjags.bin_HN[i])
}
coverageSNjags.bin_HN
coverSNjags.bin_HN <- round(mean(coverageSNjags.bin_HN)*100, 2)

################
##### AVERAGE BIAS OF tau2 #########
bias_tau2SNjags.bin_HN <- c()
for (i in 1:N.sim){
  bias_tau2SNjags.bin_HN[i] <- abs(data_SN_HN$tau2_SN_HN[i] - tau2)
}
avg_bias_tau2SNjags.bin_HN <- round(mean(bias_tau2SNjags.bin_HN),2)

pbias_tau2SNjags.bin_HN = round(mean(bias_tau2SNjags.bin_HN / tau2), 2)


###########  MSE OF tau2 ################
MSE_SNjags.tau2.bin_HN <- c()
for (i in 1:N.sim){
  MSE_SNjags.tau2.bin_HN[i] <- RMSE(estimate = data_SN_HN$tau2_SN_HN[i],
                                    parameter = tau2,
                                    type = "RMSE",
                                    MSE = TRUE           
  )
}
avg_MSE_SNjags.tau2.bin_HN <- round(mean(MSE_SNjags.tau2.bin_HN),4)

nMSE_SNjags.tau2.bin_HN <- round(mean(MSE_SNjags.tau2.bin_HN / (tau2)^2), 2)
############# COVERAGE OF tau2 ###################
interval_contains_true_mean <- function(i) { 
  ( tau2 >= data_SN_HN$LBtau2_SN_HN[i]) && ( tau2 <=  data_SN_HN$UBtau2_SN_HN[i])
}
interval_contains_true_mean(i)

coverage_tau2SNjags.bin_HN <- c()
for (i in 1:N.sim){
  if (interval_contains_true_mean(i) == TRUE){
    coverage_tau2SNjags.bin_HN[i] = 1
  }
  else  {
    coverage_tau2SNjags.bin_HN[i] = 0
  }
  print(coverage_tau2SNjags.bin_HN[i])
}
coverage_tau2SNjags.bin_HN
cover_tau2SNjags.bin_HN <- round(mean(coverage_tau2SNjags.bin_HN)*100,2)


##### AVERAGE relative BIAS OF tau #########
bias_rel_tauSNjags.bin_HN <- c()
for (i in 1:N.sim){
  bias_rel_tauSNjags.bin_HN[i] <- (data_SN_HN$tau_SN_HN[i] - tau)
}
avg_rel_bias_tauSNjags.bin_HN <- round(mean(bias_rel_tauSNjags.bin_HN),2)

prel_bias_tauSNjags.bin_HN = round(mean(bias_rel_tauSNjags.bin_HN / tau)*100, 2)

##### AVERAGE absolute BIAS OF tau #########
abs_bias_tauSNjags.bin_HN <- c()
for (i in 1:N.sim){
  abs_bias_tauSNjags.bin_HN[i] <- abs(data_SN_HN$tau_SN_HN[i] - tau)
}
avg_abs_bias_tauSNjags.bin_HN <- round(mean(abs_bias_tauSNjags.bin_HN),2)

pabs_bias_tauSNjags.bin_HN = round(mean(abs_bias_tauSNjags.bin_HN / tau)*100, 2)

###########  MSE OF tau ################
MSE_SNjags.tau.bin_HN <- c()
for (i in 1:N.sim){
  MSE_SNjags.tau.bin_HN[i] <- RMSE(estimate = data_SN_HN$tau_SN_HN[i],
                                   parameter = tau,
                                   type = "RMSE",
                                   MSE = TRUE           
  )
}
avg_MSE_SNjags.tau.bin_HN <- round(mean(MSE_SNjags.tau.bin_HN),4)

nMSE_SNjags.tau.bin_HN <- round(mean(MSE_SNjags.tau.bin_HN / tau2), 2)
############# COVERAGE OF tau ###################
interval_contains_true_mean <- function(i) { 
  ( tau >= data_SN_HN$LBtau_SN_HN[i]) && ( tau <=  data_SN_HN$UBtau_SN_HN[i])
}
interval_contains_true_mean(i)

coverage_tauSNjags.bin_HN <- c()
for (i in 1:N.sim){
  if (interval_contains_true_mean(i) == TRUE){
    coverage_tauSNjags.bin_HN[i] = 1
  }
  else  {
    coverage_tauSNjags.bin_HN[i] = 0
  }
  print(coverage_tauSNjags.bin_HN[i])
}
coverage_tauSNjags.bin_HN
cover_tauSNjags.bin_HN <- round(mean(coverage_tauSNjags.bin_HN)*100,2)


######################## ABSOLUTE BIAS OF RELATIVE EFFECTS ##################
abs_bias_relSN_HN <- list()
mean_abs_bias_relSN_HN <- list()
for(i in 1:N.sim){
  for(j in 1:nrow(pSN_HN[[i]])){
    abs_bias_relSN_HN[[i]] <- round(abs(delta_SN_HN[[i]]$median - pSN_HN[[i]]$true_eff ),2)
    mean_abs_bias_relSN_HN[[i]] <- mean(abs_bias_relSN_HN[[i]])
  }
}


################## mean bias of the mean biases of each data set for relative effects  ###################
mean_abs_bias_of_mean_abs_bias_relSN_HN <- unlist(mean_abs_bias_relSN_HN)
mean_abs_bias_of_mean_abs_bias_relSN_HN  <- round(mean(mean_abs_bias_of_mean_abs_bias_relSN_HN ),2)

######################## RELATIVE BIAS OF RELATIVE EFFECTS ##################
rel_bias_relSN_HN <- list()
mean_rel_bias_relSN_HN <- list()
for(i in 1:N.sim){
  for(j in 1:nrow(pSN_HN[[i]])){
    rel_bias_relSN_HN[[i]] <- round((delta_SN_HN[[i]]$median - pSN_HN[[i]]$true_eff ),2)
    mean_rel_bias_relSN_HN[[i]] <- mean(rel_bias_relSN_HN[[i]])
  }
}


################## mean bias of the mean biases of each data set for relative effects  ###################
mean_rel_bias_of_mean_rel_bias_relSN_HN <- unlist(mean_rel_bias_relSN_HN)
mean_rel_bias_of_mean_rel_bias_relSN_HN  <- round(mean(mean_rel_bias_of_mean_rel_bias_relSN_HN ),2)

########################### FOR THE OVERLAPPING ######################
dataSN1_HN <- list()
overlapSN1_HN <- list()
overlapSN_HN <- list()
for(i in 1:N.sim){
  dataSN1_HN[[i]] <- list(SN = delta_SN_HN[[i]]$median , True = pSN_HN[[i]]$true_eff)
  overlapSN1_HN[[i]] <- overlap(dataSN1_HN[[i]], type = "2")
  overlapSN_HN[[i]] <- overlapSN1_HN[[i]]$OV
}
overlapSN_HN <- unlist(overlapSN_HN)
avg_overlapSN_HN <- round(mean(overlapSN_HN),2)



listaSN_rel_eff_HN <- cbind(mean_abs_bias_relSN_HN,mean_abs_bias_of_mean_abs_bias_relSN_HN,
                            mean_rel_bias_relSN_HN,mean_rel_bias_of_mean_rel_bias_relSN_HN) 

len_mu_SN_HN = mean(data_SN_HN$prec_muSN_HN)
len_PI_SN_HN = mean(data_SN_HN$prec_predSN_HN)

SNormal_HN = cbind(data_SN_HN, abs_biasSNjags.bin_HN, rel_biasSNjags.bin_HN, 
                   avg_abs_biasSNjags.bin_HN, avg_rel_biasSNjags.bin_HN,
                   MSE_SNjags.bin_HN, avg_MSE_SNjags.bin_HN,
                   coverageSNjags.bin_HN, coverSNjags.bin_HN,
                   bias_tau2SNjags.bin_HN, avg_bias_tau2SNjags.bin_HN, pbias_tau2SNjags.bin_HN,
                   MSE_SNjags.tau2.bin_HN, avg_MSE_SNjags.tau2.bin_HN, nMSE_SNjags.tau2.bin_HN,
                   coverage_tau2SNjags.bin_HN, cover_tau2SNjags.bin_HN,
                   
                   abs_bias_tauSNjags.bin_HN, avg_abs_bias_tauSNjags.bin_HN, pabs_bias_tauSNjags.bin_HN,
                   
                   bias_rel_tauSNjags.bin_HN, avg_rel_bias_tauSNjags.bin_HN, prel_bias_tauSNjags.bin_HN,
                   MSE_SNjags.tau.bin_HN, avg_MSE_SNjags.tau.bin_HN, nMSE_SNjags.tau.bin_HN,
                   coverage_tauSNjags.bin_HN, cover_tauSNjags.bin_HN,
                   len_mu_SN_HN ,
                   len_PI_SN_HN,
                   RhatSN_out_HN)

Overlap_SNmodel_HN <- cbind(overlapSN_HN, avg_overlapSN_HN)

listaSN_HN = cbind(data_SN_HN$skewSN_HN, data_SN_HN$LBskewSN_HN, data_SN_HN$UBskewSN_HN)


write.csv(SNormal_HN, "C:\\Users\\Tianqi YU\\Desktop\\SIMULATION SCENARIOS FINAL\\NORMAL DATA\\SKEW NORMAL MODEL\\rev_HN_ResultsSN_Normal_bin_data_k=14_mu=0_tau2=0.12.csv",row.names=FALSE )  
write.csv(deltaSN_HN, "C:\\Users\\Tianqi YU\\Desktop\\SIMULATION SCENARIOS FINAL\\NORMAL DATA\\SKEW NORMAL MODEL\\rev_HN_Relative_eff_ResultsSN_Normal_bin_data_k=14_mu=0_tau2=0.12.csv" ) 
write.csv(listaSN_rel_eff_HN, "C:\\Users\\Tianqi YU\\Desktop\\SIMULATION SCENARIOS FINAL\\NORMAL DATA\\SKEW NORMAL MODEL\\rev_HN_Bias_Relative_eff_ResultsSN_Normal_bin_data_k=14_mu=0_tau2=0.12.csv" ) 
write.csv(listaSN_HN, "C:\\Users\\Tianqi YU\\Desktop\\SIMULATION SCENARIOS FINAL\\NORMAL DATA\\SKEW NORMAL MODEL\\rev_HN_Skew_ResultsSN_Normal_bin_data_k=14_mu=0_tau2=0.12.csv" ) 
write.csv(Overlap_SNmodel_HN, "C:\\Users\\Tianqi YU\\Desktop\\SIMULATION SCENARIOS FINAL\\NORMAL DATA\\SKEW NORMAL MODEL\\rev_HN_Overlap_SNmodel_Normal_bin_data_k=14_mu=0_tau2=0.12.csv" ) 


####### Skew normal(Unif) ##########

library(parallel)
library(cmdstanr)
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
library(bayesplot)
color_scheme_set("brightblue")
check_cmdstan_toolchain()

no_cores <- detectCores()-1
cl <- makeCluster(no_cores)
N.sim=1000
mu = 0
tau2 = 0.12
tau = sqrt(tau2)

X1=list()
for (i in 1:N.sim){ X1[[i]]=list("data"=p[[i]],"logOR"=mu,"tau_2"=tau2)}

################################################################

SN=function(X)
{
  delta = c()
  median <- c()
  LBmu <- c()
  UBmu <- c()
  pred <- c()
  LB_pred <- c()
  UB_pred <- c()
  tau <- c()
  LBtau <- c()
  UBtau <- c()
  tau2 <- c()
  LBtau2 <- c()
  UBtau2 <- c()
  skew <- c()
  UBskew <- c()
  LBskew <- c()
  Rhat_muSN <- c()
  Rhat_tau2SN <- c()
  Rhat_skewSN <- c()
  Rhat_predSN <- c()
  
  fileSN <- file.path(cmdstan_path(), "examples", "Stanmodels", "Skew_normal_model_binUpred.stan")
  modSN <- cmdstan_model(fileSN)
  
  fitSN <- modSN$sample(
    data =list(ns = nrow(X$data),
               r1 = X$data$ci,
               r2 = X$data$ti,
               n1 = X$data$nci,
               n2 = X$data$nti) , 
    seed = 2107, 
    chains = 2, 
    parallel_chains = 2,
    refresh = 500,
    iter_warmup = 10000,
    iter_sampling = 50000,
    adapt_delta = 0.99
    # print update every 500 iters
  )
  
  N.studies <- nrow(X$data)
  
  sn.res <- as.data.frame(fitSN$summary())
  
  d <- grepl("delta",sn.res$variable)
  
  delta <- sn.res[d,]
  
  
  dm <- grepl("mu",sn.res$variable)
  
  median <- sn.res[dm,]$median
  LBmu <- sn.res[dm,]$q5
  UBmu <- sn.res[dm,]$q95
  Rhat_muSN <-  sn.res[dm,]$rhat
  
  
  tau2 <- c()
  dtau2 <- grepl("tau_sqr",sn.res$variable)
  
  tau2 <- sn.res[dtau2,]$median
  LBtau2 <- sn.res[dtau2,]$q5
  UBtau2 <- sn.res[dtau2,]$q95
  Rhat_tau2SN <- sn.res[dtau2,]$rhat
  
  tau <- sqrt(sn.res[dtau2,]$median)
  LBtau <- sqrt(sn.res[dtau2,]$q5)
  UBtau <- sqrt(sn.res[dtau2,]$q95)
  
  dxn <- grepl("pred1",sn.res$variable)
  pred <- sn.res[dxn,]$median
  LBpred <- sn.res[dxn,]$q5
  UBpred <- sn.res[dxn,]$q95
  Rhat_predSN <-  sn.res[dxn,]$rhat
  
  skew <- c()
  dskew <- grepl("skew",sn.res$variable)
  
  skew <- sn.res[dskew,]$median
  LBskew <- sn.res[dskew,]$q5
  UBskew <- sn.res[dskew,]$q95
  Rhat_skewSN <- sn.res[dskew,]$rhat
  
  sn.res1_U<-data.frame(median=median, 
                        lowerCI=LBmu,
                        upperCI=UBmu,
                        pred=pred, 
                        LBpred=LBpred,
                        UBpred=UBpred,
                        tau2=tau2,
                        l_tau2 = LBtau2,
                        u_tau2 = UBtau2,
                        tau=tau,
                        l_tau = LBtau,
                        u_tau = UBtau,
                        skew = skew,
                        l_skew = LBskew,
                        u_skew = UBskew,
                        Rhat_muSN = Rhat_muSN,
                        Rhat_predSN = Rhat_predSN,
                        Rhat_tau2SN = Rhat_tau2SN,
                        Rhat_skewSN = Rhat_skewSN
  )
  
  
  return(list("res"= sn.res,
              "res1"= sn.res1_U,
              "delta_SN"=delta,
              "skew" = skew,
              "Rhat_muSN" = Rhat_muSN,
              "Rhat_tau2SN" = Rhat_tau2SN,
              "Rhat_skewSN" = Rhat_skewSN,
              "Rhat_predSN" = Rhat_predSN
              
  ))
  
}

clusterExport(cl,"X1")

clusterExport(cl,"N.sim")

clusterExport(cl, "SN")

clusterEvalQ(cl, {library(cmdstanr)}) 

clusterEvalQ(cl, {library(SimDesign)})

l2_U <- parLapply(cl,1:N.sim, function(x) SN(X1[[x]]))



Rhat_muSN_U = c()
Rhat_tau2SN_U = c()
Rhat_skewSN_U = c()
Rhat_predSN_U = c()
skewSN_U = c()
LBskewSN_U = c()
UBskewSN_U= c()
mu_SN_U = c()
LBmu_SN_U = c()
UBmu_SN_U = c()
pred_SN_U = c()
LBpred_SN_U = c()
UBpred_SN_U = c()
tau_SN_U = c()
LBtau_SN_U = c()
UBtau_SN_U = c()
tau2_SN_U = c()
LBtau2_SN_U = c()
UBtau2_SN_U = c()
prec_skewSN_U = c()
prec_muSN_U = c()
prec_predSN_U = c()
prec_tauSN_U = c()
prec_tau2SN_U = c()
delta_SN_U = list()


for(i in 1:N.sim){
  Rhat_muSN_U = c(Rhat_muSN_U, l2_U[[i]]$Rhat_muSN)
  Rhat_tau2SN_U = c(Rhat_tau2SN_U, l2_U[[i]]$Rhat_tau2SN)
  Rhat_skewSN_U = c(Rhat_skewSN_U, l2_U[[i]]$Rhat_skewSN)
  Rhat_predSN_U = c(Rhat_predSN_U, l2_U[[i]]$Rhat_predSN)
  
  skewSN_U=c(skewSN_U, l2_U[[i]]$skew)
  LBskewSN_U =c(LBskewSN_U, l2_U[[i]]$res1$l_skew)
  UBskewSN_U = c(UBskewSN_U, l2_U[[i]]$res1$u_skew)
  prec_skewSN_U = c(prec_skewSN_U,(l2_U[[i]]$res1$u_skew - l2_U[[i]]$res1$l_skew) )
  
  mu_SN_U = c(mu_SN_U,l2_U[[i]]$res1$median)
  LBmu_SN_U = c(LBmu_SN_U,l2_U[[i]]$res1$lowerCI)
  UBmu_SN_U = c(UBmu_SN_U,l2_U[[i]]$res1$upperCI)
  prec_muSN_U = c(prec_muSN_U,(l2_U[[i]]$res1$upperCI - l2_U[[i]]$res1$lowerCI) )
  
  pred_SN_U = c(pred_SN_U,l2_U[[i]]$res1$pred)
  LBpred_SN_U = c(LBpred_SN_U,l2_U[[i]]$res1$LBpred)
  UBpred_SN_U = c(UBpred_SN_U,l2_U[[i]]$res1$UBpred)
  prec_predSN_U = c(prec_predSN_U,(l2_U[[i]]$res1$UBpred - l2_U[[i]]$res1$LBpred) )
  
  tau2_SN_U = c(tau2_SN_U,l2_U[[i]]$res1$tau2)
  LBtau2_SN_U = c(LBtau2_SN_U,l2_U[[i]]$res1$l_tau2)
  UBtau2_SN_U = c(UBtau2_SN_U,l2_U[[i]]$res1$u_tau2)
  prec_tau2SN_U = c(prec_tau2SN_U,(l2_U[[i]]$res1$u_tau2 - l2_U[[i]]$res1$l_tau2) )
  
  tau_SN_U = c(tau_SN_U,sqrt(l2_U[[i]]$res1$tau2))
  LBtau_SN_U = c(LBtau_SN_U,sqrt(l2_U[[i]]$res1$l_tau2))
  UBtau_SN_U = c(UBtau_SN_U,sqrt(l2_U[[i]]$res1$u_tau2))
  prec_tauSN_U = c(prec_tauSN_U,(sqrt(l2_U[[i]]$res1$u_tau2) - sqrt(l2_U[[i]]$res1$l_tau2)) )
  
  delta_SN_U[[i]] = (l2_U[[i]]$delta_SN)
}

data_SN_U = cbind.data.frame(mu_SN_U, LBmu_SN_U, UBmu_SN_U,
                             pred_SN_U, LBpred_SN_U, UBpred_SN_U,
                             tau_SN_U, LBtau_SN_U, UBtau_SN_U,
                             tau2_SN_U, LBtau2_SN_U, UBtau2_SN_U, skewSN_U, LBskewSN_U, UBskewSN_U, Rhat_muSN_U,
                             Rhat_predSN_U,
                             Rhat_tau2SN_U , Rhat_skewSN_U ,prec_muSN_U, prec_tau2SN_U,prec_predSN_U, prec_tauSN_U, prec_skewSN_U)

############### REMOVE DATA SETS WITH RHATS > 1.05 ############

condition1_U <- which(data_SN_U$Rhat_muSN_U > 1.05)
condition2_U <- which(data_SN_U$Rhat_tau2SN_U > 1.05)
condition3_U <- which(data_SN_U$Rhat_skewSN_U > 1.05)

dist.condSN_U = c(condition1_U,condition2_U,condition3_U)
dist.condSN_U = unique(dist.condSN_U)

RhatSN_out_U = round((length(dist.condSN_U)/N.sim),4)

if (length(dist.condSN_U)== 0) {
  data_SN_U <- data_SN_U
  delta_SN_U <- delta_SN_U
  N.sim <- nrow(data_SN_U)
  pSN_U <- p
} else {
  data_SN_U = data_SN_U[-dist.condSN_U, ]
  delta_SN_U <- delta_SN_U[-dist.condSN_U]
  N.sim = nrow(data_SN_U)
  pSN_U <- p[-dist.condSN_U]
}

count_SN_U = list()
for(i in length(delta_SN_U)){
  count_SN_U[[i]] = which(delta_SN_U[[i]]$rhat > 1.05 )
  tell_meSN_U = which(count_SN_U[[i]] != 0)
}

tell_meSN_U

if(length(tell_meSN_U) == 0){
  data_SN_U <- data_SN_U
  delta_SN_U <- delta_SN_U
  N.sim <- nrow(data_SN_U)
  pSN_U <- pSN_U
  RhatSN_out_U <- RhatSN_out_U
} else {
  data_SN_U = data_SN_U[-tell_meSN_U, ]
  delta_SN_U <- delta_SN_U[-tell_meSN_U]
  N.sim = nrow(data_SN_U)
  pSN_U <- pSN_U[-tell_meSN_U]
  RhatSN_out_U <- RhatSN_out_U + ((length(tell_meSN_U))/N.sim)
}

RhatSN_out_U
############ EXTRACTION OF RELATIVE EFFECTS ########
deltaSN_U <- data.frame()
for (i in 1:N.sim) {
  temp_df_U<- data.frame(simulation = i, deltaSN_U = delta_SN_U[[i]])
  deltaSN_U <- rbind(deltaSN_U, temp_df_U)
}

###################### SIMULATION MEASURES ###############
########## AVERAGE ABSOLUTE BIAS OF mu ##########
abs_biasSNjags.bin_U <- c()
for (i in 1:N.sim){
  abs_biasSNjags.bin_U[i] <- abs( data_SN_U$mu_SN_U[i] - mu) 
}
avg_abs_biasSNjags.bin_U <- round(mean(abs_biasSNjags.bin_U), 2)

########## AVERAGE RELATIVE BIAS OF mu ##########
rel_biasSNjags.bin_U <- c()
for (i in 1:N.sim){
  rel_biasSNjags.bin_U[i] <- ( data_SN_U$mu_SN_U[i] - mu) 
}
avg_rel_biasSNjags.bin_U <- round(mean(rel_biasSNjags.bin_U), 2)

############### MSE OF mu  ##################
MSE_SNjags.bin_U <- c()
for (i in 1:N.sim){
  MSE_SNjags.bin_U[i] <- RMSE(estimate =data_SN_U$mu_SN_U[i],
                              parameter = mu,
                              type = "RMSE",
                              MSE = TRUE,           
                              percent = FALSE,
                              unname = FALSE)
}
avg_MSE_SNjags.bin_U <-  round(mean(MSE_SNjags.bin_U),4)

########## COVERAGE OF mu ##########
interval_contains_true_mean <- function(i) { 
  ( mu >= data_SN_U$LBmu_SN_U[i]) && ( mu <=  data_SN_U$UBmu_SN_U[i])
}
interval_contains_true_mean(i)

coverageSNjags.bin_U <- c()
for (i in 1:N.sim){0
  if (interval_contains_true_mean(i) == TRUE){
    coverageSNjags.bin_U[i] = 1
  }
  else {
    coverageSNjags.bin_U[i]=0
  }
  print(coverageSNjags.bin_U[i])
}
coverageSNjags.bin_U
coverSNjags.bin_U <- round(mean(coverageSNjags.bin_U)*100, 2)

################
##### AVERAGE BIAS OF tau2 #########
bias_tau2SNjags.bin_U <- c()
for (i in 1:N.sim){
  bias_tau2SNjags.bin_U[i] <- abs(data_SN_U$tau2_SN_U[i] - tau2)
}
avg_bias_tau2SNjags.bin_U <- round(mean(bias_tau2SNjags.bin_U),2)

pbias_tau2SNjags.bin_U = round(mean(bias_tau2SNjags.bin_U / tau2), 2)

###########  MSE OF tau2 ################
MSE_SNjags.tau2.bin_U <- c()
for (i in 1:N.sim){
  MSE_SNjags.tau2.bin_U[i] <- RMSE(estimate = data_SN_U$tau2_SN_U[i],
                                   parameter = tau2,
                                   type = "RMSE",
                                   MSE = TRUE           
  )
}
avg_MSE_SNjags.tau2.bin_U <- round(mean(MSE_SNjags.tau2.bin_U),4)

nMSE_SNjags.tau2.bin_U <- round(mean(MSE_SNjags.tau2.bin_U / (tau2)^2), 2)
############# COVERAGE OF tau2 ###################
interval_contains_true_mean <- function(i) { 
  ( tau2 >= data_SN_U$LBtau2_SN_U[i]) && ( tau2 <=  data_SN_U$UBtau2_SN_U[i])
}
interval_contains_true_mean(i)

coverage_tau2SNjags.bin_U <- c()
for (i in 1:N.sim){
  if (interval_contains_true_mean(i) == TRUE){
    coverage_tau2SNjags.bin_U[i] = 1
  }
  else  {
    coverage_tau2SNjags.bin_U[i] = 0
  }
  print(coverage_tau2SNjags.bin_U[i])
}
coverage_tau2SNjags.bin_U
cover_tau2SNjags.bin_U <- round(mean(coverage_tau2SNjags.bin_U)*100,2)


##### AVERAGE relative BIAS OF tau #########
bias_rel_tauSNjags.bin_U <- c()
for (i in 1:N.sim){
  bias_rel_tauSNjags.bin_U[i] <- (data_SN_U$tau_SN_U[i] - tau)
}
avg_rel_bias_tauSNjags.bin_U <- round(mean(bias_rel_tauSNjags.bin_U),2)

prel_bias_tauSNjags.bin_U = round(mean(bias_rel_tauSNjags.bin_U / tau)*100, 2)

##### AVERAGE absolute BIAS OF tau #########
abs_bias_tauSNjags.bin_U <- c()
for (i in 1:N.sim){
  abs_bias_tauSNjags.bin_U[i] <- abs(data_SN_U$tau_SN_U[i] - tau)
}
avg_abs_bias_tauSNjags.bin_U <- round(mean(abs_bias_tauSNjags.bin_U),2)

pabs_bias_tauSNjags.bin_U = round(mean(abs_bias_tauSNjags.bin_U / tau)*100, 2)

###########  MSE OF tau ################
MSE_SNjags.tau.bin_U <- c()
for (i in 1:N.sim){
  MSE_SNjags.tau.bin_U[i] <- RMSE(estimate = data_SN_U$tau_SN_U[i],
                                  parameter = tau,
                                  type = "RMSE",
                                  MSE = TRUE           
  )
}
avg_MSE_SNjags.tau.bin_U <- round(mean(MSE_SNjags.tau.bin_U),4)

nMSE_SNjags.tau.bin_U <- round(mean(MSE_SNjags.tau.bin_U / tau2), 2)
############# COVERAGE OF tau ###################
interval_contains_true_mean <- function(i) { 
  ( tau >= data_SN_U$LBtau_SN_U[i]) && ( tau <=  data_SN_U$UBtau_SN_U[i])
}
interval_contains_true_mean(i)

coverage_tauSNjags.bin_U <- c()
for (i in 1:N.sim){
  if (interval_contains_true_mean(i) == TRUE){
    coverage_tauSNjags.bin_U[i] = 1
  }
  else  {
    coverage_tauSNjags.bin_U[i] = 0
  }
  print(coverage_tauSNjags.bin_U[i])
}
coverage_tauSNjags.bin_U
cover_tauSNjags.bin_U <- round(mean(coverage_tauSNjags.bin_U)*100,2)



######################## ABSOLUTE BIAS OF RELATIVE EFFECTS ##################
abs_bias_relSN_U <- list()
mean_abs_bias_relSN_U <- list()
for(i in 1:N.sim){
  for(j in 1:nrow(pSN_U[[i]])){
    abs_bias_relSN_U[[i]] <- round(abs(delta_SN_U[[i]]$median - pSN_U[[i]]$true_eff ),2)
    mean_abs_bias_relSN_U[[i]] <- mean(abs_bias_relSN_U[[i]])
  }
}

################## mean bias of the mean biases of each data set for relative effects  ###################
mean_abs_bias_of_mean_abs_bias_relSN_U <- unlist(mean_abs_bias_relSN_U)
mean_abs_bias_of_mean_abs_bias_relSN_U  <- round(mean(mean_abs_bias_of_mean_abs_bias_relSN_U ),2)


######################## RELATIVE BIAS OF RELATIVE EFFECTS ##################
rel_bias_relSN_U <- list()
mean_rel_bias_relSN_U <- list()
for(i in 1:N.sim){
  for(j in 1:nrow(pSN_U[[i]])){
    rel_bias_relSN_U[[i]] <- round((delta_SN_U[[i]]$median - pSN_U[[i]]$true_eff ),2)
    mean_rel_bias_relSN_U[[i]] <- mean(rel_bias_relSN_U[[i]])
  }
}

################## mean bias of the mean biases of each data set for relative effects  ###################
mean_rel_bias_of_mean_rel_bias_relSN_U <- unlist(mean_rel_bias_relSN_U)
mean_rel_bias_of_mean_rel_bias_relSN_U  <- round(mean(mean_rel_bias_of_mean_rel_bias_relSN_U ),2)


########################### FOR THE OVERLAPPING ######################
dataSN1_U <- list()
overlapSN1_U <- list()
overlapSN_U <- list()
for(i in 1:N.sim){
  dataSN1_U[[i]] <- list(SN = delta_SN_U[[i]]$median , True = pSN_U[[i]]$true_eff)
  overlapSN1_U[[i]] <- overlap(dataSN1_U[[i]], type = "2")
  overlapSN_U[[i]] <- overlapSN1_U[[i]]$OV
}
overlapSN_U <- unlist(overlapSN_U)
avg_overlapSN_U <- round(mean(overlapSN_U),2)



listaSN_rel_eff_U <- cbind(mean_rel_bias_relSN_U,mean_rel_bias_of_mean_rel_bias_relSN_U,
                           mean_abs_bias_relSN_U,mean_abs_bias_of_mean_abs_bias_relSN_U) 

len_mu_SN_U = mean(data_SN_U$prec_muSN_U)
len_PI_SN_U = mean(data_SN_U$prec_predSN_U)

SNormal_U = cbind(data_SN_U, abs_biasSNjags.bin_U, rel_biasSNjags.bin_U, 
                  avg_abs_biasSNjags.bin_U, avg_rel_biasSNjags.bin_U, 
                  MSE_SNjags.bin_U, avg_MSE_SNjags.bin_U,
                  coverageSNjags.bin_U, coverSNjags.bin_U,
                  bias_tau2SNjags.bin_U, avg_bias_tau2SNjags.bin_U, pbias_tau2SNjags.bin_U,
                  MSE_SNjags.tau2.bin_U, avg_MSE_SNjags.tau2.bin_U, nMSE_SNjags.tau2.bin_U,
                  coverage_tau2SNjags.bin_U, cover_tau2SNjags.bin_U,
                  
                  abs_bias_tauSNjags.bin_U, avg_abs_bias_tauSNjags.bin_U, pabs_bias_tauSNjags.bin_U,
                  
                  bias_rel_tauSNjags.bin_U, avg_rel_bias_tauSNjags.bin_U, prel_bias_tauSNjags.bin_U,
                  MSE_SNjags.tau.bin_U, avg_MSE_SNjags.tau.bin_U, nMSE_SNjags.tau.bin_U,
                  coverage_tauSNjags.bin_U, cover_tauSNjags.bin_U,
                  len_mu_SN_U ,
                  len_PI_SN_U ,
                  RhatSN_out_U)

Overlap_SNmodel_U <- cbind(overlapSN_U, avg_overlapSN_U)

listaSN_U = cbind(data_SN_U$skewSN_U, data_SN_U$LBskewSN_U, data_SN_U$UBskewSN_U)


write.csv(SNormal_U, "C:\\Users\\Tianqi YU\\Desktop\\SIMULATION SCENARIOS FINAL\\NORMAL DATA\\SKEW NORMAL MODEL\\rev_U_ResultsSN_Normal_bin_data_k=14_mu=0_tau2=0.12.csv",row.names=FALSE )  
write.csv(deltaSN_U, "C:\\Users\\Tianqi YU\\Desktop\\SIMULATION SCENARIOS FINAL\\NORMAL DATA\\SKEW NORMAL MODEL\\rev_U_Relative_eff_ResultsSN_Normal_bin_data_k=14_mu=0_tau2=0.12.csv" ) 
write.csv(listaSN_rel_eff_U, "C:\\Users\\Tianqi YU\\Desktop\\SIMULATION SCENARIOS FINAL\\NORMAL DATA\\SKEW NORMAL MODEL\\rev_U_Bias_Relative_eff_ResultsSN_Normal_bin_data_k=14_mu=0_tau2=0.12.csv" ) 
write.csv(listaSN_U, "C:\\Users\\Tianqi YU\\Desktop\\SIMULATION SCENARIOS FINAL\\NORMAL DATA\\SKEW NORMAL MODEL\\rev_U_Skew_ResultsSN_Normal_bin_data_k=14_mu=0_tau2=0.12.csv" ) 
write.csv(Overlap_SNmodel_U, "C:\\Users\\Tianqi YU\\Desktop\\SIMULATION SCENARIOS FINAL\\NORMAL DATA\\SKEW NORMAL MODEL\\rev_U_Overlap_SNmodel_Normal_bin_data_k=14_mu=0_tau2=0.12.csv" ) 

################################# T-MODEL #############################
############# T(HN) #########
library(parallel)
library(cmdstanr)
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
library(bayesplot)
color_scheme_set("brightblue")
check_cmdstan_toolchain()

no_cores <- detectCores()-1
cl <- makeCluster(no_cores)
N.sim=1000
mu = 0
tau2 = 0.12
tau = sqrt(tau2)

X1=list()
for (i in 1:N.sim){ X1[[i]]=list("data"=p[[i]],"logOR"=mu,"tau_2"=tau2)}

################################################################
t_distr=function(X)
{
  delta = c()
  median <- c()
  LBmu <- c()
  UBmu <- c()
  pred <- c()
  LBpred <- c()
  UBpred <- c()
  tau <- c()
  LBtau <- c()
  UBtau <- c()
  tau2 <- c()
  LBtau2 <- c()
  UBtau2 <- c()
  Rhat_mut_distr <- c()
  Rhat_tau2t_distr <- c()
  Rhat_predt_distr <- c()
  
  filet_distr <- file.path(cmdstan_path(), "examples", "Stanmodels", "t_model_bin_HNpred.stan")
  modt_distr <- cmdstan_model(filet_distr)
  
  fitt_distr <- modt_distr$sample(
    data =list(ns = nrow(X$data),
               r1 = X$data$ci,
               r2 = X$data$ti,
               n1 = X$data$nci,
               n2 = X$data$nti) , 
    seed = 21015, 
    chains = 2, 
    parallel_chains = 2,
    refresh = 500,
    iter_warmup = 10000,
    iter_sampling = 50000,
    adapt_delta = 0.99
  )
  
  N.studies <- nrow(X$data)
  
  t_distr.res <- as.data.frame(fitt_distr$summary())
  
  d <- grepl("delta",t_distr.res$variable)
  
  delta <- t_distr.res[d,]
  
  
  dm <- grepl("mu",t_distr.res$variable)
  
  median <- t_distr.res[dm,]$median
  LBmu <- t_distr.res[dm,]$q5
  UBmu <- t_distr.res[dm,]$q95
  Rhat_mut_distr <-  t_distr.res[dm,]$rhat
  
  
  tau2 <- c()
  dtau2 <- grepl("tau_sqr",t_distr.res$variable)
  
  tau2 <- t_distr.res[dtau2,]$median
  LBtau2 <- t_distr.res[dtau2,]$q5
  UBtau2 <- t_distr.res[dtau2,]$q95
  Rhat_tau2t_distr <- t_distr.res[dtau2,]$rhat
  
  tau <- sqrt(t_distr.res[dtau2,]$median)
  LBtau <- sqrt(t_distr.res[dtau2,]$q5)
  UBtau <- sqrt(t_distr.res[dtau2,]$q95)
  
  dmx <- grepl("pred",t_distr.res$variable)
  
  pred <- t_distr.res[dm,]$median
  LBpred <- t_distr.res[dmx,]$q5
  UBpred <- t_distr.res[dmx,]$q95
  Rhat_predt_distr <-  t_distr.res[dmx,]$rhat
  
  t_distr.res1_HN<-data.frame(median=median, 
                              lowerCI=LBmu,
                              upperCI=UBmu,
                              pred=pred, 
                              LBpred=LBpred,
                              UBpred=UBpred,
                              tau=tau,
                              l_tau = LBtau,
                              u_tau = UBtau,
                              tau2=tau2,
                              l_tau2 = LBtau2,
                              u_tau2 = UBtau2,
                              Rhat_predt_distr = Rhat_predt_distr,
                              Rhat_mut_distr = Rhat_mut_distr,
                              Rhat_tau2t_distr = Rhat_tau2t_distr
  )
  
  
  return(list("res"= t_distr.res,
              "res1"= t_distr.res1_HN,
              "delta_t_distr"=delta,
              "Rhat_mut_distr" = Rhat_mut_distr,
              "Rhat_tau2t_distr" = Rhat_tau2t_distr,
              "Rhat_predt_distr" = Rhat_predt_distr
              
              
  ))
  
}

clusterExport(cl,"X1")

clusterExport(cl,"N.sim")

clusterExport(cl, "t_distr")

clusterEvalQ(cl, {library(cmdstanr)}) 

clusterEvalQ(cl, {library(SimDesign)})

l2t_HN <- parLapply(cl,1:N.sim, function(x) t_distr(X1[[x]]))

Rhat_mut_distr_HN = c()
Rhat_tau2t_distr_HN = c()
Rhat_predt_distr_HN = c()
mu_t_distr_HN = c()
LBmu_t_distr_HN = c()
UBmu_t_distr_HN = c()
pred_t_distr_HN = c()
LBpred_t_distr_HN = c()
UBpred_t_distr_HN = c()
tau_t_distr_HN = c()
LBtau_t_distr_HN = c()
UBtau_t_distr_HN = c()
tau2_t_distr_HN = c()
LBtau2_t_distr_HN = c()
UBtau2_t_distr_HN = c()
prec_mut_distr_HN = c()
prec_tau2t_distr_HN = c()
prec_predt_distr_HN = c()
prec_taut_distr_HN = c()
delta_t_distr_HN = list()

for(i in 1:N.sim){
  Rhat_mut_distr_HN = c(Rhat_mut_distr_HN, l2t_HN[[i]]$Rhat_mut_distr)
  Rhat_predt_distr_HN = c(Rhat_predt_distr_HN, l2t_HN[[i]]$Rhat_predt_distr)
  Rhat_tau2t_distr_HN = c(Rhat_tau2t_distr_HN, l2t_HN[[i]]$Rhat_tau2t_distr)
  
  mu_t_distr_HN = c(mu_t_distr_HN,l2t_HN[[i]]$res1$median)
  LBmu_t_distr_HN = c(LBmu_t_distr_HN,l2t_HN[[i]]$res1$lowerCI)
  UBmu_t_distr_HN = c(UBmu_t_distr_HN,l2t_HN[[i]]$res1$upperCI)
  prec_mut_distr_HN = c(prec_mut_distr_HN,(l2t_HN[[i]]$res1$upperCI - l2t_HN[[i]]$res1$lowerCI) )
  
  pred_t_distr_HN = c(pred_t_distr_HN,l2t_HN[[i]]$res1$pred)
  LBpred_t_distr_HN = c(LBpred_t_distr_HN,l2t_HN[[i]]$res1$LBpred)
  UBpred_t_distr_HN = c(UBpred_t_distr_HN,l2t_HN[[i]]$res1$UBpred)
  prec_predt_distr_HN = c(prec_predt_distr_HN,(l2t_HN[[i]]$res1$UBpred - l2t_HN[[i]]$res1$LBpred) )
  
  tau2_t_distr_HN = c(tau2_t_distr_HN,l2t_HN[[i]]$res1$tau2)
  LBtau2_t_distr_HN = c(LBtau2_t_distr_HN,l2t_HN[[i]]$res1$l_tau2)
  UBtau2_t_distr_HN = c(UBtau2_t_distr_HN,l2t_HN[[i]]$res1$u_tau2)
  prec_tau2t_distr_HN = c(prec_tau2t_distr_HN,(l2t_HN[[i]]$res1$u_tau2 - l2t_HN[[i]]$res1$l_tau2) )
  
  tau_t_distr_HN = c(tau_t_distr_HN,sqrt(l2t_HN[[i]]$res1$tau2))
  LBtau_t_distr_HN = c(LBtau_t_distr_HN,sqrt(l2t_HN[[i]]$res1$l_tau2))
  UBtau_t_distr_HN = c(UBtau_t_distr_HN,sqrt(l2t_HN[[i]]$res1$u_tau2))
  prec_taut_distr_HN = c(prec_taut_distr_HN,(sqrt(l2t_HN[[i]]$res1$u_tau2) - sqrt(l2t_HN[[i]]$res1$l_tau2)) )
  
  delta_t_distr_HN[[i]] = (l2t_HN[[i]]$delta_t_distr)
}

data_t_distr_HN = cbind.data.frame(mu_t_distr_HN, LBmu_t_distr_HN, UBmu_t_distr_HN, 
                                   pred_t_distr_HN, LBpred_t_distr_HN, UBpred_t_distr_HN,
                                   tau2_t_distr_HN, LBtau2_t_distr_HN, UBtau2_t_distr_HN, 
                                   tau_t_distr_HN, LBtau_t_distr_HN, UBtau_t_distr_HN, 
                                   Rhat_mut_distr_HN, Rhat_tau2t_distr_HN , prec_mut_distr_HN, prec_tau2t_distr_HN,
                                   Rhat_predt_distr_HN, prec_predt_distr_HN, prec_taut_distr_HN)

############REMOVE Rhats > 1.05 ##############
condition1_HN <- which(data_t_distr_HN$Rhat_mut_distr_HN > 1.05)
condition2_HN <- which(data_t_distr_HN$Rhat_tau2t_distr_HN > 1.05)

dist.condt_HN = c(condition1_HN,condition2_HN)
dist.condt_HN = unique(dist.condt_HN)

RhatT_out_HN <- round((length(dist.condt_HN)/N.sim),4)

if (length(dist.condt_HN)== 0) {
  data_t_distr_HN <- data_t_distr_HN
  delta_t_distr_HN <- delta_t_distr_HN
  N.sim = nrow(data_t_distr_HN)
  pT_HN <- p
} else {
  data_t_distr_HN = data_t_distr_HN[-dist.condt_HN, ]
  delta_t_distr_HN <- delta_t_distr_HN[-dist.condt_HN]
  N.sim = nrow(data_t_distr_HN)
  pT_HN <- p[-dist.condt_HN]
}

count_t_HN = list()
for(i in length(delta_t_distr_HN)){
  count_t_HN[[i]] = which(delta_t_distr_HN[[i]]$rhat > 1.05 )
  tell_met_HN = which(count_t_HN[[i]] != 0)
}

tell_met_HN

if(length(tell_met_HN) == 0){
  data_t_distr_HN <- data_t_distr_HN
  delta_t_distr_HN <- delta_t_distr_HN
  N.sim = nrow(data_t_distr_HN)
  pT_HN <- pT_HN
  RhatT_out_HN <- RhatT_out_HN
} else {
  data_t_distr_HN = data_t_distr_HN[-tell_met_HN, ]
  delta_t_distr_HN <- delta_t_distr_HN[-tell_met_HN]
  N.sim = nrow(data_t_distr_HN)
  pT_HN <- pT_HN[-tell_met_HN]
  RhatT_out_HN <- RhatT_out_HN + ((length(tell_met_HN))/N.sim)
}

RhatT_out_HN

########################################################
delta_tdistr_HN <- data.frame()
for (i in 1:N.sim) {
  temp_df_HN <- data.frame(simulation = i, delta_tdistr_HN = delta_t_distr_HN[[i]])
  delta_tdistr_HN <- rbind(delta_tdistr_HN, temp_df_HN)
}

########## AVERAGE ABSOLUTE BIAS OF mu ##########
abs_biast_distrjags.bin_HN <- c()
for (i in 1:N.sim){
  abs_biast_distrjags.bin_HN[i] <- abs( data_t_distr_HN$mu_t_distr_HN[i] - mu) 
}
avg_abs_biast_distrjags.bin_HN <- round(mean(abs_biast_distrjags.bin_HN), 2)

########## AVERAGE RELATIVE  BIAS OF mu ##########
rel_biast_distrjags.bin_HN <- c()
for (i in 1:N.sim){
  rel_biast_distrjags.bin_HN[i] <- ( data_t_distr_HN$mu_t_distr_HN[i] - mu) 
}
avg_rel_biast_distrjags.bin_HN <- round(mean(rel_biast_distrjags.bin_HN), 2)

############### MSE OF mu  ##################
MSE_t_distrjags.bin_HN <- c()
for (i in 1:N.sim){
  MSE_t_distrjags.bin_HN[i] <- RMSE(estimate =data_t_distr_HN$mu_t_distr_HN[i],
                                    parameter = mu,
                                    type = "RMSE",
                                    MSE = TRUE,           
                                    percent = FALSE,
                                    unname = FALSE)
}
avg_MSE_t_distrjags.bin_HN <-  round(mean(MSE_t_distrjags.bin_HN),4)

########## COVERAGE OF mu ##########
interval_contains_true_mean <- function(i) { 
  ( mu >= data_t_distr_HN$LBmu_t_distr_HN[i]) && ( mu <=  data_t_distr_HN$UBmu_t_distr_HN[i])
}
interval_contains_true_mean(i)

coveraget_distrjags.bin_HN <- c()
for (i in 1:N.sim){0
  if (interval_contains_true_mean(i) == TRUE){
    coveraget_distrjags.bin_HN[i] = 1
  }
  else {
    coveraget_distrjags.bin_HN[i]=0
  }
  print(coveraget_distrjags.bin_HN[i])
}
coveraget_distrjags.bin_HN
covert_distrjags.bin_HN <- round(mean(coveraget_distrjags.bin_HN)*100, 2)

################
##### AVERAGE BIAS OF tau2 #########
bias_tau2t_distrjags.bin_HN <- c()
for (i in 1:N.sim){
  bias_tau2t_distrjags.bin_HN[i] <- abs(data_t_distr_HN$tau2_t_distr_HN[i] - tau2)
}
avg_bias_tau2t_distrjags.bin_HN <- round(mean(bias_tau2t_distrjags.bin_HN),2)

pbias_tau2t_distrjags.bin_HN = round(mean(bias_tau2t_distrjags.bin_HN/tau2), 2)
########### MSE OF tau2 ################
MSE_t_distrjags.tau2.bin_HN <- c()
for (i in 1:N.sim){
  MSE_t_distrjags.tau2.bin_HN[i] <- RMSE(estimate = data_t_distr_HN$tau2_t_distr_HN[i],
                                         parameter = tau2,
                                         type = "RMSE",
                                         MSE = TRUE     
  )
}
avg_MSE_t_distrjags.tau2.bin_HN <- round(mean(MSE_t_distrjags.tau2.bin_HN),4)

nMSE_t_distrjags.tau2.bin_HN  <- round(mean(MSE_t_distrjags.tau2.bin_HN / (tau2)^2),4)

############# COVERAGE OF tau2 ###################
interval_contains_true_mean <- function(i) { 
  ( tau2 >= data_t_distr_HN$LBtau2_t_distr_HN[i]) && ( tau2 <=  data_t_distr_HN$UBtau2_t_distr_HN[i])
}
interval_contains_true_mean(i)

coverage_tau2t_distrjags.bin_HN <- c()
for (i in 1:N.sim){
  if (interval_contains_true_mean(i) == TRUE){
    coverage_tau2t_distrjags.bin_HN[i] = 1
  }
  else  {
    coverage_tau2t_distrjags.bin_HN[i] = 0
  }
  print(coverage_tau2t_distrjags.bin_HN[i])
}
coverage_tau2t_distrjags.bin_HN
cover_tau2t_distrjags.bin_HN <- round(mean(coverage_tau2t_distrjags.bin_HN)*100,2)


################
##### AVERAGE RELATIVE BIAS OF tau #########
rel_bias_taut_distrjags.bin_HN <- c()
for (i in 1:N.sim){
  rel_bias_taut_distrjags.bin_HN[i] <- (data_t_distr_HN$tau_t_distr_HN[i] - tau)
}
avg_rel_bias_taut_distrjags.bin_HN <- round(mean(rel_bias_taut_distrjags.bin_HN),2)

prel_bias_taut_distrjags.bin_HN = round(mean(rel_bias_taut_distrjags.bin_HN/tau)*100, 2)

##### AVERAGE ABSOLUTE BIAS OF tau #########
abs_bias_taut_distrjags.bin_HN <- c()
for (i in 1:N.sim){
  abs_bias_taut_distrjags.bin_HN[i] <- abs(data_t_distr_HN$tau_t_distr_HN[i] - tau)
}
avg_abs_bias_taut_distrjags.bin_HN <- round(mean(abs_bias_taut_distrjags.bin_HN),2)

pabs_bias_taut_distrjags.bin_HN = round(mean(abs_bias_taut_distrjags.bin_HN/tau)*100, 2)


########### MSE OF tau ################
MSE_t_distrjags.tau.bin_HN <- c()
for (i in 1:N.sim){
  MSE_t_distrjags.tau.bin_HN[i] <- RMSE(estimate = data_t_distr_HN$tau_t_distr_HN[i],
                                        parameter = tau,
                                        type = "RMSE",
                                        MSE = TRUE     
  )
}
avg_MSE_t_distrjags.tau.bin_HN <- round(mean(MSE_t_distrjags.tau.bin_HN),4)

nMSE_t_distrjags.tau.bin_HN  <- round(mean(MSE_t_distrjags.tau.bin_HN / tau2),4)

############# COVERAGE OF tau ###################
interval_contains_true_mean <- function(i) { 
  ( tau >= data_t_distr_HN$LBtau_t_distr_HN[i]) && ( tau <=  data_t_distr_HN$UBtau_t_distr_HN[i])
}
interval_contains_true_mean(i)

coverage_taut_distrjags.bin_HN <- c()
for (i in 1:N.sim){
  if (interval_contains_true_mean(i) == TRUE){
    coverage_taut_distrjags.bin_HN[i] = 1
  }
  else  {
    coverage_taut_distrjags.bin_HN[i] = 0
  }
  print(coverage_taut_distrjags.bin_HN[i])
}
coverage_taut_distrjags.bin_HN
cover_taut_distrjags.bin_HN <- round(mean(coverage_taut_distrjags.bin_HN)*100,2)




######################## ABSOLUTE BIAS OF RELATIVE EFFECTS ##################
abs_bias_relt_distr_HN <- list()
mean_abs_bias_relt_distr_HN <- list()
for(i in 1:N.sim){
  for(j in 1:nrow(pT_HN[[i]])){
    abs_bias_relt_distr_HN[[i]] <- round(abs(delta_t_distr_HN[[i]]$median - pT_HN[[i]]$true_eff),2)
    mean_abs_bias_relt_distr_HN[[i]] <- mean(abs_bias_relt_distr_HN[[i]])
  }
}
#mean_bias_relt_distr 

################## mean bias of the mean biases of each data set for relative effects  ###################
mean_abs_bias_of_mean_abs_bias_relt_distr_HN <- unlist(mean_abs_bias_relt_distr_HN)
mean_abs_bias_of_mean_abs_bias_relt_distr_HN  <- round(mean(mean_abs_bias_of_mean_abs_bias_relt_distr_HN),2)

####################### RELATIVE BIAS OF RELATIVE EFFECTS ##################
rel_bias_relt_distr_HN <- list()
mean_rel_bias_relt_distr_HN <- list()
for(i in 1:N.sim){
  for(j in 1:nrow(pT_HN[[i]])){
    rel_bias_relt_distr_HN[[i]] <- round((delta_t_distr_HN[[i]]$median - pT_HN[[i]]$true_eff),2)
    mean_rel_bias_relt_distr_HN[[i]] <- mean(rel_bias_relt_distr_HN[[i]])
  }
}
#mean_bias_relt_distr 

################## mean bias of the mean biases of each data set for relative effects  ###################
mean_rel_bias_of_mean_rel_bias_relt_distr_HN <- unlist(mean_rel_bias_relt_distr_HN)
mean_rel_bias_of_mean_rel_bias_relt_distr_HN  <- round(mean(mean_rel_bias_of_mean_rel_bias_relt_distr_HN),2)

########################### FOR THE OVERLAPPING ######################
datat_HN <- list()
overlapt1_HN <- list()
overlapt_HN <- list()
for(i in 1:N.sim){
  datat_HN[[i]] <- list(t_distr = delta_t_distr_HN[[i]]$median , True = pT_HN[[i]]$true_eff)
  overlapt1_HN[[i]] <- overlap(datat_HN[[i]], type = "2")
  overlapt_HN[[i]] <- overlapt1_HN[[i]]$OV
}
overlapt_HN <- unlist(overlapt_HN)
avg_overlapt_HN <- round(mean(overlapt_HN),2)



listat_distr_rel_eff_HN <- cbind(mean_abs_bias_relt_distr_HN,mean_abs_bias_of_mean_abs_bias_relt_distr_HN,
                                 mean_rel_bias_relt_distr_HN,mean_rel_bias_of_mean_rel_bias_relt_distr_HN) 

len_mu_t_distr_HN = mean(data_t_distr_HN$prec_mut_distr_HN)
len_PI_t_distr_HN = mean(data_t_distr_HN$prec_predt_distr_HN)

t_distr_HN = cbind(data_t_distr_HN, abs_biast_distrjags.bin_HN,rel_biast_distrjags.bin_HN,
                   avg_abs_biast_distrjags.bin_HN, avg_rel_biast_distrjags.bin_HN, 
                   MSE_t_distrjags.bin_HN, avg_MSE_t_distrjags.bin_HN, 
                   coveraget_distrjags.bin_HN, covert_distrjags.bin_HN,bias_tau2t_distrjags.bin_HN, avg_bias_tau2t_distrjags.bin_HN, pbias_tau2t_distrjags.bin_HN,
                   MSE_t_distrjags.tau2.bin_HN, avg_MSE_t_distrjags.tau2.bin_HN, nMSE_t_distrjags.tau2.bin_HN, coverage_tau2t_distrjags.bin_HN,
                   cover_tau2t_distrjags.bin_HN,
                   rel_bias_taut_distrjags.bin_HN, avg_rel_bias_taut_distrjags.bin_HN, prel_bias_taut_distrjags.bin_HN,
                   abs_bias_taut_distrjags.bin_HN, avg_abs_bias_taut_distrjags.bin_HN, pabs_bias_taut_distrjags.bin_HN,
                   MSE_t_distrjags.tau.bin_HN, avg_MSE_t_distrjags.tau.bin_HN, nMSE_t_distrjags.tau.bin_HN, coverage_taut_distrjags.bin_HN,
                   cover_taut_distrjags.bin_HN,
                   len_mu_t_distr_HN,
                   len_PI_t_distr_HN,
                   RhatT_out_HN)

Overlap_tdistr_HN = cbind(overlapt_HN, avg_overlapt_HN)

write.csv(t_distr_HN, "C:\\Users\\Lela Panag\\Desktop\\SIMULATION SCENARIOS FINAL\\NORMAL DATA\\T MODEL\\HN\\rev_HN_Results_t_distr_Normal_bin_data_k=14_mu=0_tau2=0.12.csv",row.names=FALSE )  
write.csv(delta_tdistr_HN, "C:\\Users\\Lela Panag\\Desktop\\SIMULATION SCENARIOS FINAL\\NORMAL DATA\\T MODEL\\HN\\rev_HN_Relative_eff_Results_t_distr_Normal_bin_data_k=14_mu=0_tau2=0.12.csv" ) 
write.csv(listat_distr_rel_eff_HN, "C:\\Users\\Lela Panag\\Desktop\\SIMULATION SCENARIOS FINAL\\NORMAL DATA\\T MODEL\\HN\\rev_HN_Bias_Relative_eff_Results_t_distr_Normal_bin_data_k=14_mu=0_tau2=0.12.csv" ) 
write.csv(Overlap_tdistr_HN, "C:\\Users\\Lela Panag\\Desktop\\SIMULATION SCENARIOS FINAL\\NORMAL DATA\\T MODEL\\HN\\rev_HN_Overlap_tdistr_Normal_bin_data_k=14_mu=0_tau2=0.12.csv" ) 

################## T(Unif) #########
################################ T-MODEL #############################
library(parallel)
library(cmdstanr)
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
library(bayesplot)
color_scheme_set("brightblue")
check_cmdstan_toolchain()

no_cores <- detectCores()-1
cl <- makeCluster(no_cores)
N.sim=1000
mu = 0
tau2 = 0.12
tau = sqrt(tau2)

X1=list()
for (i in 1:N.sim){ X1[[i]]=list("data"=p[[i]],"logOR"=mu,"tau_2"=tau2)}

################################################################
t_distr=function(X)
{
  delta = c()
  median <- c()
  LBmu <- c()
  UBmu <- c()
  pred <- c()
  LBpred <- c()
  UBpred <- c()
  tau <- c()
  LBtau <- c()
  UBtau <- c()
  tau2 <- c()
  LBtau2 <- c()
  UBtau2 <- c()
  Rhat_mut_distr <- c()
  Rhat_tau2t_distr <- c()
  Rhat_skewt_distr <- c()
  Rhat_predt_distr <- c()
  
  filet_distr <- file.path(cmdstan_path(), "examples", "Stanmodels", "t_model_binUpred.stan")
  modt_distr <- cmdstan_model(filet_distr)
  
  fitt_distr <- modt_distr$sample(
    data =list(ns = nrow(X$data),
               r1 = X$data$ci,
               r2 = X$data$ti,
               n1 = X$data$nci,
               n2 = X$data$nti) , 
    seed = 2108, 
    chains = 2, 
    parallel_chains = 2,
    refresh = 500,
    iter_warmup = 10000,
    iter_sampling = 50000,
    adapt_delta = 0.99
  )
  
  N.studies <- nrow(X$data)
  
  t_distr.res <- as.data.frame(fitt_distr$summary())
  
  d <- grepl("delta",t_distr.res$variable)
  
  delta <- t_distr.res[d,]
  
  
  dm <- grepl("mu",t_distr.res$variable)
  
  median <- t_distr.res[dm,]$median
  LBmu <- t_distr.res[dm,]$q5
  UBmu <- t_distr.res[dm,]$q95
  Rhat_mut_distr <-  t_distr.res[dm,]$rhat
  
  
  tau2 <- c()
  dtau2 <- grepl("tau_sqr",t_distr.res$variable)
  
  tau2 <- t_distr.res[dtau2,]$median
  LBtau2 <- t_distr.res[dtau2,]$q5
  UBtau2 <- t_distr.res[dtau2,]$q95
  Rhat_tau2t_distr <- t_distr.res[dtau2,]$rhat
  
  tau <- sqrt(t_distr.res[dtau2,]$median)
  LBtau <- sqrt(t_distr.res[dtau2,]$q5)
  UBtau <- sqrt(t_distr.res[dtau2,]$q95)
  
  dmx <- grepl("pred",t_distr.res$variable)
  
  pred <- t_distr.res[dm,]$median
  LBpred <- t_distr.res[dmx,]$q5
  UBpred <- t_distr.res[dmx,]$q95
  Rhat_predt_distr <-  t_distr.res[dmx,]$rhat
  
  t_distr.res1_U<-data.frame(median=median, 
                             lowerCI=LBmu,
                             upperCI=UBmu,
                             pred=pred, 
                             LBpred=LBpred,
                             UBpred=UBpred,
                             tau=tau,
                             l_tau = LBtau,
                             u_tau = UBtau,
                             tau2=tau2,
                             tau2=tau2,
                             l_tau2 = LBtau2,
                             u_tau2 = UBtau2,
                             Rhat_mut_distr = Rhat_mut_distr,
                             Rhat_predt_distr = Rhat_predt_distr,
                             Rhat_tau2t_distr = Rhat_tau2t_distr
  )
  
  
  return(list("res"= t_distr.res,
              "res1"= t_distr.res1_U,
              "delta_t_distr"=delta,
              "Rhat_mut_distr" = Rhat_mut_distr,
              "Rhat_tau2t_distr" = Rhat_tau2t_distr,
              "Rhat_predt_distr" = Rhat_predt_distr
              
  ))
  
}

clusterExport(cl,"X1")

clusterExport(cl,"N.sim")

clusterExport(cl, "t_distr")

clusterEvalQ(cl, {library(cmdstanr)}) 

clusterEvalQ(cl, {library(SimDesign)})

l2t_U <- parLapply(cl,1:N.sim, function(x) t_distr(X1[[x]]))


Rhat_mut_distr_U = c()
Rhat_tau2t_distr_U = c()
Rhat_predt_distr_U = c()
mu_t_distr_U = c()
LBmu_t_distr_U = c()
UBmu_t_distr_U = c()
pred_t_distr_U = c()
LBpred_t_distr_U = c()
UBpred_t_distr_U = c()
tau_t_distr_U = c()
LBtau_t_distr_U = c()
UBtau_t_distr_U = c()
tau2_t_distr_U = c()
LBtau2_t_distr_U = c()
UBtau2_t_distr_U = c()
prec_mut_distr_U = c()
prec_tau2t_distr_U = c()
prec_predt_distr_U = c()
prec_taut_distr_U = c()
delta_t_distr_U = list()

for(i in 1:N.sim){
  Rhat_mut_distr_U = c(Rhat_mut_distr_U, l2t_U[[i]]$Rhat_mut_distr)
  Rhat_tau2t_distr_U = c(Rhat_tau2t_distr_U, l2t_U[[i]]$Rhat_tau2t_distr)
  Rhat_predt_distr_U = c(Rhat_predt_distr_U, l2t_U[[i]]$Rhat_predt_distr)
  
  mu_t_distr_U = c(mu_t_distr_U,l2t_U[[i]]$res1$median)
  LBmu_t_distr_U = c(LBmu_t_distr_U,l2t_U[[i]]$res1$lowerCI)
  UBmu_t_distr_U = c(UBmu_t_distr_U,l2t_U[[i]]$res1$upperCI)
  prec_mut_distr_U = c(prec_mut_distr_U,(l2t_U[[i]]$res1$upperCI - l2t_U[[i]]$res1$lowerCI) )
  
  pred_t_distr_U = c(pred_t_distr_U,l2t_U[[i]]$res1$pred)
  LBpred_t_distr_U = c(LBpred_t_distr_U,l2t_U[[i]]$res1$LBpred)
  UBpred_t_distr_U = c(UBpred_t_distr_U,l2t_U[[i]]$res1$UBpred)
  prec_predt_distr_U = c(prec_predt_distr_U,(l2t_U[[i]]$res1$UBpred - l2t_U[[i]]$res1$LBpred) )
  
  tau2_t_distr_U = c(tau2_t_distr_U,l2t_U[[i]]$res1$tau2)
  LBtau2_t_distr_U = c(LBtau2_t_distr_U,l2t_U[[i]]$res1$l_tau2)
  UBtau2_t_distr_U = c(UBtau2_t_distr_U,l2t_U[[i]]$res1$u_tau2)
  prec_tau2t_distr_U = c(prec_tau2t_distr_U,(l2t_U[[i]]$res1$u_tau2 - l2t_U[[i]]$res1$l_tau2) )
  
  tau_t_distr_U = c(tau_t_distr_U,sqrt(l2t_U[[i]]$res1$tau2))
  LBtau_t_distr_U = c(LBtau_t_distr_U,sqrt(l2t_U[[i]]$res1$l_tau2))
  UBtau_t_distr_U = c(UBtau_t_distr_U,sqrt(l2t_U[[i]]$res1$u_tau2))
  prec_taut_distr_U = c(prec_taut_distr_U,(sqrt(l2t_U[[i]]$res1$u_tau2) - sqrt(l2t_U[[i]]$res1$l_tau2)) )
  
  delta_t_distr_U[[i]] = (l2t_U[[i]]$delta_t_distr)
}

data_t_distr_U = cbind.data.frame(mu_t_distr_U, LBmu_t_distr_U, UBmu_t_distr_U, 
                                  pred_t_distr_U, LBpred_t_distr_U, UBpred_t_distr_U,
                                  tau2_t_distr_U, LBtau2_t_distr_U, UBtau2_t_distr_U, 
                                  tau_t_distr_U, LBtau_t_distr_U, UBtau_t_distr_U, 
                                  Rhat_mut_distr_U, Rhat_tau2t_distr_U , prec_mut_distr_U, prec_tau2t_distr_U,
                                  Rhat_predt_distr_U, prec_predt_distr_U, prec_taut_distr_U)

############REMOVE Rhats > 1.05 ##############
condition1_U <- which(data_t_distr_U$Rhat_mut_distr_U > 1.05)
condition2_U <- which(data_t_distr_U$Rhat_tau2t_distr_U > 1.05)

dist.condt_U = c(condition1_U,condition2_U)
dist.condt_U = unique(dist.condt_U)

RhatT_out_U <- round((length(dist.condt_U)/N.sim),4)

if (length(dist.condt_U)== 0) {
  data_t_distr_U <- data_t_distr_U
  delta_t_distr_U <- delta_t_distr_U
  N.sim = nrow(data_t_distr_U)
  pT_U <- p
} else {
  data_t_distr_U = data_t_distr_U[-dist.condt_U, ]
  delta_t_distr_U <- delta_t_distr_U[-dist.condt_U]
  N.sim = nrow(data_t_distr_U)
  pT_U <- p[-dist.condt_U]
}

count_t_U = list()
for(i in length(delta_t_distr_U)){
  count_t_U[[i]] = which(delta_t_distr_U[[i]]$rhat > 1.05 )
  tell_met_U = which(count_t_U[[i]] != 0)
}

tell_met_U

if(length(tell_met_U) == 0){
  data_t_distr_U <- data_t_distr_U
  delta_t_distr_U <- delta_t_distr_U
  N.sim = nrow(data_t_distr_U)
  pT_U <- pT_U
  RhatT_out_U <- RhatT_out_U
} else {
  data_t_distr_U = data_t_distr_U[-tell_met_U, ]
  delta_t_distr_U <- delta_t_distr_U[-tell_met_U]
  N.sim = nrow(data_t_distr_U)
  pT_U <- pT_U[-tell_met_U]
  RhatT_out_U <- RhatT_out_U + ((length(tell_met_U))/N.sim)
}

RhatT_out_U

########################################################
delta_tdistr_U <- data.frame()
for (i in 1:N.sim) {
  temp_df_U <- data.frame(simulation = i, delta_tdistr_U = delta_t_distr_U[[i]])
  delta_tdistr_U <- rbind(delta_tdistr_U, temp_df_U)
}

########## AVERAGE ABSOLUTE BIAS OF mu ##########
abs_biast_distrjags.bin_U <- c()
for (i in 1:N.sim){
  abs_biast_distrjags.bin_U[i] <- abs( data_t_distr_U$mu_t_distr_U[i] - mu) 
}
avg_abs_biast_distrjags.bin_U <- round(mean(abs_biast_distrjags.bin_U), 2)

########## AVERAGE RELATIVE BIAS OF mu ##########
rel_biast_distrjags.bin_U <- c()
for (i in 1:N.sim){
  rel_biast_distrjags.bin_U[i] <- ( data_t_distr_U$mu_t_distr_U[i] - mu) 
}
avg_rel_biast_distrjags.bin_U <- round(mean(rel_biast_distrjags.bin_U), 2)

############### MSE OF mu  ##################
MSE_t_distrjags.bin_U <- c()
for (i in 1:N.sim){
  MSE_t_distrjags.bin_U[i] <- RMSE(estimate =data_t_distr_U$mu_t_distr_U[i],
                                   parameter = mu,
                                   type = "RMSE",
                                   MSE = TRUE,           
                                   percent = FALSE,
                                   unname = FALSE)
}
avg_MSE_t_distrjags.bin_U <-  round(mean(MSE_t_distrjags.bin_U),4)

########## COVERAGE OF mu ##########
interval_contains_true_mean <- function(i) { 
  ( mu >= data_t_distr_U$LBmu_t_distr_U[i]) && ( mu <=  data_t_distr_U$UBmu_t_distr_U[i])
}
interval_contains_true_mean(i)

coveraget_distrjags.bin_U <- c()
for (i in 1:N.sim){0
  if (interval_contains_true_mean(i) == TRUE){
    coveraget_distrjags.bin_U[i] = 1
  }
  else {
    coveraget_distrjags.bin_U[i]=0
  }
  print(coveraget_distrjags.bin_U[i])
}
coveraget_distrjags.bin_U
covert_distrjags.bin_U <- round(mean(coveraget_distrjags.bin_U)*100, 2)

################
##### AVERAGE BIAS OF tau2 #########
bias_tau2t_distrjags.bin_U <- c()
for (i in 1:N.sim){
  bias_tau2t_distrjags.bin_U[i] <- abs(data_t_distr_U$tau2_t_distr_U[i] - tau2)
}
avg_bias_tau2t_distrjags.bin_U <- round(mean(bias_tau2t_distrjags.bin_U),2)

pbias_tau2t_distrjags.bin_U = round(mean(bias_tau2t_distrjags.bin_U/tau2), 2)
########### MSE OF tau2 ################
MSE_t_distrjags.tau2.bin_U <- c()
for (i in 1:N.sim){
  MSE_t_distrjags.tau2.bin_U[i] <- RMSE(estimate = data_t_distr_U$tau2_t_distr_U[i],
                                        parameter = tau2,
                                        type = "RMSE",
                                        MSE = TRUE     
  )
}
avg_MSE_t_distrjags.tau2.bin_U <- round(mean(MSE_t_distrjags.tau2.bin_U),4)

nMSE_t_distrjags.tau2.bin_U  <- round(mean(MSE_t_distrjags.tau2.bin_U / (tau2)^2),4)

############# COVERAGE OF tau2 ###################
interval_contains_true_mean <- function(i) { 
  ( tau2 >= data_t_distr_U$LBtau2_t_distr_U[i]) && ( tau2 <=  data_t_distr_U$UBtau2_t_distr_U[i])
}
interval_contains_true_mean(i)

coverage_tau2t_distrjags.bin_U <- c()
for (i in 1:N.sim){
  if (interval_contains_true_mean(i) == TRUE){
    coverage_tau2t_distrjags.bin_U[i] = 1
  }
  else  {
    coverage_tau2t_distrjags.bin_U[i] = 0
  }
  print(coverage_tau2t_distrjags.bin_U[i])
}
coverage_tau2t_distrjags.bin_U
cover_tau2t_distrjags.bin_U <- round(mean(coverage_tau2t_distrjags.bin_U)*100,2)



################
##### AVERAGE RELATIVE BIAS OF tau #########
rel_bias_taut_distrjags.bin_U <- c()
for (i in 1:N.sim){
  rel_bias_taut_distrjags.bin_U[i] <- (data_t_distr_U$tau_t_distr_U[i] - tau)
}
avg_rel_bias_taut_distrjags.bin_U <- round(mean(rel_bias_taut_distrjags.bin_U),2)

prel_bias_taut_distrjags.bin_U = round(mean(rel_bias_taut_distrjags.bin_U/tau)*100, 2)

##### AVERAGE ABSOLUTE BIAS OF tau #########
abs_bias_taut_distrjags.bin_U <- c()
for (i in 1:N.sim){
  abs_bias_taut_distrjags.bin_U[i] <- abs(data_t_distr_U$tau_t_distr_U[i] - tau)
}
avg_abs_bias_taut_distrjags.bin_U <- round(mean(abs_bias_taut_distrjags.bin_U),2)

pabs_bias_taut_distrjags.bin_U = round(mean(abs_bias_taut_distrjags.bin_U/tau)*100, 2)


########### MSE OF tau ################
MSE_t_distrjags.tau.bin_U <- c()
for (i in 1:N.sim){
  MSE_t_distrjags.tau.bin_U[i] <- RMSE(estimate = data_t_distr_U$tau_t_distr_U[i],
                                       parameter = tau,
                                       type = "RMSE",
                                       MSE = TRUE     
  )
}
avg_MSE_t_distrjags.tau.bin_U <- round(mean(MSE_t_distrjags.tau.bin_U),4)

nMSE_t_distrjags.tau.bin_U  <- round(mean(MSE_t_distrjags.tau.bin_U / tau2),4)

############# COVERAGE OF tau ###################
interval_contains_true_mean <- function(i) { 
  ( tau >= data_t_distr_U$LBtau_t_distr_U[i]) && ( tau <=  data_t_distr_U$UBtau_t_distr_U[i])
}
interval_contains_true_mean(i)

coverage_taut_distrjags.bin_U <- c()
for (i in 1:N.sim){
  if (interval_contains_true_mean(i) == TRUE){
    coverage_taut_distrjags.bin_U[i] = 1
  }
  else  {
    coverage_taut_distrjags.bin_U[i] = 0
  }
  print(coverage_taut_distrjags.bin_U[i])
}
coverage_taut_distrjags.bin_U
cover_taut_distrjags.bin_U <- round(mean(coverage_taut_distrjags.bin_U)*100,2)


######################## ABSOLUTE BIAS OF RELATIVE EFFECTS ##################
abs_bias_relt_distr_U <- list()
mean_abs_bias_relt_distr_U <- list()
for(i in 1:N.sim){
  for(j in 1:nrow(pT_U[[i]])){
    abs_bias_relt_distr_U[[i]] <- round(abs(delta_t_distr_U[[i]]$median - pT_U[[i]]$true_eff),2)
    mean_abs_bias_relt_distr_U[[i]] <- mean(abs_bias_relt_distr_U[[i]])
  }
}
#mean_bias_relt_distr 

################## mean bias of the mean biases of each data set for relative effects  ###################
mean_abs_bias_of_mean_abs_bias_relt_distr_U <- unlist(mean_abs_bias_relt_distr_U)
mean_abs_bias_of_mean_abs_bias_relt_distr_U  <- round(mean(mean_abs_bias_of_mean_abs_bias_relt_distr_U),2)

######################## RELATIVE BIAS OF RELATIVE EFFECTS ##################
rel_bias_relt_distr_U <- list()
mean_rel_bias_relt_distr_U <- list()
for(i in 1:N.sim){
  for(j in 1:nrow(pT_U[[i]])){
    rel_bias_relt_distr_U[[i]] <- round((delta_t_distr_U[[i]]$median - pT_U[[i]]$true_eff),2)
    mean_rel_bias_relt_distr_U[[i]] <- mean(rel_bias_relt_distr_U[[i]])
  }
}
#mean_bias_relt_distr 

################## mean bias of the mean biases of each data set for relative effects  ###################
mean_rel_bias_of_mean_rel_bias_relt_distr_U <- unlist(mean_rel_bias_relt_distr_U)
mean_rel_bias_of_mean_rel_bias_relt_distr_U  <- round(mean(mean_rel_bias_of_mean_rel_bias_relt_distr_U),2)

########################### FOR THE OVERLAPPING ######################
datat_U <- list()
overlapt1_U <- list()
overlapt_U <- list()
for(i in 1:N.sim){
  datat_U[[i]] <- list(t_distr = delta_t_distr_U[[i]]$median , True = pT_U[[i]]$true_eff)
  overlapt1_U[[i]] <- overlap(datat_U[[i]], type = "2")
  overlapt_U[[i]] <- overlapt1_U[[i]]$OV
}
overlapt_U <- unlist(overlapt_U)
avg_overlapt_U <- round(mean(overlapt_U),2)



listat_distr_rel_eff_U <- cbind(mean_abs_bias_relt_distr_U,mean_abs_bias_of_mean_abs_bias_relt_distr_U,
                                mean_rel_bias_relt_distr_U,mean_rel_bias_of_mean_rel_bias_relt_distr_U) 

len_mu_t_distr_U = mean(data_t_distr_U$prec_mut_distr_U)
len_PI_t_distr_U = mean(data_t_distr_U$prec_predt_distr_U)

t_distr_U = cbind(data_t_distr_U, abs_biast_distrjags.bin_U,rel_biast_distrjags.bin_U,
                  avg_abs_biast_distrjags.bin_U,avg_rel_biast_distrjags.bin_U,
                  MSE_t_distrjags.bin_U, avg_MSE_t_distrjags.bin_U, 
                  coveraget_distrjags.bin_U, covert_distrjags.bin_U,bias_tau2t_distrjags.bin_U, avg_bias_tau2t_distrjags.bin_U, pbias_tau2t_distrjags.bin_U,
                  MSE_t_distrjags.tau2.bin_U, avg_MSE_t_distrjags.tau2.bin_U, nMSE_t_distrjags.tau2.bin_U, coverage_tau2t_distrjags.bin_U,
                  cover_tau2t_distrjags.bin_U,
                  rel_bias_taut_distrjags.bin_U, avg_rel_bias_taut_distrjags.bin_U, prel_bias_taut_distrjags.bin_U,
                  abs_bias_taut_distrjags.bin_U, avg_abs_bias_taut_distrjags.bin_U, pabs_bias_taut_distrjags.bin_U,
                  
                  MSE_t_distrjags.tau.bin_U, avg_MSE_t_distrjags.tau.bin_U, nMSE_t_distrjags.tau.bin_U, coverage_taut_distrjags.bin_U,
                  cover_taut_distrjags.bin_U,
                  len_mu_t_distr_U ,
                  len_PI_t_distr_U,
                  RhatT_out_U)

Overlap_tdistr_U = cbind(overlapt_U, avg_overlapt_U)


write.csv(t_distr_U, "C:\\Users\\Lela Panag\\Desktop\\SIMULATION SCENARIOS FINAL\\NORMAL DATA\\T MODEL\\U\\rev_U_Results_t_distr_Normal_bin_data_k=14_mu=0_tau2=0.12.csv",row.names=FALSE )  
write.csv(delta_tdistr_U, "C:\\Users\\Lela Panag\\Desktop\\SIMULATION SCENARIOS FINAL\\NORMAL DATA\\T MODEL\U\\\rev_U_Relative_eff_Results_t_distr_Normal_bin_data_k=14_mu=0_tau2=0.12.csv" ) 
write.csv(listat_distr_rel_eff_U, "C:\\Users\\Lela Panag\\Desktop\\SIMULATION SCENARIOS FINAL\\NORMAL DATA\\T MODEL\\U\\rev_U_Bias_Relative_eff_Results_t_distr_Normal_bin_data_k=14_mu=0_tau2=0.12.csv" ) 
write.csv(Overlap_tdistr_U, "C:\\Users\\Lela Panag\\Desktop\\SIMULATION SCENARIOS FINAL\\NORMAL DATA\\T MODEL\\U\\rev_U_Overlap_tdistr_Normal_bin_data_k=14_mu=0_tau2=0.12.csv" ) 


################## FREQUENTIST MODELS ####################
################## NORMAL ###########
################## FREQUENTIST MODELS ####################
################## Binomial-Normal ###########

mu = 0
tau2 = 0.12
N.sim = 1000
tau = sqrt(tau2)

######## NORMAL RANDOM EFFECTS WITH METAFOR PACKAGE ########

freq_norm1bin <- list()
freq_normbin <- list()
for(i in 1:N.sim){
  set.seed(21013)
  freq_norm1bin[[i]] <- rma.glmm(measure = "OR",
                                 ai = p[[i]]$ti, bi = p[[i]]$nti - p[[i]]$ti,
                                 ci = p[[i]]$ci, di = p[[i]]$nci - p[[i]]$ci,
                                 data = p[[i]],  model="UM.RS")
  
  freq_normbin[[i]] <- summary(freq_norm1bin[[i]])
  
}


##########
mu_fr_norm.bin = c()
tau2_fr_norm.bin = c()
LB_mu_fr_norm.bin = c()
UB_mu_fr_norm.bin = c()
prec_mu_fr_norm.bin = c()
for(i in 1:N.sim){
  mu_fr_norm.bin[i] <- freq_normbin[[i]]$beta
  LB_mu_fr_norm.bin[i] <- freq_normbin[[i]]$ci.lb
  UB_mu_fr_norm.bin[i] <- freq_normbin[[i]]$ci.ub
  prec_mu_fr_norm.bin[i] <- (UB_mu_fr_norm.bin[i] - LB_mu_fr_norm.bin[i]) 
  tau2_fr_norm.bin[i] <- freq_normbin[[i]]$tau2
  
}

tau_fr_norm.bin = sqrt(tau2_fr_norm.bin)

############## SIMULATION MEASURES ##############
########## AVERAGE ABSOLUTE BIAS OF mu ##########
abs_bias_frnormal.bin <- c()
for (i in 1:N.sim){
  abs_bias_frnormal.bin[i] <- abs( mu_fr_norm.bin[i] - mu) 
}
avg_abs_bias_frnormal.bin <- round(mean(abs_bias_frnormal.bin), 2)

########## AVERAGE RELATIVE BIAS OF mu ##########
rel_bias_frnormal.bin <- c()
for (i in 1:N.sim){
  rel_bias_frnormal.bin[i] <- ( mu_fr_norm.bin[i] - mu) 
}
avg_rel_bias_frnormal.bin <- round(mean(rel_bias_frnormal.bin), 2)

############### MSE OF mu  ##################
MSE_frnormal.bin <- c()
for (i in 1:N.sim){
  MSE_frnormal.bin[i] <- RMSE(estimate = mu_fr_norm.bin[i],
                              parameter = mu,
                              type = "RMSE",
                              MSE = TRUE,           
                              percent = FALSE,
                              unname = FALSE)
}
avg_MSE_frnormal.bin <-  round(mean(MSE_frnormal.bin),4)

########## COVERAGE OF mu ##########
interval_contains_true_mean <- function(i) { 
  ( mu >= LB_mu_fr_norm.bin[i] ) && ( mu <= UB_mu_fr_norm.bin[i])
}
interval_contains_true_mean(i)

coverage_frnormal.bin <- c()
for (i in 1:N.sim){0
  if (interval_contains_true_mean(i) == TRUE){
    coverage_frnormal.bin[i] = 1
  }
  else {
    coverage_frnormal.bin[i]=0
  }
  print(coverage_frnormal.bin[i])
}
coverage_frnormal.bin
cover_frnormal.bin <- round(mean(coverage_frnormal.bin)*100, 2)


################
##### AVERAGE BIAS OF tau2 #########
bias_tau2_frnormal.bin <- c()
for (i in 1:N.sim){
  bias_tau2_frnormal.bin[i] <- abs( tau2_fr_norm.bin[i] - tau2)
  
}
### average bias of tau2
avg_bias_tau2_frnormal.bin <- round(mean(bias_tau2_frnormal.bin),2)
## percent bias of tau2 
pbias_tau2_frnormal.bin <- round(mean(bias_tau2_frnormal.bin / tau2), 2)
#

########### MSE OF tau2 ################
MSE_frnormaltau2.bin <- c()
for (i in 1:N.sim){
  MSE_frnormaltau2.bin[i] <- RMSE(estimate = tau2_fr_norm.bin[i],
                                  parameter = tau2,
                                  type = "RMSE",
                                  MSE = TRUE           
  )
  
}
avg_MSE_frnormaltau2.bin <- round(mean(MSE_frnormaltau2.bin),4)
#### normalized MSE
nMSE_frnormaltau2.bin = round(mean(MSE_frnormaltau2.bin / (tau2)^2 ), 4)


################
##### AVERAGE RELATIVE BIAS OF tau #########
rel_bias_tau_frnormal.bin <- c()
for (i in 1:N.sim){
  rel_bias_tau_frnormal.bin[i] <- (tau_fr_norm.bin[i] - tau)
  
}
### average bias of tau
avg_rel_bias_tau_frnormal.bin <- round(mean(rel_bias_tau_frnormal.bin),2)
## percent bias of tau2 
prel_bias_tau_frnormal.bin <- round(mean(rel_bias_tau_frnormal.bin / tau)*100, 2)
#
################
##### AVERAGE ABSOLUTE BIAS OF tau #########
abs_bias_tau_frnormal.bin <- c()
for (i in 1:N.sim){
  abs_bias_tau_frnormal.bin[i] <- abs(tau_fr_norm.bin[i] - tau)
  
}
### average bias of tau
avg_abs_bias_tau_frnormal.bin <- round(mean(abs_bias_tau_frnormal.bin),2)
## percent bias of tau2 
pabs_bias_tau_frnormal.bin <- round(mean(abs_bias_tau_frnormal.bin / tau)*100, 2)
#
########### MSE OF tau ################
MSE_frnormaltau.bin <- c()
for (i in 1:N.sim){
  MSE_frnormaltau.bin[i] <- RMSE(estimate = tau_fr_norm.bin[i],
                                 parameter = tau,
                                 type = "RMSE",
                                 MSE = TRUE           
  )
  
}
avg_MSE_frnormaltau.bin <- round(mean(MSE_frnormaltau.bin),4)
#### normalized MSE
nMSE_frnormaltau.bin = round(mean(MSE_frnormaltau.bin / tau2 ), 4)


fr_normal.bin <- cbind(mu_fr_norm.bin,LB_mu_fr_norm.bin, UB_mu_fr_norm.bin, tau2_fr_norm.bin,  tau_fr_norm.bin,
                       prec_mu_fr_norm.bin, 
                       abs_bias_frnormal.bin, rel_bias_frnormal.bin,
                       avg_abs_bias_frnormal.bin ,avg_rel_bias_frnormal.bin ,
                       MSE_frnormal.bin , avg_MSE_frnormal.bin ,coverage_frnormal.bin, cover_frnormal.bin,
                       bias_tau2_frnormal.bin, avg_bias_tau2_frnormal.bin, pbias_tau2_frnormal.bin,
                       MSE_frnormaltau2.bin, avg_MSE_frnormaltau2.bin, nMSE_frnormaltau2.bin,
                       rel_bias_tau_frnormal.bin, avg_rel_bias_tau_frnormal.bin, prel_bias_tau_frnormal.bin,
                       abs_bias_tau_frnormal.bin, avg_abs_bias_tau_frnormal.bin, pabs_bias_tau_frnormal.bin,
                       MSE_frnormaltau.bin, avg_MSE_frnormaltau.bin, nMSE_frnormaltau.bin)

write.csv(fr_normal.bin, "C:\\Users\\Tianqi YU\\Desktop\\SIMULATION SCENARIOS FINAL\\NORMAL DATA\\freq_normal.bin\\rev_Results_freq_N_jags_Normal_bin_data_k=14_mu=0_tau2=0.12.csv",row.names=FALSE )  

############# FREQUENTIST MODELS ############
############ NORMAL-NORMAL ###########
mu = 0
tau2 = 0.12
N.sim = 1000
tau = sqrt(tau2)

freq_norm1 <- list()
freq_norm <- list()
fr_norm <- list()
pred <- list()

for(i in 1:N.sim){
  set.seed(2109)
  freq_norm1[[i]] <- rma(measure = "OR",
                         ai = p[[i]]$ti, bi = p[[i]]$nti - p[[i]]$ti,
                         ci = p[[i]]$ci, di = p[[i]]$nci - p[[i]]$ci,
                         data = p[[i]], method = "REML")
  freq_norm[[i]] <- summary(freq_norm1[[i]])
  
  fr_norm[[i]] = as.data.frame(confint(freq_norm1[[i]]))
  
  pred[[i]] = predict.rma(freq_norm[[i]])
  
}


# ### study specific effects ####
stud_eff = list()
for(i in 1:N.sim){
  stud_eff[[i]] = blup( freq_norm[[i]], level = 95 )
  
}

write.csv(stud_eff, "C:\\Users\\Tianqi YU\\Desktop\\SIMULATION SCENARIOS FINAL\\NORMAL DATA\\freq_normal\\rev_Rel_eff_freq_N_jags_Normal_bin_data_k=14_mu=0_tau2=0.12.csv",row.names=FALSE )  

# ######################## ABSOLUTE BIAS OF RELATIVE EFFECTS ##################
abs_bias_relnormfr_distr <- list()
mean_abs_bias_relnormfr_distr <- list()
for(i in 1:N.sim){
  for(j in 1:nrow(p[[i]])){
    abs_bias_relnormfr_distr[[i]] <- round(abs(stud_eff[[i]]$pred - p[[i]]$true_eff),2)
    mean_abs_bias_relnormfr_distr[[i]] <- mean(abs_bias_relnormfr_distr[[i]])
  }
}

#mean_bias_relt_distr 

################## mean bias of the mean biases of each data set for relative effects  ###################
mean_abs_bias_of_mean_abs_bias_relnormfr_distr <- unlist(mean_abs_bias_relnormfr_distr)
mean_abs_bias_of_mean_abs_bias_relnormfr_distr  <- round(mean(mean_abs_bias_of_mean_abs_bias_relnormfr_distr),2)

write.csv(mean_abs_bias_of_mean_abs_bias_relnormfr_distr, "C:\\Users\\Tianqi YU\\Desktop\\SIMULATION SCENARIOS FINAL\\NORMAL DATA\\freq_normal\\rev_Mean_abs_bias_Rel_eff_freq_N_jags_Normal_bin_data_k=14_mu=0_tau2=0.12.csv",row.names=FALSE )  

# ######################## ARELATIVE BIAS OF RELATIVE EFFECTS ##################
rel_bias_relnormfr_distr <- list()
mean_rel_bias_relnormfr_distr <- list()
for(i in 1:N.sim){
  for(j in 1:nrow(p[[i]])){
    rel_bias_relnormfr_distr[[i]] <- round((stud_eff[[i]]$pred - p[[i]]$true_eff),2)
    mean_rel_bias_relnormfr_distr[[i]] <- mean(rel_bias_relnormfr_distr[[i]])
  }
}

#mean_bias_relt_distr 

################## mean bias of the mean biases of each data set for relative effects  ###################
mean_rel_bias_of_mean_rel_bias_relnormfr_distr <- unlist(mean_rel_bias_relnormfr_distr)
mean_rel_bias_of_mean_rel_bias_relnormfr_distr  <- round(mean(mean_rel_bias_of_mean_rel_bias_relnormfr_distr),2)

write.csv(mean_rel_bias_of_mean_rel_bias_relnormfr_distr, "C:\\Users\\Tianqi YU\\Desktop\\SIMULATION SCENARIOS FINAL\\NORMAL DATA\\freq_normal\\rev_rel_Mean_bias_Rel_eff_freq_N_jags_Normal_bin_data_k=14_mu=0_tau2=0.12.csv",row.names=FALSE )  

# #### OVERLAP #######
datanorm_fr <- list()
overlap_norm_fr1 <- list()
overlap_norm_fr <- list()
for(i in 1:N.sim){
  datanorm_fr[[i]] <- list(norm_fr = stud_eff[[i]]$pred , True = p[[i]]$true_eff)
  overlap_norm_fr1[[i]] <- overlap(datanorm_fr[[i]], type = "2")
  overlap_norm_fr[[i]] <- overlap_norm_fr1[[i]]$OV
}
overlap_norm_fr <- unlist(overlap_norm_fr)
avg_overlap_norm_fr <- round(mean(overlap_norm_fr),2)

overlap_normfr = cbind.data.frame(overlap_norm_fr, avg_overlap_norm_fr)

write.csv(overlap_normfr, "C:\\Users\\Tianqi YU\\Desktop\\SIMULATION SCENARIOS FINAL\\NORMAL DATA\\freq_normal\\rev_Overlap_freq_N_jags_Normal_bin_data_k=14_mu=0_tau2=0.12.csv",row.names=FALSE )  

#########
mu_fr_norm = c()
predict = c()
LB_predict = c()
UB_predict = c()
tau2_fr_norm = c()
LB_mu_fr_norm = c()
UB_mu_fr_norm = c()
LB_tau2_fr_norm = c()
UB_tau2_fr_norm = c()
tau_fr_norm = c()
LB_tau_fr_norm = c()
UB_tau_fr_norm = c()
prec_mu_fr_norm = c()
prec_tau2_fr_norm = c()
prec_pred_fr_norm = c()
prec_tau_fr_norm = c()
prec_predict = c()
delta_frnorm = list()
fr_tau2 = list()
fr_tau2_n = list()
fr_tau = list()
fr_tau_n = list()
for(i in 1:N.sim){
  mu_fr_norm[i] <- freq_norm[[i]]$beta
  LB_mu_fr_norm[i] <- freq_norm[[i]]$ci.lb
  UB_mu_fr_norm[i] <- freq_norm[[i]]$ci.ub
  prec_mu_fr_norm[i] <- (UB_mu_fr_norm[i] - LB_mu_fr_norm[i]) 
  
  fr_tau2[[i]] <- grepl("tau\\^2", row.names(fr_norm[[i]]))
  fr_tau2_n[[i]] <- fr_norm[[i]][fr_tau2[[i]],]
  tau2_fr_norm[i] <- fr_tau2_n[[i]][1]
  LB_tau2_fr_norm[i] <- fr_tau2_n[[i]][2]
  UB_tau2_fr_norm[i] <- fr_tau2_n[[i]][3]
  prec_tau2_fr_norm[i] <- (UB_tau2_fr_norm[i] - LB_tau2_fr_norm[i]) 
  
  
  tau_fr_norm[i] <- fr_norm[[i]][2]
  LB_tau_fr_norm[i] <- fr_norm[[i]][6]
  UB_tau_fr_norm[i] <- fr_norm[[i]][10]
  prec_tau_fr_norm[i] <- (UB_tau_fr_norm[i] - LB_tau_fr_norm[i]) 
  
  predict[i] = pred[[i]]$pred
  LB_predict[i] = pred[[i]]$pi.lb
  UB_predict[i] = pred[[i]]$pi.ub
  prec_predict[i] = UB_predict[i] - LB_predict[i]
  
}

############## SIMULATION MEASURES ##############
########## AVERAGE ABSOLUTE BIAS OF mu ##########
abs_bias_frnormal <- c()
for (i in 1:N.sim){
  abs_bias_frnormal[i] <- abs( mu_fr_norm[i] - mu) 
}
avg_abs_bias_frnormal <- round(mean(abs_bias_frnormal), 2)

########## AVERAGE RELATIVE BIAS OF mu ##########
rel_bias_frnormal <- c()
for (i in 1:N.sim){
  rel_bias_frnormal[i] <- ( mu_fr_norm[i] - mu) 
}
avg_rel_bias_frnormal <- round(mean(rel_bias_frnormal), 2)

############### MSE OF mu  ##################
MSE_frnormal <- c()
for (i in 1:N.sim){
  MSE_frnormal[i] <- RMSE(estimate = mu_fr_norm[i],
                          parameter = mu,
                          type = "RMSE",
                          MSE = TRUE,           
                          percent = FALSE,
                          unname = FALSE)
}
avg_MSE_frnormal <-  round(mean(MSE_frnormal),4)

########## COVERAGE OF mu ##########
interval_contains_true_mean <- function(i) { 
  ( mu >= LB_mu_fr_norm[i] ) && ( mu <= UB_mu_fr_norm[i])
}
interval_contains_true_mean(i)

coverage_frnormal <- c()
for (i in 1:N.sim){0
  if (interval_contains_true_mean(i) == TRUE){
    coverage_frnormal[i] = 1
  }
  else {
    coverage_frnormal[i]=0
  }
  print(coverage_frnormal[i])
}
coverage_frnormal
cover_frnormal <- round(mean(coverage_frnormal)*100, 2)


################
##### AVERAGE BIAS OF tau2 #########
bias_tau2_frnormal <- c()
for (i in 1:N.sim){
  bias_tau2_frnormal[i] <- abs( tau2_fr_norm[i] - tau2)
}
avg_bias_tau2_frnormal <- round(mean(bias_tau2_frnormal),2)

pbias_tau2_frnormal <- round(mean(bias_tau2_frnormal/ tau2),2)

########### MSE OF tau2 ################
MSE_frnormaltau2 <- c()
for (i in 1:N.sim){
  MSE_frnormaltau2[i] <- RMSE(estimate = tau2_fr_norm[i],
                              parameter = tau2,
                              type = "RMSE",
                              MSE = TRUE           
  )
}
avg_MSE_frnormaltau2 <- round(mean(MSE_frnormaltau2),4)

nMSE_frnormaltau2 <- round(mean(MSE_frnormaltau2 / (tau2)^2),4)

########## COVERAGE OF tau2 ##########
interval_contains_true_mean <- function(i) { #
  ( tau2 >= LB_tau2_fr_norm[i] ) && ( tau2 <= UB_tau2_fr_norm[i])
}
interval_contains_true_mean(i)

coverage_frnormal_tau2 <- c()
for (i in 1:N.sim){0
  if (interval_contains_true_mean(i) == TRUE){
    coverage_frnormal_tau2[i] = 1
  }
  else {
    coverage_frnormal_tau2[i]=0
  }
  print(coverage_frnormal_tau2[i])
}
coverage_frnormal_tau2
cover_frnormal_tau2 <- round(mean(coverage_frnormal_tau2)*100, 2)



################
##### AVERAGE RELATIVE BIAS OF tau #########
rel_bias_tau_frnormal <- c()
for (i in 1:N.sim){
  rel_bias_tau_frnormal[i] <- ( tau_fr_norm[i] - tau)
}
avg_rel_bias_tau_frnormal <- round(mean(rel_bias_tau_frnormal),2)

prel_bias_tau_frnormal <- round(mean(rel_bias_tau_frnormal/ tau)* 100,2)

##### AVERAGE absolute BIAS OF tau #########
abs_bias_tau_frnormal <- c()
for (i in 1:N.sim){
  abs_bias_tau_frnormal[i] <- abs( tau_fr_norm[i] - tau)
}
avg_abs_bias_tau_frnormal <- round(mean(abs_bias_tau_frnormal),2)

pabs_bias_tau_frnormal <- round(mean(abs_bias_tau_frnormal/ tau)*100,2)

########### MSE OF tau ################
MSE_frnormaltau <- c()
for (i in 1:N.sim){
  MSE_frnormaltau[i] <- RMSE(estimate = tau_fr_norm[i],
                             parameter = tau,
                             type = "RMSE",
                             MSE = TRUE           
  )
}
avg_MSE_frnormaltau <- round(mean(MSE_frnormaltau),4)

nMSE_frnormaltau <- round(mean(MSE_frnormaltau / tau2),4)

########## COVERAGE OF tau ##########
interval_contains_true_mean <- function(i) { #
  ( tau >= LB_tau_fr_norm[i] ) && ( tau <= UB_tau_fr_norm[i])
}
interval_contains_true_mean(i)

coverage_frnormal_tau <- c()
for (i in 1:N.sim){0
  if (interval_contains_true_mean(i) == TRUE){
    coverage_frnormal_tau[i] = 1
  }
  else {
    coverage_frnormal_tau[i]=0
  }
  print(coverage_frnormal_tau[i])
}
coverage_frnormal_tau
cover_frnormal_tau <- round(mean(coverage_frnormal_tau)*100, 2)

len_mu_frnormal = mean( prec_mu_fr_norm)
len_PI_frnormal = mean(prec_predict)

fr_normal <- cbind(mu_fr_norm,LB_mu_fr_norm, UB_mu_fr_norm, tau2_fr_norm, LB_tau2_fr_norm, UB_tau2_fr_norm, prec_mu_fr_norm, prec_tau2_fr_norm,
                   predict,LB_predict, UB_predict, prec_predict,tau_fr_norm, LB_tau_fr_norm, UB_tau_fr_norm, prec_tau_fr_norm,
                   abs_bias_frnormal,  rel_bias_frnormal, 
                   avg_abs_bias_frnormal,avg_rel_bias_frnormal,
                   MSE_frnormal, avg_MSE_frnormal,coverage_frnormal, cover_frnormal,
                   bias_tau2_frnormal, avg_bias_tau2_frnormal, pbias_tau2_frnormal, MSE_frnormaltau2, avg_MSE_frnormaltau2, nMSE_frnormaltau2, 
                   coverage_frnormal_tau2, cover_frnormal_tau2,
                   rel_bias_tau_frnormal, avg_rel_bias_tau_frnormal, prel_bias_tau_frnormal,
                   abs_bias_tau_frnormal, avg_abs_bias_tau_frnormal, pabs_bias_tau_frnormal,
                   MSE_frnormaltau, avg_MSE_frnormaltau, nMSE_frnormaltau,
                   len_mu_frnormal ,
                   len_PI_frnormal,
                   coverage_frnormal_tau, cover_frnormal_tau)


write.csv(fr_normal, "C:\\Users\\Tianqi YU\\Desktop\\SIMULATION SCENARIOS FINAL\\NORMAL DATA\\freq_normal\\rev_Results_freq_N_jags_Normal_bin_data_k=14_mu=0_tau2=0.12.csv",row.names=FALSE )  

############################## METAPLUS PACKAGE ######################

library(metaplus)

######################## METAPLUS T-DISTRIBUTION ######################
N.sim = 1000
mu = 0
tau2 = 0.12
tau = sqrt(tau2)

dat1 <- list()
for(i in 1:N.sim){
  dat1[[i]] = escalc(measure="OR",ai = p[[i]]$ti, bi =p[[i]]$nti - p[[i]]$ti,
                     ci =p[[i]]$ci, di = p[[i]]$nci - p[[i]]$ci,
                     data = p[[i]])
}


metaplust1 <- list()
metaplust <- list()
outliers <- list()
det_outliers <- list()
for(i in 1:N.sim){
  set.seed(21010)
  metaplust1[[i]] <- metaplus(yi =  dat1[[i]]$yi, sei = sqrt(dat1[[i]]$vi), data = dat1[[i]], random = "t-dist")
  metaplust[[i]] <- summary(metaplust1[[i]])
  
  # outliers[[i]] <- testOutliers(metaplust1[[i]])
  # det_outliers[[i]] <- summary(outliers[[i]])
}

mu_metaplust = c()
tau2_metaplust = c()
LB_mu_metaplust = c()
UB_mu_metaplust = c()
ind_normal <- c()
for(i in 1:N.sim){
  mu_metaplust[i] <- metaplust[[i]]$results[1]
  LB_mu_metaplust[i] <- metaplust[[i]]$results[4]
  UB_mu_metaplust[i] <- metaplust[[i]]$results[7]
  
  tau2_metaplust[i] <- metaplust[[i]]$results[2]
  
  ind_normal[i] <- metaplust[[i]]$results[3]
}

tau_metaplust = sqrt(tau2_metaplust)
ind_normal = round(ind_normal, 2)

##### REMOVE DATA SETS WHERE THE CALCULATION OF LOWER OR UPPER LIMIT OF THE PARAMETER OF INTEREST WAS NOT POSSIBLE ######

condition1 <- which(is.na(LB_mu_metaplust))
condition2 <- which(is.na(UB_mu_metaplust))

remove_t = c(condition1,condition2)
remove_t = unique(remove_t)

data_out_t <- length(remove_t) / 1000

if (length(remove_t)== 0) {
  LB_mu_metaplust <- LB_mu_metaplust
  UB_mu_metaplust <- UB_mu_metaplust
  mu_metaplust <- mu_metaplust
  tau2_metaplust <- tau2_metaplust
  ind_normal <- ind_normal
  prec_mu_metaplust = UB_mu_metaplust - LB_mu_metaplust
  N.sim = N.sim
} else {
  LB_mu_metaplust <- LB_mu_metaplust[-remove_t]
  UB_mu_metaplust <- UB_mu_metaplust[-remove_t]
  mu_metaplust <- mu_metaplust[-remove_t]
  tau2_metaplust <- tau2_metaplust[-remove_t]
  ind_normal <- ind_normal[-remove_t]
  prec_mu_metaplust = UB_mu_metaplust - LB_mu_metaplust
  N.sim = length(mu_metaplust )
}


############## SIMULATION MEASURES ##############
########## AVERAGE ABSOLUTE BIAS OF mu ##########
abs_bias_metaplust <- c()
for (i in 1:N.sim){
  abs_bias_metaplust[i] <- abs( mu_metaplust[i] - mu) 
}
avg_abs_bias_metaplust <- round(mean(abs_bias_metaplust), 2)

########## AVERAGE RELATIVE BIAS OF mu ##########
rel_bias_metaplust <- c()
for (i in 1:N.sim){
  rel_bias_metaplust[i] <- ( mu_metaplust[i] - mu) 
}
avg_rel_bias_metaplust <- round(mean(rel_bias_metaplust), 2)

############### MSE OF mu  ##################
MSE_metaplust <- c()
for (i in 1:N.sim){
  MSE_metaplust[i] <- RMSE(estimate = mu_metaplust[i],
                           parameter = mu,
                           type = "RMSE",
                           MSE = TRUE,           
                           percent = FALSE,
                           unname = FALSE)
}
avg_MSE_metaplust <-  round(mean(MSE_metaplust),4)

########## COVERAGE OF mu ##########
interval_contains_true_mean <- function(i) { 
  ( mu >= LB_mu_metaplust[i] ) && ( mu <= UB_mu_metaplust[i])
}
interval_contains_true_mean(i)

coverage_metaplust <- c()
for (i in 1:N.sim){0
  if (interval_contains_true_mean(i) == TRUE){
    coverage_metaplust[i] = 1
  }
  else {
    coverage_metaplust[i]=0
  }
  print(coverage_metaplust[i])
}
coverage_metaplust
cover_metaplust <- round(mean(coverage_metaplust)*100, 2)


################
##### AVERAGE BIAS OF tau2 #########
bias_tau2_metaplust <- c()
for (i in 1:N.sim){
  bias_tau2_metaplust[i] <- abs( tau2_metaplust[i] - tau2)
}
avg_bias_tau2_metaplust <- round(mean(bias_tau2_metaplust),2)
pbias_tau2_metaplust <- round(mean(bias_tau2_metaplust / tau2),2)

########### MSE OF tau2 ################
MSE_tau2_metaplust <- c()
for (i in 1:N.sim){
  MSE_tau2_metaplust[i] <- RMSE(estimate = tau2_metaplust[i],
                                parameter = tau2,
                                type = "RMSE",
                                MSE = TRUE         
  )
}

avg_MSE_tau2_metaplust <- round(mean(MSE_tau2_metaplust),4)
nMSE_tau2_metaplust <- round(mean(MSE_tau2_metaplust / (tau2)^2),4)


##### AVERAGE RELATVE BIAS OF tau #########
rel_bias_tau_metaplust <- c()
for (i in 1:N.sim){
  rel_bias_tau_metaplust[i] <- (tau_metaplust[i] - tau)
}
avg_rel_bias_tau_metaplust <- round(mean(rel_bias_tau_metaplust),2)
prel_bias_tau_metaplust <- round(mean(rel_bias_tau_metaplust / tau)*100,2)

##### AVERAGE ABSOLUTE BIAS OF tau #########
abs_bias_tau_metaplust <- c()
for (i in 1:N.sim){
  abs_bias_tau_metaplust[i] <- abs(tau_metaplust[i] - tau)
}
avg_abs_bias_tau_metaplust <- round(mean(abs_bias_tau_metaplust),2)
pabs_bias_tau_metaplust <- round(mean(abs_bias_tau_metaplust / tau)*100,2)

########### MSE OF tau ################
MSE_tau_metaplust <- c()
for (i in 1:N.sim){
  MSE_tau_metaplust[i] <- RMSE(estimate = tau_metaplust[i],
                               parameter = tau,
                               type = "RMSE",
                               MSE = TRUE         
  )
}

avg_MSE_tau_metaplust <- round(mean(MSE_tau_metaplust),4)
nMSE_tau_metaplust <- round(mean(MSE_tau_metaplust / tau2),4)



metaplus_t <- cbind(mu_metaplust ,LB_mu_metaplust, UB_mu_metaplust, tau2_metaplust, tau_metaplust,
                    prec_mu_metaplust, abs_bias_metaplust, rel_bias_metaplust,
                    avg_abs_bias_metaplust,avg_rel_bias_metaplust,
                    MSE_metaplust, avg_MSE_metaplust,coverage_metaplust, cover_metaplust, bias_tau2_metaplust,
                    avg_bias_tau2_metaplust,pbias_tau2_metaplust,
                    MSE_tau2_metaplust, avg_MSE_tau2_metaplust, nMSE_tau2_metaplust, 
                    rel_bias_tau_metaplust,abs_bias_tau_metaplust,
                    avg_rel_bias_tau_metaplust,avg_abs_bias_tau_metaplust,
                    prel_bias_tau_metaplust, pabs_bias_tau_metaplust,
                    MSE_tau_metaplust, avg_MSE_tau_metaplust, nMSE_tau_metaplust,ind_normal, data_out_t)

write.csv(metaplus_t, "C:\\Users\\Tianqi YU\\Desktop\\SIMULATION SCENARIOS FINAL\\NORMAL DATA\\freq_tmodel\\rev_Results_metaplust_jags_Normal_bin_data_k=14_mu=0_tau2=0.12.csv",row.names=FALSE )  

####################### METAPLUS TMIXTURE ######################

N.sim = 1000
mu = 0
tau2 = 0.12
tau = sqrt(tau2)

metaplusmix1 <- list()
metaplusmix <- list()
# outliersmix <- list()
# det_outliersmix <- list()
for(i in 1:N.sim){
  set.seed(21011)
  metaplusmix1[[i]] <- metaplus(yi =  dat1[[i]]$yi, sei = sqrt(dat1[[i]]$vi), data = dat1[[i]], random = "mixture")
  metaplusmix[[i]] <- summary(metaplusmix1[[i]])
  
  # outliers[[i]] <- testOutliers(metaplusmix1[[i]])
  # det_outliers[[i]] <- summary(outliers[[i]])
}

mu_metaplusmix = c()
tau2_metaplusmix = c()
LB_mu_metaplusmix = c()
UB_mu_metaplusmix = c()
tau2_out_metaplusmix = c()
prob_out <- c()
for(i in 1:N.sim){
  mu_metaplusmix[i] <- metaplusmix[[i]]$results[1]
  LB_mu_metaplusmix[i] <- metaplusmix[[i]]$results[5]
  UB_mu_metaplusmix[i] <- metaplusmix[[i]]$results[9]
  
  tau2_metaplusmix[i] <- metaplusmix[[i]]$results[2]
  
  tau2_out_metaplusmix[i] <- metaplusmix[[i]]$results[3]
  prob_out[i] <- metaplusmix[[i]]$results[4]
}

tau_metaplusmix = sqrt(tau2_metaplusmix)
tau_out_metaplusmix = sqrt(tau2_out_metaplusmix)
##### REMOVE DATA SETS WHERE THE CALCULATION OF LOWER OR UPPER LIMIT OF THE PARAMETER OF INTEREST WAS NOT POSSIBLE ######

condition1 <- which(is.na(LB_mu_metaplusmix))
condition2 <- which(is.na(UB_mu_metaplusmix))

remove_mix = c(condition1,condition2)
remove_mix = unique(remove_mix)

data_out_mix <- length(remove_mix) / 1000

if (length(remove_mix)== 0) {
  LB_mu_metaplusmix <- LB_mu_metaplusmix
  UB_mu_metaplusmix <- UB_mu_metaplusmix
  mu_metaplusmix <- mu_metaplusmix
  tau2_metaplusmix <- tau2_metaplusmix
  prec_mu_metaplusmix <- UB_mu_metaplusmix - LB_mu_metaplusmix
  tau2_out_metaplusmix <- tau2_out_metaplusmix 
  prob_out <- prob_out
  N.sim = N.sim
} else {
  LB_mu_metaplusmix <- LB_mu_metaplusmix[-remove_mix]
  UB_mu_metaplusmix <- UB_mu_metaplusmix[-remove_mix]
  mu_metaplusmix <- mu_metaplusmix[-remove_mix]
  tau2_metaplusmix <- tau2_metaplusmix[-remove_mix]
  prec_mu_metaplusmix <- UB_mu_metaplusmix - LB_mu_metaplusmix
  tau2_out_metaplusmix <- tau2_out_metaplusmix[-remove_mix] 
  prob_out <- prob_out[-remove_mix]
  N.sim = length(mu_metaplusmix)
}


############## SIMULATION MEASURES ##############
########## AVERAGE ABSOLUTE BIAS OF mu ##########
abs_bias_metaplusmix <- c()
for (i in 1:N.sim){
  abs_bias_metaplusmix[i] <- abs( mu_metaplusmix[i] - mu) 
}
avg_abs_bias_metaplusmix <- round(mean(abs_bias_metaplusmix), 2)

########## AVERAGE RELATIVE BIAS OF mu ##########
rel_bias_metaplusmix <- c()
for (i in 1:N.sim){
  rel_bias_metaplusmix[i] <- ( mu_metaplusmix[i] - mu) 
}
avg_rel_bias_metaplusmix <- round(mean(rel_bias_metaplusmix), 2)

############### MSE OF mu  ##################
MSE_metaplusmix <- c()
for (i in 1:N.sim){
  MSE_metaplusmix[i] <- RMSE(estimate = mu_metaplusmix[i],
                             parameter = mu,
                             type = "RMSE",
                             MSE = TRUE,           
                             percent = FALSE,
                             unname = FALSE)
}
avg_MSE_metaplusmix <-  round(mean(MSE_metaplusmix),4)

########## COVERAGE OF mu ##########
interval_contains_true_mean <- function(i) { 
  ( mu >= LB_mu_metaplusmix[i] ) && ( mu <= UB_mu_metaplusmix[i])
}
interval_contains_true_mean(i)

coverage_metaplusmix <- c()
for (i in 1:N.sim){0
  if (interval_contains_true_mean(i) == TRUE){
    coverage_metaplusmix[i] = 1
  }
  else {
    coverage_metaplusmix[i]=0
  }
  print(coverage_metaplusmix[i])
}
coverage_metaplusmix
cover_metaplusmix <- round(mean(coverage_metaplusmix)*100, 2)


################
##### AVERAGE BIAS OF tau2 #########
bias_tau2_metaplusmix <- c()
for (i in 1:N.sim){
  bias_tau2_metaplusmix[i] <- abs( tau2_out_metaplusmix[i] - tau2)
}
avg_bias_tau2_metaplusmix <- round(mean(bias_tau2_metaplusmix),2)
pbias_tau2_metaplusmix <- round(mean(bias_tau2_metaplusmix / tau2),2)

########### MSE OF tau2 ################
MSE_tau2_metaplusmix <- c()
for (i in 1:N.sim){
  MSE_tau2_metaplusmix[i] <- RMSE(estimate = tau2_out_metaplusmix[i],
                                  parameter = tau2,
                                  type = "RMSE",
                                  MSE = TRUE         
  )
}
avg_MSE_tau2_metaplusmix <- round(mean(MSE_tau2_metaplusmix),4)
nMSE_tau2_metaplusmix <- round(mean(MSE_tau2_metaplusmix / (tau2)^2 ),4)


###############
##### AVERAGE RELATIVE BIAS OF tau #########
rel_bias_tau_metaplusmix <- c()
for (i in 1:N.sim){
  rel_bias_tau_metaplusmix[i] <- (tau_out_metaplusmix[i] - tau)
}
avg_rel_bias_tau_metaplusmix <- round(mean(rel_bias_tau_metaplusmix),2)
prel_bias_tau_metaplusmix <- round(mean(rel_bias_tau_metaplusmix / tau)*100,2)

##### AVERAGE ABSOLUTE BIAS OF tau #########
abs_bias_tau_metaplusmix <- c()
for (i in 1:N.sim){
  abs_bias_tau_metaplusmix[i] <- abs(tau_out_metaplusmix[i] - tau)
}
avg_abs_bias_tau_metaplusmix <- round(mean(abs_bias_tau_metaplusmix),2)
pabs_bias_tau_metaplusmix <- round(mean(abs_bias_tau_metaplusmix / tau)*100,2)

########### MSE OF tau ################
MSE_tau_metaplusmix <- c()
for (i in 1:N.sim){
  MSE_tau_metaplusmix[i] <- RMSE(estimate = tau_out_metaplusmix[i],
                                 parameter = tau,
                                 type = "RMSE",
                                 MSE = TRUE         
  )
}
avg_MSE_tau_metaplusmix <- round(mean(MSE_tau_metaplusmix),4)
nMSE_tau_metaplusmix <- round(mean(MSE_tau_metaplusmix / tau2 ),4)


###############
##### AVERAGE RELATIVE BIAS OF tau_common #########
rel_bias_tau_metaplusmix1 <- c()
for (i in 1:N.sim){
  rel_bias_tau_metaplusmix1[i] <- (tau_metaplusmix[i] - tau)
}
avg_rel_bias_tau_metaplusmix1 <- round(mean(rel_bias_tau_metaplusmix1),2)
prel_bias_tau_metaplusmix1 <- round(mean(rel_bias_tau_metaplusmix1 / tau)*100,2)

##### AVERAGE ABSOLUTE BIAS OF tau #########
abs_bias_tau_metaplusmix1 <- c()
for (i in 1:N.sim){
  abs_bias_tau_metaplusmix1[i] <- abs(tau_metaplusmix[i] - tau)
}
avg_abs_bias_tau_metaplusmix1 <- round(mean(abs_bias_tau_metaplusmix1),2)
pabs_bias_tau_metaplusmix1 <- round(mean(abs_bias_tau_metaplusmix1 / tau)*100,2)

########### MSE OF tau ################
MSE_tau_metaplusmix1 <- c()
for (i in 1:N.sim){
  MSE_tau_metaplusmix1[i] <- RMSE(estimate = tau_metaplusmix[i],
                                  parameter = tau,
                                  type = "RMSE",
                                  MSE = TRUE         
  )
}
avg_MSE_tau_metaplusmix1 <- round(mean(MSE_tau_metaplusmix1),4)
nMSE_tau_metaplusmix1 <- round(mean(MSE_tau_metaplusmix1 / tau2 ),4)




metaplus_mix <- cbind(mu_metaplusmix ,LB_mu_metaplusmix, UB_mu_metaplusmix, tau2_metaplusmix, tau2_out_metaplusmix,
                      tau_metaplusmix, tau_out_metaplusmix,
                      prec_mu_metaplusmix, abs_bias_metaplusmix,  rel_bias_metaplusmix, 
                      avg_abs_bias_metaplusmix,avg_rel_bias_metaplusmix,
                      MSE_metaplusmix, avg_MSE_metaplusmix,coverage_metaplusmix, cover_metaplusmix, bias_tau2_metaplusmix, 
                      avg_bias_tau2_metaplusmix, pbias_tau2_metaplusmix, 
                      MSE_tau2_metaplusmix, avg_MSE_tau2_metaplusmix, nMSE_tau2_metaplusmix, 
                      avg_rel_bias_tau_metaplusmix, prel_bias_tau_metaplusmix,
                      avg_abs_bias_tau_metaplusmix, pabs_bias_tau_metaplusmix,
                      rel_bias_tau_metaplusmix, abs_bias_tau_metaplusmix,
                      MSE_tau_metaplusmix, avg_MSE_tau_metaplusmix, nMSE_tau_metaplusmix, 
                      
                      avg_rel_bias_tau_metaplusmix1, prel_bias_tau_metaplusmix1,
                      avg_abs_bias_tau_metaplusmix1, pabs_bias_tau_metaplusmix1,
                      rel_bias_tau_metaplusmix1, abs_bias_tau_metaplusmix1,
                      MSE_tau_metaplusmix1, avg_MSE_tau_metaplusmix1, nMSE_tau_metaplusmix1, 
                      
                      prob_out,data_out_mix )

write.csv(metaplus_mix, "C:\\Users\\Tianqi YU\\Desktop\\SIMULATION SCENARIOS FINAL\\NORMAL DATA\\freq_mixture\\rev_Results_metaplusmix_jags_Normal_bin_data_k=14_mu=0_tau2=0.12.csv",row.names=FALSE )  



################## summary of the results #####################
r = 16
c = 17
matrix = matrix(nrow = r, ncol = c)

# Assign values to each cell in the matrix
for(i in 1:r){
  for(j in 1:c){
    matrix[1,1] <- "Models"
    
    matrix[1,2] <- "avg_abs_bias_mu" 
    matrix[1,3] <- "avg_MSE_mu"
    matrix[1,4] <- "avg_bias_tau2"
    matrix[1,5] <- "avg_MSE_tau2"
    
    matrix[1,6] <- "coverage_mu"
    matrix[1,7] <- "coverage_tau2"
    
    matrix[1,8] <- "pbias_tau2"
    matrix[1,9] <- "nMSE_tau2"
    
    matrix[1,10] <- "avg_rel_bias_tau"
    matrix[1,11] <- "MSE_tau"
    matrix[1,12] <- "coverage_tau"
    matrix[1,13] <- "prel_bias_tau"
    matrix[1,14] <- "nMSE_tau"
    
    matrix[1,15] <- "abs_bias_tau"
    matrix[1,16] <- "pabs_bias_tau"
    
    matrix[1,17] <- "avg_rel_bias_mu"
    
    matrix[2,1] <- "Binomial-Normal(Unif)"
    matrix[3,1] <- "Binomial-Normal(HN)"
    matrix[4,1] <- "Normal-Normal"
    matrix[5,1] <- "Binomial-t(Unif)"
    matrix[6,1] <- "Normal-t"
    matrix[7,1] <- "Binomial-SN(Unif)"
    matrix[8,1] <- "Common-mean-mixture"
    matrix[9,1] <- "Binomial-DP-26(HN/Unif)"
    matrix[10,1] <- "Binomial-DP-51(HN/Unif)"
    matrix[11,1] <- "Binomial-DP-26(Unif/Unif)"
    matrix[12,1] <- "Binomial-DP-n(Unif/Gamma)"
    matrix[13,1] <- "Binomial-DP-51(Unif/Unif)"
    matrix[14,1] <- "Binomial-normal-freq"
    matrix[15,1] <- "Binomial-SN(HN)"
    matrix[16,1] <- "Binomial-t(HN)"
    
    #### FOR ABSOLUTE BIAS MU 
    if(RhatN_outU > 0.05 ){
      matrix[2,2] <-"NON-CONVERGENCE"}
    else{matrix[2,2] <-avg_abs_biasNormaljags.binU}
    
    if(RhatN_outHN > 0.05 ){
      matrix[3,2] <-"NON-CONVERGENCE"}
    else{matrix[3,2] <-avg_abs_biasNormaljags.binHN}
    
    matrix[4,2] <-avg_abs_bias_frnormal
    
    if(RhatT_out_U > 0.05 ){
      matrix[5,2] <-"NON-CONVERGENCE"}
    else{matrix[5,2] <- avg_abs_biast_distrjags.bin_U}
    
    if(data_out_t > 0.05){
      matrix[6,2] <- "NON-CONVERGENCE"
    }
    else{matrix[6,2] <- avg_abs_bias_metaplust}
    
    if(RhatSN_out_U > 0.05){
      matrix[7,2] <- "NON-CONVERGENCE"
    }
    else{ matrix[7,2] <- avg_abs_biasSNjags.bin_U}
    
    if(data_out_mix > 0.05){
      matrix[8,2] <- "NON-CONVERGENCE"
    }
    else{ matrix[8,2] <- avg_abs_bias_metaplusmix}
    
    if(RhatDP_out_HN_U_26>0.05){
      matrix[9,2] <- "NON-CONVERGENCE"
    }
    else{matrix[9,2] <- avg_abs_biasDPjags.bin_HN_U_26}
    
    if(RhatDP_out_HN_U_51>0.05){
      matrix[10,2] <- "NON-CONVERGENCE"
    }
    else{matrix[10,2] <- avg_abs_biasDPjags.bin_HN_U_51}
    
    if(RhatDP_out_U_U_26>0.05){
      matrix[11,2] <- "NON-CONVERGENCE"
    }
    else{matrix[11,2] <- avg_abs_biasDPjags.bin_U_U_26}
    
    if(RhatDP_out_U_G_n>0.05){
      matrix[12,2] <- "NON-CONVERGENCE"
    }
    else{matrix[12,2] <- avg_abs_biasDPjags.bin_U_G_n}
    
    if(RhatDP_out_U_U_51>0.05){
      matrix[13,2] <- "NON-CONVERGENCE"
    }
    else{matrix[13,2] <- avg_abs_biasDPjags.bin_U_U_51}
    
    matrix[14,2] <-avg_abs_bias_frnormal.bin
    
    if(RhatSN_out_HN > 0.05){
      matrix[15,2] <- "NON-CONVERGENCE"
    }
    else{ matrix[15,2] <- avg_abs_biasSNjags.bin_HN}
    
    if(RhatT_out_HN > 0.05 ){
      matrix[16,2] <-"NON-CONVERGENCE"}
    else{matrix[16,2] <- avg_abs_biast_distrjags.bin_HN}
    
    ######### FOR MSE MU 
    
    if(RhatN_outU > 0.05 ){
      matrix[2,3] <-"NON-CONVERGENCE"}
    else{matrix[2,3] <-avg_MSE_Normaljags.binU}
    
    if(RhatN_outHN > 0.05 ){
      matrix[3,3] <-"NON-CONVERGENCE"}
    else{matrix[3,3] <-avg_MSE_Normaljags.binHN}
    
    matrix[4,3] <- avg_MSE_frnormal
    
    if(RhatT_out_U > 0.05 ){
      matrix[5,3] <-"NON-CONVERGENCE"}
    else{matrix[5,3] <- avg_MSE_t_distrjags.bin_U}
    
    if(data_out_t > 0.05){
      matrix[6,3] <- "NON-CONVERGENCE"
    }
    else{matrix[6,3] <- avg_MSE_metaplust}
    
    if(RhatSN_out_U > 0.05){
      matrix[7,3] <- "NON-CONVERGENCE"
    }
    else{ matrix[7,3] <- avg_MSE_SNjags.bin_U}
    
    if(data_out_mix > 0.05){
      matrix[8,3] <- "NON-CONVERGENCE"
    }
    else{ matrix[8,3] <- avg_MSE_metaplusmix}
    
    if(RhatDP_out_HN_U_26>0.05){
      matrix[9,3] <- "NON-CONVERGENCE"
    }
    else{matrix[9,3] <-  avg_MSE_DPjags.bin_HN_U_26}
    
    if(RhatDP_out_HN_U_51>0.05){
      matrix[10,3] <- "NON-CONVERGENCE"
    }
    else{matrix[10,3] <- avg_MSE_DPjags.bin_HN_U_51}
    
    if(RhatDP_out_U_U_26>0.05){
      matrix[11,3] <- "NON-CONVERGENCE"
    }
    else{matrix[11,3] <- avg_MSE_DPjags.bin_U_U_26}
    
    if(RhatDP_out_U_G_n>0.05){
      matrix[12,3] <- "NON-CONVERGENCE"
    }
    else{matrix[12,3] <- avg_MSE_DPjags.bin_U_G_n}
    
    if(RhatDP_out_U_U_51 >0.05){
      matrix[13,3] <- "NON-CONVERGENCE"
    }
    else{matrix[13,3] <- avg_MSE_DPjags.bin_U_U_51}
    
    matrix[14,3] <- avg_MSE_frnormal.bin
    
    if(RhatSN_out_HN > 0.05){
      matrix[15,3] <- "NON-CONVERGENCE"
    }
    else{ matrix[15,3] <- avg_MSE_SNjags.bin_HN}
    
    if(RhatT_out_HN > 0.05 ){
      matrix[16,3] <-"NON-CONVERGENCE"}
    else{matrix[16,3] <- avg_MSE_t_distrjags.bin_HN}
    
    ######### FOR ABSOLUTE BIAS TAU2
    
    if(RhatN_outU > 0.05 ){
      matrix[2,4] <-"NON-CONVERGENCE"}
    else{matrix[2,4] <-avg_bias_tau2Normaljags.binU}
    
    if(RhatN_outHN > 0.05 ){
      matrix[3,4] <-"NON-CONVERGENCE"}
    else{matrix[3,4] <-avg_bias_tau2Normaljags.binHN}
    
    matrix[4,4] <- avg_bias_tau2_frnormal
    
    if(RhatT_out_U > 0.05 ){
      matrix[5,4] <-"NON-CONVERGENCE"}
    else{matrix[5,4] <- avg_bias_tau2t_distrjags.bin_U}
    
    if(data_out_t > 0.05){
      matrix[6,4] <- "NON-CONVERGENCE"
    }
    else{matrix[6,4] <- avg_bias_tau2_metaplust}
    
    if(RhatSN_out_U > 0.05){
      matrix[7,4] <- "NON-CONVERGENCE"
    }
    else{ matrix[7,4] <- avg_bias_tau2SNjags.bin_U}
    
    if(data_out_mix > 0.05){
      matrix[8,4] <- "NON-CONVERGENCE"
    }
    else{ matrix[8,4] <- avg_bias_tau2_metaplusmix}
    
    if(RhatDP_out_HN_U_26>0.05){
      matrix[9,4] <- "NON-CONVERGENCE"
    }
    else{matrix[9,4] <- avg_bias_tau2DPjags.bin_HN_U_26}
    
    if(RhatDP_out_HN_U_51>0.05){
      matrix[10,4] <- "NON-CONVERGENCE"
    }
    else{matrix[10,4] <- avg_bias_tau2DPjags.bin_HN_U_51}
    
    if(RhatDP_out_U_U_26>0.05){
      matrix[11,4] <- "NON-CONVERGENCE"
    }
    else{matrix[11,4] <- avg_bias_tau2DPjags.bin_U_U_26}
    
    if(RhatDP_out_U_G_n>0.05){
      matrix[12,4] <- "NON-CONVERGENCE"
    }
    else{matrix[12,4] <- avg_bias_tau2DPjags.bin_U_G_n}
    
    if(RhatDP_out_U_U_51>0.05){
      matrix[13,4] <- "NON-CONVERGENCE"
    }
    else{matrix[13,4] <- avg_bias_tau2DPjags.bin_U_U_51}
    
    matrix[14,4] <- avg_bias_tau2_frnormal.bin
    
    if(RhatSN_out_HN > 0.05){
      matrix[15,4] <- "NON-CONVERGENCE"
    }
    else{ matrix[15,4] <- avg_bias_tau2SNjags.bin_HN}
    
    if(RhatT_out_HN > 0.05 ){
      matrix[16,4] <-"NON-CONVERGENCE"}
    else{matrix[16,4] <- avg_bias_tau2t_distrjags.bin_HN}
    
    
    ######### FOR MSE TAU2
    
    if(RhatN_outU > 0.05 ){
      matrix[2,5] <-"NON-CONVERGENCE"}
    else{matrix[2,5] <-avg_MSE_Normaljags.tau2.binU}
    
    if(RhatN_outHN > 0.05 ){
      matrix[3,5] <-"NON-CONVERGENCE"}
    else{matrix[3,5] <-avg_MSE_Normaljags.tau2.binHN}
    
    matrix[4,5] <- avg_MSE_frnormaltau2
    
    if(RhatT_out_U > 0.05 ){
      matrix[5,5] <-"NON-CONVERGENCE"}
    else{matrix[5,5] <- avg_MSE_t_distrjags.tau2.bin_U}
    
    if(data_out_t > 0.05){
      matrix[6,5] <- "NON-CONVERGENCE"
    }
    else{matrix[6,5] <- avg_MSE_tau2_metaplust}
    
    if(RhatSN_out_U > 0.05){
      matrix[7,5] <- "NON-CONVERGENCE"
    }
    else{ matrix[7,5] <- avg_MSE_SNjags.tau2.bin_U}
    
    if(data_out_mix > 0.05){
      matrix[8,5] <- "NON-CONVERGENCE"
    }
    else{ matrix[8,5] <- avg_MSE_tau2_metaplusmix}
    
    if(RhatDP_out_HN_U_26>0.05){
      matrix[9,5] <- "NON-CONVERGENCE"
    }
    else{matrix[9,5] <- avg_MSE_DPjags.bin.tau2_HN_U_26}
    
    if(RhatDP_out_HN_U_51>0.05){
      matrix[10,5] <- "NON-CONVERGENCE"
    }
    else{matrix[10,5] <- avg_MSE_DPjags.bin.tau2_HN_U_51}
    
    if(RhatDP_out_U_U_26>0.05){
      matrix[11,5] <- "NON-CONVERGENCE"
    }
    else{matrix[11,5] <- avg_MSE_DPjags.bin.tau2_U_U_26}
    
    if(RhatDP_out_U_G_n>0.05){
      matrix[12,5] <- "NON-CONVERGENCE"
    }
    else{matrix[12,5] <- avg_MSE_DPjags.bin.tau2_U_G_n}
    
    if(RhatDP_out_U_U_51>0.05){
      matrix[13,5] <- "NON-CONVERGENCE"
    }
    else{matrix[13,5] <- avg_MSE_DPjags.bin.tau2_U_U_51}
    
    matrix[14,5] <- avg_MSE_frnormaltau2.bin
    
    if(RhatSN_out_HN > 0.05){
      matrix[15,5] <- "NON-CONVERGENCE"
    }
    else{ matrix[15,5] <- avg_MSE_SNjags.tau2.bin_HN}
    
    if(RhatT_out_HN > 0.05 ){
      matrix[16,5] <-"NON-CONVERGENCE"}
    else{matrix[16,5] <- avg_MSE_t_distrjags.tau2.bin_HN}
    
    ######### FOR COVERAGE MU
    
    if(RhatN_outU > 0.05 ){
      matrix[2,6] <-"NON-CONVERGENCE"}
    else{matrix[2,6] <-coverNormaljags.binU}
    
    if(RhatN_outHN > 0.05 ){
      matrix[3,6] <-"NON-CONVERGENCE"}
    else{matrix[3,6] <-coverNormaljags.binHN}
    
    matrix[4,6] <- cover_frnormal
    
    if(RhatT_out_U > 0.05 ){
      matrix[5,6] <-"NON-CONVERGENCE"}
    else{matrix[5,6] <- covert_distrjags.bin_U}
    
    if(data_out_t > 0.05){
      matrix[6,6] <- "NON-CONVERGENCE"
    }
    else{matrix[6,6] <- cover_metaplust}
    
    if(RhatSN_out_U > 0.05){
      matrix[7,6] <- "NON-CONVERGENCE"
    }
    else{ matrix[7,6] <- coverSNjags.bin_U}
    
    if(data_out_mix > 0.05){
      matrix[8,6] <- "NON-CONVERGENCE"
    }
    else{ matrix[8,6] <- cover_metaplusmix}
    
    if(RhatDP_out_HN_U_26>0.05){
      matrix[9,6] <- "NON-CONVERGENCE"
    }
    else{matrix[9,6] <- coverDPjags.bin_HN_U_26}
    
    if(RhatDP_out_HN_U_51>0.05){
      matrix[10,6] <- "NON-CONVERGENCE"
    }
    else{matrix[10,6] <- coverDPjags.bin_HN_U_51}
    
    if(RhatDP_out_U_U_26>0.05){
      matrix[11,6] <- "NON-CONVERGENCE"
    }
    else{matrix[11,6] <- coverDPjags.bin_U_U_26}
    
    if(RhatDP_out_U_G_n>0.05){
      matrix[12,6] <- "NON-CONVERGENCE"
    }
    else{matrix[12,6] <- coverDPjags.bin_U_G_n}
    
    if(RhatDP_out_U_U_51>0.05){
      matrix[13,6] <- "NON-CONVERGENCE"
    }
    else{matrix[13,6] <- coverDPjags.bin_U_U_51}
    
    matrix[14,6] <- cover_frnormal.bin
    
    if(RhatSN_out_HN > 0.05){
      matrix[15,6] <- "NON-CONVERGENCE"
    }
    else{ matrix[15,6] <- coverSNjags.bin_HN}
    
    if(RhatT_out_HN > 0.05 ){
      matrix[16,6] <-"NON-CONVERGENCE"}
    else{matrix[16,6] <- covert_distrjags.bin_HN}
    
    ######### FOR COVERAGE TAU2
    
    if(RhatN_outU > 0.05 ){
      matrix[2,7] <-"NON-CONVERGENCE"}
    else{matrix[2,7] <-cover_tau2Normaljags.binU}
    
    if(RhatN_outHN > 0.05 ){
      matrix[3,7] <-"NON-CONVERGENCE"}
    else{matrix[3,7] <-cover_tau2Normaljags.binHN}
    
    matrix[4,7] <- cover_frnormal_tau2
    
    if(RhatT_out_U > 0.05 ){
      matrix[5,7] <-"NON-CONVERGENCE"}
    else{matrix[5,7] <- cover_tau2t_distrjags.bin_U}
    
    matrix[6,7] <- "NA"
    
    if(RhatSN_out_U > 0.05){
      matrix[7,7] <- "NON-CONVERGENCE"
    }
    else{ matrix[7,7] <- cover_tau2SNjags.bin_U}
    
    matrix[8,7] <- "NA"
    
    if(RhatDP_out_HN_U_26>0.05){
      matrix[9,7] <- "NON-CONVERGENCE"
    }
    else{matrix[9,7] <- cover_tau2DPjags.bin_HN_U_26}
    
    if(RhatDP_out_HN_U_51>0.05){
      matrix[10,7] <- "NON-CONVERGENCE"
    }
    else{matrix[10,7] <- cover_tau2DPjags.bin_HN_U_51}
    
    if(RhatDP_out_U_U_26>0.05){
      matrix[11,7] <- "NON-CONVERGENCE"
    }
    else{matrix[11,7] <- cover_tau2DPjags.bin_U_U_26}
    
    if(RhatDP_out_U_G_n>0.05){
      matrix[12,7] <- "NON-CONVERGENCE"
    }
    else{matrix[12,7] <- cover_tau2DPjags.bin_U_G_n}
    
    if(RhatDP_out_U_U_51>0.05){
      matrix[13,7] <- "NON-CONVERGENCE"
    }
    else{matrix[13,7] <- cover_tau2DPjags.bin_U_U_51}
    
    matrix[14,7] <- "NA"
    
    if(RhatSN_out_HN > 0.05){
      matrix[15,7] <- "NON-CONVERGENCE"
    }
    else{ matrix[15,7] <- cover_tau2SNjags.bin_HN}
    
    if(RhatT_out_HN > 0.05 ){
      matrix[16,7] <-"NON-CONVERGENCE"}
    else{matrix[16,7] <- cover_tau2t_distrjags.bin_HN}
    
    ######### FOR pABSOLUTE_BIAS TAU2
    
    if(RhatN_outU > 0.05 ){
      matrix[2,8] <-"NON-CONVERGENCE"}
    else{matrix[2,8] <- pbias_tau2Normaljags.binU}
    
    if(RhatN_outHN > 0.05 ){
      matrix[3,8] <-"NON-CONVERGENCE"}
    else{matrix[3,8] <-pbias_tau2Normaljags.binHN}
    
    matrix[4,8] <- pbias_tau2_frnormal
    
    if(RhatT_out_U > 0.05 ){
      matrix[5,8] <-"NON-CONVERGENCE"}
    else{matrix[5,8] <- pbias_tau2t_distrjags.bin_U}
    
    if(data_out_t > 0.05){
      matrix[6,8] <- "NON-CONVERGENCE"
    }
    else{matrix[6,8] <- pbias_tau2_metaplust}
    
    if(RhatSN_out_U > 0.05){
      matrix[7,8] <- "NON-CONVERGENCE"
    }
    else{ matrix[7,8] <- pbias_tau2SNjags.bin_U}
    
    if(data_out_mix > 0.05){
      matrix[8,8] <- "NON-CONVERGENCE"
    }
    else{ matrix[8,8] <- pbias_tau2_metaplusmix}
    
    if(RhatDP_out_HN_U_26>0.05){
      matrix[9,8] <- "NON-CONVERGENCE"
    }
    else{matrix[9,8] <- pbias_tau2DPjags.bin_HN_U_26}
    
    if(RhatDP_out_HN_U_51>0.05){
      matrix[10,8] <- "NON-CONVERGENCE"
    }
    else{matrix[10,8] <- pbias_tau2DPjags.bin_HN_U_51}
    
    if(RhatDP_out_U_U_26>0.05){
      matrix[11,8] <- "NON-CONVERGENCE"
    }
    else{matrix[11,8] <- pbias_tau2DPjags.bin_U_U_26}
    
    if(RhatDP_out_U_G_n>0.05){
      matrix[12,8] <- "NON-CONVERGENCE"
    }
    else{matrix[12,8] <- pbias_tau2DPjags.bin_U_G_n}
    
    if(RhatDP_out_U_U_51>0.05){
      matrix[13,8] <- "NON-CONVERGENCE"
    }
    else{matrix[13,8] <- pbias_tau2DPjags.bin_U_U_51}
    
    matrix[14,8] <- pbias_tau2_frnormal.bin
    
    if(RhatSN_out_HN > 0.05){
      matrix[15,8] <- "NON-CONVERGENCE"
    }
    else{ matrix[15,8] <- pbias_tau2SNjags.bin_HN}
    
    if(RhatT_out_HN > 0.05 ){
      matrix[16,8] <-"NON-CONVERGENCE"}
    else{matrix[16,8] <- pbias_tau2t_distrjags.bin_HN}
    
    ######### FOR nMSE of TAU2
    
    if(RhatN_outU > 0.05 ){
      matrix[2,9] <-"NON-CONVERGENCE"}
    else{matrix[2,9] <-nMSE_Normaljags.tau2.binU}
    
    if(RhatN_outHN > 0.05 ){
      matrix[3,9] <-"NON-CONVERGENCE"}
    else{matrix[3,9] <-nMSE_Normaljags.tau2.binHN}
    
    matrix[4,9] <- nMSE_frnormaltau2
    
    if(RhatT_out_U > 0.05 ){
      matrix[5,9] <-"NON-CONVERGENCE"}
    else{matrix[5,9] <- nMSE_t_distrjags.tau2.bin_U}
    
    if(data_out_t > 0.05){
      matrix[6,9] <- "NON-CONVERGENCE"
    }
    else{matrix[6,9] <- nMSE_tau2_metaplust}
    
    if(RhatSN_out_U > 0.05){
      matrix[7,9] <- "NON-CONVERGENCE"
    }
    else{ matrix[7,9] <- nMSE_SNjags.tau2.bin_U}
    
    if(data_out_mix > 0.05){
      matrix[8,9] <- "NON-CONVERGENCE"
    }
    else{ matrix[8,9] <- nMSE_tau2_metaplusmix}
    
    if(RhatDP_out_HN_U_26>0.05){
      matrix[9,9] <- "NON-CONVERGENCE"
    }
    else{matrix[9,9] <- nMSE_DPjags.bin.tau2_HN_U_26}
    
    if(RhatDP_out_HN_U_51>0.05){
      matrix[10,9] <- "NON-CONVERGENCE"
    }
    else{matrix[10,9] <- nMSE_DPjags.bin.tau2_HN_U_51}
    
    if(RhatDP_out_U_U_26>0.05){
      matrix[11,9] <- "NON-CONVERGENCE"
    }
    else{matrix[11,9] <- nMSE_DPjags.bin.tau2_U_U_26}
    
    if(RhatDP_out_U_G_n>0.05){
      matrix[12,9] <- "NON-CONVERGENCE"
    }
    else{matrix[12,9] <- nMSE_DPjags.bin.tau2_U_G_n}
    
    if(RhatDP_out_U_U_51>0.05){
      matrix[13,9] <- "NON-CONVERGENCE"
    }
    else{matrix[13,9] <- nMSE_DPjags.bin.tau2_U_U_51}
    
    matrix[14,9] <- nMSE_frnormaltau2.bin
    
    if(RhatSN_out_HN > 0.05){
      matrix[15,9] <- "NON-CONVERGENCE"
    }
    else{ matrix[15,9] <- nMSE_SNjags.tau2.bin_HN}
    
    if(RhatT_out_HN > 0.05 ){
      matrix[16,9] <-"NON-CONVERGENCE"}
    else{matrix[16,9] <- nMSE_t_distrjags.tau2.bin_HN}
    
    ######### FOR RELATIVE BIAS TAU
    
    if(RhatN_outU > 0.05 ){
      matrix[2,10] <-"NON-CONVERGENCE"}
    else{matrix[2,10] <-avg_rel_bias_tauNormaljags.binU}
    
    if(RhatN_outHN > 0.05 ){
      matrix[3,10] <-"NON-CONVERGENCE"}
    else{matrix[3,10] <-avg_rel_bias_tauNormaljags.binHN}
    
    matrix[4,10] <- avg_rel_bias_tau_frnormal
    
    if(RhatT_out_U > 0.05 ){
      matrix[5,10] <-"NON-CONVERGENCE"}
    else{matrix[5,10] <- avg_rel_bias_taut_distrjags.bin_U}
    
    if(data_out_t > 0.05){
      matrix[6,10] <- "NON-CONVERGENCE"
    }
    else{matrix[6,10] <- avg_rel_bias_tau_metaplust}
    
    if(RhatSN_out_U > 0.05){
      matrix[7,10] <- "NON-CONVERGENCE"
    }
    else{ matrix[7,10] <- avg_rel_bias_tauSNjags.bin_U}
    
    if(data_out_mix > 0.05){
      matrix[8,10] <- "NON-CONVERGENCE"
    }
    else{ matrix[8,10] <- avg_rel_bias_tau_metaplusmix}
    
    if(RhatDP_out_HN_U_26>0.05){
      matrix[9,10] <- "NON-CONVERGENCE"
    }
    else{matrix[9,10] <- avg_rel_bias_tauDPjags.bin_HN_U_26}
    
    if(RhatDP_out_HN_U_51>0.05){
      matrix[10,10] <- "NON-CONVERGENCE"
    }
    else{matrix[10,10] <- avg_rel_bias_tauDPjags.bin_HN_U_51}
    
    if(RhatDP_out_U_U_26>0.05){
      matrix[11,10] <- "NON-CONVERGENCE"
    }
    else{matrix[11,10] <- avg_rel_bias_tauDPjags.bin_U_U_26}
    
    if(RhatDP_out_U_G_n>0.05){
      matrix[12,10] <- "NON-CONVERGENCE"
    }
    else{matrix[12,10] <- avg_rel_bias_tauDPjags.bin_U_G_n}
    
    if(RhatDP_out_U_U_51>0.05){
      matrix[13,10] <- "NON-CONVERGENCE"
    }
    else{matrix[13,10] <- avg_rel_bias_tauDPjags.bin_U_U_51}
    
    matrix[14,10] <- avg_rel_bias_tau_frnormal.bin
    
    if(RhatSN_out_HN > 0.05){
      matrix[15,10] <- "NON-CONVERGENCE"
    }
    else{ matrix[15,10] <- avg_rel_bias_tauSNjags.bin_HN}
    
    if(RhatT_out_HN > 0.05 ){
      matrix[16,10] <-"NON-CONVERGENCE"}
    else{matrix[16,10] <- avg_rel_bias_taut_distrjags.bin_HN}
    
    ######### FOR MSE TAU
    
    if(RhatN_outU > 0.05 ){
      matrix[2,11] <-"NON-CONVERGENCE"}
    else{matrix[2,11] <-avg_MSE_Normaljags.tau.binU}
    
    if(RhatN_outHN > 0.05 ){
      matrix[3,11] <-"NON-CONVERGENCE"}
    else{matrix[3,11] <-avg_MSE_Normaljags.tau.binHN}
    
    matrix[4,11] <- avg_MSE_frnormaltau
    
    if(RhatT_out_U > 0.05 ){
      matrix[5,11] <-"NON-CONVERGENCE"}
    else{matrix[5,11] <- avg_MSE_t_distrjags.tau.bin_U}
    
    if(data_out_t > 0.05){
      matrix[6,11] <- "NON-CONVERGENCE"
    }
    else{matrix[6,11] <- avg_MSE_tau_metaplust}
    
    if(RhatSN_out_U > 0.05){
      matrix[7,11] <- "NON-CONVERGENCE"
    }
    else{ matrix[7,11] <- avg_MSE_SNjags.tau.bin_U}
    
    if(data_out_mix > 0.05){
      matrix[8,11] <- "NON-CONVERGENCE"
    }
    else{ matrix[8,11] <- avg_MSE_tau_metaplusmix}
    
    if(RhatDP_out_HN_U_26>0.05){
      matrix[9,11] <- "NON-CONVERGENCE"
    }
    else{matrix[9,11] <- avg_MSE_DPjags.bin.tau_HN_U_26}
    
    if(RhatDP_out_HN_U_51>0.05){
      matrix[10,11] <- "NON-CONVERGENCE"
    }
    else{matrix[10,11] <- avg_MSE_DPjags.bin.tau_HN_U_51}
    
    if(RhatDP_out_U_U_26>0.05){
      matrix[11,11] <- "NON-CONVERGENCE"
    }
    else{matrix[11,11] <- avg_MSE_DPjags.bin.tau_U_U_26}
    
    if(RhatDP_out_U_G_n>0.05){
      matrix[12,11] <- "NON-CONVERGENCE"
    }
    else{matrix[12,11] <- avg_MSE_DPjags.bin.tau_U_G_n}
    
    if(RhatDP_out_U_U_51>0.05){
      matrix[13,11] <- "NON-CONVERGENCE"
    }
    else{matrix[13,11] <- avg_MSE_DPjags.bin.tau_U_U_51}
    
    matrix[14,11] <- avg_MSE_frnormaltau.bin
    
    if(RhatSN_out_HN > 0.05){
      matrix[15,11] <- "NON-CONVERGENCE"
    }
    else{ matrix[15,11] <- avg_MSE_SNjags.tau.bin_HN}
    
    if(RhatT_out_HN > 0.05 ){
      matrix[16,11] <-"NON-CONVERGENCE"}
    else{matrix[16,11] <- avg_MSE_t_distrjags.tau.bin_HN}
    
    ######### FOR COVERAGE TAU
    
    if(RhatN_outU > 0.05 ){
      matrix[2,12] <-"NON-CONVERGENCE"}
    else{matrix[2,12] <-cover_tauNormaljags.binU}
    
    if(RhatN_outHN > 0.05 ){
      matrix[3,12] <-"NON-CONVERGENCE"}
    else{matrix[3,12] <-cover_tauNormaljags.binHN}
    
    matrix[4,12] <- cover_frnormal_tau
    
    if(RhatT_out_U > 0.05 ){
      matrix[5,12] <-"NON-CONVERGENCE"}
    else{matrix[5,12] <- cover_taut_distrjags.bin_U}
    
    matrix[6,12] <- "NA"
    
    if(RhatSN_out_U > 0.05){
      matrix[7,12] <- "NON-CONVERGENCE"
    }
    else{ matrix[7,12] <- cover_tauSNjags.bin_U}
    
    matrix[8,12] <- "NA"
    
    if(RhatDP_out_HN_U_26>0.05){
      matrix[9,12] <- "NON-CONVERGENCE"
    }
    else{matrix[9,12] <- cover_tauDPjags.bin_HN_U_26}
    
    if(RhatDP_out_HN_U_51>0.05){
      matrix[10,12] <- "NON-CONVERGENCE"
    }
    else{matrix[10,12] <- cover_tauDPjags.bin_HN_U_51}
    
    if(RhatDP_out_U_U_26>0.05){
      matrix[11,12] <- "NON-CONVERGENCE"
    }
    else{matrix[11,12] <- cover_tauDPjags.bin_U_U_26}
    
    if(RhatDP_out_U_G_n>0.05){
      matrix[12,12] <- "NON-CONVERGENCE"
    }
    else{matrix[12,12] <- cover_tauDPjags.bin_U_G_n}
    
    if(RhatDP_out_U_U_51>0.05){
      matrix[13,12] <- "NON-CONVERGENCE"
    }
    else{matrix[13,12] <- cover_tauDPjags.bin_U_U_51}
    
    matrix[14,12] <- "NA"
    
    if(RhatSN_out_HN > 0.05){
      matrix[15,12] <- "NON-CONVERGENCE"
    }
    else{ matrix[15,12] <- cover_tauSNjags.bin_HN}
    
    if(RhatT_out_HN > 0.05 ){
      matrix[16,12] <-"NON-CONVERGENCE"}
    else{matrix[16,12] <- cover_taut_distrjags.bin_HN}
    
    ######### FOR prel_bias tau
    
    if(RhatN_outU > 0.05 ){
      matrix[2,13] <-"NON-CONVERGENCE"}
    else{matrix[2,13] <- prel_bias_tauNormaljags.binU}
    
    if(RhatN_outHN > 0.05 ){
      matrix[3,13] <-"NON-CONVERGENCE"}
    else{matrix[3,13] <-prel_bias_tauNormaljags.binHN}
    
    matrix[4,13] <- prel_bias_tau_frnormal
    
    if(RhatT_out_U > 0.05 ){
      matrix[5,13] <-"NON-CONVERGENCE"}
    else{matrix[5,13] <- prel_bias_taut_distrjags.bin_U}
    
    if(data_out_t > 0.05){
      matrix[6,13] <- "NON-CONVERGENCE"
    }
    else{matrix[6,13] <- prel_bias_tau_metaplust}
    
    if(RhatSN_out_U > 0.05){
      matrix[7,13] <- "NON-CONVERGENCE"
    }
    else{ matrix[7,13] <- prel_bias_tauSNjags.bin_U}
    
    if(data_out_mix > 0.05){
      matrix[8,13] <- "NON-CONVERGENCE"
    }
    else{ matrix[8,13] <- prel_bias_tau_metaplusmix}
    
    if(RhatDP_out_HN_U_26>0.05){
      matrix[9,13] <- "NON-CONVERGENCE"
    }
    else{matrix[9,13] <- prel_bias_tauDPjags.bin_HN_U_26}
    
    if(RhatDP_out_HN_U_51>0.05){
      matrix[10,13] <- "NON-CONVERGENCE"
    }
    else{matrix[10,13] <- prel_bias_tauDPjags.bin_HN_U_51}
    
    if(RhatDP_out_U_U_26>0.05){
      matrix[11,13] <- "NON-CONVERGENCE"
    }
    else{matrix[11,13] <- prel_bias_tauDPjags.bin_U_U_26}
    
    if(RhatDP_out_U_G_n>0.05){
      matrix[12,13] <- "NON-CONVERGENCE"
    }
    else{matrix[12,13] <- prel_bias_tauDPjags.bin_U_G_n}
    
    if(RhatDP_out_U_U_51>0.05){
      matrix[13,13] <- "NON-CONVERGENCE"
    }
    else{matrix[13,13] <- prel_bias_tauDPjags.bin_U_U_51}
    
    matrix[14,13] <- prel_bias_tau_frnormal.bin
    
    if(RhatSN_out_HN > 0.05){
      matrix[15,13] <- "NON-CONVERGENCE"
    }
    else{ matrix[15,13] <- prel_bias_tauSNjags.bin_HN}
    
    if(RhatT_out_HN > 0.05 ){
      matrix[16,13] <-"NON-CONVERGENCE"}
    else{matrix[16,13] <- prel_bias_taut_distrjags.bin_HN}
    
    ######### FOR nMSE of tau
    
    if(RhatN_outU > 0.05 ){
      matrix[2,14] <-"NON-CONVERGENCE"}
    else{matrix[2,14] <-nMSE_Normaljags.tau.binU}
    
    if(RhatN_outHN > 0.05 ){
      matrix[3,14] <-"NON-CONVERGENCE"}
    else{matrix[3,14] <-nMSE_Normaljags.tau.binHN}
    
    matrix[4,14] <- nMSE_frnormaltau
    
    if(RhatT_out_U > 0.05 ){
      matrix[5,14] <-"NON-CONVERGENCE"}
    else{matrix[5,14] <- nMSE_t_distrjags.tau.bin_U}
    
    if(data_out_t > 0.05){
      matrix[6,14] <- "NON-CONVERGENCE"
    }
    else{matrix[6,14] <- nMSE_tau_metaplust}
    
    if(RhatSN_out_U > 0.05){
      matrix[7,14] <- "NON-CONVERGENCE"
    }
    else{ matrix[7,14] <- nMSE_SNjags.tau.bin_U}
    
    if(data_out_mix > 0.05){
      matrix[8,14] <- "NON-CONVERGENCE"
    }
    else{ matrix[8,14] <- nMSE_tau_metaplusmix}
    
    if(RhatDP_out_HN_U_26>0.05){
      matrix[9,14] <- "NON-CONVERGENCE"
    }
    else{matrix[9,14] <- nMSE_DPjags.bin.tau_HN_U_26}
    
    if(RhatDP_out_HN_U_51>0.05){
      matrix[10,14] <- "NON-CONVERGENCE"
    }
    else{matrix[10,14] <- nMSE_DPjags.bin.tau_HN_U_51}
    
    if(RhatDP_out_U_U_26>0.05){
      matrix[11,14] <- "NON-CONVERGENCE"
    }
    else{matrix[11,14] <- nMSE_DPjags.bin.tau_U_U_26}
    
    if(RhatDP_out_U_G_n>0.05){
      matrix[12,14] <- "NON-CONVERGENCE"
    }
    else{matrix[12,14] <- nMSE_DPjags.bin.tau_U_G_n}
    
    if(RhatDP_out_U_U_51>0.05){
      matrix[13,14] <- "NON-CONVERGENCE"
    }
    else{matrix[13,14] <- nMSE_DPjags.bin.tau_U_U_51}
    
    matrix[14,14] <- nMSE_frnormaltau.bin
    
    if(RhatSN_out_HN > 0.05){
      matrix[15,14] <- "NON-CONVERGENCE"
    }
    else{ matrix[15,14] <- nMSE_SNjags.tau.bin_HN}
    
    if(RhatT_out_HN > 0.05 ){
      matrix[16,14] <-"NON-CONVERGENCE"}
    else{matrix[16,14] <- nMSE_t_distrjags.tau.bin_HN}
    
    ######## FOR ABSOLUTE BIAS TAU
    
    if(RhatN_outU > 0.05 ){
      matrix[2,15] <-"NON-CONVERGENCE"}
    else{matrix[2,15] <-avg_abs_bias_tauNormaljags.binU}
    
    if(RhatN_outHN > 0.05 ){
      matrix[3,15] <-"NON-CONVERGENCE"}
    else{matrix[3,15] <-avg_abs_bias_tauNormaljags.binHN}
    
    matrix[4,15] <- avg_abs_bias_tau_frnormal
    
    if(RhatT_out_U > 0.05 ){
      matrix[5,15] <-"NON-CONVERGENCE"}
    else{matrix[5,15] <- avg_abs_bias_taut_distrjags.bin_U}
    
    if(data_out_t > 0.05){
      matrix[6,15] <- "NON-CONVERGENCE"
    }
    else{matrix[6,15] <- avg_abs_bias_tau_metaplust}
    
    if(RhatSN_out_U > 0.05){
      matrix[7,15] <- "NON-CONVERGENCE"
    }
    else{ matrix[7,15] <- avg_abs_bias_tauSNjags.bin_U}
    
    if(data_out_mix > 0.05){
      matrix[8,15] <- "NON-CONVERGENCE"
    }
    else{ matrix[8,15] <- avg_abs_bias_tau_metaplusmix}
    
    if(RhatDP_out_HN_U_26>0.05){
      matrix[9,15] <- "NON-CONVERGENCE"
    }
    else{matrix[9,15] <- avg_abs_bias_tauDPjags.bin_HN_U_26}
    
    if(RhatDP_out_HN_U_51>0.05){
      matrix[10,15] <- "NON-CONVERGENCE"
    }
    else{matrix[10,15] <- avg_abs_bias_tauDPjags.bin_HN_U_51}
    
    if(RhatDP_out_U_U_26>0.05){
      matrix[11,15] <- "NON-CONVERGENCE"
    }
    else{matrix[11,15] <- avg_abs_bias_tauDPjags.bin_U_U_26}
    
    if(RhatDP_out_U_G_n>0.05){
      matrix[12,15] <- "NON-CONVERGENCE"
    }
    else{matrix[12,15] <- avg_abs_bias_tauDPjags.bin_U_G_n}
    
    if(RhatDP_out_U_U_51>0.05){
      matrix[13,15] <- "NON-CONVERGENCE"
    }
    else{matrix[13,15] <- avg_abs_bias_tauDPjags.bin_U_U_51}
    
    matrix[14,15] <- avg_abs_bias_tau_frnormal.bin
    
    if(RhatSN_out_HN > 0.05){
      matrix[15,15] <- "NON-CONVERGENCE"
    }
    else{ matrix[15,15] <- avg_abs_bias_tauSNjags.bin_HN}
    
    if(RhatT_out_HN > 0.05 ){
      matrix[16,15] <-"NON-CONVERGENCE"}
    else{matrix[16,15] <- avg_abs_bias_taut_distrjags.bin_HN}
    
    ######### FOR pabs_bias tau
    
    if(RhatN_outU > 0.05 ){
      matrix[2,16] <-"NON-CONVERGENCE"}
    else{matrix[2,16] <- pabs_bias_tauNormaljags.binU}
    
    if(RhatN_outHN > 0.05 ){
      matrix[3,16] <-"NON-CONVERGENCE"}
    else{matrix[3,16] <-pabs_bias_tauNormaljags.binHN}
    
    matrix[4,16] <- pabs_bias_tau_frnormal
    
    if(RhatT_out_U > 0.05 ){
      matrix[5,16] <-"NON-CONVERGENCE"}
    else{matrix[5,16] <- pabs_bias_taut_distrjags.bin_U}
    
    if(data_out_t > 0.05){
      matrix[6,16] <- "NON-CONVERGENCE"
    }
    else{matrix[6,16] <- pabs_bias_tau_metaplust}
    
    if(RhatSN_out_U > 0.05){
      matrix[7,16] <- "NON-CONVERGENCE"
    }
    else{ matrix[7,16] <- pabs_bias_tauSNjags.bin_U}
    
    if(data_out_mix > 0.05){
      matrix[8,16] <- "NON-CONVERGENCE"
    }
    else{ matrix[8,16] <- pabs_bias_tau_metaplusmix}
    
    if(RhatDP_out_HN_U_26>0.05){
      matrix[9,16] <- "NON-CONVERGENCE"
    }
    else{matrix[9,16] <- pabs_bias_tauDPjags.bin_HN_U_26}
    
    if(RhatDP_out_HN_U_51>0.05){
      matrix[10,16] <- "NON-CONVERGENCE"
    }
    else{matrix[10,16] <- pabs_bias_tauDPjags.bin_HN_U_51}
    
    if(RhatDP_out_U_U_26>0.05){
      matrix[11,16] <- "NON-CONVERGENCE"
    }
    else{matrix[11,16] <- pabs_bias_tauDPjags.bin_U_U_26}
    
    if(RhatDP_out_U_G_n>0.05){
      matrix[12,16] <- "NON-CONVERGENCE"
    }
    else{matrix[12,16] <- pabs_bias_tauDPjags.bin_U_G_n}
    
    if(RhatDP_out_U_U_51>0.05){
      matrix[13,16] <- "NON-CONVERGENCE"
    }
    else{matrix[13,16] <- pabs_bias_tauDPjags.bin_U_U_51}
    
    matrix[14,16] <- pabs_bias_tau_frnormal.bin
    
    if(RhatSN_out_HN > 0.05){
      matrix[15,16] <- "NON-CONVERGENCE"
    }
    else{ matrix[15,16] <- pabs_bias_tauSNjags.bin_HN}
    
    if(RhatT_out_HN > 0.05 ){
      matrix[16,16] <-"NON-CONVERGENCE"}
    else{matrix[16,16] <- pabs_bias_taut_distrjags.bin_HN}
    
    #### FOR RELATIVE BIAS MU 
    if(RhatN_outU > 0.05 ){
      matrix[2,17] <-"NON-CONVERGENCE"}
    else{matrix[2,17] <-avg_rel_biasNormaljags.binU}
    
    if(RhatN_outHN > 0.05 ){
      matrix[3,17] <-"NON-CONVERGENCE"}
    else{matrix[3,17] <-avg_rel_biasNormaljags.binHN}
    
    matrix[4,17] <-avg_rel_bias_frnormal
    
    if(RhatT_out_U > 0.05 ){
      matrix[5,17] <-"NON-CONVERGENCE"}
    else{matrix[5,17] <- avg_rel_biast_distrjags.bin_U}
    
    if(data_out_t > 0.05){
      matrix[6,17] <- "NON-CONVERGENCE"
    }
    else{matrix[6,17] <- avg_rel_bias_metaplust}
    
    if(RhatSN_out_U > 0.05){
      matrix[7,17] <- "NON-CONVERGENCE"
    }
    else{ matrix[7,17] <- avg_rel_biasSNjags.bin_U}
    
    if(data_out_mix > 0.05){
      matrix[8,17] <- "NON-CONVERGENCE"
    }
    else{ matrix[8,17] <- avg_rel_bias_metaplusmix}
    
    if(RhatDP_out_HN_U_26>0.05){
      matrix[9,17] <- "NON-CONVERGENCE"
    }
    else{matrix[9,17] <- avg_rel_biasDPjags.bin_HN_U_26}
    
    if(RhatDP_out_HN_U_51>0.05){
      matrix[10,17] <- "NON-CONVERGENCE"
    }
    else{matrix[10,17] <- avg_rel_biasDPjags.bin_HN_U_51}
    
    if(RhatDP_out_U_U_26>0.05){
      matrix[11,17] <- "NON-CONVERGENCE"
    }
    else{matrix[11,17] <- avg_rel_biasDPjags.bin_U_U_26}
    
    if(RhatDP_out_U_G_n>0.05){
      matrix[12,17] <- "NON-CONVERGENCE"
    }
    else{matrix[12,17] <- avg_rel_biasDPjags.bin_U_G_n}
    
    if(RhatDP_out_U_U_51>0.05){
      matrix[13,17] <- "NON-CONVERGENCE"
    }
    else{matrix[13,17] <- avg_rel_biasDPjags.bin_U_U_51}
    
    matrix[14,17] <-avg_rel_bias_frnormal.bin
    
    if(RhatSN_out_HN > 0.05){
      matrix[15,17] <- "NON-CONVERGENCE"
    }
    else{ matrix[15,17] <- avg_rel_biasSNjags.bin_HN}
    
    if(RhatT_out_HN > 0.05 ){
      matrix[16,17] <-"NON-CONVERGENCE"}
    else{matrix[16,17] <- avg_rel_biast_distrjags.bin_HN}
    
  }
}

write.csv(matrix, "C:\\Users\\Tianqi YU\\Desktop\\SIMULATION SCENARIOS FINAL\\NORMAL DATA\\ALL MODELS RESULTS\\rev_Matrix_All_models_Normal_bin_data_k=14_mu=0_tau2=0.12.csv" ) 


################### summary of the results #####################
r = 13
c = 3
matrix1 = matrix(nrow = r, ncol = c)

# Assign values to each cell in the matrix
for(i in 1:r){
  for(j in 1:c){
    matrix1[1,1] <- "Models"
    
    matrix1[1,2] <- "mean_of_mean_abs_bias_mu" 
    
    matrix1[1,3] <- "mean_of_mean_rel_bias_mu" 
    
    matrix1[2,1] <- "Binomial-Normal(Unif)"
    matrix1[3,1] <- "Binomial-Normal(HN)"
    matrix1[4,1] <- "Normal-Normal"
    matrix1[5,1] <- "Binomial-t(Unif)"
    matrix1[6,1] <- "Binomial-SN(Unif)"
    matrix1[7,1] <- "Binomial-DP-26(HN/Unif)"
    matrix1[8,1] <- "Binomial-DP-51(HN/Unif)"
    matrix1[9,1] <- "Binomial-DP-26(Unif/Unif)"
    matrix1[10,1] <- "Binomial-DP-n(Unif/Gamma)"
    matrix1[11,1] <- "Binomial-DP-51(Unif/Unif)"
    matrix1[12,1] <- "Binomial-SN(HN)"
    matrix1[13,1] <- "Binomial-t(HN)"
    #### FOR BIAS OF REL_EFF ###
    
    if(RhatN_outU > 0.05 ){
      matrix1[2,2] <-"NON-CONVERGENCE"}
    else{matrix1[2,2] <-mean_abs_bias_of_mean_abs_bias_relNU}
    
    if(RhatN_outHN > 0.05 ){
      matrix1[3,2] <-"NON-CONVERGENCE"}
    else{matrix1[3,2] <-mean_abs_bias_of_mean_abs_bias_relNHN}
    
    matrix1[4,2] <-mean_abs_bias_of_mean_abs_bias_relnormfr_distr
    
    if(RhatT_out_U > 0.05 ){
      matrix1[5,2] <-"NON-CONVERGENCE"}
    else{matrix1[5,2] <-mean_abs_bias_of_mean_abs_bias_relt_distr_U}
    
    if(RhatSN_out_U > 0.05){
      matrix1[6,2] <- "NON-CONVERGENCE"
    }
    else{ matrix1[6,2] <- mean_abs_bias_of_mean_abs_bias_relSN_U}
    
    if(RhatDP_out_HN_U_26>0.05){
      matrix1[7,2] <- "NON-CONVERGENCE"
    }
    else{matrix1[7,2] <- mean_abs_bias_of_mean_abs_bias_relDP_HN_U_26}
    
    if(RhatDP_out_HN_U_51>0.05){
      matrix1[8,2] <- "NON-CONVERGENCE"
    }
    else{matrix1[8,2] <- mean_abs_bias_of_mean_abs_bias_relDP_HN_U_51}
    
    if(RhatDP_out_U_U_26>0.05){
      matrix1[9,2] <- "NON-CONVERGENCE"
    }
    else{matrix1[9,2] <- mean_abs_bias_of_mean_abs_bias_relDP_U_U_26}
    
    if(RhatDP_out_U_G_n>0.05){
      matrix1[10,2] <- "NON-CONVERGENCE"
    }
    else{matrix1[10,2] <- mean_abs_bias_of_mean_abs_bias_relDP_U_G_n}
    
    if(RhatDP_out_U_U_51>0.05){
      matrix1[11,2] <- "NON-CONVERGENCE"
    }
    else{matrix1[11,2] <- mean_abs_bias_of_mean_abs_bias_relDP_U_U_51}
    
    
    if(RhatSN_out_HN > 0.05){
      matrix1[12,2] <- "NON-CONVERGENCE"
    }
    else{ matrix1[12,2] <- mean_abs_bias_of_mean_abs_bias_relSN_HN}
    
    if(RhatT_out_HN > 0.05 ){
      matrix1[13,2] <-"NON-CONVERGENCE"}
    else{matrix1[13,2] <-mean_abs_bias_of_mean_abs_bias_relt_distr_HN}
    
    #### FOR BIAS OF REL_EFF ###
    
    if(RhatN_outU > 0.05 ){
      matrix1[2,3] <-"NON-CONVERGENCE"}
    else{matrix1[2,3] <-mean_rel_bias_of_mean_rel_bias_relNU}
    
    if(RhatN_outHN > 0.05 ){
      matrix1[3,3] <-"NON-CONVERGENCE"}
    else{matrix1[3,3] <-mean_rel_bias_of_mean_rel_bias_relNHN}
    
    matrix1[4,3] <-mean_rel_bias_of_mean_rel_bias_relnormfr_distr
    
    if(RhatT_out_U > 0.05 ){
      matrix1[5,3] <-"NON-CONVERGENCE"}
    else{matrix1[5,3] <-mean_rel_bias_of_mean_rel_bias_relt_distr_U}
    
    if(RhatSN_out_U > 0.05){
      matrix1[6,3] <- "NON-CONVERGENCE"
    }
    else{ matrix1[6,3] <- mean_rel_bias_of_mean_rel_bias_relSN_U}
    
    if(RhatDP_out_HN_U_26>0.05){
      matrix1[7,3] <- "NON-CONVERGENCE"
    }
    else{matrix1[7,3] <- mean_rel_bias_of_mean_rel_bias_relDP_HN_U_26}
    
    if(RhatDP_out_HN_U_51>0.05){
      matrix1[8,3] <- "NON-CONVERGENCE"
    }
    else{matrix1[8,3] <- mean_rel_bias_of_mean_rel_bias_relDP_HN_U_51}
    
    if(RhatDP_out_U_U_26>0.05){
      matrix1[9,3] <- "NON-CONVERGENCE"
    }
    else{matrix1[9,3] <- mean_rel_bias_of_mean_rel_bias_relDP_U_U_26}
    
    if(RhatDP_out_U_G_n>0.05){
      matrix1[10,3] <- "NON-CONVERGENCE"
    }
    else{matrix1[10,3] <- mean_rel_bias_of_mean_rel_bias_relDP_U_G_n}
    
    if(RhatDP_out_U_U_51>0.05){
      matrix1[11,3] <- "NON-CONVERGENCE"
    }
    else{matrix1[11,3] <- mean_rel_bias_of_mean_rel_bias_relDP_U_U_51}
    
    
    if(RhatSN_out_HN > 0.05){
      matrix1[12,3] <- "NON-CONVERGENCE"
    }
    else{ matrix1[12,3] <- mean_rel_bias_of_mean_rel_bias_relSN_HN}
    
    if(RhatT_out_HN > 0.05 ){
      matrix1[13,3] <-"NON-CONVERGENCE"}
    else{matrix1[13,3] <-mean_rel_bias_of_mean_rel_bias_relt_distr_HN}
    
  }
}

write.csv(matrix1, "C:\\Users\\Tianqi YU\\Desktop\\SIMULATION SCENARIOS FINAL\\NORMAL DATA\\ALL MODELS RESULTS\\rev_Matrix_Bias_Rel_eff_Normal_data_k=14_mu=0_tau2=0.12.csv" ) 

