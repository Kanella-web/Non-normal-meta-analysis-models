############################# Scenario 6:Normal data k=26, mu=0, tau2=2.63 ######################
### DATA GENERATING PROCESS ####
N.sim = 1000

df <- list()

dich_data <- function(k, tau2, mu){
  
  mati <- matrix(1:k,nrow = k, ncol=2) 
  
  for (i in 1:N.sim){
    
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
set.seed(25)
m <- dich_data(26,2.63,0) 
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
  colnames(p[[i]]) = c("ci", "ti","nci","nti","true_eff")
  
}

#### CHECK IF THERE ARE STUDIES WITH 0 IN TWO ARMS ####

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
#d

count <- c()
for(i in 1:N.sim){
  count[i] <- sum(d[[i]])
}
count
studies_zero_2arms <- which(count == 1)

# ### DENSITY FUNCTION OF THE 26 TRUE STUDY-SPECIFIC EFFECTS ####
### SELECTED DATA SET ###
plot(density(p[[4]]$true_eff),main = "Normal distribution for true treatment effects", lwd= 3, xlim = c(-6,6))
write.csv(p[[4]]$true_eff , "C:\\Users\\Lela Panag\\Desktop\\Github\\GITHUB_FINAL\\Forest plots of 3 selected simulated datasets\\true_eff_normal.csv",row.names=FALSE )  

################ BAYESIAN M-A ################
library(meta)
library(metafor)
library(rjags)
library(R2jags)

############ Binomial-Normal(HN) model ##############
set.seed(1234600)
writeLines("
  model{
  for(i in 1:ns){ 
  m[i] ~ dnorm(0,.0001) 
  
  r1[i] ~ dbin(p1[i], n1[i])
  logit(p1[i]) <- m[i] 
  
  r2[i] ~ dbin(p2[i], n2[i])
  logit(p2[i]) <- m[i] + delta12[i]
  
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

dat<-list(ns = nrow(p[[4]]),
          r1 = p[[4]]$ci,
          r2 = p[[4]]$ti,
          n1 = p[[4]]$nci,
          n2 = p[[4]]$nti )

run.modelN = jags(
  data = dat,
  inits = NULL,
  parameters.to.save = c(
    "mu",  ## mean treatment effect
    "tau", ## between-study sd
    "delta12", ## study-specific effects
    "delta_new" ## prediction
  ),
  n.chains = 2,
  n.iter = 50000,
  
  n.burnin = 10000,
  DIC = T,
  model.file = modfile
)

resultsN <- as.data.frame(run.modelN$BUGSoutput$summary) 


#### The density plot of pred
plot(density(run.modelN$BUGSoutput$sims.matrix[  ,"delta_new"]))

### Pr of pred <- 0
Pr_mu = mean(run.modelN$BUGSoutput$sims.matrix[  ,"mu"] < 0 )
Pr_pred = mean(run.modelN$BUGSoutput$sims.matrix[  ,"delta_new"] < 0 )

f1 <- grepl("mu", row.names(resultsN))
m1 <- resultsN[f1,]
mu_Normaljags.bin <- m1$`50%` 
LBmu.jags.bin <- m1$`2.5%`
UBmu.jags.bin <- m1$`97.5%`
Rhat_muN <- m1$Rhat
precN_mu <-  UBmu.jags.bin - LBmu.jags.bin

fdn <- grepl("delta_new", row.names(resultsN))
mdn <- resultsN[fdn,]
delta_new_Normaljags.bin <- mdn$`50%` 
LBdelta_new.jags.bin <- mdn$`2.5%`
UBdelta_new.jags.bin <- mdn$`97.5%`
Rhat_delta_newN <- mdn$Rhat
precN_delta_new <-  UBdelta_new.jags.bin - LBdelta_new.jags.bin

f2 <- grepl("tau", row.names(resultsN))
m2 <- resultsN[f2,]
tau_Normaljags.bin <- m2$`50%`
LBtau.jags.bin <- m2$`2.5%`
UBtau.jags.bin <- m2$`97.5%`
tau2_Normaljags.bin <- (m2$`50%`)^2
LBtau2.jags.bin <- (m2$`2.5%`)^2
UBtau2.jags.bin <- (m2$`97.5%`)^2
Rhat_tauN <- m2$Rhat
precN_tau2 <-  UBtau2.jags.bin - LBtau2.jags.bin
precN_tau <-  UBtau.jags.bin - LBtau.jags.bin


fd <- grepl("delta12", row.names(resultsN))
md <- resultsN[fd,]
rel_eff <- md$`50%`
LB_rel_eff <- md$`2.5%`
UB_rel_eff <- md$`97.5%`
sd_rel_eff <- md$sd
Rhat_deltaN <- md$Rhat
lista <- cbind.data.frame(rel_eff,sd_rel_eff,LB_rel_eff, UB_rel_eff, Rhat_deltaN)


RhatsN <- cbind.data.frame(mu_Normaljags.bin, LBmu.jags.bin, UBmu.jags.bin, tau2_Normaljags.bin,
                           LBtau2.jags.bin, UBtau2.jags.bin,tau_Normaljags.bin,
                           LBtau.jags.bin, UBtau.jags.bin,
                           Rhat_muN, Rhat_tauN,
                           precN_mu, precN_tau2,precN_tau, delta_new_Normaljags.bin, 
                           LBdelta_new.jags.bin , UBdelta_new.jags.bin, Pr_mu, Pr_pred)


write.csv(RhatsN , "C:\\Users\\kanel\\OneDrive\\Υπολογιστής\\Simulated data application\\Normal_data\\NORMAL\\res_Normal_data_Bayesian_Normal_model_HN.csv",row.names=FALSE )  
write.csv(lista , "C:\\Users\\kanel\\OneDrive\\Υπολογιστής\\Simulated data application\\Normal_data\\NORMAL\\rel_effN_Normal_data_Bayesian_Normal_model_HN.csv",row.names=FALSE )  


############ Binomial-Normal(Unif) model ##############
set.seed(1234601)
library(R2jags)
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

dat<-list(ns = nrow(p[[4]]),
          r1 = p[[4]]$ci,
          r2 = p[[4]]$ti,
          n1 = p[[4]]$nci,
          n2 = p[[4]]$nti )

run.modelN_U = jags(
  data = dat,
  inits = NULL,
  parameters.to.save = c(
    "mu", ## mean treatment effect
    "tau", ## between-study sd
    "delta12", ##study-specific effect
    "delta_new" ## prediction
  ),
  n.chains = 2,
  n.iter = 50000,
  
  n.burnin = 10000,
  DIC = T,
  model.file = modfile
)

resultsNU <- as.data.frame(run.modelN_U$BUGSoutput$summary) 


#### The density plot of pred
plot(density(run.modelN_U$BUGSoutput$sims.matrix[  ,"delta_new"]))

### Pr of pred <- 0
Pr_muU = mean(run.modelN_U$BUGSoutput$sims.matrix[  ,"mu"] < 0 )
Pr_predU = mean(run.modelN_U$BUGSoutput$sims.matrix[  ,"delta_new"] < 0 )


f1 <- grepl("mu", row.names(resultsNU))
m1 <- resultsNU[f1,]
mu_Normaljags.binU <- m1$`50%` 
LBmu.jags.binU <- m1$`2.5%`
UBmu.jags.binU <- m1$`97.5%`
Rhat_muNU <- m1$Rhat
precN_muU <-  UBmu.jags.binU - LBmu.jags.binU

fdn <- grepl("delta_new", row.names(resultsNU))
mdn <- resultsNU[fdn,]
delta_new_Normaljags.binU <- mdn$`50%` 
LBdelta_new.jags.binU <- mdn$`2.5%`
UBdelta_new.jags.binU <- mdn$`97.5%`
Rhat_delta_newNU <- mdn$Rhat
precN_delta_newU <-  UBdelta_new.jags.binU - LBdelta_new.jags.binU

f2 <- grepl("tau", row.names(resultsNU))
m2 <- resultsNU[f2,]
tau_Normaljags.binU <- m2$`50%`
LBtau.jags.binU <- m2$`2.5%`
UBtau.jags.binU <- m2$`97.5%`
tau2_Normaljags.binU <- (m2$`50%`)^2
LBtau2.jags.binU <- (m2$`2.5%`)^2
UBtau2.jags.binU <- (m2$`97.5%`)^2
Rhat_tauNU <- m2$Rhat
precN_tau2U <-  UBtau2.jags.binU - LBtau2.jags.binU
precN_tauU <-  UBtau.jags.binU - LBtau.jags.binU

fd <- grepl("delta12", row.names(resultsNU))
md <- resultsNU[fd,]
rel_effU <- md$`50%`
LB_rel_effU <- md$`2.5%`
UB_rel_effU <- md$`97.5%`
sd_rel_effU <- md$sd
Rhat_deltaNU <- md$Rhat
listaU <- cbind.data.frame(rel_effU,sd_rel_effU,LB_rel_effU, UB_rel_effU, Rhat_deltaNU)


RhatsNU <- cbind.data.frame(mu_Normaljags.binU, LBmu.jags.binU, UBmu.jags.binU, tau2_Normaljags.binU,
                            LBtau2.jags.binU, UBtau2.jags.binU,tau_Normaljags.binU,
                            LBtau.jags.binU, UBtau.jags.binU,
                            Rhat_muNU, Rhat_tauNU,
                            precN_muU, precN_tau2U,precN_tauU, delta_new_Normaljags.binU, 
                            LBdelta_new.jags.binU , UBdelta_new.jags.binU, Pr_muU, Pr_predU)


write.csv(RhatsNU , "C:\\Users\\Lela Panag\\Desktop\\Github\\GITHUB_FINAL\\Forest plots of 3 selected simulated datasets\\res_Normal_data_Bayesian_Normal_model_Unif.csv",row.names=FALSE )  
write.csv(listaU , "C:\\Users\\Lela Panag\\Desktop\\Github\\GITHUB_FINAL\\Forest plots of 3 selected simulated datasets\\rel_effN_Normal_data_Bayesian_Normal_model_Unif.csv",row.names=FALSE )  

####################### Binomial-DP-n(Unif/Gamma) model ######################
set.seed(12345)

cat("model{
  # Random effects logistic regression part of model
  for( i in 1: ns ) {

    m[i] ~ dnorm(0,.0001)
    
    logit(p1[i]) <- m[i]
    r1[i] ~ dbin(p1[i],n1[i])
    
    logit(p2[i]) <- m[i] + delta12[i]
    r2[i] ~ dbin(p2[i],n2[i])
    
    delta12[i] <- theta[Z[i]]
    
    Z[i] ~ dcat(p[]) #Z is an integer variable
  }
  # Constructive DPP
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
  
   for(i in 1:ns){
   rho[i] <- (exp(delta12[i]) / (1+exp(delta12[i])))
  }
  
   for(i in 1:ns){
   for(j in 1:ns){
    equalsmatrix[i,j]<-equals(rho[i],rho[j])
   }
  equalsres[i]<-sum(equalsmatrix[i,])
  }
  
}",file="DPmodel2.bin_U_G_n.txt")
modfile = 'DPmodel2.bin_U_G_n.txt'


run.modelDP = jags(
  data =list(ns = nrow(p[[4]]),
             r1 = p[[4]]$ci,
             r2 = p[[4]]$ti,
             n1 = p[[4]]$nci,
             n2 = p[[4]]$nti,
             N= 26 ## TRUNCATION POINT
  ) ,  
  inits = NULL,
  parameters.to.save = c(
    "basemu",    ## mu of the base Normal distribution
    "basetau",   ## sd of the base Normal distribution
    "poptrue",   ## mean treatment effect
    "var.true", ## between-study variance
    "delta12", #### study-specific effects
    "K",       ### the total number of clusters 
    "p",      ## the weights of the process 
    "alpha",  ## the concentration parameter
    "SC"     ## the probability of each cluster assignment 
  ),   
  
  n.chains = 2,
  n.iter = 50000,
  
  
  n.burnin = 10000,
  DIC = T,
  model.file = modfile
  
)

DPresults <- as.data.frame( run.modelDP$BUGSoutput$summary) 


Pr_mu = mean(run.modelDP$BUGSoutput$sims.matrix[  ,"poptrue"] < 0 )


DP_results <-DPresults


DPresults$ind <- row.names(DPresults)

############# extraction of the parameters of interest ###########

f55sc <- grepl("SC", row.names(DPresults))
m55sc <- DPresults[f55sc,]
SC <- m55sc$mean

#### BASED ON  THE SC[i,j] ####
##### Distribution of data points to clusters ###

extra_col = rownames(m55sc)
m55sc$names = extra_col

prob_list <- data.frame(Column10 = m55sc[, 10], Column1 = m55sc[, 1])

split_dataframe <- function(df, chunk_size) {
  split(df, rep(1:ceiling(nrow(df) / chunk_size), each = chunk_size, length.out = nrow(df)))
}

split_prob_list <- split_dataframe(prob_list,nrow(p[[4]]) ) 

split_prob_list

# Function to find the maximum probabilities for data points to be distributed into clusters ####
find_max_values_and_indices <- function(split_list) {
  num_rows <- nrow(split_list[[1]])
  num_sublists <- length(split_list)
  
  max_values <- numeric(num_rows)
  max_indices <- integer(num_rows)
  
  for (i in 1:num_rows) {
    values <- sapply(split_list, function(x) x[i, 2])
    max_value <- max(values, na.rm = TRUE)
    max_index <- which(values == max_value)[1]  
    max_values[i] <- max_value
    max_indices[i] <- max_index
  }
  
  return(list(max_values = max_values, max_indices = max_indices))
}

result <- find_max_values_and_indices(split_prob_list)
max_values <- result$max_values #### the probabilities of cluster assignment
max_indices <- result$max_indices #### the cluster assignments 

### clusters with data points ###
k <- unique(max_indices)

positions <- vector("list", length(k))
names(positions) <- k

for(value in k) {
  positions[[as.character(value)]] <- which(max_indices == value)
}


### clusters with data points ###

cluster1 = positions[["1"]]
cluster2 = positions[["2"]]
cluster3 = positions[["3"]] 
cluster7 = positions[["7"]] 


f11 <- grepl("poptrue", row.names(DP_results))
m11 <- DP_results[f11,]
mu_DP<- m11$`50%`
LB_mu_DP <- m11$`2.5%`
UB_mu_DP <- m11$`97.5%`
Rhat_muDP <- m11$Rhat
precDP_mu <- UB_mu_DP - LB_mu_DP

f22 <- grepl("var.true", row.names(DP_results))
m22 <- DP_results[f22,]
tau2_DP <- m22$`50%`
LB_tau2_DP <- m22$`2.5%`
UB_tau2_DP <- m22$`97.5%`
Rhat_tauDP <- m22$Rhat

tau_DP <- sqrt(tau2_DP)
LB_tau_DP <- sqrt(LB_tau2_DP)
UB_tau_DP <- sqrt(UB_tau2_DP)
precDP_tau2 <- UB_tau2_DP -  LB_tau2_DP
precDP_tau <- UB_tau_DP -  LB_tau_DP


fdd <- grepl("delta12", row.names(DP_results))
mdd <- DP_results[fdd,]
rel_effDP <- mdd$`50%`
LB_rel_effDP <- mdd$`2.5%`
UB_rel_effDP <- mdd$`97.5%`
sd_rel_effDP <- mdd$sd
Rhat_deltaDP <- mdd$Rhat

#### Mean of each cluster ######
cl1 = rel_effDP[cluster1]
mean_cl1 = mean(cl1)
cl2 = rel_effDP[cluster2]
mean_cl2 = mean(cl2)
cl3 = rel_effDP[cluster3]
mean_cl3 = mean(cl3)
cl7 = rel_effDP[cluster7]
mean_cl7 = mean(cl7)

cluster_mean = cbind.data.frame(mean_cl1, mean_cl2, mean_cl3,mean_cl7 )


f33 <- grepl("basemu", row.names(DP_results))
m33 <- DP_results[f33,]
base_mu <- m33$`50%`
LB_base_mu <- m33$`2.5%`
UB_base_mu <- m33$`97.5%`
Rhat_basemu <- m33$Rhat
prec_basemu <- UB_base_mu - LB_base_mu

f44 <- grepl("basetau", row.names(DP_results))
m44 <- DP_results[f44,]
base_tau <- m44$`50%`
base_tau2 <- (m44$`50%`)^2
LB_base_tau2 <- (m44$`2.5%`)^2
UB_base_tau2 <- (m44$`97.5%`)^2
Rhat_basetau <- (m44$Rhat)^2
prec_basetau2 <- UB_base_tau2 - LB_base_tau2

f55 <- grepl("alpha", row.names(DP_results))
m55 <- DP_results[f55,]
alpha <- m55$`50%`
LB_alpha <- m55$`2.5%`
UB_alpha <- m55$`97.5%`
Rhat_alpha <- m55$Rhat
prec_alpha <- UB_alpha - LB_alpha

fcl <- grepl("K", row.names(DP_results))   
mcl <- DP_results[fcl,]
median_K <- mcl$`50%` 
LB_K <- mcl$`2.5%`
UB_K <- mcl$`97.5%`

fp <- grepl("p", row.names(DP_results))
mp <- DP_results[fp,]
mp <- mp[!grepl("pop", row.names(mp)),]
mp <- mp[!grepl("alpha", row.names(mp)),]
pi <- mp$mean


listaDP1 <- cbind.data.frame(rel_effDP,sd_rel_effDP, LB_rel_effDP,UB_rel_effDP,Rhat_deltaDP)


numclus <- unlist(median_K)
LB_K <- unlist(LB_K)
UB_K <- unlist(UB_K)

RhatsDP <- cbind.data.frame(base_mu,LB_base_mu,UB_base_mu,base_tau,base_tau2, LB_base_tau2,UB_base_tau2,
                            mu_DP, LB_mu_DP, UB_mu_DP, tau2_DP, LB_tau2_DP, UB_tau2_DP,
                            tau_DP, LB_tau_DP, UB_tau_DP,
                            precDP_mu, precDP_tau2,  precDP_tau,
                            prec_basemu, prec_basetau2,prec_alpha,  alpha , LB_alpha, UB_alpha,
                            Rhat_muDP, Rhat_tauDP, Rhat_basemu, Rhat_basetau, Rhat_alpha, numclus, LB_K, UB_K, Pr_mu)

extra_col = row.names(DPresults)
DPresults$extra_col = extra_col
prob = round(max_values,2)

write.csv(listaDP1 , "C:\\Users\\Lela Panag\\Desktop\\Github\\GITHUB_FINAL\\Forest plots of 3 selected simulated datasets\\UGDP_rel_eff_Normal_data.csv",row.names=FALSE )  
write.csv(DPresults , "C:\\Users\\Lela Panag\\Desktop\\Github\\GITHUB_FINAL\\Forest plots of 3 selected simulated datasets\\all_UGDP_res_Normal_data.csv",row.names=FALSE )  
write.csv(RhatsDP, "C:\\Users\\Lela Panag\\Desktop\\Github\\GITHUB_FINAL\\Forest plots of 3 selected simulated datasets\\UGDP_res_Normal_data.csv",row.names=FALSE )  
write.csv(prob , "C:\\Users\\Lela Panag\\Desktop\\Github\\GITHUB_FINAL\\Forest plots of 3 selected simulated datasets\\UGDP_prob_Normal_data.csv",row.names=FALSE )  
write.csv( cluster_mean, "C:\\Users\\Lela Panag\\Desktop\\Github\\GITHUB_FINAL\\Forest plots of 3 selected simulated datasets\\UGDP_cluster_mean_Normal_data.csv",row.names=FALSE )  


############# BINOMIAL-DP-26(HN/UNIF) #############
set.seed(123451)

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
  
}",file="DPmodel2.binHN_U_26.txt")
modfile = 'DPmodel2.binHN_U_26.txt'

run.modelDP_HN_U_26 = jags(
  data =list(ns = nrow(p[[4]]),
             r1 = p[[4]]$ci,
             r2 = p[[4]]$ti,
             n1 = p[[4]]$nci,
             n2 = p[[4]]$nti,
             N= 26 ## TRUNCATION POINT
  ) ,  
  inits = NULL,
  parameters.to.save = c(
    "basemu",    ## mu of the base Normal distribution
    "basetau",   ## sd of the base Normal distribution
    "poptrue",   ## mean treatment effect
    "var.true", ## between-study variance 
    "delta12", #### study-specific effects
    "K",       ### the total number of clusters 
    "p",      ## the weights of the process 
    "alpha",  ## the concentration parameter
    "SC"     ## the probability of each cluster assignment 
  ),   
  
  n.chains = 2,
  n.iter = 50000,
  
  
  n.burnin = 10000,
  DIC = T,
  model.file = modfile
  
)

DPresults_HN_U_26 <- as.data.frame( run.modelDP_HN_U_26$BUGSoutput$summary) 


Pr_mu_HN_U_26 = mean(run.modelDP_HN_U_26$BUGSoutput$sims.matrix[  ,"poptrue"] < 0 )


DP_results_HN_U_26 <-DPresults_HN_U_26


DPresults_HN_U_26$ind <- row.names(DPresults_HN_U_26)

############# extraction of the parameters of interest ###########

f55sc <- grepl("SC", row.names(DPresults_HN_U_26))
m55sc <- DPresults_HN_U_26[f55sc,]
SC <- m55sc$mean

#### BASED ON  THE SC[i,j] ####
##### Distribution of data points to clusters ###

extra_col = rownames(m55sc)
m55sc$names = extra_col

prob_list <- data.frame(Column10 = m55sc[, 10], Column1 = m55sc[, 1])

split_dataframe <- function(df, chunk_size) {
  split(df, rep(1:ceiling(nrow(df) / chunk_size), each = chunk_size, length.out = nrow(df)))
}

split_prob_list <- split_dataframe(prob_list, nrow(p[[4]]))  

split_prob_list

# Function to find the maximum probabilities for data points to be distributed into clusters ####
find_max_values_and_indices <- function(split_list) {
  num_rows <- nrow(split_list[[1]])
  num_sublists <- length(split_list)
  
  max_values <- numeric(num_rows)
  max_indices <- integer(num_rows)
  
  for (i in 1:num_rows) {
    values <- sapply(split_list, function(x) x[i, 2])
    max_value <- max(values, na.rm = TRUE)
    max_index <- which(values == max_value)[1]  
    max_values[i] <- max_value
    max_indices[i] <- max_index
  }
  
  return(list(max_values = max_values, max_indices = max_indices))
}

result <- find_max_values_and_indices(split_prob_list)
max_values <- result$max_values #### the probabilities of cluster assignment
max_indices <- result$max_indices #### the cluster assignments 

### clusters with data points ###
k <- unique(max_indices)

positions <- vector("list", length(k))
names(positions) <- k

for(value in k) {
  positions[[as.character(value)]] <- which(max_indices == value)
}


### clusters with data points ###

cluster1_HN_U_26 = positions[["1"]]
cluster2_HN_U_26 = positions[["2"]]
cluster3_HN_U_26 = positions[["3"]] 
cluster4_HN_U_26 = positions[["4"]] 
cluster9_HN_U_26 = positions[["9"]] 


f11 <- grepl("poptrue", row.names(DP_results_HN_U_26))
m11 <- DP_results_HN_U_26[f11,]
mu_DP_HN_U_26<- m11$`50%`
LB_mu_DP_HN_U_26 <- m11$`2.5%`
UB_mu_DP_HN_U_26 <- m11$`97.5%`
Rhat_muDP_HN_U_26 <- m11$Rhat
precDP_mu_HN_U_26 <- UB_mu_DP_HN_U_26 - LB_mu_DP_HN_U_26

f22 <- grepl("var.true", row.names(DP_results_HN_U_26))
m22 <- DP_results_HN_U_26[f22,]
tau2_DP_HN_U_26 <- m22$`50%`
LB_tau2_DP_HN_U_26 <- m22$`2.5%`
UB_tau2_DP_HN_U_26 <- m22$`97.5%`
Rhat_tauDP_HN_U_26 <- m22$Rhat
precDP_tau2_HN_U_26 <- UB_tau2_DP_HN_U_26 -  LB_tau2_DP_HN_U_26

tau_DP_HN_U_26 <- sqrt(tau2_DP_HN_U_26)
LB_tau_DP_HN_U_26 <- sqrt(LB_tau2_DP_HN_U_26)
UB_tau_DP_HN_U_26 <- sqrt(UB_tau2_DP_HN_U_26)
precDP_tau_HN_U_26 <- UB_tau_DP_HN_U_26 -  LB_tau_DP_HN_U_26


fdd <- grepl("delta12", row.names(DP_results_HN_U_26))
mdd <- DP_results_HN_U_26[fdd,]
rel_effDP_HN_U_26 <- mdd$`50%`
LB_rel_effDP_HN_U_26 <- mdd$`2.5%`
UB_rel_effDP_HN_U_26 <- mdd$`97.5%`
sd_rel_effDP_HN_U_26 <- mdd$sd
Rhat_deltaDP_HN_U_26 <- mdd$Rhat

#### Mean of each cluster ######
cl1_HN_U_26 = rel_effDP_HN_U_26[cluster1_HN_U_26]
mean_cl1_HN_U_26 = mean(cl1_HN_U_26)
cl2_HN_U_26 = rel_effDP_HN_U_26[cluster2_HN_U_26]
mean_cl2_HN_U_26 = mean(cl2_HN_U_26)
cl3_HN_U_26 = rel_effDP_HN_U_26[cluster3_HN_U_26]
mean_cl3_HN_U_26 = mean(cl3_HN_U_26)
cl4_HN_U_26 = rel_effDP_HN_U_26[cluster4_HN_U_26]
mean_cl4_HN_U_26 = mean(cl4_HN_U_26)
cl9_HN_U_26 = rel_effDP_HN_U_26[cluster9_HN_U_26]
mean_cl9_HN_U_26 = mean(cl9_HN_U_26)


cluster_mean_HN_U_26 = cbind.data.frame(mean_cl1_HN_U_26, mean_cl2_HN_U_26, mean_cl3_HN_U_26,mean_cl4_HN_U_26, mean_cl9_HN_U_26 )


f33 <- grepl("basemu", row.names(DP_results_HN_U_26))
m33 <- DP_results_HN_U_26[f33,]
base_mu_HN_U_26 <- m33$`50%`
LB_base_mu_HN_U_26 <- m33$`2.5%`
UB_base_mu_HN_U_26 <- m33$`97.5%`
Rhat_basemu_HN_U_26 <- m33$Rhat
prec_basemu_HN_U_26 <- UB_base_mu_HN_U_26 - LB_base_mu_HN_U_26

f44 <- grepl("basetau", row.names(DP_results_HN_U_26))
m44 <- DP_results_HN_U_26[f44,]
base_tau_HN_U_26 <- m44$`50%`
base_tau2_HN_U_26 <- (m44$`50%`)^2
LB_base_tau2_HN_U_26 <- (m44$`2.5%`)^2
UB_base_tau2_HN_U_26 <- (m44$`97.5%`)^2
Rhat_basetau_HN_U_26 <- (m44$Rhat)^2
prec_basetau2_HN_U_26 <- UB_base_tau2_HN_U_26 - LB_base_tau2_HN_U_26

f55 <- grepl("alpha", row.names(DP_results_HN_U_26))
m55 <- DP_results_HN_U_26[f55,]
alpha_HN_U_26 <- m55$`50%`
LB_alpha_HN_U_26 <- m55$`2.5%`
UB_alpha_HN_U_26 <- m55$`97.5%`
Rhat_alpha_HN_U_26 <- m55$Rhat
prec_alpha_HN_U_26 <- UB_alpha_HN_U_26 - LB_alpha_HN_U_26

fcl <- grepl("K", row.names(DP_results_HN_U_26))   
mcl <- DP_results_HN_U_26[fcl,]
median_K_HN_U_26 <- mcl$`50%` 
LB_K_HN_U_26 <- mcl$`2.5%`
UB_K_HN_U_26 <- mcl$`97.5%`

fp <- grepl("p", row.names(DP_results_HN_U_26))
mp <- DP_results_HN_U_26[fp,]
mp <- mp[!grepl("pop", row.names(mp)),]
mp <- mp[!grepl("alpha", row.names(mp)),]
pi_HN_U_26 <- mp$mean


listaDP1_HN_U_26 <- cbind.data.frame(rel_effDP_HN_U_26,sd_rel_effDP_HN_U_26, LB_rel_effDP_HN_U_26,UB_rel_effDP_HN_U_26,Rhat_deltaDP_HN_U_26)


numclus_HN_U_26 <- unlist(median_K_HN_U_26)
LB_K_HN_U_26 <- unlist(LB_K_HN_U_26)
UB_K_HN_U_26 <- unlist(UB_K_HN_U_26)

RhatsDP_HN_U_26 <- cbind.data.frame(base_mu_HN_U_26,LB_base_mu_HN_U_26,UB_base_mu_HN_U_26,
                                    base_tau_HN_U_26,base_tau2_HN_U_26, LB_base_tau2_HN_U_26,UB_base_tau2_HN_U_26,
                                    mu_DP_HN_U_26, LB_mu_DP_HN_U_26, UB_mu_DP_HN_U_26, 
                                    tau2_DP_HN_U_26, LB_tau2_DP_HN_U_26, UB_tau2_DP_HN_U_26,
                                    tau_DP_HN_U_26, LB_tau_DP_HN_U_26, UB_tau_DP_HN_U_26,
                                    precDP_mu_HN_U_26, precDP_tau2_HN_U_26,  precDP_tau_HN_U_26,
                                    prec_basemu_HN_U_26, prec_basetau2_HN_U_26,
                                    prec_alpha_HN_U_26,  alpha_HN_U_26 , LB_alpha_HN_U_26, UB_alpha_HN_U_26,
                                    Rhat_muDP_HN_U_26, Rhat_tauDP_HN_U_26, Rhat_basemu_HN_U_26, Rhat_basetau_HN_U_26, 
                                    Rhat_alpha_HN_U_26, numclus_HN_U_26, LB_K_HN_U_26, UB_K_HN_U_26, Pr_mu_HN_U_26)

extra_col_HN_U_26 = row.names(DP_results_HN_U_26)
DPresults_HN_U_26$extra_col = extra_col_HN_U_26
prob_HN_U_26 = round(max_values,2)

write.csv(listaDP1_HN_U_26 , "C:\\Users\\kanel\\OneDrive\\Υπολογιστής\\Simulated data application\\Normal_data\\DP\\DP_HN_U_26_rel_eff_Normal_data.csv",row.names=FALSE )  
write.csv(DPresults_HN_U_26 , "C:\\Users\\kanel\\OneDrive\\Υπολογιστής\\Simulated data application\\Normal_data\\DP\\all_DP_HN_U_26_res_Normal_data.csv",row.names=FALSE )  
write.csv(RhatsDP_HN_U_26, "C:\\Users\\kanel\\OneDrive\\Υπολογιστής\\Simulated data application\\Normal_data\\DP\\DP_HN_U_26_res_Normal_data.csv",row.names=FALSE )  
write.csv(prob_HN_U_26 , "C:\\Users\\kanel\\OneDrive\\Υπολογιστής\\Simulated data application\\Normal_data\\DP\\DP_HN_U_26_prob_Normal_data.csv",row.names=FALSE )  
write.csv( cluster_mean_HN_U_26, "C:\\Users\\kanel\\OneDrive\\Υπολογιστής\\Simulated data application\\Normal_data\\DP\\DP_HN_U_26_cluster_mean_Normal_data.csv",row.names=FALSE )  



############# BINOMIAL-DP-51(HN/UNIF) #############
set.seed(123452)

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
  
}",file="DPmodel2.binHN_U_51.txt")
modfile = 'DPmodel2.binHN_U_51.txt'

run.modelDP_HN_U_51 = jags(
  data =list(ns = nrow(p[[4]]),
             r1 = p[[4]]$ci,
             r2 = p[[4]]$ti,
             n1 = p[[4]]$nci,
             n2 = p[[4]]$nti,
             N= 51 ## TRUNCATION POINT
  ) ,  
  inits = NULL,
  parameters.to.save = c(
    "basemu",    ## mu of the base Normal distribution
    "basetau",   ## sd of the base Normal distribution
    "poptrue",   ## mean treatment effect
    "var.true", ## between study variance 
    "delta12", #### study-specific effects
    "K",       ### the total number of clusters 
    "p",      ## the weights of the process 
    "alpha",  ## the concentration parameter
    "SC"     ## the probability of each cluster assignment 
  ),   
  
  n.chains = 2,
  n.iter = 50000,
  
  
  n.burnin = 10000,
  DIC = T,
  model.file = modfile
  
)

DPresults_HN_U_51 <- as.data.frame( run.modelDP_HN_U_51$BUGSoutput$summary) 


Pr_mu_HN_U_51 = mean(run.modelDP_HN_U_51$BUGSoutput$sims.matrix[  ,"poptrue"] < 0 )


DP_results_HN_U_51 <-DPresults_HN_U_51


DPresults_HN_U_51$ind <- row.names(DPresults_HN_U_51)

############# extraction of the parameters of interest ###########

f55sc <- grepl("SC", row.names(DPresults_HN_U_51))
m55sc <- DPresults_HN_U_51[f55sc,]
SC <- m55sc$mean

#### BASED ON  THE SC[i,j] ####
##### Distribution of data points to clusters ###

extra_col = rownames(m55sc)
m55sc$names = extra_col

prob_list <- data.frame(Column10 = m55sc[, 10], Column1 = m55sc[, 1])

split_dataframe <- function(df, chunk_size) {
  split(df, rep(1:ceiling(nrow(df) / chunk_size), each = chunk_size, length.out = nrow(df)))
}

split_prob_list <- split_dataframe(prob_list, nrow(p[[4]]))  

split_prob_list

# Function to find the maximum probabilities for data points to be distributed into clusters ####
find_max_values_and_indices <- function(split_list) {
  num_rows <- nrow(split_list[[1]])
  num_sublists <- length(split_list)
  
  max_values <- numeric(num_rows)
  max_indices <- integer(num_rows)
  
  for (i in 1:num_rows) {
    values <- sapply(split_list, function(x) x[i, 2])
    max_value <- max(values, na.rm = TRUE)
    max_index <- which(values == max_value)[1]  
    max_values[i] <- max_value
    max_indices[i] <- max_index
  }
  
  return(list(max_values = max_values, max_indices = max_indices))
}

result <- find_max_values_and_indices(split_prob_list)
max_values <- result$max_values #### the probabilities of cluster assignment
max_indices <- result$max_indices #### the cluster assignments 

### clusters with data points ###
# Find unique indices
k <- unique(max_indices)

# Initialize the positions list
positions <- vector("list", length(k))
names(positions) <- k

# Fill the positions list with indices
for(value in k) {
  positions[[as.character(value)]] <- which(max_indices == value)
}


### clusters with data points ###

cluster1_HN_U_51 = positions[["1"]]
cluster2_HN_U_51 = positions[["2"]]
cluster3_HN_U_51 = positions[["3"]] 
cluster4_HN_U_51 = positions[["4"]] 
cluster5_HN_U_51 = positions[["5"]] 


f11 <- grepl("poptrue", row.names(DP_results_HN_U_51))
m11 <- DP_results_HN_U_51[f11,]
mu_DP_HN_U_51<- m11$`50%`
LB_mu_DP_HN_U_51 <- m11$`2.5%`
UB_mu_DP_HN_U_51 <- m11$`97.5%`
Rhat_muDP_HN_U_51 <- m11$Rhat
precDP_mu_HN_U_51 <- UB_mu_DP_HN_U_51 - LB_mu_DP_HN_U_51

f22 <- grepl("var.true", row.names(DP_results_HN_U_51))
m22 <- DP_results_HN_U_51[f22,]
tau2_DP_HN_U_51 <- m22$`50%`
LB_tau2_DP_HN_U_51 <- m22$`2.5%`
UB_tau2_DP_HN_U_51 <- m22$`97.5%`
Rhat_tauDP_HN_U_51 <- m22$Rhat
precDP_tau2_HN_U_51 <- UB_tau2_DP_HN_U_51 -  LB_tau2_DP_HN_U_51

tau_DP_HN_U_51 <- sqrt(tau2_DP_HN_U_51)
LB_tau_DP_HN_U_51 <-sqrt(LB_tau2_DP_HN_U_51)
UB_tau_DP_HN_U_51 <-sqrt(UB_tau2_DP_HN_U_51)
precDP_tau_HN_U_51 <- UB_tau_DP_HN_U_51 -  LB_tau_DP_HN_U_51

fdd <- grepl("delta12", row.names(DP_results_HN_U_51))
mdd <- DP_results_HN_U_51[fdd,]
rel_effDP_HN_U_51 <- mdd$`50%`
LB_rel_effDP_HN_U_51 <- mdd$`2.5%`
UB_rel_effDP_HN_U_51 <- mdd$`97.5%`
sd_rel_effDP_HN_U_51 <- mdd$sd
Rhat_deltaDP_HN_U_51 <- mdd$Rhat

#### Mean of each cluster ######
cl1_HN_U_51 = rel_effDP_HN_U_51[cluster1_HN_U_51]
mean_cl1_HN_U_51 = mean(cl1_HN_U_51)
cl2_HN_U_51 = rel_effDP_HN_U_51[cluster2_HN_U_51]
mean_cl2_HN_U_51 = mean(cl2_HN_U_51)
cl3_HN_U_51 = rel_effDP_HN_U_51[cluster3_HN_U_51]
mean_cl3_HN_U_51 = mean(cl3_HN_U_51)
cl4_HN_U_51 = rel_effDP_HN_U_51[cluster4_HN_U_51]
mean_cl4_HN_U_51 = mean(cl4_HN_U_51)
cl5_HN_U_51 = rel_effDP_HN_U_51[cluster5_HN_U_51]
mean_cl5_HN_U_51 = mean(cl5_HN_U_51)


cluster_mean_HN_U_51 = cbind.data.frame(mean_cl1_HN_U_51, mean_cl2_HN_U_51, mean_cl3_HN_U_51,mean_cl4_HN_U_51, mean_cl5_HN_U_51 )


f33 <- grepl("basemu", row.names(DP_results_HN_U_51))
m33 <- DP_results_HN_U_51[f33,]
base_mu_HN_U_51 <- m33$`50%`
LB_base_mu_HN_U_51 <- m33$`2.5%`
UB_base_mu_HN_U_51 <- m33$`97.5%`
Rhat_basemu_HN_U_51 <- m33$Rhat
prec_basemu_HN_U_51 <- UB_base_mu_HN_U_51 - LB_base_mu_HN_U_51

f44 <- grepl("basetau", row.names(DP_results_HN_U_51))
m44 <- DP_results_HN_U_51[f44,]
base_tau_HN_U_51 <- m44$`50%`
base_tau2_HN_U_51 <- (m44$`50%`)^2
LB_base_tau2_HN_U_51 <- (m44$`2.5%`)^2
UB_base_tau2_HN_U_51 <- (m44$`97.5%`)^2
Rhat_basetau_HN_U_51 <- (m44$Rhat)^2
prec_basetau2_HN_U_51 <- UB_base_tau2_HN_U_51 - LB_base_tau2_HN_U_51

f55 <- grepl("alpha", row.names(DP_results_HN_U_51))
m55 <- DP_results_HN_U_51[f55,]
alpha_HN_U_51 <- m55$`50%`
LB_alpha_HN_U_51 <- m55$`2.5%`
UB_alpha_HN_U_51 <- m55$`97.5%`
Rhat_alpha_HN_U_51 <- m55$Rhat
prec_alpha_HN_U_51 <- UB_alpha_HN_U_51 - LB_alpha_HN_U_51

fcl <- grepl("K", row.names(DP_results_HN_U_51))   
mcl <- DP_results_HN_U_51[fcl,]
median_K_HN_U_51 <- mcl$`50%` 
LB_K_HN_U_51 <- mcl$`2.5%`
UB_K_HN_U_51 <- mcl$`97.5%`

fp <- grepl("p", row.names(DP_results_HN_U_51))
mp <- DP_results_HN_U_51[fp,]
mp <- mp[!grepl("pop", row.names(mp)),]
mp <- mp[!grepl("alpha", row.names(mp)),]
pi_HN_U_51 <- mp$mean


listaDP1_HN_U_51 <- cbind.data.frame(rel_effDP_HN_U_51,sd_rel_effDP_HN_U_51, LB_rel_effDP_HN_U_51,UB_rel_effDP_HN_U_51,Rhat_deltaDP_HN_U_51)


numclus_HN_U_51 <- unlist(median_K_HN_U_51)
LB_K_HN_U_51 <- unlist(LB_K_HN_U_51)
UB_K_HN_U_51 <- unlist(UB_K_HN_U_51)

RhatsDP_HN_U_51 <- cbind.data.frame(base_mu_HN_U_51,LB_base_mu_HN_U_51,UB_base_mu_HN_U_51,
                                    base_tau_HN_U_51,base_tau2_HN_U_51, LB_base_tau2_HN_U_51,UB_base_tau2_HN_U_51,
                                    mu_DP_HN_U_51, LB_mu_DP_HN_U_51, UB_mu_DP_HN_U_51, 
                                    tau2_DP_HN_U_51, LB_tau2_DP_HN_U_51, UB_tau2_DP_HN_U_51,
                                    tau_DP_HN_U_51, LB_tau_DP_HN_U_51, UB_tau_DP_HN_U_51,
                                    precDP_mu_HN_U_51, precDP_tau2_HN_U_51, precDP_tau_HN_U_51,
                                    prec_basemu_HN_U_51, prec_basetau2_HN_U_51,
                                    prec_alpha_HN_U_51,  alpha_HN_U_51 , LB_alpha_HN_U_51, UB_alpha_HN_U_51,
                                    Rhat_muDP_HN_U_51, Rhat_tauDP_HN_U_51, Rhat_basemu_HN_U_51, Rhat_basetau_HN_U_51, 
                                    Rhat_alpha_HN_U_51, numclus_HN_U_51, LB_K_HN_U_51, UB_K_HN_U_51, Pr_mu_HN_U_51)

extra_col_HN_U_51 = row.names(DP_results_HN_U_51)
DPresults_HN_U_51$extra_col = extra_col_HN_U_51
prob_HN_U_51 = round(max_values,2)

write.csv(listaDP1_HN_U_51 , "C:\\Users\\kanel\\OneDrive\\Υπολογιστής\\Simulated data application\\Normal_data\\DP\\DP_HN_U_51_rel_eff_Normal_data.csv",row.names=FALSE )  
write.csv(DPresults_HN_U_51 , "C:\\Users\\kanel\\OneDrive\\Υπολογιστής\\Simulated data application\\Normal_data\\DP\\all_DP_HN_U_51_res_Normal_data.csv",row.names=FALSE )  
write.csv(RhatsDP_HN_U_51, "C:\\Users\\kanel\\OneDrive\\Υπολογιστής\\Simulated data application\\Normal_data\\DP\\DP_HN_U_51_res_Normal_data.csv",row.names=FALSE )  
write.csv(prob_HN_U_51 , "C:\\Users\\kanel\\OneDrive\\Υπολογιστής\\Simulated data application\\Normal_data\\DP\\DP_HN_U_51_prob_Normal_data.csv",row.names=FALSE )  
write.csv( cluster_mean_HN_U_51, "C:\\Users\\kanel\\OneDrive\\Υπολογιστής\\Simulated data application\\Normal_data\\DP\\DP_HN_U_51_cluster_mean_Normal_data.csv",row.names=FALSE )  



############# BINOMIAL-DP-26(Unif/Unif) #############
set.seed(123453)
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
  
}",file="DPmodel2.binU_U_26.txt")
modfile = 'DPmodel2.binU_U_26.txt'

run.modelDP_U_U_26 = jags(
  data =list(ns = nrow(p[[4]]),
             r1 = p[[4]]$ci,
             r2 = p[[4]]$ti,
             n1 = p[[4]]$nci,
             n2 = p[[4]]$nti,
             N= 26 ## TRUNCATION POINT
  ) ,  
  inits = NULL,
  parameters.to.save = c(
    "basemu",    ## mu of the base Normal distribution
    "basetau",   ## sd of the base Normal distribution
    "poptrue",   ## mean treatment effect
    "var.true", ## between study variance 
    "delta12", #### study-specific effects
    "K",       ### the total number of clusters 
    "p",      ## the weights of the process 
    "alpha",  ## the concentration parameter
    "SC"     ## the probability of each cluster assignment 
  ),   
  
  n.chains = 2,
  n.iter = 50000,
  
  
  n.burnin = 10000,
  DIC = T,
  model.file = modfile
  
)

DPresults_U_U_26 <- as.data.frame( run.modelDP_U_U_26$BUGSoutput$summary) 


Pr_mu_U_U_26 = mean(run.modelDP_U_U_26$BUGSoutput$sims.matrix[  ,"poptrue"] < 0 )


DP_results_U_U_26 <-DPresults_U_U_26 


DPresults_U_U_26$ind <- row.names(DPresults_U_U_26)

############# extraction of the parameters of interest ###########

f55sc <- grepl("SC", row.names(DPresults_U_U_26))
m55sc <- DPresults_U_U_26[f55sc,]
SC <- m55sc$mean

#### BASED ON  THE SC[i,j] ####
##### Distribution of data points to clusters ###

extra_col = rownames(m55sc)
m55sc$names = extra_col

prob_list <- data.frame(Column10 = m55sc[, 10], Column1 = m55sc[, 1])

split_dataframe <- function(df, chunk_size) {
  split(df, rep(1:ceiling(nrow(df) / chunk_size), each = chunk_size, length.out = nrow(df)))
}

split_prob_list <- split_dataframe(prob_list, nrow(p[[4]]))  

split_prob_list

# Function to find the maximum probabilities for data points to be distributed into clusters ####
find_max_values_and_indices <- function(split_list) {
  num_rows <- nrow(split_list[[1]])
  num_sublists <- length(split_list)
  
  max_values <- numeric(num_rows)
  max_indices <- integer(num_rows)
  
  for (i in 1:num_rows) {
    values <- sapply(split_list, function(x) x[i, 2])
    max_value <- max(values, na.rm = TRUE)
    max_index <- which(values == max_value)[1]  
    max_values[i] <- max_value
    max_indices[i] <- max_index
  }
  
  return(list(max_values = max_values, max_indices = max_indices))
}

result <- find_max_values_and_indices(split_prob_list)
max_values <- result$max_values #### the probabilities of cluster assignment
max_indices <- result$max_indices #### the cluster assignments 

### clusters with data points ###
k <- unique(max_indices)

positions <- vector("list", length(k))
names(positions) <- k

for(value in k) {
  positions[[as.character(value)]] <- which(max_indices == value)
}


### clusters with data points ###

cluster1_U_U_26 = positions[["1"]]
cluster2_U_U_26 = positions[["2"]]
cluster3_U_U_26 = positions[["3"]] 
cluster4_U_U_26 = positions[["4"]] 
cluster5_U_U_26 = positions[["5"]] 


f11 <- grepl("poptrue", row.names(DP_results_U_U_26))
m11 <- DP_results_U_U_26[f11,]
mu_DP_U_U_26<- m11$`50%`
LB_mu_DP_U_U_26 <- m11$`2.5%`
UB_mu_DP_U_U_26 <- m11$`97.5%`
Rhat_muDP_U_U_26 <- m11$Rhat
precDP_mu_U_U_26 <- UB_mu_DP_U_U_26 - LB_mu_DP_U_U_26

f22 <- grepl("var.true", row.names(DP_results_U_U_26))
m22 <- DP_results_U_U_26[f22,]
tau2_DP_U_U_26 <- m22$`50%`
LB_tau2_DP_U_U_26 <- m22$`2.5%`
UB_tau2_DP_U_U_26 <- m22$`97.5%`
Rhat_tauDP_U_U_26 <- m22$Rhat
precDP_tau2_U_U_26 <- UB_tau2_DP_U_U_26 -  LB_tau2_DP_U_U_26

tau_DP_U_U_26 <- sqrt(tau2_DP_U_U_26)
LB_tau_DP_U_U_26 <- sqrt(LB_tau2_DP_U_U_26)
UB_tau_DP_U_U_26 <- sqrt(UB_tau2_DP_U_U_26)
precDP_tau_U_U_26 <- UB_tau_DP_U_U_26 -  LB_tau_DP_U_U_26


fdd <- grepl("delta12", row.names(DP_results_U_U_26))
mdd <- DP_results_U_U_26[fdd,]
rel_effDP_U_U_26 <- mdd$`50%`
LB_rel_effDP_U_U_26 <- mdd$`2.5%`
UB_rel_effDP_U_U_26 <- mdd$`97.5%`
sd_rel_effDP_U_U_26 <- mdd$sd
Rhat_deltaDP_U_U_26 <- mdd$Rhat

#### Mean of each cluster ######
cl1_U_U_26 = rel_effDP_U_U_26[cluster1_U_U_26]
mean_cl1_U_U_26 = mean(cl1_U_U_26)
cl2_U_U_26 = rel_effDP_U_U_26[cluster2_U_U_26]
mean_cl2_U_U_26 = mean(cl2_U_U_26)
cl3_U_U_26 = rel_effDP_U_U_26[cluster3_U_U_26]
mean_cl3_U_U_26 = mean(cl3_U_U_26)
cl4_U_U_26 = rel_effDP_U_U_26[cluster4_U_U_26]
mean_cl4_U_U_26 = mean(cl4_U_U_26)
cl5_U_U_26 = rel_effDP_U_U_26[cluster5_U_U_26]
mean_cl5_U_U_26 = mean(cl5_U_U_26)


cluster_mean_U_U_26 = cbind.data.frame(mean_cl1_U_U_26, mean_cl2_U_U_26, mean_cl3_U_U_26,mean_cl4_U_U_26, mean_cl5_U_U_26 )


f33 <- grepl("basemu", row.names(DP_results_U_U_26))
m33 <- DP_results_U_U_26[f33,]
base_mu_U_U_26 <- m33$`50%`
LB_base_mu_U_U_26 <- m33$`2.5%`
UB_base_mu_U_U_26 <- m33$`97.5%`
Rhat_basemu_U_U_26 <- m33$Rhat
prec_basemu_U_U_26 <- UB_base_mu_U_U_26 - LB_base_mu_U_U_26

f44 <- grepl("basetau", row.names(DP_results_U_U_26))
m44 <- DP_results_U_U_26[f44,]
base_tau_U_U_26 <- m44$`50%`
base_tau2_U_U_26 <- (m44$`50%`)^2
LB_base_tau2_U_U_26 <- (m44$`2.5%`)^2
UB_base_tau2_U_U_26 <- (m44$`97.5%`)^2
Rhat_basetau_U_U_26 <- (m44$Rhat)^2
prec_basetau2_U_U_26 <- UB_base_tau2_U_U_26 - LB_base_tau2_U_U_26

f55 <- grepl("alpha", row.names(DP_results_U_U_26))
m55 <- DP_results_U_U_26[f55,]
alpha_U_U_26 <- m55$`50%`
LB_alpha_U_U_26 <- m55$`2.5%`
UB_alpha_U_U_26 <- m55$`97.5%`
Rhat_alpha_U_U_26 <- m55$Rhat
prec_alpha_U_U_26 <- UB_alpha_U_U_26 - LB_alpha_U_U_26

fcl <- grepl("K", row.names(DP_results_U_U_26))   
mcl <- DP_results_U_U_26[fcl,]
median_K_U_U_26 <- mcl$`50%` 
LB_K_U_U_26 <- mcl$`2.5%`
UB_K_U_U_26 <- mcl$`97.5%`

fp <- grepl("p", row.names(DP_results_U_U_26))
mp <- DP_results_U_U_26[fp,]
mp <- mp[!grepl("pop", row.names(mp)),]
mp <- mp[!grepl("alpha", row.names(mp)),]
pi_U_U_26 <- mp$mean


listaDP1_U_U_26 <- cbind.data.frame(rel_effDP_U_U_26,sd_rel_effDP_U_U_26, LB_rel_effDP_U_U_26,UB_rel_effDP_U_U_26,Rhat_deltaDP_U_U_26)


numclus_U_U_26 <- unlist(median_K_U_U_26)
LB_K_U_U_26 <- unlist(LB_K_U_U_26)
UB_K_U_U_26 <- unlist(UB_K_U_U_26)

RhatsDP_U_U_26 <- cbind.data.frame(base_mu_U_U_26,LB_base_mu_U_U_26,UB_base_mu_U_U_26,
                                   base_tau_U_U_26,base_tau2_U_U_26, LB_base_tau2_U_U_26,UB_base_tau2_U_U_26,
                                   mu_DP_U_U_26, LB_mu_DP_U_U_26, UB_mu_DP_U_U_26, 
                                   tau2_DP_U_U_26, LB_tau2_DP_U_U_26, UB_tau2_DP_U_U_26,
                                   tau_DP_U_U_26, LB_tau_DP_U_U_26, UB_tau_DP_U_U_26,
                                   precDP_mu_U_U_26, precDP_tau2_U_U_26,precDP_tau_U_U_26,
                                   prec_basemu_U_U_26, prec_basetau2_U_U_26,
                                   prec_alpha_U_U_26,  alpha_U_U_26 , LB_alpha_U_U_26, UB_alpha_U_U_26,
                                   Rhat_muDP_U_U_26, Rhat_tauDP_U_U_26, Rhat_basemu_U_U_26, Rhat_basetau_U_U_26, 
                                   Rhat_alpha_U_U_26, numclus_U_U_26, LB_K_U_U_26, UB_K_U_U_26, Pr_mu_U_U_26)

extra_col_U_U_26 = row.names(DP_results_U_U_26)
DPresults_U_U_26$extra_col = extra_col_U_U_26
prob_U_U_26 = round(max_values,2)

write.csv(listaDP1_U_U_26 , "C:\\Users\\kanel\\OneDrive\\Υπολογιστής\\Simulated data application\\Normal_data\\DP\\DP_U_U_26_rel_eff_Normal_data.csv",row.names=FALSE )  
write.csv(DPresults_U_U_26 , "C:\\Users\\kanel\\OneDrive\\Υπολογιστής\\Simulated data application\\Normal_data\\DP\\all_DP_U_U_26_res_Normal_data.csv",row.names=FALSE )  
write.csv(RhatsDP_U_U_26, "C:\\Users\\kanel\\OneDrive\\Υπολογιστής\\Simulated data application\\Normal_data\\DP\\DP_U_U_26_res_Normal_data.csv",row.names=FALSE )  
write.csv(prob_U_U_26 , "C:\\Users\\kanel\\OneDrive\\Υπολογιστής\\Simulated data application\\Normal_data\\DP\\DP_U_U_26_prob_Normal_data.csv",row.names=FALSE )  
write.csv( cluster_mean_U_U_26, "C:\\Users\\kanel\\OneDrive\\Υπολογιστής\\Simulated data application\\Normal_data\\DP\\DP_U_U_26_cluster_mean_Normal_data.csv",row.names=FALSE )  

############# BINOMIAL-DP-51(Unif/Unif) #############
set.seed(123454)
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
  
}",file="DPmodel2.binU_U_51.txt")
modfile = 'DPmodel2.binU_U_51.txt'

run.modelDPU_U_51 = jags(
  data =list(ns = nrow(p[[4]]),
             r1 = p[[4]]$ci,
             r2 = p[[4]]$ti,
             n1 = p[[4]]$nci,
             n2 = p[[4]]$nti,
             N= 51 ## TRUNCATION POINT
  ) ,  
  inits = NULL,
  parameters.to.save = c(
    "basemu",    ## mu of the base Normal distribution
    "basetau",   ## sd of the base Normal distribution
    "poptrue",   ## mean treatment effect
    "var.true", ## between study variance 
    "delta12", #### study-specific effects
    "K",       ### the total number of clusters 
    "p",      ## the weights of the process 
    "alpha",  ## the concentration parameter
    "SC"     ## the probability of each cluster assignment 
  ),   
  
  n.chains = 2,
  n.iter = 50000,
  
  
  n.burnin = 10000,
  DIC = T,
  model.file = modfile
  
)

DPresultsU_U_51 <- as.data.frame( run.modelDPU_U_51$BUGSoutput$summary) 


Pr_muU_U_51 = mean(run.modelDPU_U_51$BUGSoutput$sims.matrix[  ,"poptrue"] < 0 )


DP_resultsU_U_51 <-DPresultsU_U_51


DPresultsU_U_51$ind <- row.names(DPresultsU_U_51)

############# extraction of the parameters of interest ###########

f55sc <- grepl("SC", row.names(DPresultsU_U_51))
m55sc <- DPresultsU_U_51[f55sc,]
SC <- m55sc$mean

#### BASED ON  THE SC[i,j] ####
##### Distribution of data points to clusters ###

extra_col = rownames(m55sc)
m55sc$names = extra_col

prob_list <- data.frame(Column10 = m55sc[, 10], Column1 = m55sc[, 1])

split_dataframe <- function(df, chunk_size) {
  split(df, rep(1:ceiling(nrow(df) / chunk_size), each = chunk_size, length.out = nrow(df)))
}

split_prob_list <- split_dataframe(prob_list, nrow(p[[4]]))  

split_prob_list

# Function to find the maximum probabilities for data points to be distributed into clusters ####
find_max_values_and_indices <- function(split_list) {
  num_rows <- nrow(split_list[[1]])
  num_sublists <- length(split_list)
  
  max_values <- numeric(num_rows)
  max_indices <- integer(num_rows)
  
  for (i in 1:num_rows) {
    values <- sapply(split_list, function(x) x[i, 2])
    max_value <- max(values, na.rm = TRUE)
    max_index <- which(values == max_value)[1]  
    max_values[i] <- max_value
    max_indices[i] <- max_index
  }
  
  return(list(max_values = max_values, max_indices = max_indices))
}

result <- find_max_values_and_indices(split_prob_list)
max_values <- result$max_values #### the probabilities of cluster assignment
max_indices <- result$max_indices #### the cluster assignments 

### clusters with data points ###
k <- unique(max_indices)

positions <- vector("list", length(k))
names(positions) <- k

for(value in k) {
  positions[[as.character(value)]] <- which(max_indices == value)
}


### clusters with data points ###

cluster1U_U_51 = positions[["1"]]
cluster2U_U_51 = positions[["2"]]
cluster3U_U_51 = positions[["3"]] 
cluster4U_U_51 = positions[["4"]] 


f11 <- grepl("poptrue", row.names(DP_resultsU_U_51))
m11 <- DP_resultsU_U_51[f11,]
mu_DPU_U_51<- m11$`50%`
LB_mu_DPU_U_51 <- m11$`2.5%`
UB_mu_DPU_U_51 <- m11$`97.5%`
Rhat_muDPU_U_51 <- m11$Rhat
precDP_muU_U_51 <- UB_mu_DPU_U_51 - LB_mu_DPU_U_51

f22 <- grepl("var.true", row.names(DP_resultsU_U_51))
m22 <- DP_resultsU_U_51[f22,]
tau2_DPU_U_51 <- m22$`50%`
LB_tau2_DPU_U_51 <- m22$`2.5%`
UB_tau2_DPU_U_51 <- m22$`97.5%`
Rhat_tauDPU_U_51 <- m22$Rhat
precDP_tau2U_U_51 <- UB_tau2_DPU_U_51 -  LB_tau2_DPU_U_51

tau_DPU_U_51 <- sqrt(tau2_DPU_U_51)
LB_tau_DPU_U_51 <- sqrt(LB_tau2_DPU_U_51)
UB_tau_DPU_U_51 <- sqrt(UB_tau2_DPU_U_51)
precDP_tauU_U_51 <- UB_tau_DPU_U_51 -  LB_tau_DPU_U_51


fdd <- grepl("delta12", row.names(DP_resultsU_U_51))
mdd <- DP_resultsU_U_51[fdd,]
rel_effDPU_U_51 <- mdd$`50%`
LB_rel_effDPU_U_51 <- mdd$`2.5%`
UB_rel_effDPU_U_51 <- mdd$`97.5%`
sd_rel_effDPU_U_51 <- mdd$sd
Rhat_deltaDPU_U_51 <- mdd$Rhat

#### Mean of each cluster ######
cl1U_U_51 = rel_effDPU_U_51[cluster1U_U_51]
mean_cl1U_U_51 = mean(cl1U_U_51)
cl2U_U_51 = rel_effDPU_U_51[cluster2U_U_51]
mean_cl2U_U_51 = mean(cl2U_U_51)
cl3U_U_51 = rel_effDPU_U_51[cluster3U_U_51]
mean_cl3U_U_51 = mean(cl3U_U_51)
cl4U_U_51 = rel_effDPU_U_51[cluster4U_U_51]
mean_cl4U_U_51 = mean(cl4U_U_51)


cluster_meanU_U_51 = cbind.data.frame(mean_cl1U_U_51, mean_cl2U_U_51, mean_cl3U_U_51,mean_cl4U_U_51 )


f33 <- grepl("basemu", row.names(DP_resultsU_U_51))
m33 <- DP_resultsU_U_51[f33,]
base_muU_U_51 <- m33$`50%`
LB_base_muU_U_51 <- m33$`2.5%`
UB_base_muU_U_51 <- m33$`97.5%`
Rhat_basemuU_U_51 <- m33$Rhat
prec_basemuU_U_51 <- UB_base_muU_U_51 - LB_base_muU_U_51

f44 <- grepl("basetau", row.names(DP_resultsU_U_51))
m44 <- DP_resultsU_U_51[f44,]
base_tauU_U_51 <- m44$`50%`
base_tau2U_U_51 <- (m44$`50%`)^2
LB_base_tau2U_U_51 <- (m44$`2.5%`)^2
UB_base_tau2U_U_51 <- (m44$`97.5%`)^2
Rhat_basetauU_U_51 <- (m44$Rhat)^2
prec_basetau2U_U_51 <- UB_base_tau2U_U_51 - LB_base_tau2U_U_51

f55 <- grepl("alpha", row.names(DP_resultsU_U_51))
m55 <- DP_resultsU_U_51[f55,]
alphaU_U_51 <- m55$`50%`
LB_alphaU_U_51 <- m55$`2.5%`
UB_alphaU_U_51 <- m55$`97.5%`
Rhat_alphaU_U_51 <- m55$Rhat
prec_alphaU_U_51 <- UB_alphaU_U_51 - LB_alphaU_U_51

fcl <- grepl("K", row.names(DP_resultsU_U_51))   
mcl <- DP_resultsU_U_51[fcl,]
median_KU_U_51 <- mcl$`50%` 
LB_KU_U_51 <- mcl$`2.5%`
UB_KU_U_51 <- mcl$`97.5%`

fp <- grepl("p", row.names(DP_resultsU_U_51))
mp <- DP_resultsU_U_51[fp,]
mp <- mp[!grepl("pop", row.names(mp)),]
mp <- mp[!grepl("alpha", row.names(mp)),]
piU_U_51 <- mp$mean


listaDP1U_U_51 <- cbind.data.frame(rel_effDPU_U_51,sd_rel_effDPU_U_51, LB_rel_effDPU_U_51,UB_rel_effDPU_U_51,Rhat_deltaDPU_U_51)


numclusU_U_51 <- unlist(median_KU_U_51)
LB_KU_U_51 <- unlist(LB_KU_U_51)
UB_KU_U_51 <- unlist(UB_KU_U_51)

RhatsDPU_U_51 <- cbind.data.frame(base_muU_U_51,LB_base_muU_U_51,UB_base_muU_U_51,
                                  base_tauU_U_51,base_tau2U_U_51, LB_base_tau2U_U_51,UB_base_tau2U_U_51,
                                  mu_DPU_U_51, LB_mu_DPU_U_51, UB_mu_DPU_U_51, tau2_DPU_U_51, LB_tau2_DPU_U_51, UB_tau2_DPU_U_51,
                                  tau_DPU_U_51, LB_tau_DPU_U_51, UB_tau_DPU_U_51,
                                  precDP_muU_U_51, precDP_tau2U_U_51, precDP_tauU_U_51, 
                                  prec_basemuU_U_51, prec_basetau2U_U_51,
                                  prec_alphaU_U_51,  alphaU_U_51 , LB_alphaU_U_51, UB_alphaU_U_51,
                                  Rhat_muDPU_U_51, Rhat_tauDPU_U_51, Rhat_basemuU_U_51, Rhat_basetauU_U_51, 
                                  Rhat_alphaU_U_51, numclusU_U_51, LB_KU_U_51, UB_KU_U_51, Pr_muU_U_51)

extra_colU_U_51 = row.names(DP_resultsU_U_51)
DPresultsU_U_51$extra_col = extra_colU_U_51
probU_U_51 = round(max_values,2)

write.csv(listaDP1U_U_51 , "C:\\Users\\kanel\\OneDrive\\Υπολογιστής\\Simulated data application\\Normal_data\\DP\\DPU_U_51_rel_eff_Normal_data.csv",row.names=FALSE )  
write.csv(DPresultsU_U_51 , "C:\\Users\\kanel\\OneDrive\\Υπολογιστής\\Simulated data application\\Normal_data\\DP\\all_DPU_U_51_res_Normal_data.csv",row.names=FALSE )  
write.csv(RhatsDPU_U_51, "C:\\Users\\kanel\\OneDrive\\Υπολογιστής\\Simulated data application\\Normal_data\\DP\\DPU_U_51_res_Normal_data.csv",row.names=FALSE )  
write.csv(probU_U_51 , "C:\\Users\\kanel\\OneDrive\\Υπολογιστής\\Simulated data application\\Normal_data\\DP\\DPU_U_51_prob_Normal_data.csv",row.names=FALSE )  
write.csv( cluster_meanU_U_51, "C:\\Users\\kanel\\OneDrive\\Υπολογιστής\\Simulated data application\\Normal_data\\DP\\DPU_U_51_cluster_mean_Normal_data.csv",row.names=FALSE )  

####################### Binomial-t(HN) model ########################

library(cmdstanr)
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
check_cmdstan_toolchain()

filet_distr <- file.path(cmdstan_path(), "examples", "Stanmodels", "t_model_bin_HNpred.stan") 

modt_distr <- cmdstan_model(filet_distr)

fitt_distr <- modt_distr$sample(
  data = list(ns = nrow(p[[4]]),
              r1 = p[[4]]$ci,
              r2 = p[[4]]$ti,
              n1 = p[[4]]$nci,
              n2 = p[[4]]$nti), 
  seed = 12347, 
  chains = 2, 
  parallel_chains = 2,
  refresh = 500,
  iter_warmup = 10000,
  iter_sampling = 50000,
  adapt_delta = 0.99
)


t_distr.res <- as.data.frame(fitt_distr$summary())

#### To visualize the predictive distribution
bayesplot::color_scheme_set("brightblue")
bayesplot::mcmc_dens(fitt_distr$draws(c("mu", "pred")))

### Pr of mu <- 0
Pr_mu = mean(fitt_distr$draws("mu") < 0 )

### Pr of pred <- 0
Pr_pred = mean(fitt_distr$draws("pred") < 0 )

### the study specific effects
d <- grepl("delta",t_distr.res$variable)

deltat <- t_distr.res[d,]

## the mean treatment effect
dm <- grepl("mu",t_distr.res$variable)

median <- t_distr.res[dm,]$median
LBmu <- t_distr.res[dm,]$q5
UBmu <- t_distr.res[dm,]$q95
Rhat_mut_distr <-  t_distr.res[dm,]$rhat

## prediction
ddn <- grepl("pred",t_distr.res$variable)

delta_new <- t_distr.res[ddn,]$median
LBdelta_new <- t_distr.res[ddn,]$q5
UBdelta_new <- t_distr.res[ddn,]$q95
Rhat_delta_newt_distr <-  t_distr.res[ddn,]$rhat

## between-study variance 
tau2 <- c()
dtau2 <- grepl("tau_sqr",t_distr.res$variable)

tau2 <- t_distr.res[dtau2,]$median
LBtau2 <- t_distr.res[dtau2,]$q5
UBtau2 <- t_distr.res[dtau2,]$q95
Rhat_tau2t_distr <- t_distr.res[dtau2,]$rhat

## between-study sd
tau <- sqrt(t_distr.res[dtau2,]$median)
LBtau <- sqrt(t_distr.res[dtau2,]$q5)
UBtau <- sqrt(t_distr.res[dtau2,]$q95)

t_distr.res1<-data.frame(median=median, 
                         lowerCI=LBmu,
                         upperCI=UBmu,
                         tau2=tau2,
                         l_tau2 = LBtau2,
                         u_tau2 = UBtau2,
                         tau=tau,
                         l_tau = LBtau,
                         u_tau = UBtau,
                         Rhat_mut_distr = Rhat_mut_distr,
                         Rhat_tau2t_distr = Rhat_tau2t_distr,
                         delta_new = delta_new,
                         LBdelta_new = LBdelta_new,
                         UBdelta_new = UBdelta_new,
                         Rhat_delta_newt_distr = Rhat_delta_newt_distr,
                         Pr_mu = Pr_mu,
                         Pr_pred = Pr_pred
)


write.csv(t_distr.res1 , "C:\\Users\\kanel\\OneDrive\\Υπολογιστής\\Simulated data application\\Normal_data\\T\\pred_res_Normal_data_Bayesian_t_model_HN.csv",row.names=FALSE )  
write.csv(deltat, "C:\\Users\\kanel\\OneDrive\\Υπολογιστής\\Simulated data application\\Normal_data\\T\\pred_rel_effT_Normal_data_Bayesian_t_model_HN.csv",row.names=FALSE )  
# 
# 

####### Binomial-t(Unif)
library(cmdstanr)
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
check_cmdstan_toolchain()

filet_distr_U <- file.path(cmdstan_path(), "examples", "Stanmodels", "t_model_binUpred.stan") 

modt_distr_U <- cmdstan_model(filet_distr_U)

fitt_distr_U <- modt_distr_U$sample(
  data = list(ns = nrow(p[[4]]),
              r1 = p[[4]]$ci,
              r2 = p[[4]]$ti,
              n1 = p[[4]]$nci,
              n2 = p[[4]]$nti), 
  seed = 123470, 
  chains = 2, 
  parallel_chains = 2,
  refresh = 500,
  iter_warmup = 10000,
  iter_sampling = 50000,
  adapt_delta = 0.99
)


t_distr.res_U <- as.data.frame(fitt_distr_U$summary())

#### To visualize the predictive distribution
bayesplot::color_scheme_set("brightblue")
bayesplot::mcmc_dens(fitt_distr_U$draws(c("mu", "pred")))

### Pr of mu <- 0
Pr_mu_U = mean(fitt_distr_U$draws("mu") < 0 )

### Pr of pred <- 0
Pr_pred_U = mean(fitt_distr_U$draws("pred") < 0 )

## the study specific effects
d_U <- grepl("delta",t_distr.res_U$variable)

deltat_U <- t_distr.res_U[d_U,]

## the mean treatment effect
dm_U <- grepl("mu",t_distr.res_U$variable)

median_U <- t_distr.res_U[dm_U,]$median
LBmu_U <- t_distr.res_U[dm_U,]$q5
UBmu_U <- t_distr.res_U[dm_U,]$q95
Rhat_mut_distr_U <-  t_distr.res_U[dm_U,]$rhat

## prediction
ddn_U <- grepl("pred",t_distr.res_U$variable)

delta_new_U <- t_distr.res_U[ddn_U,]$median
LBdelta_new_U <- t_distr.res_U[ddn_U,]$q5
UBdelta_new_U <- t_distr.res_U[ddn_U,]$q95
Rhat_delta_newt_distr_U <-  t_distr.res_U[ddn_U,]$rhat

## the between-study variance
tau2_U <- c()
dtau2_U <- grepl("tau_sqr",t_distr.res_U$variable)

tau2_U <- t_distr.res_U[dtau2_U,]$median
LBtau2_U <- t_distr.res_U[dtau2_U,]$q5
UBtau2_U <- t_distr.res_U[dtau2_U,]$q95
Rhat_tau2t_distr_U <- t_distr.res_U[dtau2_U,]$rhat

## the between-study sd
tau_U <- sqrt(t_distr.res_U[dtau2_U,]$median)
LBtau_U <- sqrt(t_distr.res_U[dtau2_U,]$q5)
UBtau_U <- sqrt(t_distr.res_U[dtau2_U,]$q95)

t_distr.res_U1<-data.frame(median=median_U, 
                           lowerCI=LBmu_U,
                           upperCI=UBmu_U,
                           tau2=tau2_U,
                           l_tau2 = LBtau2_U,
                           u_tau2 = UBtau2_U,
                           tau=tau_U,
                           l_tau = LBtau_U,
                           u_tau = UBtau_U,
                           Rhat_mut_distr = Rhat_mut_distr_U,
                           Rhat_tau2t_distr = Rhat_tau2t_distr_U,
                           delta_new = delta_new_U,
                           LBdelta_new = LBdelta_new_U,
                           UBdelta_new = UBdelta_new_U,
                           Rhat_delta_newt_distr = Rhat_delta_newt_distr_U,
                           Pr_mu = Pr_mu_U,
                           Pr_pred = Pr_pred_U
)


write.csv(t_distr.res_U1 , "C:\\Users\\Lela Panag\\Desktop\\Github\\GITHUB_FINAL\\Forest plots of 3 selected simulated datasets\\pred_res_Normal_data_Bayesian_t_model_U.csv",row.names=FALSE )  
write.csv(deltat_U, "C:\\Users\\Lela Panag\\Desktop\\Github\\GITHUB_FINAL\\Forest plots of 3 selected simulated datasets\\pred_rel_effT_Normal_data_Bayesian_t_model_U.csv",row.names=FALSE )  
# 
# 
# ################################## Binomial-SN(HN) model #############################

library(cmdstanr)
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
check_cmdstan_toolchain()

fileSN <- file.path(cmdstan_path(), "examples", "Stanmodels", "Skew_normal_model_bin_HNpred.stan")
modSN <- cmdstan_model(fileSN)

fitSN <- modSN$sample(
  data =list(ns = nrow(p[[4]]),
             r1 = p[[4]]$ci,
             r2 = p[[4]]$ti,
             n1 = p[[4]]$nci,
             n2 = p[[4]]$nti), 
  seed = 12348, 
  chains = 2, 
  parallel_chains = 2,
  refresh = 500,
  iter_warmup = 10000,
  iter_sampling = 50000,
  adapt_delta = 0.99
)

sn.res <- as.data.frame(fitSN$summary())

#### To visualize the predictive distribution
bayesplot::color_scheme_set("brightblue")
bayesplot::mcmc_dens(fitSN$draws(c("mu", "pred1")))

### Pr of mu <- 0
Pr_mu = mean(fitSN$draws("mu") < 0 )

### Pr of prediction <- 0
Pr_pred = mean(fitSN$draws("pred1") < 0 )

## the lerative effects
d <- grepl("delta",sn.res$variable)

deltaSN <- sn.res[d,]

## the mean treatment effect
dm <- grepl("mu",sn.res$variable)

median <- sn.res[dm,]$median
LBmu <- sn.res[dm,]$q5
UBmu <- sn.res[dm,]$q95
Rhat_muSN <-  sn.res[dm,]$rhat

### the between-study variance
dtau2 <- grepl("tau_sqr",sn.res$variable)

tau2 <- sn.res[dtau2,]$median
LBtau2 <- sn.res[dtau2,]$q5
UBtau2 <- sn.res[dtau2,]$q95
Rhat_tau2SN <- sn.res[dtau2,]$rhat

## the between-study sd
tau <- sqrt(sn.res[dtau2,]$median)
LBtau <- sqrt(sn.res[dtau2,]$q5)
UBtau <- sqrt(sn.res[dtau2,]$q95)

## prediction
dxn <- grepl("pred1",sn.res$variable)
pred <- sn.res[dxn,]$median
LBpred <- sn.res[dxn,]$q5
UBpred <- sn.res[dxn,]$q95
Rhat_predSN <-  sn.res[dxn,]$rhat


sn.res1<-data.frame(median=median, 
                    lowerCI=LBmu,
                    upperCI=UBmu,
                    tau2=tau2,
                    l_tau2 = LBtau2,
                    u_tau2 = UBtau2,
                    tau=tau,
                    l_tau = LBtau,
                    u_tau = UBtau,
                    Rhat_muSN = Rhat_muSN,
                    Rhat_tau2SN = Rhat_tau2SN,
                    pred = pred,
                    LBpred = LBpred,
                    UBpred = UBpred,
                    Rhat_predSN = Rhat_predSN,
                    Pr_mu = Pr_mu,
                    Pr_pred1 = Pr_pred
)


write.csv(sn.res1 , "C:\\Users\\kanel\\OneDrive\\Υπολογιστής\\Simulated data application\\Normal_data\\SN\\pred_res_Normal_data_Bayesian_SN_model_HN.csv",row.names=FALSE )  
write.csv(deltaSN, "C:\\Users\\kanel\\OneDrive\\Υπολογιστής\\Simulated data application\\Normal_data\\SN\\pred_rel_effSN_Normal_data_Bayesian_SN_model_HN.csv",row.names=FALSE )  


# ################################## Binomial-SN(Unif) model #############################

library(cmdstanr)
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
check_cmdstan_toolchain()

fileSN_U <- file.path(cmdstan_path(), "examples", "Stanmodels", "Skew_normal_model_binUpred.stan")
modSN_U <- cmdstan_model(fileSN_U)

fitSN_U <- modSN_U$sample(
  data =list(ns = nrow(p[[4]]),
             r1 = p[[4]]$ci,
             r2 = p[[4]]$ti,
             n1 = p[[4]]$nci,
             n2 = p[[4]]$nti), 
  seed = 123480, 
  chains = 2, 
  parallel_chains = 2,
  refresh = 500,
  iter_warmup = 10000,
  iter_sampling = 50000,
  adapt_delta = 0.99
)

sn.res_U <- as.data.frame(fitSN_U$summary())

#### To visualize the predictive distribution
bayesplot::color_scheme_set("brightblue")
bayesplot::mcmc_dens(fitSN_U$draws(c("mu", "pred1")))

### Pr of mu <- 0
Pr_mu_U = mean(fitSN_U$draws("mu") < 0 )

### Pr of prediction <- 0
Pr_pred_U = mean(fitSN_U$draws("pred1") < 0 )

## the study specific effects
d_U <- grepl("delta",sn.res_U$variable)

deltaSN_U <- sn.res_U[d_U,]

## the mean treatment effect
dm_U <- grepl("mu",sn.res_U$variable)

median_U <- sn.res_U[dm_U,]$median
LBmu_U <- sn.res_U[dm_U,]$q5
UBmu_U <- sn.res_U[dm_U,]$q95
Rhat_muSN_U <-  sn.res_U[dm_U,]$rhat

## the between-study variance 
dtau2_U <- grepl("tau_sqr",sn.res_U$variable)

tau2_U <- sn.res_U[dtau2_U,]$median
LBtau2_U <- sn.res_U[dtau2_U,]$q5
UBtau2_U <- sn.res_U[dtau2_U,]$q95
Rhat_tau2SN_U <- sn.res_U[dtau2_U,]$rhat

### the between-study sd
tau_U <- sqrt(sn.res_U[dtau2_U,]$median)
LBtau_U <- sqrt(sn.res_U[dtau2_U,]$q5)
UBtau_U <- sqrt(sn.res_U[dtau2_U,]$q95)

dxn_U <- grepl("pred1",sn.res_U$variable)
pred_U <- sn.res_U[dxn_U,]$median
LBpred_U <- sn.res_U[dxn_U,]$q5
UBpred_U <- sn.res_U[dxn_U,]$q95
Rhat_predSN_U <-  sn.res_U[dxn_U,]$rhat


sn.res_U1<-data.frame(median=median_U, 
                      lowerCI=LBmu_U,
                      upperCI=UBmu_U,
                      tau2=tau2_U,
                      l_tau2 = LBtau2_U,
                      u_tau2 = UBtau2_U,
                      tau=tau_U,
                      l_tau = LBtau_U,
                      u_tau = UBtau_U,
                      Rhat_muSN = Rhat_muSN_U,
                      Rhat_tau2SN = Rhat_tau2SN_U,
                      pred = pred_U,
                      LBpred = LBpred_U,
                      UBpred = UBpred_U,
                      Rhat_predSN = Rhat_predSN_U,
                      Pr_mu = Pr_mu_U,
                      Pr_pred1 = Pr_pred_U
)


write.csv(sn.res_U1 , "C:\\Users\\Lela Panag\\Desktop\\Github\\GITHUB_FINAL\\Forest plots of 3 selected simulated datasets\\pred_res_Normal_data_Bayesian_SN_model_U.csv",row.names=FALSE )  
write.csv(deltaSN_U, "C:\\Users\\Lela Panag\\Desktop\\Github\\GITHUB_FINAL\\Forest plots of 3 selected simulated datasets\\pred_rel_effSN_Normal_data_Bayesian_SN_model_U.csv",row.names=FALSE )  


################## FREQUENTIST MODELS ####################
################## Binomial-Normal ###########

######## NORMAL RANDOM EFFECTS WITH METAFOR PACKAGE ########

set.seed(12349)
freq_norm1bin <- rma.glmm(measure = "OR",
                          ai = p[[4]]$ti, bi = p[[4]]$nti - p[[4]]$ti,
                          ci = p[[4]]$ci, di = p[[4]]$nci - p[[4]]$ci,
                          data = p[[4]],  model="UM.RS")

freq_normbin <- summary(freq_norm1bin)

mu_fr_norm.bin <- freq_normbin$beta   ## mean treatment effect
LB_mu_fr_norm.bin <- freq_normbin$ci.lb
UB_mu_fr_norm.bin <- freq_normbin$ci.ub
prec_mu_fr_norm.bin <- (UB_mu_fr_norm.bin - LB_mu_fr_norm.bin) 
tau2_fr_norm.bin <- freq_normbin$tau2  ## between-study variance
tau_fr_norm.bin <- sqrt(freq_normbin$tau)


fr_normal.bin <- cbind(mu_fr_norm.bin,LB_mu_fr_norm.bin, UB_mu_fr_norm.bin, tau2_fr_norm.bin, tau_fr_norm.bin,   prec_mu_fr_norm.bin)

write.csv(fr_normal.bin, "C:\\Users\\kanel\\OneDrive\\Υπολογιστής\\Simulated data application\\Normal_data\\norm_fr_bin\\Results_freq_Normal_bin_model_normal_data.csv",row.names=FALSE )  

################ Frequentist models######
###### Normal-normal model#############
######## NORMAL RANDOM EFFECTS WITH METAFOR PACKAGE ########
library(metafor)
set.seed(123499)
freq_norm1 <- rma(measure = "OR",
                  ai = p[[4]]$ti, bi = p[[4]]$nti - p[[4]]$ti,
                  ci = p[[4]]$ci, di = p[[4]]$nci - p[[4]]$ci,
                  data = p[[4]], method = "REML")
freq_norm <- summary(freq_norm1)
fr_norm = as.data.frame(confint(freq_norm1))

### prediction intervals for mean treatment effect ###
pred = predict.rma(freq_norm1)
pred_mu = pred$pred
LB_pred_mu = pred$pi.lb
UB_pred_mu = pred$pi.ub
prec_pred_mu = UB_pred_mu - LB_pred_mu

### study specific effects ####
stud_eff = blup( freq_norm, level = 95 )


write.csv(stud_eff ,   "C:\\Users\\kanel\\OneDrive\\Υπολογιστής\\Simulated data application\\Normal_data\\norm_fr\\Rel_eff_freq_Normal_normal_model_normal_data.csv",row.names=FALSE )  

mu_fr_norm <- freq_norm$beta  ### mean treatment effect
LB_mu_fr_norm <- freq_norm$ci.lb
UB_mu_fr_norm <- freq_norm$ci.ub
prec_mu_fr_norm <- (UB_mu_fr_norm - LB_mu_fr_norm) 


fr_tau2 <- fr_norm[1,]
tau2_fr_norm <- fr_tau2[1] ## between-study variance
LB_tau2_fr_norm <- fr_tau2[2]
UB_tau2_fr_norm <- fr_tau2[3]
prec_tau2_fr_norm <- (UB_tau2_fr_norm - LB_tau2_fr_norm) 

fr_tau <- fr_norm[2,]
tau_fr_norm <- fr_tau[1]
LB_tau_fr_norm <- fr_tau[2]
UB_tau_fr_norm <- fr_tau[3]
prec_tau_fr_norm <- (UB_tau_fr_norm - LB_tau_fr_norm) 


fr_normal <- cbind(mu_fr_norm,LB_mu_fr_norm, UB_mu_fr_norm, tau2_fr_norm, LB_tau2_fr_norm, UB_tau2_fr_norm, tau_fr_norm, LB_tau_fr_norm, UB_tau_fr_norm, 
                   prec_mu_fr_norm, prec_tau2_fr_norm, prec_tau_fr_norm,
                   pred_mu, LB_pred_mu, UB_pred_mu)

write.csv(fr_normal, "C:\\Users\\kanel\\OneDrive\\Υπολογιστής\\Simulated data application\\Normal_data\\norm_fr\\Results_freq_Normal_normal_model_normal_data.csv",row.names=FALSE )  


############################## METAPLUS PACKAGE ######################
library(metaplus)

######################## METAPLUS T-DISTRIBUTION ######################

dat1 = escalc(measure="OR",ai = p[[4]]$ti, bi =p[[4]]$nti - p[[4]]$ti,
              ci =p[[4]]$ci, di = p[[4]]$nci - p[[4]]$ci,
              data = p[[4]])


set.seed(1234999)
metaplust1 <- metaplus(yi =  dat1$yi, sei = sqrt(dat1$vi), data = dat1, random = "t-dist")
metaplust <- summary(metaplust1)

mu_metaplust <- metaplust$results[1]  ## mean treatment effect
LB_mu_metaplust <- metaplust$results[4]
UB_mu_metaplust <- metaplust$results[7]
prec_mu_metaplust <- UB_mu_metaplust  - LB_mu_metaplust

tau2_metaplust <- metaplust$results[2]  ## between-study variance
tau_metaplust <- sqrt(tau2_metaplust)

ind_normal <- metaplust$results[3] 



metaplus_t <- cbind(mu_metaplust ,LB_mu_metaplust, UB_mu_metaplust, tau2_metaplust, tau_metaplust,
                    prec_mu_metaplust, ind_normal)

write.csv(metaplus_t, "C:\\Users\\kanel\\OneDrive\\Υπολογιστής\\Simulated data application\\Normal_data\\METAPLUS_T\\Results_metaplust_Normal_data.csv",row.names=FALSE )  

######################## METAPLUS TMIXTURE ######################

set.seed(12349990)
metaplusmix1 <- metaplus(yi =  dat1$yi, sei = sqrt(dat1$vi), data = dat1, random = "mixture")
metaplusmix <- summary(metaplusmix1)


mu_metaplusmix <- metaplusmix$results[1]  ## mean treatment effect
LB_mu_metaplusmix <- metaplusmix$results[5]
UB_mu_metaplusmix <- metaplusmix$results[9]
prec_mu_metaplusmix <- UB_mu_metaplusmix - LB_mu_metaplusmix 

tau2_metaplusmix <- metaplusmix$results[2] ### homogeneous studies variance
tau_metaplusmix <- sqrt(tau2_metaplusmix)

tau2_out_metaplusmix <- metaplusmix$results[3] ### outlying studies variance
tau_out_metaplusmix <- sqrt(tau2_out_metaplusmix)

prob_out <- metaplusmix$results[4]

metaplus_mix <- cbind(mu_metaplusmix ,LB_mu_metaplusmix, UB_mu_metaplusmix, tau2_metaplusmix, tau_metaplusmix,
                      prec_mu_metaplusmix, tau2_out_metaplusmix , tau_out_metaplusmix ,prob_out )

write.csv(metaplus_mix, "C:\\Users\\kanel\\OneDrive\\Υπολογιστής\\Simulated data application\\Normal_data\\METAPLUS_MIXTURE\\Results_metaplust_Normal_data.csv",row.names=FALSE )  



############### SAVE THE RESULTS ############

NHN = read.csv( "C:\\Users\\kanel\\OneDrive\\Υπολογιστής\\Simulated data application\\Normal_data\\NORMAL\\res_Normal_data_Bayesian_Normal_model_HN.csv")
NU = read.csv( "C:\\Users\\kanel\\OneDrive\\Υπολογιστής\\Simulated data application\\Normal_data\\NORMAL\\res_Normal_data_Bayesian_Normal_model_Unif.csv" )
DP_UGn = read.csv( "C:\\Users\\kanel\\OneDrive\\Υπολογιστής\\Simulated data application\\Normal_data\\DP\\UGDP_res_Normal_data.csv")
DP_26_HN_U = read.csv( "C:\\Users\\kanel\\OneDrive\\Υπολογιστής\\Simulated data application\\Normal_data\\DP\\DP_HN_U_26_res_Normal_data.csv" )
DP_HN_U_51 = read.csv( "C:\\Users\\kanel\\OneDrive\\Υπολογιστής\\Simulated data application\\Normal_data\\DP\\DP_HN_U_51_res_Normal_data.csv")
DP_U_U_26 = read.csv( "C:\\Users\\kanel\\OneDrive\\Υπολογιστής\\Simulated data application\\Normal_data\\DP\\DP_U_U_26_res_Normal_data.csv")
DPU_U_51 = read.csv( "C:\\Users\\kanel\\OneDrive\\Υπολογιστής\\Simulated data application\\Normal_data\\DP\\DPU_U_51_res_Normal_data.csv" )
t_distrHN = read.csv( "C:\\Users\\kanel\\OneDrive\\Υπολογιστής\\Simulated data application\\Normal_data\\T\\pred_res_Normal_data_Bayesian_t_model_HN.csv")
t_distrU = read.csv( "C:\\Users\\kanel\\OneDrive\\Υπολογιστής\\Simulated data application\\Normal_data\\T\\pred_res_Normal_data_Bayesian_t_model_U.csv" )
SNHN = read.csv("C:\\Users\\kanel\\OneDrive\\Υπολογιστής\\Simulated data application\\Normal_data\\SN\\pred_res_Normal_data_Bayesian_SN_model_HN.csv" )
SNU = read.csv( "C:\\Users\\kanel\\OneDrive\\Υπολογιστής\\Simulated data application\\Normal_data\\SN\\pred_res_Normal_data_Bayesian_SN_model_U.csv")
fr_normal.bin = read.csv("C:\\Users\\kanel\\OneDrive\\Υπολογιστής\\Simulated data application\\Normal_data\\norm_fr_bin\\Results_freq_Normal_bin_model_normal_data.csv")
fr_normal = read.csv( "C:\\Users\\kanel\\OneDrive\\Υπολογιστής\\Simulated data application\\Normal_data\\norm_fr\\Results_freq_Normal_normal_model_normal_data.csv" )
metaplus_t = read.csv( "C:\\Users\\kanel\\OneDrive\\Υπολογιστής\\Simulated data application\\Normal_data\\METAPLUS_T\\Results_metaplust_Normal_data.csv" )
metaplus_mix = read.csv( "C:\\Users\\kanel\\OneDrive\\Υπολογιστής\\Simulated data application\\Normal_data\\METAPLUS_MIXTURE\\Results_metaplust_Normal_data.csv")



################## summary of the results #####################
r = 16
c = 10
matrix = matrix(nrow = r, ncol = c)
# Assign values to each cell in the matrix
for(i in 1:r){
  for(j in 1:c){
    matrix[1,1] <- "Models"
    
    matrix[1,2] <- "mu" 
    matrix[1,3] <- "lb_mu"
    matrix[1,4] <- "ub_mu"
    
    matrix[1,5] <- "tau2"
    matrix[1,6] <- "lb_tau2"
    matrix[1,7] <- "ub_tau2"
    
    matrix[1,8] <- "pred"
    matrix[1,9] <- "lb_pred"
    matrix[1,10] <- "ub_pred"
    
    matrix[2,1] <- "Binomial-Normal(HN)"
    matrix[3,1] <- "Binomial-Normal(Unif)"
    matrix[4,1] <- "Binomial-t(HN)"
    matrix[5,1] <- "Binomial-t(Unif)"
    matrix[6,1] <- "Binomial-SN(HN)"
    matrix[7,1] <- "Binomial-SN(Unif)"
    matrix[8,1] <- "Binomial-DP-26(HN/Unif)"
    matrix[9,1] <- "Binomial-DP-51(HN/Unif)"
    matrix[10,1] <- "Binomial-DP-26(Unif/Unif)"
    matrix[11,1] <- "Binomial-DP-51(Unif/Unif)"
    matrix[12,1] <- "Binomial-DP-n(Unif/Gamma)"
    matrix[13,1] <- "Binomial-Normal(ML)"
    matrix[14,1] <- "Normal-Normal(REML)"
    matrix[15,1] <- "Normal-t"
    matrix[16,1] <- "Common-mean-mixture"
    
    #### MU 
    matrix[2,2] <- NHN$mu_Normaljags.bin
    matrix[3,2] <- NU$mu_Normaljags.binU
    matrix[4,2] <- t_distrHN$median
    matrix[5,2] <- t_distrU$median
    matrix[6,2] <- SNHN$median
    matrix[7,2] <- SNU$median
    matrix[8,2] <- DP_26_HN_U$mu_DP_HN_U_26
    matrix[9,2] <- DP_HN_U_51$mu_DP_HN_U_51
    matrix[10,2] <- DP_U_U_26$mu_DP_U_U_26
    matrix[11,2] <- DPU_U_51$mu_DPU_U_51
    matrix[12,2] <- DP_UGn$mu_DP
    matrix[13,2] <- fr_normal.bin$X
    matrix[14,2] <- fr_normal$X
    matrix[15,2] <- metaplus_t$mu_metaplust
    matrix[16,2] <- metaplus_mix$mu_metaplusmix
    
    #### LB MU 
    matrix[2,3] <- NHN$LBmu.jags.bin
    matrix[3,3] <- NU$LBmu.jags.binU
    matrix[4,3] <- t_distrHN$lowerCI
    matrix[5,3] <- t_distrU$lowerCI
    matrix[6,3] <- SNHN$lowerCI
    matrix[7,3] <- SNU$lowerCI
    matrix[8,3] <- DP_26_HN_U$LB_mu_DP_HN_U_26
    matrix[9,3] <- DP_HN_U_51$LB_mu_DP_HN_U_51
    matrix[10,3] <- DP_U_U_26$LB_mu_DP_U_U_26
    matrix[11,3] <- DPU_U_51$LB_mu_DPU_U_51
    matrix[12,3] <- DP_UGn$LB_mu_DP
    matrix[13,3] <- fr_normal.bin$LB_mu_fr_norm.bin
    matrix[14,3] <- fr_normal$LB_mu_fr_norm
    matrix[15,3] <- metaplus_t$LB_mu_metaplust
    matrix[16,3] <-  metaplus_mix$LB_mu_metaplusmix
    
    #### UB MU 
    matrix[2,4] <- NHN$UBmu.jags.bin
    matrix[3,4] <- NU$UBmu.jags.binU
    matrix[4,4] <- t_distrHN$upperCI
    matrix[5,4] <- t_distrU$upperCI
    matrix[6,4] <- SNHN$upperCI
    matrix[7,4] <- SNU$upperCI
    matrix[8,4] <- DP_26_HN_U$UB_mu_DP_HN_U_26
    matrix[9,4] <- DP_HN_U_51$UB_mu_DP_HN_U_51
    matrix[10,4] <- DP_U_U_26$UB_mu_DP_U_U_26
    matrix[11,4] <- DPU_U_51$UB_mu_DPU_U_51
    matrix[12,4] <- DP_UGn$UB_mu_DP
    matrix[13,4] <- fr_normal.bin$UB_mu_fr_norm.bin
    matrix[14,4] <- fr_normal$UB_mu_fr_norm
    matrix[15,4] <- metaplus_t$UB_mu_metaplust
    matrix[16,4] <-  metaplus_mix$UB_mu_metaplusmix
    
    #### tau2 
    matrix[2,5] <- NHN$tau2_Normaljags.bin
    matrix[3,5] <- NU$tau2_Normaljags.binU
    matrix[4,5] <- t_distrHN$tau2
    matrix[5,5] <- t_distrU$tau2
    matrix[6,5] <- SNHN$tau2
    matrix[7,5] <- SNU$tau2
    matrix[8,5] <- DP_26_HN_U$tau2_DP_HN_U_26
    matrix[9,5] <- DP_HN_U_51$tau2_DP_HN_U_51
    matrix[10,5] <- DP_U_U_26$tau2_DP_U_U_26
    matrix[11,5] <- DPU_U_51$tau2_DPU_U_51
    matrix[12,5] <- DP_UGn$tau2_DP
    matrix[13,5] <- fr_normal.bin$tau2_fr_norm.bin
    matrix[14,5] <- fr_normal$tau2_fr_norm
    matrix[15,5] <- metaplus_t$tau2_metaplust
    matrix[16,5] <- metaplus_mix$tau2_out_metaplusmix
    
    #### LB tau2 
    matrix[2,6] <- NHN$LBtau2.jags.bin
    matrix[3,6] <- NU$LBtau2.jags.binU
    matrix[4,6] <- t_distrHN$l_tau2
    matrix[5,6] <- t_distrU$l_tau2
    matrix[6,6] <- SNHN$l_tau2
    matrix[7,6] <- SNU$l_tau2
    matrix[8,6] <- DP_26_HN_U$LB_tau2_DP_HN_U_26
    matrix[9,6] <- DP_HN_U_51$LB_tau2_DP_HN_U_51
    matrix[10,6] <- DP_U_U_26$LB_tau2_DP_U_U_26
    matrix[11,6] <- DPU_U_51$LB_tau2_DPU_U_51
    matrix[12,6] <- DP_UGn$LB_tau2_DP
    matrix[13,6] <- NA
    matrix[14,6] <- fr_normal$LB_tau2_fr_norm
    matrix[15,6] <- NA
    matrix[16,6] <- NA
    
    #### UB tau2 
    matrix[2,7] <- NHN$UBtau2.jags.bin
    matrix[3,7] <- NU$UBtau2.jags.binU
    matrix[4,7] <- t_distrHN$u_tau2
    matrix[5,7] <- t_distrU$u_tau2
    matrix[6,7] <- SNHN$u_tau2
    matrix[7,7] <- SNU$u_tau2
    matrix[8,7] <- DP_26_HN_U$UB_tau2_DP_HN_U_26
    matrix[9,7] <- DP_HN_U_51$UB_tau2_DP_HN_U_51
    matrix[10,7] <- DP_U_U_26$UB_tau2_DP_U_U_26
    matrix[11,7] <- DPU_U_51$UB_tau2_DPU_U_51
    matrix[12,7] <- DP_UGn$UB_tau2_DP
    matrix[13,7] <- NA
    matrix[14,7] <- fr_normal$UB_tau2_fr_norm
    matrix[15,7] <- NA
    matrix[16,7] <- NA
    
    #### pred 
    matrix[2,8] <- NHN$delta_new_Normaljags.bin
    matrix[3,8] <- NU$delta_new_Normaljags.binU
    matrix[4,8] <- t_distrHN$delta_new
    matrix[5,8] <- t_distrU$delta_new
    matrix[6,8] <- SNHN$pred
    matrix[7,8] <- SNU$pred
    matrix[8,8] <- NA
    matrix[9,8] <- NA
    matrix[10,8] <- NA
    matrix[11,8] <- NA
    matrix[12,8] <- NA
    matrix[13,8] <- NA
    matrix[14,8] <- fr_normal$pred_mu
    matrix[15,8] <- NA
    matrix[16,8] <- NA
    
    #### lb pred 
    matrix[2,9] <- NHN$LBdelta_new.jags.bin
    matrix[3,9] <- NU$LBdelta_new.jags.binU
    matrix[4,9] <- t_distrHN$LBdelta_new
    matrix[5,9] <- t_distrU$LBdelta_new
    matrix[6,9] <- SNHN$LBpred
    matrix[7,9] <- SNU$LBpred
    matrix[8,9] <- NA
    matrix[9,9] <-NA
    matrix[10,9] <- NA
    matrix[11,9] <- NA
    matrix[12,9] <- NA
    matrix[13,9] <- NA
    matrix[14,9] <- fr_normal$LB_pred_mu
    matrix[15,9] <- NA
    matrix[16,9] <- NA  
    
    #### UB pred 
    matrix[2,10] <- NHN$UBdelta_new.jags.bin
    matrix[3,10] <- NU$UBdelta_new.jags.binU
    matrix[4,10] <- t_distrHN$UBdelta_new
    matrix[5,10] <- t_distrU$UBdelta_new
    matrix[6,10] <- SNHN$UBpred
    matrix[7,10] <- SNU$UBpred
    matrix[8,10] <- NA
    matrix[9,10] <-NA
    matrix[10,10] <- NA
    matrix[11,10] <- NA
    matrix[12,10] <- NA
    matrix[13,10] <- NA
    matrix[14,10] <- fr_normal$UB_pred_mu
    matrix[15,10] <- NA
    matrix[16,10] <- NA  
  }
}


write.csv(matrix, "C:\\Users\\kanel\\OneDrive\\Υπολογιστής\\Simulated data application\\Normal_data\\Results_matrix_normal_data.csv",row.names=FALSE )  

