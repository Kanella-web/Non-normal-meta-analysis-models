####################### Scenario 13: Skew Normal data, k=26, mu=0, tau2= 0.12 #######################
### DATA GENERATING PROCESS ####
library(sn)

N.sim = 1000
df <- list()
dich_data <- function(k, params){
  
  mati <- matrix(1:k,nrow = k, ncol=2) 
  
  for (i in 1:N.sim){
    
    thetai <- rsn(k, dp=params)
    
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
set.seed(610)
### In order to convert the mu, tau2 in location and scale parameters
params <- cp2dp(c(0 , sqrt(0.12), 0.785), "SN")
m <- dich_data(26,params)
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
d

count <- c()
for(i in 1:N.sim){
  count[i] <- sum(d[[i]])
}
count
studies_zero_2arms <- which(count == 1)
studies_zero_2arms


mu = 0
tau2 = 0.12

# ### DENSITY FUNCTION OF THE TRUE STUDY-SPECIFIC EFFECTS ####
### SELECTED DATA SET ###
plot(density(p[[1]]$true_eff), main = "Skewed Normal distribution for true treatment effects", lwd= 3, xlim = c(-6,6))
write.csv(p[[1]]$true_eff , "C:\\Users\\Lela Panag\\Desktop\\SIMULATION SCENARIOS FINAL\\Datasets application\\Simulated data application\\Skew_Normal_data\\true_eff\\true_eff_skew_normal.csv",row.names=FALSE )  

################ BAYESIAN M-A ################
library(meta)
library(metafor)
library(rjags)
library(R2jags)

############ Binomial-Normal(HN) model ##############
set.seed(11221)
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

dat<-list(ns = nrow(p[[1]]),
          r1 = p[[1]]$ci,
          r2 = p[[1]]$ti,
          n1 = p[[1]]$nci,
          n2 = p[[1]]$nti )

run.modelN = jags(
  data = dat,
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
tau2_Normaljags.bin <- (m2$`50%`)^2
LBtau2.jags.bin <- (m2$`2.5%`)^2
UBtau2.jags.bin <- (m2$`97.5%`)^2
Rhat_tauN <- m2$Rhat
precN_tau2 <-  UBtau2.jags.bin - LBtau2.jags.bin

fd <- grepl("delta12", row.names(resultsN))
md <- resultsN[fd,]
rel_eff <- md$`50%`
LB_rel_eff <- md$`2.5%`
UB_rel_eff <- md$`97.5%`
sd_rel_eff <- md$sd
Rhat_deltaN <- md$Rhat

lista <- cbind.data.frame(rel_eff,sd_rel_eff,LB_rel_eff, UB_rel_eff, Rhat_deltaN)


RhatsN <- cbind.data.frame(mu_Normaljags.bin, LBmu.jags.bin, UBmu.jags.bin, tau2_Normaljags.bin,
                           LBtau2.jags.bin, UBtau2.jags.bin,Rhat_muN, Rhat_tauN,
                           precN_mu, precN_tau2, delta_new_Normaljags.bin, 
                           LBdelta_new.jags.bin , UBdelta_new.jags.bin, Pr_mu, Pr_pred)

write.csv(RhatsN , "C:\\Users\\Lela Panag\\Desktop\\SIMULATION SCENARIOS FINAL\\Datasets application\\Simulated data application\\Skew_normal_data\\NORMAL\\res_SN_data_Bayesian_NORMAL_HN_model.csv",row.names=FALSE )  
write.csv(lista , "C:\\Users\\Lela Panag\\Desktop\\SIMULATION SCENARIOS FINAL\\Datasets application\\Simulated data application\\Skew_normal_data\\NORMAL\\rel_eff_SN_data_Bayesian_NORMAL_HN_model.csv",row.names=FALSE )  

####################### Binomial-DP-26(HN/Unif) model ######################
set.seed(11222)
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
  
  for(i in 1:ns){
   rho[i] <- (exp(delta12[i]) / (1+exp(delta12[i])))
  }
  
  for(i in 1:ns){
   for(j in 1:ns){
    equalsmatrix[i,j]<-equals(rho[i],rho[j])
   }
  equalsres[i]<-sum(equalsmatrix[i,])
  }
  
}",file="DPmodel2.bin_HN_U_26.txt")
modfile = 'DPmodel2.bin_HN_U_26.txt'


run.modelDP = jags(
  data =list(ns = nrow(p[[1]]),
             r1 = p[[1]]$ci,
             r2 = p[[1]]$ti,
             n1 = p[[1]]$nci,
             n2 = p[[1]]$nti,
             N= 26
  ) ,  
  inits = NULL,
  parameters.to.save = c(
    "basemu",    ## tau2 of the base Normal distribution
    "basetau",   ## mu of the base Normal distribution
    "poptrue",   ## overall mu
    "var.true", ## the heterogeneity 
    "delta12", #### the random effects
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
ncol(run.modelDP$BUGSoutput$sims.matrix)

Pr_mu = mean(run.modelDP$BUGSoutput$sims.matrix[  ,"poptrue"] < 0 )

DP_results <-round(DPresults, digits=2) 


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

split_prob_list <- split_dataframe(prob_list, 26)

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

cluster1 = positions[["1"]]
cluster2 = positions[["2"]]
cluster3 = positions[["4"]] 


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
precDP_tau2 <- UB_tau2_DP -  LB_tau2_DP


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

cluster_mean = cbind.data.frame(mean_cl1, mean_cl2, mean_cl3)

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
                            precDP_mu, precDP_tau2, prec_basemu, prec_basetau2,prec_alpha,  alpha , LB_alpha, UB_alpha,
                            Rhat_muDP, Rhat_tauDP, Rhat_basemu, Rhat_basetau, Rhat_alpha, numclus, LB_K, UB_K, Pr_mu)


extra_col = row.names(DPresults)
DPresults$extra_col = extra_col
prob = round(max_values,2)

write.csv(listaDP1 , "C:\\Users\\Lela Panag\\Desktop\\SIMULATION SCENARIOS FINAL\\Datasets application\\Simulated data application\\Skew_normal_data\\DP\\rel_eff_HN_U26_SN_data_Bayesian_DP_model.csv",row.names=FALSE )  
write.csv( DPresults, "C:\\Users\\Lela Panag\\Desktop\\SIMULATION SCENARIOS FINAL\\Datasets application\\Simulated data application\\Skew_normal_data\\DP\\res_all_HN_U26_SN_data_Bayesian_DP_model.csv",row.names=FALSE )  
write.csv( RhatsDP , "C:\\Users\\Lela Panag\\Desktop\\SIMULATION SCENARIOS FINAL\\Datasets application\\Simulated data application\\Skew_normal_data\\DP\\res_HN_U26_SN_data_Bayesian_DP_model.csv",row.names=FALSE )  
write.csv( prob, "C:\\Users\\Lela Panag\\Desktop\\SIMULATION SCENARIOS FINAL\\Datasets application\\Simulated data application\\Skew_normal_data\\DP\\prob_HN_U26_res_SN_data_Bayesian_DP_model.csv",row.names=FALSE )  
write.csv( cluster_mean, "C:\\Users\\Lela Panag\\Desktop\\SIMULATION SCENARIOS FINAL\\Datasets application\\Simulated data application\\Skew_normal_data\\DP\\cluster_mean_HN_U26_res_SN_data_Bayesian_DP_model.csv",row.names=FALSE )  


# ####################### Binomial-t(Unif) model ########################

library(cmdstanr)
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
check_cmdstan_toolchain()

filet_distr <- file.path(cmdstan_path(), "examples", "Stanmodels",  "t_model_binUpred.stan")

modt_distr <- cmdstan_model(filet_distr)

fitt_distr <- modt_distr$sample(
  data = list(ns = nrow(p[[1]]),
              r1 = p[[1]]$ci,
              r2 = p[[1]]$ti,
              n1 = p[[1]]$nci,
              n2 = p[[1]]$nti), 
  seed = 11223, 
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

### Pr of pred1 <- 0
Pr_pred = mean(fitt_distr$draws("pred") < 0 )


d <- grepl("delta",t_distr.res$variable)

deltat <- t_distr.res[d,]


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

ddn <- grepl("pred",t_distr.res$variable)

delta_new <- t_distr.res[ddn,]$median
LBdelta_new <- t_distr.res[ddn,]$q5
UBdelta_new <- t_distr.res[ddn,]$q95
Rhat_delta_newt_distr <-  t_distr.res[ddn,]$rhat

t_distr.res1<-data.frame(median=median, 
                         lowerCI=LBmu,
                         upperCI=UBmu,
                         tau2=tau2,
                         l_tau2 = LBtau2,
                         u_tau2 = UBtau2,
                         Rhat_mut_distr = Rhat_mut_distr,
                         Rhat_tau2t_distr = Rhat_tau2t_distr,
                         delta_new = delta_new,
                         LBdelta_new = LBdelta_new,
                         UBdelta_new = UBdelta_new,
                         Rhat_delta_newt_distr = Rhat_delta_newt_distr,
                         Pr_mu = Pr_mu,
                         Pr_pred = Pr_pred
)

# write.csv(t_distr.res1 , "C:\\Users\\Lela Panag\\Desktop\\SIMULATION SCENARIOS FINAL\\Datasets application\\Simulated data application\\Skew_normal_data\\T\\SN_data_Bayesian_T_model_U.csv",row.names=FALSE )  
# write.csv(deltat , "C:\\Users\\Lela Panag\\Desktop\\SIMULATION SCENARIOS FINAL\\Datasets application\\Simulated data application\\Skew_normal_data\\T\\rel_effT_SN_data_Bayesian_T_model_U.csv",row.names=FALSE )  
# 
write.csv(t_distr.res1 , "C:\\Users\\Lela Panag\\Desktop\\SIMULATION SCENARIOS FINAL\\Datasets application\\Simulated data application\\Skew_normal_data\\T\\pred_SN_data_Bayesian_T_model_U.csv",row.names=FALSE )  
write.csv(deltat , "C:\\Users\\Lela Panag\\Desktop\\SIMULATION SCENARIOS FINAL\\Datasets application\\Simulated data application\\Skew_normal_data\\T\\pred_rel_effT_SN_data_Bayesian_T_model_U.csv",row.names=FALSE )  



################################## Binomial-SN(Unif) model #############################

library(cmdstanr)
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
check_cmdstan_toolchain()

fileSN <- file.path(cmdstan_path(), "examples", "Stanmodels", "Skew_normal_model_binUpred.stan")
modSN <- cmdstan_model(fileSN)

fitSN <- modSN$sample(
  data =list(ns = nrow(p[[1]]),
             r1 = p[[1]]$ci,
             r2 = p[[1]]$ti,
             n1 = p[[1]]$nci,
             n2 = p[[1]]$nti), 
  seed = 11224, 
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

### Pr of pred1 <- 0
Pr_pred = mean(fitSN$draws("pred1") < 0 )

d <- grepl("delta",sn.res$variable)

deltaSN <- sn.res[d,]


dm <- grepl("mu",sn.res$variable)

median <- sn.res[dm,]$median
LBmu <- sn.res[dm,]$q5
UBmu <- sn.res[dm,]$q95
Rhat_muSN <-  sn.res[dm,]$rhat


dtau2 <- grepl("tau_sqr",sn.res$variable)

tau2 <- sn.res[dtau2,]$median
LBtau2 <- sn.res[dtau2,]$q5
UBtau2 <- sn.res[dtau2,]$q95
Rhat_tau2SN <- sn.res[dtau2,]$rhat

dx <- grepl("xi",sn.res$variable)

xi <- sn.res[dx, ]$median
LBxi <- sn.res[dx,]$q5
UBxi <- sn.res[dx,]$q95
Rhat_SN <-  sn.res[dx,]$rhat

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

sn.res1<-data.frame(median=median, 
                    lowerCI=LBmu,
                    upperCI=UBmu,
                    tau2=tau2,
                    l_tau2 = LBtau2,
                    u_tau2 = UBtau2,
                    skew = skew,
                    l_skew = LBskew,
                    u_skew = UBskew,
                    Rhat_muSN = Rhat_muSN,
                    Rhat_tau2SN = Rhat_tau2SN,
                    Rhat_skewSN = Rhat_skewSN,
                    pred = pred,
                    LBpred = LBpred,
                    UBpred = UBpred,
                    Rhat_predSN = Rhat_predSN,
                    Pr_mu = Pr_mu,
                    Pr_pred = Pr_pred
)

# write.csv(sn.res1 , "C:\\Users\\Lela Panag\\Desktop\\SIMULATION SCENARIOS FINAL\\Datasets application\\Simulated data application\\Skew_normal_data\\SN\\SN_data_Bayesian_SN_model_U.csv",row.names=FALSE )  
# write.csv(deltaSN , "C:\\Users\\Lela Panag\\Desktop\\SIMULATION SCENARIOS FINAL\\Datasets application\\Simulated data application\\Skew_normal_data\\SN\\rel_effSN_SN_data_Bayesian_SN_model_U.csv",row.names=FALSE )  

write.csv(sn.res1 , "C:\\Users\\Lela Panag\\Desktop\\SIMULATION SCENARIOS FINAL\\Datasets application\\Simulated data application\\Skew_normal_data\\SN\\pred_SN_data_Bayesian_SN_model_U.csv",row.names=FALSE )  
write.csv(deltaSN , "C:\\Users\\Lela Panag\\Desktop\\SIMULATION SCENARIOS FINAL\\Datasets application\\Simulated data application\\Skew_normal_data\\SN\\pred_rel_effSN_SN_data_Bayesian_SN_model_U.csv",row.names=FALSE )  

################################# Binomial-Normal(HN) model with STAN #############################

library(cmdstanr)
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
check_cmdstan_toolchain()

fileN <- file.path(cmdstan_path(), "examples", "Stanmodels", "normal_model_bin_HN.stan")
modN <- cmdstan_model(fileN)

fitN <- modN$sample(
  data =list(ns = nrow(p[[1]]),
             r1 = p[[1]]$ci,
             r2 = p[[1]]$ti,
             n1 = p[[1]]$nci,
             n2 = p[[1]]$nti), 
  seed = 11225, 
  chains = 2, 
  parallel_chains = 2,
  refresh = 500,
  iter_warmup = 10000,
  iter_sampling = 50000,
  adapt_delta = 0.99
)

n.res <- as.data.frame(fitN$summary())

d <- grepl("delta",n.res$variable)

deltaN <- n.res[d,]


dm <- grepl("mu",n.res$variable)

median <- n.res[dm,]$median
LBmu <- n.res[dm,]$q5
UBmu <- n.res[dm,]$q95
Rhat_muN <-  n.res[dm,]$rhat

ddn <- grepl("pred",n.res$variable)

delta_new <- n.res[ddn,]$median
LBdelta_new <- n.res[ddn,]$q5
UBdelta_new <- n.res[ddn,]$q95
Rhat_delta_newN <-  n.res[ddn,]$rhat

tau2 <- c()
dtau2 <- grepl("tau_sqr",n.res$variable)

tau2 <- n.res[dtau2,]$median
LBtau2 <- n.res[dtau2,]$q5
UBtau2 <- n.res[dtau2,]$q95
Rhat_tau2N <- n.res[dtau2,]$rhat


n.res1<-data.frame(median=median, 
                   lowerCI=LBmu,
                   upperCI=UBmu,
                   tau2=tau2,
                   l_tau2 = LBtau2,
                   u_tau2 = UBtau2,
                   Rhat_muSN = Rhat_muN,
                   Rhat_tau2SN = Rhat_tau2N,
                   delta_new = delta_new,
                   LBdelta_new = LBdelta_new,
                   UBdelta_new = UBdelta_new,
                   Rhat_delta_newN = Rhat_delta_newN
)

write.csv(n.res1 , "C:\\Users\\Lela Panag\\Desktop\\SIMULATION SCENARIOS FINAL\\Datasets application\\Simulated data application\\Skew_Normal_data\\norm_stan\\res_Normal_data_Bayesian_N_model_HN.csv",row.names=FALSE )  
write.csv(deltaN, "C:\\Users\\Lela Panag\\Desktop\\SIMULATION SCENARIOS FINAL\\Datasets application\\Simulated data application\\Skew_Normal_data\\norm_stan\\rel_effN_Normal_data_Bayesian_N_model_HN.csv",row.names=FALSE )  
