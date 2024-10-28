
############ FOREST PLOTS OF THE NORMAL DATA APPLICATION############
############# FOREST PLOT FOR BINOMIAL-NORMAL(HN) MODEL ##########
NormalHN_mod_rel_eff = read.csv("C:\\Users\\Lela Panag\\Desktop\\SIMULATION SCENARIOS FINAL\\Datasets application\\Simulated data application\\Normal_data\\NORMAL\\rel_effN_Normal_data_Bayesian_Normal_model_HN.csv")
NormalHN_mod_es = read.csv("C:\\Users\\Lela Panag\\Desktop\\SIMULATION SCENARIOS FINAL\\Datasets application\\Simulated data application\\Normal_data\\NORMAL\\res_Normal_data_Bayesian_Normal_model_HN.csv")

extra_dataNHN = NormalHN_mod_rel_eff$rel_eff
lb_NHN = NormalHN_mod_rel_eff$LB_rel_eff
ub_NHN = NormalHN_mod_rel_eff$UB_rel_eff


p_true = read.csv("C:\\Users\\Lela Panag\\Desktop\\SIMULATION SCENARIOS FINAL\\Datasets application\\Simulated data application\\Normal_data\\true_eff\\true_eff_normal.csv")
p_true = round(p_true, 2)
p_true = as.vector(p_true)

add_data = cbind.data.frame(p_true) 

forest(x = extra_dataNHN, 
       ci.lb = lb_NHN, 
       ci.ub = ub_NHN, 
       xlim = c(-7.2,11), 
       psize = 2,
       cex = 1.3,
       ilab = add_data, 
       ilab.xpos = 8, 
       lwd = 3,
       ylim =c(-1,29))

text(c(-6.7, 8, 9.9), 28, c("Studies" ,"true","Estimate[95% CI]"),   font=2, cex=1.5)
text(0, 29, c("Binomial-normal(HN)"),   font=2, cex=1.5)

addpoly(x= NormalHN_mod_es$mu_Normaljags.bin, ci.lb = NormalHN_mod_es$LBmu.jags.bin,ci.ub = NormalHN_mod_es$UBmu.jags.bin, rows=-1)
abline(h=0, lwd=0.1, col="black", lty=1)

# Open the TIFF device
tiff("forest_Normal_NHN.tiff", units="in", width=15, height=9, res=300)

forest(x = extra_dataNHN, 
       ci.lb = lb_NHN, 
       ci.ub = ub_NHN, 
       xlim = c(-7.2,11), 
       psize = 2,
       cex = 1.3,
       ilab = add_data, 
       ilab.xpos = 8, 
       lwd = 3,
       ylim =c(-1,29))

text(c(-6.6, 8, 9.7), 28, c("Studies" ,"true","Estimate[95% CI]"),   font=2, cex=1.5)
text(0, 29, c("Binomial-normal(HN)"),   font=2, cex=1.5)

addpoly(x= NormalHN_mod_es$mu_Normaljags.bin, ci.lb = NormalHN_mod_es$LBmu.jags.bin,ci.ub = NormalHN_mod_es$UBmu.jags.bin, rows=-1)
abline(h=0, lwd=0.1, col="black", lty=1)


# Close the TIFF device
dev.off()

######################## FOREST PLOT OF DP ################

DPUG_mod_rel_eff = read.csv("C:\\Users\\Lela Panag\\Desktop\\SIMULATION SCENARIOS FINAL\\Datasets application\\Simulated data application\\Normal_data\\DP\\UGDP_rel_eff_Normal_data.csv")
DPUG_mod_es = read.csv("C:\\Users\\Lela Panag\\Desktop\\SIMULATION SCENARIOS FINAL\\Datasets application\\Simulated data application\\Normal_data\\DP\\UGDP_res_Normal_data.csv")

extra_dataDPUG = DPUG_mod_rel_eff$rel_effDP
lb_DPUG = DPUG_mod_rel_eff$LB_rel_effDP
ub_DPUG = DPUG_mod_rel_eff$UB_rel_effDP
# 
 

mean_cluster = read.csv("C:\\Users\\Lela Panag\\Desktop\\SIMULATION SCENARIOS FINAL\\Datasets application\\Simulated data application\\Normal_data\\DP\\UGDP_cluster_mean_Normal_data.csv")  
mean_cluster1 = mean_cluster$ mean_cl1
mean_cluster2 = mean_cluster$ mean_cl2
mean_cluster3 = mean_cluster$ mean_cl3
mean_cluster4 = mean_cluster$ mean_cl7

 p_true = read.csv("C:\\Users\\Lela Panag\\Desktop\\SIMULATION SCENARIOS FINAL\\Datasets application\\Simulated data application\\Normal_data\\true_eff\\true_eff_normal.csv")
 p_true = round(p_true, 2)
 p_true = as.vector(p_true)
 
 # clusters = c("7/8/2/4","1","1","3","1/5/6","1","2","3/1","3","3","1/3","1","1","2","3/2/4","2/1/4",
 #              "1","3/1","2/1","2","2",
 #              "2","1","1/2/5","2","1/7")
 # 
 clusters = c("4/8/2/7","1","1","3","1/5/6","1","2","3/1","3","3","1/3","1","1","2","3/2/7","2/1/7",
              "1","3/1","2/1","2","2",
              "2","1","1/2/5","2","1/4")
 
 pr = c("0.07","0.13","0.13","0.10","0.09/0.08/0.08","0.12","0.13","0.09","0.10","0.09","0.09","0.14",
        "0.13","0.12","0.09/0.09/0.08","0.12/0.11","0.10","0.11/0.10","0.12","0.14","0.13","0.13","0.11","0.10/0.09/0.09","13","0.08")
 
 
 add_data = cbind.data.frame(clusters, pr,p_true) 
 
 forest(x = extra_dataDPUG, 
        ci.lb = lb_DPUG, 
        ci.ub = ub_DPUG, 
        xlim = c(-7.2,11), 
        psize = 2,
        cex = 1.3,
        ilab = add_data, 
        ilab.xpos = c(5.3,6.9,8.3), 
        lwd = 3,
        col = c("black", "red","red","black","red","red","black","black","black","black","red","red","red","black","black","black",
                "red", "black","black","black","black","black",
                "red","red","black","red"),
        pch = c(16,15,15,16,15,15,16,16,16,16,15,15,15,16,16,16,15,16,16,16,16,16,15,15,16,15),
        ylim =c(-1,29))
 
 text(c(-6.7, 5.2,7,8.3, 9.9), 28, c("Studies" ,"clusters", "probabilities", "true","Estimate[95% CI]"),   font=2, cex=1.5)
 text(0, 29, "Binomial-DP-n(Unif/Gamma)",   font=2, cex=1.5)
 addpoly(x= DPUG_mod_es$mu_DP, ci.lb = DPUG_mod_es$LB_mu_DP ,ci.ub = DPUG_mod_es$UB_mu_DP , rows=-1)
 abline(h=0, lwd=0.1, col="black", lty=1)

 line_positions <- c(mean_cluster1, mean_cluster2, mean_cluster3, mean_cluster4)
 line_colors <- c("red", "black",  "black", "black")
 
 for (i in seq_along(line_positions)) {
   segments(x0 = line_positions[i], y0 = -4, x1 = line_positions[i], y1 = length(extra_dataDPUG) + 1, 
            lwd = 3, col = line_colors[i], lty = 4)
 }
 
 legend(x = -1.5, y= 29, legend = c("Mean cluster1","Mean cluster2", "Mean cluster3","Mean cluster4"),
        col = c("red", "black", "black", "black"), lwd = c(3,3,3,3),
        lty = c(2,2,2,2),
        cex = 1.3 , y.intersp = 1, box.lwd = 1, box.lty = 1, box. = 10, xjust = 0)
 
 # Open the TIFF device
 tiff("forest_Normal_DP_UGn.tiff", units="in", width=15, height=9, res=300)
 
 forest(x = extra_dataDPUG, 
        ci.lb = lb_DPUG, 
        ci.ub = ub_DPUG, 
        xlim = c(-7.2,11), 
        psize = 2,
        cex = 1.3,
        ilab = add_data, 
        ilab.xpos = c(5.3,6.9,8.3), 
        lwd = 3,
        col = c("black", "red","red","black","red","red","black","black","black","black","red","red","red","black","black","black",
                "red", "black","black","black","black","black",
                "red","red","black","red"),
        pch = c(16,15,15,16,15,15,16,16,16,16,15,15,15,16,16,16,15,16,16,16,16,16,15,15,16,15),
        ylim =c(-1,29))
 
 text(c(-6.6, 5.2,7,8.3, 9.7), 28, c("Studies" ,"clusters", "probabilities", "true","Estimate[95% CI]"),   font=2, cex=1.3)
 text(0, 29, "Binomial-DP-n(Unif/Gamma)",   font=2, cex=1.5)
 addpoly(x= DPUG_mod_es$mu_DP, ci.lb = DPUG_mod_es$LB_mu_DP ,ci.ub = DPUG_mod_es$UB_mu_DP , rows=-1)
 abline(h=0, lwd=0.1, col="black", lty=1)
 
 line_positions <- c(mean_cluster1, mean_cluster2, mean_cluster3, mean_cluster4)
 line_colors <- c("red", "black",  "black", "black")
 
 for (i in seq_along(line_positions)) {
   segments(x0 = line_positions[i], y0 = -4, x1 = line_positions[i], y1 = length(extra_dataDPUG) + 1, 
            lwd = 3, col = line_colors[i], lty = 4)
 }
 
 # Close the TIFF device
 dev.off()
############ FOREST PLOT OF BINOMIAL-T(HN) MODEL ########### 
# t_mod_rel_eff = read.csv("C:\\Users\\Lela Panag\\Desktop\\SIMULATION SCENARIOS FINAL\\Datasets application\\Simulated data application\\Normal_data\\T\\rel_effT_Normal_data_Bayesian_t_model_HN.csv")
# t_mod_es = read.csv("C:\\Users\\Lela Panag\\Desktop\\SIMULATION SCENARIOS FINAL\\Datasets application\\Simulated data application\\Normal_data\\T\\res_Normal_data_Bayesian_t_model_HN.csv")

t_mod_rel_eff = read.csv("C:\\Users\\Lela Panag\\Desktop\\SIMULATION SCENARIOS FINAL\\Datasets application\\Simulated data application\\Normal_data\\T\\pred_rel_effT_Normal_data_Bayesian_t_model_HN.csv")
t_mod_es = read.csv("C:\\Users\\Lela Panag\\Desktop\\SIMULATION SCENARIOS FINAL\\Datasets application\\Simulated data application\\Normal_data\\T\\pred_res_Normal_data_Bayesian_t_model_HN.csv")


extra_data2 = t_mod_rel_eff$median
lb_t = t_mod_rel_eff$q5
ub_t = t_mod_rel_eff$q95

p_true = read.csv("C:\\Users\\Lela Panag\\Desktop\\SIMULATION SCENARIOS FINAL\\Datasets application\\Simulated data application\\Normal_data\\true_eff\\true_eff_normal.csv")
p_true = round(p_true, 2)
p_true = as.vector(p_true)

add_data = cbind.data.frame(p_true) 

forest(x = extra_data2, 
       ci.lb = lb_t, 
       ci.ub = ub_t, 
       psize = 2,
       xlim = c(-7.2,11), 
       cex = 1.3,
       ilab = add_data, 
       ilab.xpos = 8, 
       lwd = 3,
       ylim =c(-1,29))

text(c(-6.7, 8, 9.9), 28, c("Studies" ,"true","Estimate[95% CI]"),   font=2, cex=1.5)
text(0, 29, c("Binomial-t(HN)"),   font=2, cex=1.5)
addpoly(x= t_mod_es$median , ci.lb =t_mod_es$lowerCI, ci.ub = t_mod_es$upperCI,rows=-1)

abline(h=0, lwd=0.1, col="black", lty=1)

# Open the TIFF device
tiff("forest_Normal_t.tiff", units="in", width=15, height=9, res=300)


forest(x = extra_data2, 
       ci.lb = lb_t, 
       ci.ub = ub_t, 
       psize = 2,
       xlim = c(-7.2,11), 
       cex = 1.3,
       ilab = add_data, 
       ilab.xpos = 8, 
       lwd = 3,
       ylim =c(-1,29))

text(c(-6.6, 8, 9.7), 28, c("Studies" ,"true","Estimate[95% CI]"),   font=2, cex=1.5)
text(0, 29, c("Binomial-t(HN)"),   font=2, cex=1.5)
addpoly(x= t_mod_es$median , ci.lb =t_mod_es$lowerCI, ci.ub = t_mod_es$upperCI,rows=-1)

abline(h=0, lwd=0.1, col="black", lty=1)

# Close the TIFF device
dev.off()

################## FOR THE BINOMIAL-SKEW NORMAL(HN) MODEL #############
# sn_mod_rel_eff = read.csv("C:\\Users\\Lela Panag\\Desktop\\SIMULATION SCENARIOS FINAL\\Datasets application\\Simulated data application\\Normal_data\\SN\\rel_effSN_Normal_data_Bayesian_SN_model_HN.csv")
# sn_mod_es = read.csv("C:\\Users\\Lela Panag\\Desktop\\SIMULATION SCENARIOS FINAL\\Datasets application\\Simulated data application\\Normal_data\\SN\\res_Normal_data_Bayesian_SN_model_HN.csv")

sn_mod_rel_eff = read.csv("C:\\Users\\Lela Panag\\Desktop\\SIMULATION SCENARIOS FINAL\\Datasets application\\Simulated data application\\Normal_data\\SN\\pred_rel_effSN_Normal_data_Bayesian_SN_model_HN.csv")
sn_mod_es = read.csv("C:\\Users\\Lela Panag\\Desktop\\SIMULATION SCENARIOS FINAL\\Datasets application\\Simulated data application\\Normal_data\\SN\\pred_res_Normal_data_Bayesian_SN_model_HN.csv")

extra_data3 = sn_mod_rel_eff$median
lb_sn = sn_mod_rel_eff$q5
ub_sn = sn_mod_rel_eff$q95

p_true = read.csv("C:\\Users\\Lela Panag\\Desktop\\SIMULATION SCENARIOS FINAL\\Datasets application\\Simulated data application\\Normal_data\\true_eff\\true_eff_normal.csv")
p_true = round(p_true, 2)
p_true = as.vector(p_true)

add_data = cbind.data.frame(p_true) 

forest(x = extra_data3, 
       ci.lb = lb_sn, 
       ci.ub = ub_sn, 
       xlim = c(-7.2,11), 
       psize = 2,
       cex = 1.3,
       ilab = add_data, 
       ilab.xpos = 8, 
       lwd = 3,
       ylim =c(-1,29))

text(c(-6.7, 8, 9.9), 28, c("Studies" ,"true","Estimate[95% CI]"),   font=2, cex=1.5)
text(0, 29, c("Binomial-SN(HN)"),   font=2, cex=1.5)
addpoly(x= sn_mod_es$median, ci.lb =sn_mod_es$lowerCI,ci.ub = sn_mod_es$upperCI, rows=-1)

abline(h=0, lwd=0.1, col="black", lty=1)

# Open the TIFF device
tiff("forest_Normal_SN.tiff", units="in", width=15, height=9, res=300)

forest(x = extra_data3, 
       ci.lb = lb_sn, 
       ci.ub = ub_sn, 
       xlim = c(-7.2,11), 
       psize = 2,
       cex = 1.3,
       ilab = add_data, 
       ilab.xpos = 8, 
       lwd = 3,
       ylim =c(-1,29))

text(c(-6.6, 8, 9.7), 28, c("Studies" ,"true","Estimate[95% CI]"),   font=2, cex=1.5)
text(0, 29, c("Binomial-SN(HN)"),   font=2, cex=1.5)
addpoly(x= sn_mod_es$median, ci.lb =sn_mod_es$lowerCI,ci.ub = sn_mod_es$upperCI, rows=-1)

abline(h=0, lwd=0.1, col="black", lty=1)

# Close the TIFF device
dev.off()

############## THE OVERLAP PLOTS ###########
library(overlapping)

scN <- list(True_study_effects_posterior = p_true$x, Binomial_Normal_HN = extra_dataNHN)
scT <- list(True_study_effects_posterior = p_true$x,  Binomial_t_HN =extra_data2)
scSN <- list(True_study_effects_posterior = p_true$x,  Binomial_SN_HN=extra_data3)
scDP <- list(True_study_effects_posterior = p_true$x,  Binomial_DP_n = extra_dataDPUG)

outN <- overlap(scN,type="2",  plot=TRUE, boundaries  = c(-6,6))
outT <- overlap(scT,type="2", plot=TRUE, boundaries  = c(-6,6))
outSN <- overlap(scSN,type="2", plot=TRUE, boundaries  = c(-6,6))
outDP <- overlap(scDP,type="2", plot=TRUE, boundaries  = c(-6,6))


# Open the TIFF device
tiff("Overlap_normal_data_NHN.tiff", units="in", width=15, height=9, res=300)

outN <- overlap(scN,type="2",  plot=TRUE, boundaries  = c(-6,6))

# Close the TIFF device
dev.off()

########## SKEW NORMAL DATA ###########
############### BAYESIAN NORMAL MODEL ############
########### BINOMIAL-NORMAL(HN)
NormalHN_mod_rel_eff = read.csv("C:\\Users\\Lela Panag\\Desktop\\SIMULATION SCENARIOS FINAL\\Datasets application\\Simulated data application\\Skew_Normal_data\\NORMAL\\rel_eff_SN_data_Bayesian_NORMAL_HN_model.csv")
NormalHN_mod_es = read.csv("C:\\Users\\Lela Panag\\Desktop\\SIMULATION SCENARIOS FINAL\\Datasets application\\Simulated data application\\Skew_Normal_data\\NORMAL\\res_SN_data_Bayesian_NORMAL_HN_model.csv")

extra_dataNHN = NormalHN_mod_rel_eff$rel_eff
lb_NHN = NormalHN_mod_rel_eff$LB_rel_eff
ub_NHN = NormalHN_mod_rel_eff$UB_rel_eff

p_true = read.csv("C:\\Users\\Lela Panag\\Desktop\\SIMULATION SCENARIOS FINAL\\Datasets application\\Simulated data application\\Skew_Normal_data\\true_eff\\true_eff_skew_normal.csv")
p_true = round(p_true, 2)
p_true = as.vector(p_true)

add_data = cbind.data.frame(p_true) 

forest(x = extra_dataNHN, 
       ci.lb = lb_NHN, 
       ci.ub = ub_NHN, 
       xlim = c(-3,4), 
       psize = 2,
       cex = 1.3,
       ilab = add_data, 
       ilab.xpos = 2.7, 
       lwd = 3,
       ylim =c(-1,29))

text(c(-2.8, 2.7, 3.6), 28, c("Studies" ,"true","Estimate[95% CI]"),   font=2, cex=1.5)
text(0, 29, "Binomial-Normal(HN)",   font=2, cex=1.5)
addpoly(x= NormalHN_mod_es$mu_Normaljags.bin, ci.lb = NormalHN_mod_es$LBmu.jags.bin,ci.ub = NormalHN_mod_es$UBmu.jags.bin,rows=-1)
abline(h=0, lwd=0.1, col="black", lty=1)

# Open the TIFF device
tiff("forest_SNormal_NHN.tiff", units="in", width=15, height=9, res=300)

forest(x = extra_dataNHN, 
       ci.lb = lb_NHN, 
       ci.ub = ub_NHN, 
       xlim = c(-3,4), 
       psize = 2,
       cex = 1.3,
       ilab = add_data, 
       ilab.xpos = 2.7, 
       lwd = 3,
       ylim =c(-1,29))

text(c(-2.7, 2.7, 3.5), 28, c("Studies" ,"true","Estimate[95% CI]"),   font=2, cex=1.5)
text(0, 29, "Binomial-Normal(HN)",   font=2, cex=1.5)
addpoly(x= NormalHN_mod_es$mu_Normaljags.bin, ci.lb = NormalHN_mod_es$LBmu.jags.bin,ci.ub = NormalHN_mod_es$UBmu.jags.bin,rows=-1)
abline(h=0, lwd=0.1, col="black", lty=1)

# Close the TIFF device
dev.off()

############### FOREST PLOT OF THE DP MODEL #########
DPHNU_mod_rel_eff = read.csv("C:\\Users\\Lela Panag\\Desktop\\SIMULATION SCENARIOS FINAL\\Datasets application\\Simulated data application\\Skew_Normal_data\\DP\\rel_eff_HN_U26_SN_data_Bayesian_DP_model.csv")
DPHNU_mod_es = read.csv("C:\\Users\\Lela Panag\\Desktop\\SIMULATION SCENARIOS FINAL\\Datasets application\\Simulated data application\\Skew_Normal_data\\DP\\res_HN_U26_SN_data_Bayesian_DP_model.csv")

extra_dataDPHNU = DPHNU_mod_rel_eff$rel_effDP
sd_extra_dataDPHNU = DPHNU_mod_rel_eff$sd_rel_effDP
lb_DPHNU = DPHNU_mod_rel_eff$LB_rel_effDP
ub_DPHNU = DPHNU_mod_rel_eff$UB_rel_effDP

prob_DP = read.csv("C:\\Users\\Lela Panag\\Desktop\\SIMULATION SCENARIOS FINAL\\Datasets application\\Simulated data application\\Skew_Normal_data\\DP\\prob_HN_U26_res_SN_data_Bayesian_DP_model.csv")
cluster_mean_DP = read.csv("C:\\Users\\Lela Panag\\Desktop\\SIMULATION SCENARIOS FINAL\\Datasets application\\Simulated data application\\Skew_Normal_data\\DP\\cluster_mean_HN_U26_res_SN_data_Bayesian_DP_model.csv")

p_true = read.csv("C:\\Users\\Lela Panag\\Desktop\\SIMULATION SCENARIOS FINAL\\Datasets application\\Simulated data application\\Skew_Normal_data\\true_eff\\true_eff_skew_normal.csv")
p_true = round(p_true, 2)
p_true = as.vector(p_true)

# ####### clusters means ######
mean_cl1 = cluster_mean_DP$mean_cl1
mean_cl2 = cluster_mean_DP$mean_cl2
mean_cl3 = cluster_mean_DP$mean_cl3

clusters = c("1","1","1","1","1","1","1","1","2","1","1","1","1","2","1","2","1","1","1","1","1","3","1","1","1","1")
pr = prob_DP
### study 9 has 0.1590 prob of being in cluster 2 and 0.1560 of being in cluster 3 

add_data = cbind.data.frame(clusters, pr, p_true) 

forest(x = extra_dataDPHNU, 
       ci.lb = lb_DPHNU, 
       ci.ub = ub_DPHNU, 
       psize = 2,
       xlim = c(-3,4), 
       cex = 1.3,
       ilab = add_data, 
       ilab.xpos = c(2,2.5,3), 
       lwd = 3,
       ylim =c(-1,29),
       col = c("black","black","black", "black","black","black","black","black", "red", "black","black","black","black","red","black","red","black","black","black","black",
       "black","blue","black","black","black","black"),
       pch = c(15,15,15,15,15,15,15,15,16,15,15,15,15,16,15,16,15,15,15,15,15,14,15,15,15,15))


text(c(-2.8,2,2.55,3, 3.6), 28, c("Studies" ,"clusters", "probabilities", "true","Estimate[95% CI]"),   font=2, cex=1.5)
text(0, 29, "Binomial-DP-26(HN/Unif)",   font=2, cex=1.5)
addpoly(x= DPHNU_mod_es$mu_DP, ci.lb = DPHNU_mod_es$LB_mu_DP ,ci.ub = DPHNU_mod_es$UB_mu_DP ,rows=-1)
abline(h=0, lwd=0.1, col="black", lty=1)

line_positions <- c(mean_cl1, mean_cl2, mean_cl3)
line_colors <- c("black", "red", "blue")

for (i in seq_along(line_positions)) {
  segments(x0 = line_positions[i], y0 = -4, x1 = line_positions[i], y1 = length(extra_dataDPHNU) + 1, 
           lwd = 3, col = line_colors[i], lty = 4)
}

legend(x = -1.5, y= 29, legend = c("Mean cluster1","Mean cluster2","Mean cluster 3"),
       col = c("black","red", "blue"), lwd = c(3,3),
       lty = c(2,2),
       cex = 1.3 , y.intersp = 1, box.lwd = 1, box.lty = 1, box. = 10, xjust = 0)

# Open the TIFF device
tiff("forest_SNormal_DP.tiff", units="in", width=15, height=9, res=300)

forest(x = extra_dataDPHNU, 
       ci.lb = lb_DPHNU, 
       ci.ub = ub_DPHNU, 
       psize = 2,
       xlim = c(-3,4), 
       cex = 1.3,
       ilab = add_data, 
       ilab.xpos = c(1.8,2.3,2.8), 
       lwd = 3,
       ylim =c(-1,29),
       col = c("black","black","black", "black","black","black","black","black", "red", "black","black","black","black","red","black","red","black","black","black","black",
               "black","blue","black","black","black","black"),
       pch = c(15,15,15,15,15,15,15,15,16,15,15,15,15,16,15,16,15,15,15,15,15,14,15,15,15,15))


text(c(-2.7,1.7,2.3,2.8, 3.5), 28, c("Studies" ,"clusters", "probabilities", "true","Estimate[95% CI]"),   font=2, cex=1.3)
text(0, 29, "Binomial-DP-26(HN/Unif)",   font=2, cex=1.5)
addpoly(x= DPHNU_mod_es$mu_DP, ci.lb = DPHNU_mod_es$LB_mu_DP ,ci.ub = DPHNU_mod_es$UB_mu_DP ,rows=-1)
abline(h=0, lwd=0.1, col="black", lty=1)

line_positions <- c(mean_cl1, mean_cl2, mean_cl3)
line_colors <- c("black", "red", "blue")

for (i in seq_along(line_positions)) {
  segments(x0 = line_positions[i], y0 = -4, x1 = line_positions[i], y1 = length(extra_dataDPHNU) + 1, 
           lwd = 3, col = line_colors[i], lty = 4)
}
# Close the TIFF device
dev.off()
################ Binomial-t(Unif) MODEL FOREST PLOT ###########
# t_mod_rel_eff = read.csv("C:\\Users\\Lela Panag\\Desktop\\SIMULATION SCENARIOS FINAL\\Datasets application\\Simulated data application\\Skew_Normal_data\\T\\rel_effT_SN_data_Bayesian_T_model_U.csv")
# t_mod_es = read.csv("C:\\Users\\Lela Panag\\Desktop\\SIMULATION SCENARIOS FINAL\\Datasets application\\Simulated data application\\Skew_Normal_data\\T\\SN_data_Bayesian_T_model_U.csv")

t_mod_rel_eff = read.csv("C:\\Users\\Lela Panag\\Desktop\\SIMULATION SCENARIOS FINAL\\Datasets application\\Simulated data application\\Skew_Normal_data\\T\\pred_rel_effT_SN_data_Bayesian_T_model_U.csv")
t_mod_es = read.csv("C:\\Users\\Lela Panag\\Desktop\\SIMULATION SCENARIOS FINAL\\Datasets application\\Simulated data application\\Skew_Normal_data\\T\\pred_SN_data_Bayesian_T_model_U.csv")

extra_data2 = t_mod_rel_eff$median
lb_t = t_mod_rel_eff$q5
ub_t = t_mod_rel_eff$q95

p_true = read.csv("C:\\Users\\Lela Panag\\Desktop\\SIMULATION SCENARIOS FINAL\\Datasets application\\Simulated data application\\Skew_Normal_data\\true_eff\\true_eff_skew_normal.csv")
p_true = round(p_true, 2)
p_true = as.vector(p_true)

add_data = cbind.data.frame(p_true) 

forest(x = extra_data2, 
       ci.lb = lb_t, 
       ci.ub = ub_t, 
       xlim = c(-3,4), 
       psize = 2,
       cex = 1.3,
       ilab = add_data, 
       ilab.xpos = 2.7, 
       lwd = 3,
       ylim =c(-1,29))


text(c(-2.8, 2.7, 3.6), 28, c("Studies" ,"true","Estimate[95% CI]"),   font=2, cex=1.5)
text(0, 29, "Binomial-t(Unif)",   font=2, cex=1.5)
addpoly(x= t_mod_es$median , ci.lb =t_mod_es$lowerCI, ci.ub = t_mod_es$upperCI,rows=-1)
abline(h=0, lwd=0.1, col="black", lty=1)

# Open the TIFF device
tiff("forest_SNormal_t.tiff", units="in", width=15, height=9, res=300)

forest(x = extra_data2, 
       ci.lb = lb_t, 
       ci.ub = ub_t, 
       xlim = c(-3,4), 
       psize = 2,
       cex = 1.3,
       ilab = add_data, 
       ilab.xpos = 2.7, 
       lwd = 3,
       ylim =c(-1,29))


text(c(-2.7, 2.7, 3.5), 28, c("Studies" ,"true","Estimate[95% CI]"),   font=2, cex=1.5)
text(0, 29, "Binomial-t(Unif)",   font=2, cex=1.5)
addpoly(x= t_mod_es$median , ci.lb =t_mod_es$lowerCI, ci.ub = t_mod_es$upperCI,rows=-1)
abline(h=0, lwd=0.1, col="black", lty=1)

# Close the TIFF device
dev.off()
################# SKEW NORMAL MODEL ################
# sn_mod_rel_eff = read.csv("C:\\Users\\Lela Panag\\Desktop\\SIMULATION SCENARIOS FINAL\\Datasets application\\Simulated data application\\Skew_Normal_data\\SN\\rel_effSN_SN_data_Bayesian_SN_model_U.csv")
# sn_mod_es = read.csv("C:\\Users\\Lela Panag\\Desktop\\SIMULATION SCENARIOS FINAL\\Datasets application\\Simulated data application\\Skew_Normal_data\\SN\\SN_data_Bayesian_SN_model_U.csv")

sn_mod_rel_eff = read.csv("C:\\Users\\Lela Panag\\Desktop\\SIMULATION SCENARIOS FINAL\\Datasets application\\Simulated data application\\Skew_Normal_data\\SN\\pred_rel_effSN_SN_data_Bayesian_SN_model_U.csv")
sn_mod_es = read.csv("C:\\Users\\Lela Panag\\Desktop\\SIMULATION SCENARIOS FINAL\\Datasets application\\Simulated data application\\Skew_Normal_data\\SN\\pred_SN_data_Bayesian_SN_model_U.csv")

extra_data3 = sn_mod_rel_eff$median
lb_sn = sn_mod_rel_eff$q5
ub_sn = sn_mod_rel_eff$q95

p_true = read.csv("C:\\Users\\Lela Panag\\Desktop\\SIMULATION SCENARIOS FINAL\\Datasets application\\Simulated data application\\Skew_Normal_data\\true_eff\\true_eff_skew_normal.csv")
p_true = round(p_true, 2)
p_true = as.vector(p_true)

add_data = cbind.data.frame(p_true) 

forest(x = extra_data3, 
       ci.lb = lb_sn, 
       ci.ub = ub_sn, 
       xlim = c(-3,4),
       psize = 2,
       cex = 1.3,
       ilab = add_data, 
       ilab.xpos = 2.7, 
       lwd = 3,
       ylim =c(-1,29))


text(c(-2.8, 2.7, 3.6), 28, c("Studies" ,"true","Estimate[95% CI]"),   font=2, cex=1.5)
text(0, 29, "Binomial-SN(Unif)",   font=2, cex=1.5)
addpoly(x= sn_mod_es$median, ci.lb =sn_mod_es$lowerCI,ci.ub = sn_mod_es$upperCI, rows=-1)
abline(h=0, lwd=0.1, col="black", lty=1)

# Open the TIFF device
tiff("forest_SNormal_SN.tiff", units="in", width=15, height=9, res=300)

forest(x = extra_data3, 
       ci.lb = lb_sn, 
       ci.ub = ub_sn, 
       xlim = c(-3,4),
       psize = 2,
       cex = 1.3,
       ilab = add_data, 
       ilab.xpos = 2.7, 
       lwd = 3,
       ylim =c(-1,29))


text(c(-2.7, 2.7, 3.5), 28, c("Studies" ,"true","Estimate[95% CI]"),   font=2, cex=1.5)
text(0, 29, "Binomial-SN(Unif)",   font=2, cex=1.5)
addpoly(x= sn_mod_es$median, ci.lb =sn_mod_es$lowerCI,ci.ub = sn_mod_es$upperCI, rows=-1)
abline(h=0, lwd=0.1, col="black", lty=1)

# Close the TIFF device
dev.off()
############## THE OVERLAP PLOTS ###########
library(overlapping)

scN <- list(True = p_true$x, Binomial_Normal= extra_dataNHN)
scT <- list(True = p_true$x,  Binomial_t_U =extra_data2)
scSN <- list(True = p_true$x,  Binomial_SN_U=extra_data3)
scDP <- list(True = p_true$x,  Binomial_DP12 = extra_dataDPHNU)

outN <- overlap(scN,type="2",  plot=TRUE, boundaries  = c(-6,6))
outT <- overlap(scT,type="2", plot=TRUE, boundaries  = c(-6,6))
outSN <- overlap(scSN,type="2", plot=TRUE, boundaries  = c(-6,6))
outDP <- overlap(scDP,type="2", plot=TRUE, boundaries  = c(-6,6))


################ MIXTURE DATA FOREST PLOTS #############
######### NORMAL(HN) FOREST PLOT ##############

NormalHN_mod_rel_eff = read.csv("C:\\Users\\Lela Panag\\Desktop\\SIMULATION SCENARIOS FINAL\\Datasets application\\Simulated data application\\Mix_data\\NORMAL\\HN_rel_effN_Mix_data_Bayesian_Normal_model.csv")
NormalHN_mod_es = read.csv("C:\\Users\\Lela Panag\\Desktop\\SIMULATION SCENARIOS FINAL\\Datasets application\\Simulated data application\\Mix_data\\NORMAL\\HN_Mix_data_Bayesian_Normal_model.csv")

extra_dataNHN = NormalHN_mod_rel_eff$rel_eff
lb_NHN = NormalHN_mod_rel_eff$LB_rel_eff
ub_NHN = NormalHN_mod_rel_eff$UB_rel_eff

p_true = read.csv("C:\\Users\\Lela Panag\\Desktop\\SIMULATION SCENARIOS FINAL\\Datasets application\\Simulated data application\\Mix_data\\true_eff\\true_eff_mix.csv")
p_true = round(p_true, 2)
p_true = as.vector(p_true)

add_data = cbind.data.frame(p_true) 

forest(x = extra_dataNHN, 
       ci.lb = lb_NHN, 
       ci.ub = ub_NHN, 
       xlim = c(-2,4), 
       psize = 2,
       cex = 1.3,
       ilab = add_data, 
       ilab.xpos = 3, 
       lwd = 3,
       ylim =c(-1,29))

text(c(-1.8, 3, 3.6), 28, c("Studies" ,"true","Estimate[95% CI]"),   font=2, cex=1.5)
text(0, 29, "Binomial-Normal(HN)",   font=2, cex=1.5)
addpoly(x= NormalHN_mod_es$mu_Normaljags.bin, ci.lb = NormalHN_mod_es$LBmu.jags.bin,ci.ub = NormalHN_mod_es$UBmu.jags.bin, rows=-1)
abline(h=0, lwd=0.1, col="black", lty=1)

# Open the TIFF device
tiff("forest_Mix_NHN.tiff", units="in", width=15, height=9, res=300)

forest(x = extra_dataNHN, 
       ci.lb = lb_NHN, 
       ci.ub = ub_NHN, 
       xlim = c(-2,4), 
       psize = 2,
       cex = 1.3,
       ilab = add_data, 
       ilab.xpos = 3, 
       lwd = 3,
       ylim =c(-1,29))

text(c(-1.8, 3, 3.55), 28, c("Studies" ,"true","Estimate[95% CI]"),   font=2, cex=1.5)
text(0, 29, "Binomial-Normal(HN)",   font=2, cex=1.5)
addpoly(x= NormalHN_mod_es$mu_Normaljags.bin, ci.lb = NormalHN_mod_es$LBmu.jags.bin,ci.ub = NormalHN_mod_es$UBmu.jags.bin, rows=-1)
abline(h=0, lwd=0.1, col="black", lty=1)

# Close the TIFF device
dev.off()
################# 

############## DP MODEL FOREST PLOT ############
DPHNU51_mod_rel_eff = read.csv("C:\\Users\\Lela Panag\\Desktop\\SIMULATION SCENARIOS FINAL\\Datasets application\\Simulated data application\\Mix_data\\DP\\HNU51_rel_effDP_Mix_data_Bayesian_DP_model.csv")
DPHNU51_mod_es = read.csv("C:\\Users\\Lela Panag\\Desktop\\SIMULATION SCENARIOS FINAL\\Datasets application\\Simulated data application\\Mix_data\\DP\\HNU51_Mix_data_Bayesian_DP_model.csv")

extra_dataDPHNU51 = DPHNU51_mod_rel_eff$rel_effDP
sd_extra_dataDPHNU51 = DPHNU51_mod_rel_eff$sd_rel_effDP
lb_DPHNU51 = DPHNU51_mod_rel_eff$LB_rel_effDP
ub_DPHNU51 = DPHNU51_mod_rel_eff$UB_rel_effDP
# 
cluster_mean_DP = read.csv("C:\\Users\\Lela Panag\\Desktop\\SIMULATION SCENARIOS FINAL\\Datasets application\\Simulated data application\\Mix_data\\DP\\HNU51_cluster_mean_Mix_data_Bayesian_DP_model.csv")
# ####### clusters means ######
mean_cl1 = cluster_mean_DP$mean_cl1
mean_cl2 = cluster_mean_DP$mean_cl2
mean_cl3 = cluster_mean_DP$mean_cl3

prob = read.csv("C:\\Users\\Lela Panag\\Desktop\\SIMULATION SCENARIOS FINAL\\Datasets application\\Simulated data application\\Mix_data\\DP\\prob_HNU51_Mix_data_Bayesian_DP_model.csv")

clusters = c("3","2","3","3","3","1/2","3","3","3","3","3","2","3","2","3","3/2","3","3","3","3","3","3","2","2","3","3")


pr = prob

p_true = read.csv("C:\\Users\\Lela Panag\\Desktop\\SIMULATION SCENARIOS FINAL\\Datasets application\\Simulated data application\\Mix_data\\true_eff\\true_eff_mix.csv")
p_true = round(p_true, 2)
p_true = as.vector(p_true)

add_data = cbind.data.frame(clusters, pr, p_true) 

forest(x = extra_dataDPHNU51, 
       ci.lb = lb_DPHNU51, 
       ci.ub = ub_DPHNU51, 
       xlim = c(-2,4), 
       psize = 2,
       cex = 1.3,
       ilab = add_data, 
       ilab.xpos = c(2,2.5,3), 
       lwd = 3,
       ylim =c(-1,29),
       col = c("black","red","black","black","black","blue","black","black","black", "black","black","red","black",
       "red","black","black","black","black","black","black","black","black","red","red","black","black"),
       pch = c(15,16,15,15,15,14,15,15,15,15,15,16, 15,16,15,15,15,15,15,15,15,15,16,16,15,15))
 

text(c(-1.8,2,2.5,3, 3.6), 28, c("Studies" ,"clusters", "probabilities", "true","Estimate[95% CI]"),   font=2, cex=1.5)
text(0, 29, "Binomial-DP-51(HN/Unif)",   font=2, cex=1.5)
addpoly(x= DPHNU51_mod_es$mu_DP, ci.lb = DPHNU51_mod_es$LB_mu_DP ,ci.ub = DPHNU51_mod_es$UB_mu_DP , rows=-1)
abline(h=0, lwd=0.1, col="black", lty=1)

line_positions <- c(mean_cl3, mean_cl2, mean_cl1)
line_colors <- c("red", "blue", "black")

for (i in seq_along(line_positions)) {
  segments(x0 = line_positions[i], y0 = -4, x1 = line_positions[i], y1 = length(extra_dataDPHNU51) + 1, 
           lwd = 3, col = line_colors[i], lty = 2)
}

legend(x = -1.5, y= 29, legend = c("Mean cluster1","Mean cluster2", "Mean cluster3"),
       col = c("blue","red","black"), lwd = c(3,3,3),
       lty = c(2,2,2),
       cex = 1.3 , y.intersp = 1, box.lwd = 1, box.lty = 1, box. = 10, xjust = 0)


# Open the TIFF device
tiff("forest_Mix_DP.tiff", units="in", width=15, height=9, res=300)

forest(x = extra_dataDPHNU51, 
       ci.lb = lb_DPHNU51, 
       ci.ub = ub_DPHNU51, 
       xlim = c(-2,4), 
       psize = 2,
       cex = 1.3,
       ilab = add_data, 
       ilab.xpos = c(2,2.5,3), 
       lwd = 3,
       ylim =c(-1,29),
       col = c("black","red","black","black","black","blue","black","black","black", "black","black","red","black",
               "red","black","black","black","black","black","black","black","black","red","red","black","black"),
       pch = c(15,16,15,15,15,14,15,15,15,15,15,16, 15,16,15,15,15,15,15,15,15,15,16,16,15,15))


text(c(-1.8,2,2.5,3, 3.55), 28, c("Studies" ,"clusters", "probabilities", "true","Estimate[95% CI]"),   font=2, cex=1.3)
text(0, 29, "Binomial-DP-51(HN/Unif)",   font=2, cex=1.5)
addpoly(x= DPHNU51_mod_es$mu_DP, ci.lb = DPHNU51_mod_es$LB_mu_DP ,ci.ub = DPHNU51_mod_es$UB_mu_DP , rows=-1)
abline(h=0, lwd=0.1, col="black", lty=1)

line_positions <- c(mean_cl3, mean_cl2, mean_cl1)
line_colors <- c("red", "blue", "black")

for (i in seq_along(line_positions)) {
  segments(x0 = line_positions[i], y0 = -4, x1 = line_positions[i], y1 = length(extra_dataDPHNU51) + 1, 
           lwd = 3, col = line_colors[i], lty = 2)
}
# Close the TIFF device
dev.off()

################## T(HN) MODEL FOREST PLOT ##############
# 
# t_mod_rel_eff = read.csv("C:\\Users\\Lela Panag\\Desktop\\SIMULATION SCENARIOS FINAL\\Datasets application\\Simulated data application\\Mix_data\\T\\rel_effT_Mix_data_Bayesian_T_model_HN.csv")
# t_mod_es = read.csv("C:\\Users\\Lela Panag\\Desktop\\SIMULATION SCENARIOS FINAL\\Datasets application\\Simulated data application\\Mix_data\\T\\Mix_data_Bayesian_T_model_HN.csv")

t_mod_rel_eff = read.csv("C:\\Users\\Lela Panag\\Desktop\\SIMULATION SCENARIOS FINAL\\Datasets application\\Simulated data application\\Mix_data\\T\\pred_rel_effT_Mix_data_Bayesian_T_model_HN.csv")
t_mod_es = read.csv("C:\\Users\\Lela Panag\\Desktop\\SIMULATION SCENARIOS FINAL\\Datasets application\\Simulated data application\\Mix_data\\T\\pred_Mix_data_Bayesian_T_model_HN.csv")


extra_data2 = t_mod_rel_eff$median
lb_t = t_mod_rel_eff$q5
ub_t = t_mod_rel_eff$q95

p_true = read.csv("C:\\Users\\Lela Panag\\Desktop\\SIMULATION SCENARIOS FINAL\\Datasets application\\Simulated data application\\Mix_data\\true_eff\\true_eff_mix.csv")
p_true = round(p_true, 2)
p_true = as.vector(p_true)

add_data = cbind.data.frame(p_true) 

forest(x = extra_data2, 
       ci.lb = lb_t, 
       ci.ub = ub_t, 
       xlim = c(-2,4), 
       psize = 2,
       cex = 1.3,
       ilab = add_data, 
       ilab.xpos = 3, 
       lwd = 3,
       ylim =c(-1,29))

text(c(-1.8, 3, 3.6), 28, c("Studies" ,"true","Estimate[95% CI]"),   font=2, cex=1.5)
text(0, 29, "Binomial-t(HN)",   font=2, cex=1.5)
addpoly(x= t_mod_es$median , ci.lb =t_mod_es$lowerCI, ci.ub = t_mod_es$upperCI, rows=-1)
abline(h=0, lwd=0.1, col="black", lty=1)

# Open the TIFF device
tiff("forest_Mix_t.tiff", units="in", width=15, height=9, res=300)

forest(x = extra_data2, 
       ci.lb = lb_t, 
       ci.ub = ub_t, 
       xlim = c(-2,4), 
       psize = 2,
       cex = 1.3,
       ilab = add_data, 
       ilab.xpos = 3, 
       lwd = 3,
       ylim =c(-1,29))

text(c(-1.8, 3, 3.55), 28, c("Studies" ,"true","Estimate[95% CI]"),   font=2, cex=1.5)
text(0, 29, "Binomial-t(HN)",   font=2, cex=1.5)
addpoly(x= t_mod_es$median , ci.lb =t_mod_es$lowerCI, ci.ub = t_mod_es$upperCI, rows=-1)
abline(h=0, lwd=0.1, col="black", lty=1)

# Close the TIFF device
dev.off()
################# SN(HN) MODEL FORESTPLOT ############
# sn_mod_rel_eff = read.csv("C:\\Users\\Lela Panag\\Desktop\\SIMULATION SCENARIOS FINAL\\Datasets application\\Simulated data application\\Mix_data\\SN\\rel_effSN_Mix_data_Bayesian_SN_model_HN.csv")
# sn_mod_es = read.csv("C:\\Users\\Lela Panag\\Desktop\\SIMULATION SCENARIOS FINAL\\Datasets application\\Simulated data application\\Mix_data\\SN\\Mix_data_Bayesian_SN_model_HN.csv")

sn_mod_rel_eff = read.csv("C:\\Users\\Lela Panag\\Desktop\\SIMULATION SCENARIOS FINAL\\Datasets application\\Simulated data application\\Mix_data\\SN\\pred_rel_effSN_Mix_data_Bayesian_SN_model_HN.csv")
sn_mod_es = read.csv("C:\\Users\\Lela Panag\\Desktop\\SIMULATION SCENARIOS FINAL\\Datasets application\\Simulated data application\\Mix_data\\SN\\pred_Mix_data_Bayesian_SN_model_HN.csv")

extra_data3 = sn_mod_rel_eff$median
lb_sn = sn_mod_rel_eff$q5
ub_sn = sn_mod_rel_eff$q95

p_true = read.csv("C:\\Users\\Lela Panag\\Desktop\\SIMULATION SCENARIOS FINAL\\Datasets application\\Simulated data application\\Mix_data\\true_eff\\true_eff_mix.csv")
p_true = round(p_true, 2)
p_true = as.vector(p_true)

add_data = cbind.data.frame(p_true) 

forest(x = extra_data3, 
       ci.lb = lb_sn, 
       ci.ub = ub_sn, 
       xlim = c(-2,4), 
       psize = 2,
       cex = 1.3,
       ilab = add_data, 
       ilab.xpos = 3, 
       lwd = 3,
       ylim =c(-1,29))

text(c(-1.8, 3, 3.6), 28, c("Studies" ,"true","Estimate[95% CI]"),   font=2, cex=1.5)
text(0, 29, "Binomial-SN(HN)",   font=2, cex=1.5)
addpoly(x= sn_mod_es$median, ci.lb =sn_mod_es$lowerCI,ci.ub = sn_mod_es$upperCI, rows=-1)
abline(h=0, lwd=0.1, col="black", lty=1)

# Open the TIFF device
tiff("forest_Mix_SN.tiff", units="in", width=15, height=9, res=300)

forest(x = extra_data3, 
       ci.lb = lb_sn, 
       ci.ub = ub_sn, 
       xlim = c(-2,4), 
       psize = 2,
       cex = 1.3,
       ilab = add_data, 
       ilab.xpos = 3, 
       lwd = 3,
       ylim =c(-1,29))

text(c(-1.8, 3, 3.55), 28, c("Studies" ,"true","Estimate[95% CI]"),   font=2, cex=1.5)
text(0, 29, "Binomial-SN(HN)",   font=2, cex=1.5)
addpoly(x= sn_mod_es$median, ci.lb =sn_mod_es$lowerCI,ci.ub = sn_mod_es$upperCI, rows=-1)
abline(h=0, lwd=0.1, col="black", lty=1)

# Close the TIFF device
dev.off()
############# THE OVERLAP PLOTS ###########
library(overlapping)

scN <- list(True = p_true$x, Binomial_Normal= extra_dataNHN)
scT <- list(True = p_true$x,  Binomial_t =extra_data2)
scSN <- list(True = p_true$x,  Binomial_SN=extra_data3)
scDP <- list(True = p_true$x,  Binomial_DP51 = extra_dataDPHNU51)

outN <- overlap(scN,type="2",  plot=TRUE, boundaries  = c(-6,6))
outT <- overlap(scT,type="2", plot=TRUE, boundaries  = c(-6,6))
outSN <- overlap(scSN,type="2", plot=TRUE, boundaries  = c(-6,6))
outDP <- overlap(scDP,type="2", plot=TRUE, boundaries  = c(-6,6))



