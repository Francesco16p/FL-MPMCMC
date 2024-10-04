# This code reproduces the simulation results reported in Section 4 of the paper
# on the application of different multiple proposals MCMC methods to samples
# from a high-dimensional logistic regression.

# Preparation of the working environment. 
# Note that the folder containing the code MUST BE SET as the working directory
# before proceeding.

rm(list = ls())
source("Log_densities.R") # File containing the target density
source("Log_derivatives.R") # File containing the first derivative of target density
source("MALA algorithm.R") # Single-proposal MALA algorithm
source("MTM algorithm.R") # Multiple-try Metropolis with random walk proposals
source("MTM MALA algorithm.R") #  Multiple-try Metropolis with MALA proposals
source("Generalised Metropolis.R") # Generalised Metropolis star-shaped proposal
source("Summary stats.R") # File containing functions to compute the ESJD

### Model specification and data simulation ###

set.seed(1) # set seed for reproducibility
n <- 50 # sample size
d <- 50 # number of explanatory variables
hyper_param <- sqrt(25/d) # hyper-parameter variance prior Bayesian model

# Regression setting
theta_0 <- rnorm(d,sd = 1/4 ) # True parameter
z <- cbind( rep(1,n) , matrix(rnorm(n*(d-1)), nrow = n )) # Explanatory variables
eta_0 <- z%*%theta_0 # linear predictors under the true parameter
pr0 <- 1/(1+exp(-eta_0)) #  probabilities under the true parameter
y <- rbinom(n,1,pr0) # simulated response variable

# Save the explanatory variables, the simulated response vector and 
# the variance of the prior distribution for future use
pars <- list(z = z, y = y, hyper_param = hyper_param)

### Sampling using different multiple-proposals algorithms ###

# We know draw a sample of dimension 50 000 for each one of the 
# MP-MCMC considered in Section 4 of the main article. For each simulation
# method the results are saved in a list.

nsim <-  5*10^4
# alpha_s <- 0.3 
# Random walk methods
sims_mtm <- list(); sims_mtm_anti <-list(); 
sims_mtm_lb <- list();  sims_mtm_anti_lb <- list();
sims_gm <- list()
# MALA methods
sims_mtm_mala <- list(); sims_mtm_lb_mala <-list(); 

# Vector containing different values of K at which the algorithms are evaluated
K_vec<-round(seq(from = 2, to =d, length.out = 8)) 
j <- 1
for(K in K_vec){
   sims_mtm[[j]] <- MTM(nsim=nsim, K=K,log_den = log_post_logistic, pars = pars,  prop_method = "ind", weight_fun = "identity",
                  sd_0 = 0.1, alpha_s = 0.3, d = d)
   sims_mtm_anti[[j]] <- MTM(nsim=nsim, K=K,log_den = log_post_logistic, pars = pars,  prop_method = "anti", weight_fun = "identity",
                  sd_0 = 0.1, alpha_s = 0.3, d = d)
   sims_mtm_lb[[j]] <- MTM(nsim=nsim, K=K,log_den = log_post_logistic, pars = pars,  prop_method = "ind", weight_fun = "sqrt",
                  sd_0 = 0.1, alpha_s = 0.3, d = d)
   sims_mtm_anti_lb[[j]] <- MTM(nsim=nsim, K=K,log_den = log_post_logistic, pars = pars,  prop_method = "anti", weight_fun = "sqrt",
                       sd_0 = 0.1, alpha_s = 0.3, d = d)
   sims_gm[[j]] <- GMH(nsim = nsim, K = K, log_den = log_post_logistic,pars = pars,sd_0 = 0.1,alpha_s = 0.3,d = d)
   sims_mtm_mala[[j]] <- MTM_MALA(nsim=nsim, K=K,log_den = log_post_logistic, dlog_den = dlog_logistic ,pars = pars, weight_fun = "identity",
                             sd_0 = 0.1, alpha_s = 0.4, d = d)
   sims_mtm_lb_mala[[j]] <- MTM_MALA(nsim=nsim, K=K,log_den = log_post_logistic,dlog_den = dlog_logistic, pars = pars, 
                                weight_fun = "sqrt", sd_0 = 0.01, alpha_s = 0.4, d = d)
   j = j+1
   print(K)  
 }
 
# Simulation K = 1 
sim_mh <- GMH(nsim = nsim, K = 1, log_den = log_post_logistic,pars = pars,sd_0 = 0.1, alpha_s = 0.3, d = d) #Standard MH
sim_mala <- MALA(nsim = nsim, log_den = log_post_logistic, dlog_den = dlog_logistic,
               pars = pars,sd_0 = 0.1,alpha_s = 0.4, d = d) # MALA with K=1

# Save the results
save(sim_mh,sims_mtm,sims_mtm_anti,sims_mtm_lb,sims_mtm_anti_lb,sims_gm,K_vec, file = "MP_logistic_d50.RData" )
save(sim_mala,sims_mtm_mala,sims_mtm_lb_mala,K_vec, file = "MP_logistic_MALA_d50.RData" )

load("MP_logistic_d50.RData")
load("MP_logistic_MALA_d50.RData")

# EJSD K = 1 
esjd_k1 <- mean(summary_stats(sim_mh$theta_sim)$ESJD_sum) #Random walk
esjd_k1_mala <- mean(summary_stats(sim_mala$theta_sim)$ESJD_sum) # MALA


# Initialise vectors for ESJD
esjd_vecMTManti <- esjd_vecMTM <- esjd_vecMTManti_lb <- esjd_vecMTM_lb <- esjd_vecGM <- c(esjd_k1,rep(NA,length(K_vec)))
esjd_vecMTM_mala <- esjd_vecMTM_lb_mala <- c(esjd_k1_mala,rep(NA,length(K_vec)))
# Summary statistics

for(K in K_vec){
  sim_mtm <- sims_mtm[[which(K==K_vec)]]
  sim_mtm_anti <- sims_mtm_anti[[which(K==K_vec)]]
  sim_mtm_lb <- sims_mtm_lb[[which(K==K_vec)]]
  sim_mtm_anti_lb <- sims_mtm_anti_lb[[which(K==K_vec)]]
  sim_gm <- sims_gm[[which(K==K_vec)]]
  sim_mtm_mala <- sims_mtm_mala[[which(K==K_vec)]]
  sim_mtm_lb_mala <- sims_mtm_lb_mala[[which(K==K_vec)]]
  
  # ESJD
  esjd_vecMTM[which(K==K_vec)+1]<-mean(summary_stats(sim_mtm$theta_sim)$ESJD_sum)
  esjd_vecMTManti[which(K==K_vec)+1]<-mean(summary_stats(sim_mtm_anti$theta_sim)$ESJD_sum)
  esjd_vecMTM_lb[which(K==K_vec)+1]<-mean(summary_stats(sim_mtm_lb$theta_sim)$ESJD_sum)
  esjd_vecMTManti_lb[which(K==K_vec)+1]<-mean(summary_stats(sim_mtm_anti_lb$theta_sim)$ESJD_sum)
  esjd_vecGM[which(K==K_vec)+1]<-mean(summary_stats(sim_gm$theta_sim)$ESJD_sum)
  esjd_vecMTM_mala[which(K==K_vec)+1]<-mean(summary_stats(sim_mtm_mala$theta_sim)$ESJD_sum)
  esjd_vecMTM_lb_mala[which(K==K_vec)+1]<-mean(summary_stats(sim_mtm_lb_mala$theta_sim)$ESJD_sum)
}

# Save ESJD
save(esjd_vecMTM,esjd_vecMTM_lb,esjd_vecGM,esjd_vecMTManti,
     esjd_vecMTManti_lb, file = "data_Gplot_logistic.RData")
save(esjd_vecMTM_mala,esjd_vecMTM_lb_mala, file = "data_Gplot_logisticMALA.RData")


# Load ESJD
load("data_Gplot_logistic.RData")
load("data_Gplot_logisticMALA.RData")

# Create two data-frame to plot the results using ggplot.

# Data-frame for random walk results
df_esjd_rw <- data.frame(ESJD = c(esjd_vecMTM,esjd_vecMTManti,esjd_vecMTM_lb,
                esjd_vecMTManti_lb,esjd_vecGM) , K = rep(c(1,K_vec),5), Method = c( rep("MTM",9),
                rep("antiMTM",9), rep("LB-MTM",9), rep("LB-antiMTM",9),rep("GM",9)))
df_esjd_rw$Method <- as.factor(df_esjd_rw$Method)

# Data frame for MALA results
df_esjd_mala <- data.frame(ESJD = c(esjd_vecMTM_mala,esjd_vecMTM_lb_mala),
                           K = rep(c(1,K_vec),2), Method = c( rep("MTM MALA",9),
                                                              rep("LB-MTM MALA",9)))
df_esjd_mala$Method <- as.factor(df_esjd_mala$Method)

# Save the two dataframes for future use
save(df_esjd_rw, file = "df_esjd_rw.RData")
save(df_esjd_mala, file = "df_esjd_mala.RData")

### Final plots ###
library(ggplot2)
library(patchwork)
# Load simulations results
load("df_esjd_rw.RData")
load("df_esjd_mala.RData")

plot_rw <- ggplot(df_esjd_rw, aes(x=K, y=ESJD, group=Method)) +
  geom_line(aes(linetype=Method, color = Method), size =1.5)+
  geom_point(aes(color = Method, shape=Method), size = 5) +
  scale_shape_manual(values = c(21, 22, 23, 24, 25))+
  theme_light()+
  scale_y_continuous(breaks = c(0.010,0.024))+
  scale_x_continuous(limits = c(1, 50),breaks = c(1,10,20,30,40,50))+
  geom_function(fun = function(x) 0.007 * log(x+1), color = "black") +
  annotate("text", x = 30, y =  0.0024 + 0.007 * log(30), label = "0.007*log(K+1)", color = "black", vjust = -1, size = 7)+
  labs(y = "Average ESJD") + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = c(0.786, 0.193),
        legend.background = element_rect(color = "gray"),
        text = element_text(size = 26),  
        axis.title = element_text(size = 26), 
        axis.text = element_text(size = 22),  
        legend.text = element_text(size = 22), 
        legend.title = element_text(size = 24))+
  ggtitle("A)")

plot_mala <- ggplot(df_esjd_mala, aes(x=K, y=ESJD, group=Method)) +
  geom_line(aes(linetype=Method, color = Method), size =1.5)+
  geom_point(aes(color = Method, shape=Method), size = 5) +
  scale_shape_manual(values = c(21, 22, 23, 24, 25))+
  theme_light()+
  scale_y_continuous(breaks = c(0.055,0.083))+
  scale_x_continuous(limits = c(1, 50),breaks = c(1,10,20,30,40,50))+
  geom_function(fun = function(x) 0.05 + 0.01 * log(x), color = "black") +
  labs(y = "")+
  annotate("text", x = 30, y =  0.051 + 0.009 * log(51), label = "0.05 + 0.01*log(K)",
           color = "black", vjust = -1,size =7)+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = c(0.747, 0.112),
        legend.background = element_rect(color = "gray"),
        text = element_text(size = 26),  
        axis.title = element_text(size = 26), 
        axis.text = element_text(size = 22),  
        legend.text = element_text(size = 22), 
        legend.title = element_text(size = 24))+
  ggtitle("B)")

#X11()
plot_rw+plot_mala
