#### Header ####

# Figure showing relationship between theta (proportion of infetions before symptom onset) and R_s and R_q
# X axis is theta, ranging from 0 to 0.2 (or whatever the max is). Y axis is R_effective under Q or SM

# Run SMC
# Run Figure_Fraser to calculate theta
# Up-weight the high-theta values to have more data points out there
# Run forward simulations under Q and SM using the same input person so we hold theta constant.

#### Load data ####
desired_root <- "20200202_2019nCoV" # Paste the desired root here "YYYYMMDD_DISEASE"

load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "/", desired_root, "_HR.Theta.RData", sep=""))

load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "/", desired_root, "_LR.Theta.RData", sep=""))

#### Inspect original data ####
nrow(data[data$theta > 0,]) /
nrow(data)

summary(data[data$theta > 0,"theta"])

#### Plot R as a function of T_Lat_Offset ####

data.hr.melt <- melt(data.hr[,c("theta", "R_s", "R_q",  "T_lat_offset")], id.vars = c("theta", "T_lat_offset"))

# Plot 
plot_Tlat_hr <- ggplot(data.hr.melt, aes(x = T_lat_offset, y = value, color = variable)) +
  geom_point() +
  stat_smooth() +
  scale_color_manual(values = c("R_0_input" = "grey", "R_s" = "yellow", "R_q" = "blue"), name = "Intervention", labels = c("Active Monitoring" , "Quarantine", "None (Input)")) +
  theme_bw() +
  geom_hline(yintercept = 1, color = "grey") +
  xlab("Proportion of Infections\nbefore Symptom Onset") +
  ylab("Reproductive Number")
plot_Tlat_hr

pdf(file=paste(root, "_Plot_TLatOffset_HR.pdf", sep=""), width = 6, height = 4)
plot(plot_Tlat_hr)
dev.off()

#### High Resource Setting ####

# Interventions
background_intervention = "u"

prob_CT <- 0.9

gamma <- 0.9

gamma_q <- 0.9

parms_epsilon = list("uniform", 0, 1, 999, "independent", "independent")
names(parms_epsilon) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_CT_delay = list("uniform", 0, 1, 999, "independent", "independent")
names(parms_CT_delay) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

parms_R_0 = list("uniform", 2.2, 2.2, 999, "independent", "independent", 1) # mean = 2.249. 95% CI [1.46, 3.3]. Inspired by Li NEJM and Riou 2020 EuroSurveillance. May want to run sensitivity analysis with a wider CI, but can't get a big enough right skew to match Li or Riou's upper bound.
names(parms_R_0) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

# Settings
n_pop = 1000
num_generations <- 5
times <- 2000
names <- c("R_0", "R_hsb", "R_s", "R_q", "Abs_Benefit","Rel_Benefit","NNQ","obs_to_iso_q","Abs_Benefit_per_Qday", "ks", "T_lat_offset", "d_inf", "pi_t_triangle_center")
data.hr.theta <- data.frame(matrix(rep(NA, length(names)*times), nrow=times))
names(data.hr.theta) <- names

# sample from joint posterior distribution
data <- function_calculate_theta(data)
sample <- sample(x = row.names(data), size = times, replace = TRUE, prob = (data$theta + 0.001*max(data.hr$theta))) # Updated to up-sample the high theta values
hist(data[sample, "theta"])
params.set <- cbind(
  T_lat_offset = data[sample, "T_lat_offset"] * runif(n = times, min = 0.25, max = 4), # permute between 25% to 400%
  d_inf = data[sample, "d_inf"] * runif(n = times, min = 0.25, max = 1), # permute between 25% to 100%
  pi_t_triangle_center = data[sample, "pi_t_triangle_center"] * runif(n = times, min = 0, max = 1), # permute between 0% to 100%
  dispersion = dispersion,
  R_0 = Draw_Dist_fcn(sample, parms_R_0$dist, parms_R_0$parm1, parms_R_0$parm2))

for (i in 1:times){
  cat(".")
  if (i%%10 == 0){cat("|")}
  if (i%%100 == 0){cat("\n")}
  
  data.hr.theta[i, "T_lat_offset"] <- parms_T_lat$anchor_value <- params.set[i,"T_lat_offset"]
  data.hr.theta[i, "d_inf"] <- parms_d_inf$parm2 <- max(1, params.set[i,"d_inf"]) # Make sure it's at least 1
  data.hr.theta[i, "pi_t_triangle_center"] <- parms_pi_t$triangle_center <- params.set[i,"pi_t_triangle_center"]
  
  data.hr.theta[i, "R_0_input"] <- params.set[i,"R_0"]
  parms_R_0_input <- list("uniform", as.numeric(params.set[i,"R_0"]), as.numeric(params.set[i,"R_0"]), 999, "independent", "independent", 1)
  names(parms_R_0_input) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")
  
  dispersion <- params.set[i, "dispersion"]
  
  for (subseq_interventions in c("s","q")){      
    if (subseq_interventions == background_intervention & data.hr.theta[i, "R_0_input"] > 1){
      n_pop_input <- 200
    } else if (subseq_interventions == "hsb" & data.hr.theta[i, "R_0_input"] * (1-gamma) > 1){ 
      n_pop_input <- 200
    } else if ((subseq_interventions == "s" | subseq_interventions == "q") & data.hr.theta[i, "R_0_input"] * (1-gamma*prob_CT) > 1.1){
      n_pop_input <- 200
    } else {n_pop_input <- n_pop}
    In_Out <- repeat_call_fcn(n_pop=n_pop_input, 
                              parms_T_inc, 
                              parms_T_lat, 
                              parms_d_inf, 
                              parms_d_symp, 
                              parms_R_0_input, 
                              parms_epsilon, 
                              parms_pi_t,
                              num_generations,
                              background_intervention,
                              subseq_interventions,
                              gamma,
                              gamma_q,
                              prob_CT,
                              parms_CT_delay,
                              parms_serial_interval,
                              dispersion = dispersion,
                              printing = printing)
    if (subseq_interventions == background_intervention){
      data.hr.theta[i,"R_0_u"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"R"], w=In_Out$output[3:nrow(In_Out$output),"n"])
      data.hr.theta[i,"ks"]  <- weighted.mean(x=In_Out$output[2:nrow(In_Out$output),"ks"], w=In_Out$output[2:nrow(In_Out$output),"n"])
    }
    if (subseq_interventions == "hsb"){
      data.hr.theta[i,"R_hsb"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"R"], w=In_Out$output[3:nrow(In_Out$output),"n"])
    }
    if (subseq_interventions == "s"){
      data.hr.theta[i,"R_s"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"R"], w=In_Out$output[3:nrow(In_Out$output),"n"])
      data.hr.theta[i,"obs_to_iso_s"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"obs_to_iso"], w=In_Out$output[3:nrow(In_Out$output),"n"]) / 24
    }
    if (subseq_interventions == "q"){
      data.hr.theta[i,"R_q"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"R"], w=In_Out$output[3:nrow(In_Out$output),"n"])
      data.hr.theta[i,"obs_to_iso_q"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"obs_to_iso"], w=In_Out$output[3:nrow(In_Out$output),"n"]) / 24
    }
  }
}

data.hr.theta <- function_calculate_theta(data.hr.theta)
# View(data.hr.theta)
hist(data.hr.theta$theta)
plot(data.hr.theta$theta, data.hr.theta$R_q)
plot(data.hr.theta$theta, data.hr.theta$R_s)

data.hr.theta.melt <- melt(data.hr.theta[,c("theta", "R_s", "R_q", "R_0_input", "T_lat_offset")], id.vars = c("theta", "T_lat_offset"))

plot_theta_hr <- ggplot(data.hr.theta.melt, aes(x = theta, y = value, color = variable)) +
  geom_point() +
  stat_smooth() +
  scale_color_manual(values = c("R_0_input" = "grey", "R_s" = "yellow", "R_q" = "blue"), name = "Intervention", labels = c("Active Monitoring" , "Quarantine", "None (Input)")) +
  theme_bw() +
  geom_hline(yintercept = 1, color = "grey") +
  xlab("Proportion of Infections\nbefore Symptom Onset") +
  ylab("Reproductive Number")
plot_theta_hr

pdf(file=paste(root, "_Plot_Theta_HR.pdf", sep=""), width = 6, height = 4)
plot(plot_theta_hr)
dev.off()

data.hr.theta[,"Abs_Benefit"] <- data.hr.theta[,"R_s"] - data.hr.theta[,"R_q"]
data.hr.theta[,"Rel_Benefit"] <- data.hr.theta[,"Abs_Benefit"] / data.hr.theta[,"R_s"]
data.hr.theta[,"NNQ"] <- 1 / data.hr.theta[,"Abs_Benefit"]
data.hr.theta[data.hr.theta$NNQ < 1,"NNQ"] <- 1
data.hr.theta[data.hr.theta$NNQ > 9999,"NNQ"] <- 9999
data.hr.theta[data.hr.theta$NNQ == Inf,"NNQ"] <- 9999
data.hr.theta[,"Abs_Benefit_per_Qday"] <- data.hr.theta[,"Abs_Benefit"] / data.hr.theta[,"obs_to_iso_q"]

summary(data.hr.theta$R_0_u)
summary(data.hr.theta$R_hsb)
summary(data.hr.theta$R_s)
summary(data.hr.theta$R_q)
summary(data.hr.theta$Abs_Benefit)
summary(data.hr.theta$NNQ)

save.image(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", root, "_HR.Theta.RData", sep=""))


# What is the maximum theta that symptom monitoring could tolerate?
summary(data.hr.theta[data.hr.theta$R_s <= 1,"theta"])
summary(data.hr.theta[data.hr.theta$theta >= 0.014 & data.hr.theta$theta < 0.015,"R_s"])

summary(data.hr.theta[data.hr.theta$theta >= 0.375 & data.hr.theta$theta < 0.385,"R_q"])


#### Low Resource Setting ####

# Interventions
background_intervention <- "u"

prob_CT <- 0.5

gamma <- 0.5

gamma_q <- 0.5

parms_epsilon = list("uniform", 0, 4, 999, "independent", "independent")
names(parms_epsilon) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_CT_delay = list("uniform", 0, 4, 999, "independent", "independent")
names(parms_CT_delay) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

parms_R_0 = list("uniform", 2.2, 2.2, 999, "independent", "independent", 1) # mean = 2.249. 95% CI [1.46, 3.3]. Inspired by Li NEJM and Riou 2020 EuroSurveillance. May want to run sensitivity analysis with a wider CI, but can't get a big enough right skew to match Li or Riou's upper bound.
names(parms_R_0) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

# Settings
n_pop = 500
num_generations <- 5
times <- 2000
names <- c("R_0", "R_hsb", "R_s", "R_q", "Abs_Benefit","Rel_Benefit","NNQ","obs_to_iso_q","Abs_Benefit_per_Qday", "ks", "T_lat_offset", "d_inf", "pi_t_triangle_center")
data.lr.theta <- data.frame(matrix(rep(NA, length(names)*times), nrow=times))
names(data.lr.theta) <- names

# sample from joint posterior distribution
data <- function_calculate_theta(data)
sample <- sample(x = row.names(data), size = times, replace = TRUE, prob = (data$theta + 0.001*max(data.lr$theta))) # Updated to up-sample the high theta values
hist(data[sample, "theta"])
params.set <- cbind(
  T_lat_offset = data[sample, "T_lat_offset"] * runif(n = times, min = 0.25, max = 4), # permute between 25% to 400%
  d_inf = data[sample, "d_inf"] * runif(n = times, min = 0.25, max = 1), # permute between 25% to 100%
  pi_t_triangle_center = data[sample, "pi_t_triangle_center"] * runif(n = times, min = 0, max = 1), # permute between 0% to 100%
  dispersion = dispersion,
  R_0 = Draw_Dist_fcn(sample, parms_R_0$dist, parms_R_0$parm1, parms_R_0$parm2))

for (i in 1:times){
  cat(".")
  if (i%%10 == 0){cat("|")}
  if (i%%100 == 0){cat("\n")}
  
  data.lr.theta[i, "T_lat_offset"] <- parms_T_lat$anchor_value <- params.set[i,"T_lat_offset"]
  data.lr.theta[i, "d_inf"] <- parms_d_inf$parm2 <- max(1, params.set[i,"d_inf"]) # Make sure it's at least 1
  data.lr.theta[i, "pi_t_triangle_center"] <- parms_pi_t$triangle_center <- params.set[i,"pi_t_triangle_center"]
  
  data.lr.theta[i, "R_0_input"] <- params.set[i,"R_0"]
  parms_R_0_input <- list("uniform", as.numeric(params.set[i,"R_0"]), as.numeric(params.set[i,"R_0"]), 999, "independent", "independent", 1)
  names(parms_R_0_input) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")
  
  dispersion <- params.set[i, "dispersion"]
  
  for (subseq_interventions in c("s","q")){      
    if (subseq_interventions == background_intervention & data.lr.theta[i, "R_0_input"] > 1){
      n_pop_input <- 200
    } else if (subseq_interventions == "hsb" & data.lr.theta[i, "R_0_input"] * (1-gamma) > 1){ 
      n_pop_input <- 200
    } else if ((subseq_interventions == "s" | subseq_interventions == "q") & data.lr.theta[i, "R_0_input"] * (1-gamma*prob_CT) > 1.1){
      n_pop_input <- 200
    } else {n_pop_input <- n_pop}
    In_Out <- repeat_call_fcn(n_pop=n_pop_input, 
                              parms_T_inc, 
                              parms_T_lat, 
                              parms_d_inf, 
                              parms_d_symp, 
                              parms_R_0_input, 
                              parms_epsilon, 
                              parms_pi_t,
                              num_generations,
                              background_intervention,
                              subseq_interventions,
                              gamma,
                              gamma_q,
                              prob_CT,
                              parms_CT_delay,
                              parms_serial_interval,
                              dispersion = dispersion,
                              printing = printing)
    if (subseq_interventions == background_intervention){
      data.lr.theta[i,"R_0_u"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"R"], w=In_Out$output[3:nrow(In_Out$output),"n"])
      data.lr.theta[i,"ks"]  <- weighted.mean(x=In_Out$output[2:nrow(In_Out$output),"ks"], w=In_Out$output[2:nrow(In_Out$output),"n"])
    }
    if (subseq_interventions == "hsb"){
      data.lr.theta[i,"R_hsb"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"R"], w=In_Out$output[3:nrow(In_Out$output),"n"])
    }
    if (subseq_interventions == "s"){
      data.lr.theta[i,"R_s"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"R"], w=In_Out$output[3:nrow(In_Out$output),"n"])
      data.lr.theta[i,"obs_to_iso_s"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"obs_to_iso"], w=In_Out$output[3:nrow(In_Out$output),"n"]) / 24
    }
    if (subseq_interventions == "q"){
      data.lr.theta[i,"R_q"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"R"], w=In_Out$output[3:nrow(In_Out$output),"n"])
      data.lr.theta[i,"obs_to_iso_q"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"obs_to_iso"], w=In_Out$output[3:nrow(In_Out$output),"n"]) / 24
    }
  }
}

data.lr.theta <- function_calculate_theta(data.lr.theta)
# View(data.lr.theta)
hist(data.lr.theta$theta)
plot(data.lr.theta$theta, data.lr.theta$R_q)
plot(data.lr.theta$theta, data.lr.theta$R_s)

data.lr.theta.melt <- melt(data.lr.theta[,c("theta", "R_s", "R_q", "R_0_input")], id.vars = c("theta"))

plot_theta_lr <- ggplot(data.lr.theta.melt, aes(x = theta, y = value, color = variable)) +
  geom_point() +
  stat_smooth() +
  scale_color_manual(values = c("R_0_input" = "grey", "R_s" = "yellow", "R_q" = "blue"), name = "Intervention", labels = c("Active Monitoring" , "Quarantine", "None (Input)")) +
  theme_bw() +
  xlab("Proportion of Infections\nbefore Symptom Onset") +
  ylab("Reproductive Number")
plot_theta_lr

pdf(file=paste(root, "_Plot_Theta_LR.pdf", sep=""), width = 6, height = 4)
plot(plot_theta_lr)
dev.off()

data.lr.theta[,"Abs_Benefit"] <- data.lr.theta[,"R_s"] - data.lr.theta[,"R_q"]
data.lr.theta[,"Rel_Benefit"] <- data.lr.theta[,"Abs_Benefit"] / data.lr.theta[,"R_s"]
data.lr.theta[,"NNQ"] <- 1 / data.lr.theta[,"Abs_Benefit"]
data.lr.theta[data.lr.theta$NNQ < 1,"NNQ"] <- 1
data.lr.theta[data.lr.theta$NNQ > 9999,"NNQ"] <- 9999
data.lr.theta[data.lr.theta$NNQ == Inf,"NNQ"] <- 9999
data.lr.theta[,"Abs_Benefit_per_Qday"] <- data.lr.theta[,"Abs_Benefit"] / data.lr.theta[,"obs_to_iso_q"]

summary(data.lr.theta$R_0_u)
summary(data.lr.theta$R_hsb)
summary(data.lr.theta$R_s)
summary(data.lr.theta$R_q)
summary(data.lr.theta$Abs_Benefit)
summary(data.lr.theta$NNQ)

save.image(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", root, "_LR.Theta.RData", sep=""))

#### Plot ####
min = 0.017
max = 0.018
mean(data.hr.theta[data.hr.theta$theta >= min & data.hr.theta$theta < max, "R_s"]) - mean(data.hr.theta[data.hr.theta$theta >= min & data.hr.theta$theta < max, "R_q"])


#### Try reducing gamma_q ####

# Interventions
background_intervention = "u"

prob_CT <- 0.9

gamma <- 0.9

gamma_q <- 0.9

parms_epsilon = list("uniform", 0, 1, 999, "independent", "independent")
names(parms_epsilon) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_CT_delay = list("uniform", 0, 1, 999, "independent", "independent")
names(parms_CT_delay) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

parms_R_0 = list("lognormal", 2.2, 1.233, 999, "independent", "independent", 1) # mean = 2.249. 95% CI [1.46, 3.3]. Inspired by Li NEJM and Riou 2020 EuroSurveillance. May want to run sensitivity analysis with a wider CI, but can't get a big enough right skew to match Li or Riou's upper bound.
names(parms_R_0) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

# Settings
n_pop = 500
num_generations <- 5
times <- 6000
names <- c("R_0", "R_hsb", "R_s", "R_q", "Abs_Benefit","Rel_Benefit","NNQ","obs_to_iso_q","Abs_Benefit_per_Qday", "ks", "T_lat_offset", "d_inf", "pi_t_triangle_center", "gamma_q")
data.hr.theta.gamma_q <- data.frame(matrix(rep(NA, length(names)*times), nrow=times))
names(data.hr.theta.gamma_q) <- names

# sample from joint posterior distribution
data <- function_calculate_theta(data)
sample <- sample(x = row.names(data), size = times, replace = TRUE, prob = (data$theta + 0.001*max(data.hr$theta))) # Updated to up-sample the high theta values
hist(data[sample, "theta"])
params.set <- cbind(
  T_lat_offset = data[sample, "T_lat_offset"] * runif(n = times, min = 0.25, max = 4), # permute between 25% to 400%
  d_inf = data[sample, "d_inf"] * runif(n = times, min = 0.25, max = 1), # permute between 25% to 100%
  pi_t_triangle_center = data[sample, "pi_t_triangle_center"] * runif(n = times, min = 0, max = 1), # permute between 0% to 100%
  dispersion = dispersion,
  R_0 = Draw_Dist_fcn(sample, parms_R_0$dist, parms_R_0$parm1, parms_R_0$parm2),
  gamma_q = sample(x = c(0, .5, .9), size = times, replace = TRUE))

for (i in 1:times){
  cat(".")
  if (i%%10 == 0){cat("|")}
  if (i%%100 == 0){cat("\n")}
  
  data.hr.theta.gamma_q[i, "T_lat_offset"] <- parms_T_lat$anchor_value <- params.set[i,"T_lat_offset"]
  data.hr.theta.gamma_q[i, "d_inf"] <- parms_d_inf$parm2 <- max(1, params.set[i,"d_inf"]) # Make sure it's at least 1
  data.hr.theta.gamma_q[i, "pi_t_triangle_center"] <- parms_pi_t$triangle_center <- params.set[i,"pi_t_triangle_center"]
  
  data.hr.theta.gamma_q[i, "R_0_input"] <- params.set[i,"R_0"]
  parms_R_0_input <- list("uniform", as.numeric(params.set[i,"R_0"]), as.numeric(params.set[i,"R_0"]), 999, "independent", "independent", 1)
  names(parms_R_0_input) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")
  
  dispersion <- params.set[i, "dispersion"]
  
  data.hr.theta.gamma_q[i, "gamma_q"] <- gamma_q <- params.set[i, "gamma_q"]
  
  for (subseq_interventions in c("q")){      
    if (subseq_interventions == background_intervention & data.hr.theta.gamma_q[i, "R_0_input"] > 1){
      n_pop_input <- 200
    } else if (subseq_interventions == "hsb" & data.hr.theta.gamma_q[i, "R_0_input"] * (1-gamma) > 1){ 
      n_pop_input <- 200
    } else if ((subseq_interventions == "s" | subseq_interventions == "q") & data.hr.theta.gamma_q[i, "R_0_input"] * (1-gamma*prob_CT) > 1.1){
      n_pop_input <- 200
    } else {n_pop_input <- n_pop}
    In_Out <- repeat_call_fcn(n_pop=n_pop_input, 
                              parms_T_inc, 
                              parms_T_lat, 
                              parms_d_inf, 
                              parms_d_symp, 
                              parms_R_0_input, 
                              parms_epsilon, 
                              parms_pi_t,
                              num_generations,
                              background_intervention,
                              subseq_interventions,
                              gamma,
                              gamma_q,
                              prob_CT,
                              parms_CT_delay,
                              parms_serial_interval,
                              dispersion = dispersion,
                              printing = printing)
    if (subseq_interventions == background_intervention){
      data.hr.theta.gamma_q[i,"R_0_u"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"R"], w=In_Out$output[3:nrow(In_Out$output),"n"])
      data.hr.theta.gamma_q[i,"ks"]  <- weighted.mean(x=In_Out$output[2:nrow(In_Out$output),"ks"], w=In_Out$output[2:nrow(In_Out$output),"n"])
    }
    if (subseq_interventions == "hsb"){
      data.hr.theta.gamma_q[i,"R_hsb"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"R"], w=In_Out$output[3:nrow(In_Out$output),"n"])
    }
    if (subseq_interventions == "s"){
      data.hr.theta.gamma_q[i,"R_s"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"R"], w=In_Out$output[3:nrow(In_Out$output),"n"])
      data.hr.theta.gamma_q[i,"obs_to_iso_s"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"obs_to_iso"], w=In_Out$output[3:nrow(In_Out$output),"n"]) / 24
    }
    if (subseq_interventions == "q"){
      data.hr.theta.gamma_q[i,"R_q"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"R"], w=In_Out$output[3:nrow(In_Out$output),"n"])
      data.hr.theta.gamma_q[i,"obs_to_iso_q"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"obs_to_iso"], w=In_Out$output[3:nrow(In_Out$output),"n"]) / 24
    }
  }
}

data.hr.theta.gamma_q <- function_calculate_theta(data.hr.theta.gamma_q)
# View(data.hr.theta.gamma_q)
hist(data.hr.theta.gamma_q$theta)
plot(data.hr.theta.gamma_q$theta, data.hr.theta.gamma_q$R_q)
plot(data.hr.theta.gamma_q$theta, data.hr.theta.gamma_q$R_s)

data.hr.theta.gamma_q.melt <- melt(data.hr.theta.gamma_q[,c("theta", "R_s", "R_q", "R_0_input", "gamma_q")], id.vars = c("theta", "gamma_q"))

plot_theta_hr_gamma_q <- ggplot(data.hr.theta.gamma_q.melt, aes(x = theta, y = value, color = variable)) +
  facet_grid(gamma_q ~ .) +
  geom_point() +
  stat_smooth() +
  scale_color_manual(values = c("R_0_input" = "grey", "R_s" = "yellow", "R_q" = "blue"), name = "Intervention", labels = c("Symptom Monitoring" , "Quarantine", "None (Input)")) +
  theme_bw() +
  geom_hline(yintercept = 1, color = "grey") +
  xlab("Proportion of Infections\nbefore Symptom Onset") +
  ylab("Reproductive Number")
plot_theta_hr_gamma_q

pdf(file=paste(root, "_PlotTheta_HR_gamma_q.pdf", sep=""), width = 6, height = 4)
plot(plot_theta_hr)
dev.off()

data.hr.theta.gamma_q[,"Abs_Benefit"] <- data.hr.theta.gamma_q[,"R_s"] - data.hr.theta.gamma_q[,"R_q"]
data.hr.theta.gamma_q[,"Rel_Benefit"] <- data.hr.theta.gamma_q[,"Abs_Benefit"] / data.hr.theta.gamma_q[,"R_s"]
data.hr.theta.gamma_q[,"NNQ"] <- 1 / data.hr.theta.gamma_q[,"Abs_Benefit"]
data.hr.theta.gamma_q[data.hr.theta.gamma_q$NNQ < 1,"NNQ"] <- 1
data.hr.theta.gamma_q[data.hr.theta.gamma_q$NNQ > 9999,"NNQ"] <- 9999
data.hr.theta.gamma_q[data.hr.theta.gamma_q$NNQ == Inf,"NNQ"] <- 9999
data.hr.theta.gamma_q[,"Abs_Benefit_per_Qday"] <- data.hr.theta.gamma_q[,"Abs_Benefit"] / data.hr.theta.gamma_q[,"obs_to_iso_q"]

summary(data.hr.theta.gamma_q$R_0_u)
summary(data.hr.theta.gamma_q$R_hsb)
summary(data.hr.theta.gamma_q$R_s)
summary(data.hr.theta.gamma_q$R_q)
summary(data.hr.theta.gamma_q$Abs_Benefit)
summary(data.hr.theta.gamma_q$NNQ)

save.image(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", root, "_HR.Theta.RData", sep=""))

#### Explore obs to iso q for data.hr.theta.gamma_q ####

hist(data.hr.theta.gamma_q[data.hr.theta.gamma_q$gamma_q ==0,]$obs_to_iso_q, xlim = c(0, max(data.hr.theta.gamma_q$obs_to_iso_q)), main = "Gamma_Q = 0.0", xlab = "", breaks = 10, freq = FALSE)
hist(data.hr.theta.gamma_q[data.hr.theta.gamma_q$gamma_q ==0.5,]$obs_to_iso_q, xlim = c(0, max(data.hr.theta.gamma_q$obs_to_iso_q)), main = "Gamma_Q = 0.5", xlab = "", breaks = 10, freq = FALSE)
hist(data.hr.theta.gamma_q[data.hr.theta.gamma_q$gamma_q ==0.9,]$obs_to_iso_q, xlim = c(0, max(data.hr.theta.gamma_q$obs_to_iso_q)), main = "Gamma_Q = 0.9", xlab = "Days under Quarantine", breaks = 10, freq = FALSE)
