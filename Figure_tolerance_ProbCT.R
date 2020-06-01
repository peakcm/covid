#### Header ####
# Create a figure that assesses the tolerance for prob_CT
# Assuming that people will be less enthusiastic about Q than for SM, we may get a lower fraction of contacts traced

# Plot Prob_CT on the x axis and R on the y axis. Plot horizontal threshold for R_0 and lines for each of R_SM and R_Q.
# Note R_SM when Prob_CT = 1. Find the value of Prob_CT that causes R_Q to equal R_SM. This is the threshold such that if we aren't able to trace at least X% of contacts for quarantine, we should just do symptom monitoring if that is tolerated such that all contacts can be traced.

# Define a range of Prob_CT

# Fix the R_0 input to 2.2

#### High resource setting ####
# Interventions
background_intervention = "u"

gamma <- 0.9
gamma_q <- 0.9

parms_epsilon = list("uniform", 0, 1, 999, "independent", "independent")
names(parms_epsilon) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_CT_delay = list("uniform", 0, 1, 999, "independent", "independent")
names(parms_CT_delay) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

parms_R_0 = list("uniform", 1.9, 1.9, 999, "independent", "independent", 1) # mean = 2.249. 95% CI [1.46, 3.3]. Inspired by Li NEJM and Riou 2020 EuroSurveillance. May want to run sensitivity analysis with a wider CI, but can't get a big enough right skew to match Li or Riou's upper bound.
names(parms_R_0) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

# Settings
n_pop = 100
num_generations <- 5
times <- 200
names <- c("R_0", "R_hsb", "R_s", "R_q", "Abs_Benefit","Rel_Benefit","NNQ","obs_to_iso_q","Abs_Benefit_per_Qday", "ks", "T_lat_offset", "d_inf", "pi_t_triangle_center", "Prob_CT")
data.hr.prob_CT <- data.frame(matrix(rep(NA, length(names)*times), nrow=times))
names(data.hr.prob_CT) <- names

# sample from joint posterior distribution
sample <- sample(x = row.names(data), size = times, replace = TRUE) # Updated to up-sample the high theta values
params.set <- cbind(
  T_lat_offset = data[sample, "T_lat_offset"], # permute between 25% to 400%
  d_inf = data[sample, "d_inf"], # permute between 25% to 100%
  pi_t_triangle_center = data[sample, "pi_t_triangle_center"], # permute between 0% to 100%
  dispersion = dispersion,
  R_0 = Draw_Dist_fcn(sample, parms_R_0$dist, parms_R_0$parm1, parms_R_0$parm2),
  Prob_CT = runif(n = times, 0,1))

# Store raw outputs
individuals = TRUE
data.hr.prob_CT.individuals <- list()

for (i in 1:times){
  cat(".")
  if (i%%10 == 0){cat("|")}
  if (i%%100 == 0){cat("\n")}
  
  data.hr.prob_CT[i, "T_lat_offset"] <- parms_T_lat$anchor_value <- params.set[i,"T_lat_offset"]
  data.hr.prob_CT[i, "d_inf"] <- parms_d_inf$parm2 <- max(1, params.set[i,"d_inf"]) # Make sure it's at least 1
  data.hr.prob_CT[i, "pi_t_triangle_center"] <- parms_pi_t$triangle_center <- params.set[i,"pi_t_triangle_center"]
  
  data.hr.prob_CT[i, "R_0_input"] <- params.set[i,"R_0"]
  parms_R_0_input <- list("uniform", as.numeric(params.set[i,"R_0"]), as.numeric(params.set[i,"R_0"]), 999, "independent", "independent", 1)
  names(parms_R_0_input) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")
  
  dispersion <- params.set[i, "dispersion"]
  
  data.hr.prob_CT[i, "Prob_CT"] <- prob_CT <- params.set[i, "Prob_CT"]
  
  for (subseq_interventions in c("s","q")){      
    if (subseq_interventions == background_intervention & data.hr.prob_CT[i, "R_0_input"] > 1){
      n_pop_input <- 200
    } else if (subseq_interventions == "hsb" & data.hr.prob_CT[i, "R_0_input"] * (1-gamma) > 1){ 
      n_pop_input <- 200
    } else if ((subseq_interventions == "s" | subseq_interventions == "q") & data.hr.prob_CT[i, "R_0_input"] * (1-gamma*prob_CT) > 1.1){
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
                              printing = printing,
                              individuals = individuals)
    if (subseq_interventions == background_intervention){
      data.hr.prob_CT[i,"R_0_u"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"R"], w=In_Out$output[3:nrow(In_Out$output),"n"])
      data.hr.prob_CT[i,"ks"]  <- weighted.mean(x=In_Out$output[2:nrow(In_Out$output),"ks"], w=In_Out$output[2:nrow(In_Out$output),"n"])
    }
    if (subseq_interventions == "hsb"){
      data.hr.prob_CT[i,"R_hsb"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"R"], w=In_Out$output[3:nrow(In_Out$output),"n"])
    }
    if (subseq_interventions == "s"){
      data.hr.prob_CT[i,"R_s"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"R"], w=In_Out$output[3:nrow(In_Out$output),"n"])
      data.hr.prob_CT[i,"obs_to_iso_s"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"obs_to_iso"], w=In_Out$output[3:nrow(In_Out$output),"n"]) / 24
    }
    if (subseq_interventions == "q"){
      data.hr.prob_CT[i,"R_q"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"R"], w=In_Out$output[3:nrow(In_Out$output),"n"])
      data.hr.prob_CT[i,"obs_to_iso_q"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"obs_to_iso"], w=In_Out$output[3:nrow(In_Out$output),"n"]) / 24
    }
    if (individuals == TRUE){
      data.hr.prob_CT.individuals[[length(data.hr.prob_CT.individuals)+1]] <- cbind(In_Out$output_individuals, "simulation" = i, "intervention" = subseq_interventions, "prob_CT" = prob_CT)
    }
  }
}

# View(data.hr.prob_CT)
hist(data.hr.prob_CT$Prob_CT)
plot(data.hr.prob_CT$Prob_CT, data.hr.prob_CT$R_q)
plot(data.hr.prob_CT$Prob_CT, data.hr.prob_CT$R_s)

data.hr.prob_CT.melt <- melt(data.hr.prob_CT[,c("Prob_CT", "R_s", "R_q", "R_0_input")], id.vars = c("Prob_CT"))

plot_prob_CT_hr <- ggplot(data.hr.prob_CT.melt, aes(x = Prob_CT, y = value, color = variable)) +
  geom_point() +
  stat_smooth(method = "lm") +
  scale_color_manual(values = c("R_0_input" = "grey", "R_s" = "yellow", "R_q" = "blue"), name = "Intervention", labels = c("Active Monitoring" , "Individual Quarantine", "None (Input)")) +
  theme_bw() +
  geom_hline(yintercept = 1, color = "grey") +
  xlab("Proportion of Contacts Traced") +
  ylab("Reproductive Number")
plot_prob_CT_hr

pdf(file=paste(root, "_Plot_Prob_CT_HR.pdf", sep=""), width = 6, height = 4)
plot(plot_prob_CT_hr)
dev.off()

save.image(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", root, "_HR.Prob_CT.RData", sep=""))

#### Calculate Re for each individual actually contact traced in each simulation, to see if Prob_CT influences individual performance, or just the average ####
data.hr.prob_CT.individuals <- bind_rows(data.hr.prob_CT.individuals)
head(data.hr.prob_CT.individuals)

# Calculate individual R_0
data.hr.prob_CT.individuals$R_individual <- NA

for (sim in 1:max(data.hr.prob_CT.individuals$simulation)){
  cat(".")
  if (sim%%10 == 0){cat("|")}
  if (sim%%100 == 0){cat("\n")}
  
  temp <- data.hr.prob_CT.individuals[data.hr.prob_CT.individuals$simulation == sim & data.hr.prob_CT.individuals$intervention == "s" & data.hr.prob_CT.individuals$generation < num_generations,]
  data.hr.prob_CT.individuals[data.hr.prob_CT.individuals$simulation == sim & data.hr.prob_CT.individuals$intervention == "s" & data.hr.prob_CT.individuals$generation < num_generations,"R_individual"] <-
    apply(temp, 1, function(x) length(temp[temp$generation == (as.numeric(x["generation"]) + 1) & temp$infector == as.numeric(x["ID"]),"ID"]))
  
  temp <- data.hr.prob_CT.individuals[data.hr.prob_CT.individuals$simulation == sim & data.hr.prob_CT.individuals$intervention == "q" & data.hr.prob_CT.individuals$generation < num_generations,]
  data.hr.prob_CT.individuals[data.hr.prob_CT.individuals$simulation == sim & data.hr.prob_CT.individuals$intervention == "q" & data.hr.prob_CT.individuals$generation < num_generations,"R_individual"] <-
    apply(temp, 1, function(x) length(temp[temp$generation == (as.numeric(x["generation"]) + 1) & temp$infector == as.numeric(x["ID"]),"ID"]))
}

# Indicate individuals who were traced before end of their disease cycle
data.hr.prob_CT.individuals$traced <- NA
data.hr.prob_CT.individuals$traced <- apply(data.hr.prob_CT.individuals, 1, function(x) (as.numeric(x["t_obs"]) < max(as.numeric(x["T_inc"]) + as.numeric(x["d_symp"]), as.numeric(x["T_lat"]) + as.numeric(x["d_inf"]))))

data.hr.prob_CT.individuals[data.hr.prob_CT.individuals$traced == TRUE,"traced"] <- "Contact Traced"
data.hr.prob_CT.individuals[data.hr.prob_CT.individuals$traced == FALSE,"traced"] <- "Contact Not Traced"

# Only use the first 50 from each generation after the first to avoid truncation
summary(data.hr.prob_CT.individuals$R_individual)
summary(data.hr.prob_CT.individuals[data.hr.prob_CT.individuals$ID <= 25, "R_individual"])
summary(data.hr.prob_CT.individuals[data.hr.prob_CT.individuals$ID <= 50 & 
                                      data.hr.prob_CT.individuals$generation > 1, "R_individual"])
summary(data.hr.prob_CT.individuals[data.hr.prob_CT.individuals$ID <= 50 & 
                                      data.hr.prob_CT.individuals$generation > 1 &
                                      data.hr.prob_CT.individuals$intervention == "q", "R_individual"])
summary(data.hr.prob_CT.individuals[data.hr.prob_CT.individuals$ID <= 50 & 
                                      data.hr.prob_CT.individuals$generation > 1 &
                                      data.hr.prob_CT.individuals$intervention == "s", "R_individual"])
summary(data.hr.prob_CT.individuals[data.hr.prob_CT.individuals$ID <= 50 & 
                                      data.hr.prob_CT.individuals$generation > 1 &
                                      data.hr.prob_CT.individuals$intervention == "q" & 
                                      data.hr.prob_CT.individuals$traced == "Contact Traced", "R_individual"])
summary(data.hr.prob_CT.individuals[data.hr.prob_CT.individuals$ID <= 50 & 
                                      data.hr.prob_CT.individuals$generation > 1 &
                                      data.hr.prob_CT.individuals$intervention == "s" & 
                                      data.hr.prob_CT.individuals$traced == "Contact Traced", "R_individual"])
summary(data.hr.prob_CT.individuals[data.hr.prob_CT.individuals$ID <= 50 & 
                                      data.hr.prob_CT.individuals$generation > 1 &
                                      data.hr.prob_CT.individuals$traced == "Contact Not Traced", "R_individual"])

figure_ProbCT_Individuals <- ggplot(data.hr.prob_CT.individuals[data.hr.prob_CT.individuals$ID <= 50 & data.hr.prob_CT.individuals$generation > 1,], aes(x = prob_CT, y = R_individual, color = intervention)) +
  theme_bw() +
  facet_grid(traced ~.) +
  ylab("Effective Reproductive Number") +
  scale_x_continuous(name = "Fraction of Contacts Traced") +
  scale_color_manual(values = c( q = "blue", s = "yellow"), name = "Intervention", labels = c("Individual Quarantine", "Active Monitoring")) +
  stat_smooth(method = "lm")
figure_ProbCT_Individuals

pdf(file=paste(root, "_Plot_Prob_CT_HR_Individuals.pdf", sep=""), width = 6, height = 4)
plot(figure_ProbCT_Individuals)
dev.off()



# Remove first two generations where intervention isn't settled in yet
data.hr.prob_CT.individuals <- data.hr.prob_CT.individuals[data.hr.prob_CT.individuals$generation > 2,] 
nrow(data.hr.prob_CT.individuals[data.hr.prob_CT.individuals$keep == FALSE,]) / nrow(data.hr.prob_CT.individuals) # Should be around the 1 - prob_CT
nrow(data.hr.prob_CT.individuals[data.hr.prob_CT.individuals$keep_2 == FALSE,]) / nrow(data.hr.prob_CT.individuals) 

data.hr.prob_CT.individuals <- data.hr.prob_CT.individuals[data.hr.prob_CT.individuals$keep == TRUE,]
data.hr.prob_CT.individuals <- data.hr.prob_CT.individuals[data.hr.prob_CT.individuals$keep_2 == TRUE,]




#### High resource setting, with larger theta values ####
# Interventions
background_intervention = "u"

gamma <- 0.9
gamma_q <- 0.9

parms_epsilon = list("uniform", 0, 1, 999, "independent", "independent")
names(parms_epsilon) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_CT_delay = list("uniform", 0, 1, 999, "independent", "independent")
names(parms_CT_delay) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

parms_R_0 = list("uniform", 2.2, 2.2, 999, "independent", "independent", 1) # mean = 2.249. 95% CI [1.46, 3.3]. Inspired by Li NEJM and Riou 2020 EuroSurveillance. May want to run sensitivity analysis with a wider CI, but can't get a big enough right skew to match Li or Riou's upper bound.
names(parms_R_0) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

# Settings
n_pop = 100
num_generations <- 5
times <- 200
names <- c("R_0", "R_hsb", "R_s", "R_q", "Abs_Benefit","Rel_Benefit","NNQ","obs_to_iso_q","Abs_Benefit_per_Qday", "ks", "T_lat_offset", "d_inf", "pi_t_triangle_center", "Prob_CT")
data.hr.prob_CT <- data.frame(matrix(rep(NA, length(names)*times), nrow=times))
names(data.hr.prob_CT) <- names

# sample from joint posterior distribution
sample <- sample(x = row.names(data), size = times, replace = TRUE, prob = (data$theta + 0.001*max(data.hr$theta))) # Updated to up-sample the high theta values
hist(data[sample, "theta"])
params.set <- cbind(
  T_lat_offset = data[sample, "T_lat_offset"] * runif(n = times, min = 0.25, max = 4), # permute between 25% to 400%
  d_inf = data[sample, "d_inf"] * runif(n = times, min = 0.25, max = 1), # permute between 25% to 100%
  pi_t_triangle_center = data[sample, "pi_t_triangle_center"] * runif(n = times, min = 0, max = 1), # permute between 0% to 100%
  dispersion = dispersion,
  R_0 = Draw_Dist_fcn(sample, parms_R_0$dist, parms_R_0$parm1, parms_R_0$parm2),
  Prob_CT = runif(n = times, 0,1))

for (i in 1:times){
  cat(".")
  if (i%%10 == 0){cat("|")}
  if (i%%100 == 0){cat("\n")}
  
  data.hr.prob_CT[i, "T_lat_offset"] <- parms_T_lat$anchor_value <- params.set[i,"T_lat_offset"]
  data.hr.prob_CT[i, "d_inf"] <- parms_d_inf$parm2 <- max(1, params.set[i,"d_inf"]) # Make sure it's at least 1
  data.hr.prob_CT[i, "pi_t_triangle_center"] <- parms_pi_t$triangle_center <- params.set[i,"pi_t_triangle_center"]
  
  data.hr.prob_CT[i, "R_0_input"] <- params.set[i,"R_0"]
  parms_R_0_input <- list("uniform", as.numeric(params.set[i,"R_0"]), as.numeric(params.set[i,"R_0"]), 999, "independent", "independent", 1)
  names(parms_R_0_input) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")
  
  dispersion <- params.set[i, "dispersion"]
  
  data.hr.prob_CT[i, "Prob_CT"] <- prob_CT <- params.set[i, "Prob_CT"]
  
  for (subseq_interventions in c("s","q")){      
    if (subseq_interventions == background_intervention & data.hr.prob_CT[i, "R_0_input"] > 1){
      n_pop_input <- 200
    } else if (subseq_interventions == "hsb" & data.hr.prob_CT[i, "R_0_input"] * (1-gamma) > 1){ 
      n_pop_input <- 200
    } else if ((subseq_interventions == "s" | subseq_interventions == "q") & data.hr.prob_CT[i, "R_0_input"] * (1-gamma*prob_CT) > 1.1){
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
      data.hr.prob_CT[i,"R_0_u"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"R"], w=In_Out$output[3:nrow(In_Out$output),"n"])
      data.hr.prob_CT[i,"ks"]  <- weighted.mean(x=In_Out$output[2:nrow(In_Out$output),"ks"], w=In_Out$output[2:nrow(In_Out$output),"n"])
    }
    if (subseq_interventions == "hsb"){
      data.hr.prob_CT[i,"R_hsb"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"R"], w=In_Out$output[3:nrow(In_Out$output),"n"])
    }
    if (subseq_interventions == "s"){
      data.hr.prob_CT[i,"R_s"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"R"], w=In_Out$output[3:nrow(In_Out$output),"n"])
      data.hr.prob_CT[i,"obs_to_iso_s"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"obs_to_iso"], w=In_Out$output[3:nrow(In_Out$output),"n"]) / 24
    }
    if (subseq_interventions == "q"){
      data.hr.prob_CT[i,"R_q"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"R"], w=In_Out$output[3:nrow(In_Out$output),"n"])
      data.hr.prob_CT[i,"obs_to_iso_q"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"obs_to_iso"], w=In_Out$output[3:nrow(In_Out$output),"n"]) / 24
    }
  }
}
data.hr.prob_CT <- function_calculate_theta(data.hr.prob_CT)

# Group theta values 
data.hr.prob_CT_grouped <- data.hr.prob_CT

# View(data.hr.prob_CT)
hist(data.hr.prob_CT$Prob_CT)
plot(data.hr.prob_CT$Prob_CT, data.hr.prob_CT$R_q)
plot(data.hr.prob_CT$Prob_CT, data.hr.prob_CT$R_s)

data.hr.prob_CT.melt <- melt(data.hr.prob_CT[data.hr.prob_CT$theta > 0. & data.hr.prob_CT$theta < 0.2,c("Prob_CT", "R_s", "R_q", "R_0_input", "theta")], id.vars = c("Prob_CT", "theta"))

plot_prob_CT_hr <- ggplot(data.hr.prob_CT.melt, aes(x = Prob_CT, y = value, color = variable)) +
  # facet_grid(theta~.) +
  geom_point() +
  stat_smooth() +
  scale_color_manual(values = c("R_0_input" = "grey", "R_s" = "yellow", "R_q" = "blue"), name = "Intervention", labels = c("Active Monitoring" , "Quarantine", "None (Input)")) +
  theme_bw() +
  geom_hline(yintercept = 1, color = "grey") +
  xlab("Proportion of Contacts Traced") +
  ylab("Reproductive Number")
plot_prob_CT_hr

pdf(file=paste(root, "_Plot_Prob_CT_HR.pdf", sep=""), width = 6, height = 4)
plot(plot_prob_CT_hr)
dev.off()



save.image(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", root, "_HR.Prob_CT.RData", sep=""))


