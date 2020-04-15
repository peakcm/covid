#### Header ####

# Figure showing relationship between T_lat_offset and R_s and R_q

# Hold R_0 constant at 2.2.

#### Load data ####
desired_root <- "20200202_2019nCoV" # Paste the desired root here "YYYYMMDD_DISEASE"

load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "/", desired_root, "_HR.Theta.RData", sep=""))

load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "/", desired_root, "_LR.Theta.RData", sep=""))

#### HR Setting, holding R constant at 2.2 ####

# Interventions
background_intervention = "u"

prob_CT <- 0.9

gamma <- 0.9

gamma_q <- 0.9

dispersion = 1/0.54

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
names <- c("R_0_input", "R_0_u", "R_hsb", "R_s", "R_q", "Abs_Benefit","Rel_Benefit","NNQ","obs_to_iso_q","obs_to_iso_s","Abs_Benefit_per_Qday", "ks", "T_lat_offset", "d_inf", "pi_t_triangle_center")
data.hr.constantR <- data.frame(matrix(rep(NA, length(names)*times), nrow=times))
names(data.hr.constantR) <- names

# sample from joint posterior distribution
sample <- sample(x = row.names(data), size = times, replace = TRUE) # Normal
params.set <- cbind(
  T_lat_offset = data[sample, "T_lat_offset"],
  d_inf = data[sample, "d_inf"],
  pi_t_triangle_center = data[sample, "pi_t_triangle_center"],
  dispersion = dispersion,
  R_0 = Draw_Dist_fcn(sample, parms_R_0$dist, parms_R_0$parm1, parms_R_0$parm2))

for (i in 1:times){
  cat(".")
  if (i%%10 == 0){cat("|")}
  if (i%%100 == 0){cat("\n")}
  
  data.hr.constantR[i, "T_lat_offset"] <- parms_T_lat$anchor_value <- params.set[i,"T_lat_offset"]
  data.hr.constantR[i, "d_inf"] <- parms_d_inf$parm2 <- params.set[i,"d_inf"]
  data.hr.constantR[i, "pi_t_triangle_center"] <- parms_pi_t$triangle_center <- params.set[i,"pi_t_triangle_center"]
  
  data.hr.constantR[i, "R_0_input"] <- params.set[i,"R_0"]
  parms_R_0_input <- list("uniform", as.numeric(params.set[i,"R_0"]), as.numeric(params.set[i,"R_0"]), 999, "independent", "independent", 1)
  names(parms_R_0_input) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")
  
  dispersion <- params.set[i, "dispersion"]
  
  for (subseq_interventions in c("s","q")){      
    if (subseq_interventions == background_intervention & data.hr.constantR[i, "R_0_input"] > 1){
      n_pop_input <- 200
    } else if (subseq_interventions == "hsb" & data.hr.constantR[i, "R_0_input"] * (1-gamma) > 1){ 
      n_pop_input <- 200
    } else if ((subseq_interventions == "s" | subseq_interventions == "q") & data.hr.constantR[i, "R_0_input"] * (1-gamma*prob_CT) > 1.1){
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
      data.hr.constantR[i,"R_0_u"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"R"], w=In_Out$output[3:nrow(In_Out$output),"n"])
      data.hr.constantR[i,"ks"]  <- weighted.mean(x=In_Out$output[2:nrow(In_Out$output),"ks"], w=In_Out$output[2:nrow(In_Out$output),"n"])
    }
    if (subseq_interventions == "hsb"){
      data.hr.constantR[i,"R_hsb"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"R"], w=In_Out$output[3:nrow(In_Out$output),"n"])
    }
    if (subseq_interventions == "s"){
      data.hr.constantR[i,"R_s"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"R"], w=In_Out$output[3:nrow(In_Out$output),"n"])
      data.hr.constantR[i,"obs_to_iso_s"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"obs_to_iso"], w=In_Out$output[3:nrow(In_Out$output),"n"]) / 24
    }
    if (subseq_interventions == "q"){
      data.hr.constantR[i,"R_q"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"R"], w=In_Out$output[3:nrow(In_Out$output),"n"])
      data.hr.constantR[i,"obs_to_iso_q"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"obs_to_iso"], w=In_Out$output[3:nrow(In_Out$output),"n"]) / 24
    }
  }
}

data.hr.constantR <- function_calculate_theta(data.hr.constantR)

data.hr.constantR[,"Abs_Benefit"] <- data.hr.constantR[,"R_s"] - data.hr.constantR[,"R_q"]
data.hr.constantR[,"Rel_Benefit"] <- data.hr.constantR[,"Abs_Benefit"] / data.hr.constantR[,"R_s"]
data.hr.constantR[,"NNQ"] <- 1 / data.hr.constantR[,"Abs_Benefit"]
data.hr.constantR[data.hr.constantR$NNQ < 1,"NNQ"] <- 1
data.hr.constantR[data.hr.constantR$NNQ > 9999,"NNQ"] <- 9999
data.hr.constantR[data.hr.constantR$NNQ == Inf,"NNQ"] <- 9999
data.hr.constantR[,"Abs_Benefit_per_Qday"] <- data.hr.constantR[,"Abs_Benefit"] / data.hr.constantR[,"obs_to_iso_q"]

# # Plot each of the covariate - outcome scatterplots
# for (covariate in names(data.mr)[11:15]){
#   panel_plot_fcn(data = data.mr, covariate = covariate)
#   cat ("Press [enter] to continue")
#   line <- readline()
# }

summary(data.hr.constantR$R_0_input)
summary(data.hr.constantR$R_hsb)
summary(data.hr.constantR$R_s)
summary(data.hr.constantR$R_q)
summary(data.hr.constantR$Abs_Benefit)
summary(data.hr.constantR$NNQ)
summary(data.hr.constantR$T_lat_offset)


summary(data.hr.constantR$obs_to_iso_q)
summary(data.hr.constantR$obs_to_iso_s)

median(data.hr.constantR$Abs_Benefit / data.hr.constantR$R_s) # Percentage improvement over SM
1/median(data.hr.constantR$Abs_Benefit) # NNQ if all contacts infected
1/median(data.hr.constantR$Abs_Benefit) * 1/0.5 # NNQ if 50% contacts infected
1/median(data.hr.constantR$Abs_Benefit) * 1/.05
1/median(data.hr.constantR$Abs_Benefit) * 1/(24/150000)

1/median(data.hr.constantR$Abs_Benefit) * mean(data.hr.constantR$obs_to_iso_q) # Number Q Days if all contacts infected
1/median(data.hr.constantR$Abs_Benefit) * mean(data.hr.constantR$obs_to_iso_q) + (1/median(data.hr.constantR$Abs_Benefit))*14 # Number of Q days if 50% contacts infected



#### Plot R as a function of T_Lat_Offset ####

# Trim a little to get rid of extremes
data.hr.constantR.melt <- melt(data.hr.constantR[data.hr.constantR$T_lat_offset < 1 & data.hr.constantR$T_lat_offset > -3.5,c("theta", "R_s", "R_q",  "T_lat_offset")], id.vars = c("theta", "T_lat_offset"))
# data.hr.constantR.melt <- melt(data.hr.constantR[,c("theta", "R_s", "R_q",  "T_lat_offset")], id.vars = c("theta", "T_lat_offset"))

# Plot 
plot_Tlat_hr <- ggplot(data.hr.constantR.melt, aes(x = T_lat_offset, y = value, color = variable)) +
  geom_point( size = 1, alpha = 0.3) +
  stat_smooth() +
  scale_color_manual(values = c("R_0_input" = "grey", "R_s" = "yellow", "R_q" = "blue"), name = "Intervention", labels = c("Active Monitoring" , "Quarantine", "None (Input)")) +
  theme_bw() +
  geom_hline(yintercept = 1, color = "grey") +
  geom_hline(yintercept = 2.2, color = "grey", size = 2) +
  ggtitle("Mean Serial Interval = 4.8 days") +
  xlab("Days between Onset of Infectiousness and Symptoms") +
  ylab("Reproductive Number")
plot_Tlat_hr

pdf(file=paste(root, "_Plot_TLatOffset_HR.pdf", sep=""), width = 6, height = 4)
plot(plot_Tlat_hr)
dev.off()
