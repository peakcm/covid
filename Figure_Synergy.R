#### Header ####
# Figure to show the impact of quarantine and monitoring above and beyond social distancing-type interventions
# X-axis is (R_0)*(1 - efficacy of social distancing)
# Y-axis is the effective reproductive number (R_s, R_q, R_s - R_q)
# Line color is Prob_CT (90%, 50%, 10%)
  # Repeat this, showing only the R for individuals who are being traced.


# Interventions
background_intervention = "u"
prob_CT <- 0.9
gamma <- 0.9
gamma_q <- 0.75

parms_epsilon = list("uniform", 0, 1, 999, "independent", "independent")
names(parms_epsilon) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_CT_delay = list("uniform", 0, 1, 999, "independent", "independent")
names(parms_CT_delay) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

parms_d_inf = list("uniform", 1, 8, 999, "independent", "independent")
names(parms_d_inf) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

# Settings
n_pop = 500
num_generations <- 5
times <- 270
names <- c("R_0", "R_hsb", "R_s", "R_q", "Abs_Benefit","Rel_Benefit","NNQ","obs_to_iso_q","Abs_Benefit_per_Qday", "Prob_CT")
data.synergy <- data.frame(matrix(rep(NA, length(names)*times), nrow=times))
names(data.synergy) <- names

R_options <- c(0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25)
Prob_CT_options <- c(0.1, 0.5, 0.9)

data.synergy$R_0 <- c(rep(R_options, each = times/(length(R_options))))
data.synergy$Prob_CT <- c(rep(Prob_CT_options, times = times/length(Prob_CT_options)))

# sample from joint posterior distribution
sample <- sample(x = row.names(data), size = times, replace = FALSE)
params.set <- cbind(
  T_lat_offset = data[sample, "T_lat_offset"],
  d_inf = data[sample, "d_inf"],
  pi_t_triangle_center = data[sample, "pi_t_triangle_center"],
  dispersion = runif(n=times, min = 1, max = 1),
  R_0 = data.synergy$R_0,
  Prob_CT = data.synergy$Prob_CT) # note this is changed

for (i in 1:times){
  cat(".")
  if (i%%10 == 0){cat("|")}
  if (i%%100 == 0){cat("\n")}
  
  parms_T_lat$anchor_value <- params.set[i,"T_lat_offset"]
  parms_d_inf$parm2 <- params.set[i,"d_inf"]
  parms_pi_t$triangle_center <- params.set[i,"pi_t_triangle_center"]
  parms_R_0$parm1 <- params.set[i,"R_0"]
  parms_R_0$parm2 <- params.set[i,"R_0"]
  dispersion <- params.set[i, "dispersion"]
  
  prob_CT <- params.set[i, "Prob_CT"]
  
  for (subseq_interventions in c("s","q")){      
    if (subseq_interventions == background_intervention & parms_R_0$parm1 > 1){
      n_pop_input <- 200
    } else if (subseq_interventions == "hsb" & parms_R_0$parm1 * (1-gamma) > 1){ 
      n_pop_input <- 200
    } else if ((subseq_interventions == "s" | subseq_interventions == "q") & parms_R_0$parm1 * (1-gamma*prob_CT) > 1.1){
      n_pop_input <- 200
    } else {n_pop_input <- n_pop}
    In_Out <- repeat_call_fcn(n_pop=n_pop_input, 
                              parms_T_inc, 
                              parms_T_lat, 
                              parms_d_inf, 
                              parms_d_symp, 
                              parms_R_0, 
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
      data.synergy[i,"R_0"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"R"], w=In_Out$output[3:nrow(In_Out$output),"n"])
    }
    if (subseq_interventions == "s"){
      data.synergy[i,"R_s"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"R"], w=In_Out$output[3:nrow(In_Out$output),"n"])
    }
    if (subseq_interventions == "q"){
      data.synergy[i,"R_q"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"R"], w=In_Out$output[3:nrow(In_Out$output),"n"])
      data.synergy[i,"obs_to_iso_q"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"obs_to_iso"], w=In_Out$output[3:nrow(In_Out$output),"n"]) / 24
    }
  }
}

data.synergy[,"Abs_Benefit"] <- data.synergy[,"R_s"] - data.synergy[,"R_q"]
data.synergy[,"Rel_Benefit"] <- data.synergy[,"Abs_Benefit"]/data.synergy[,"R_s"]
data.synergy[,"Rel_Impact_S"] <- (data.synergy[,"R_0"] - data.synergy[,"R_s"])/data.synergy[,"R_0"]
data.synergy[,"Rel_Impact_Q"] <- (data.synergy[,"R_0"] - data.synergy[,"R_q"])/data.synergy[,"R_0"]

data.synergy %>%
  group_by(Prob_CT) %>%
  summary()

summary(data.synergy[data.synergy$R_0 == 1.25 & data.synergy$Prob_CT %in% c("50% of\nContacts Traced"),])

summary(data.synergy[data.synergy$Prob_CT %in% c("10% of\nContacts Traced"), ])
summary(data.synergy[data.synergy$Prob_CT %in% c("50% of\nContacts Traced"), ])
summary(data.synergy[data.synergy$Prob_CT %in% c("90% of\nContacts Traced"), ])


data.synergy[data.synergy$Prob_CT == 0.1,"Prob_CT"] <- "10% of\nContacts Traced"
data.synergy[data.synergy$Prob_CT == 0.5,"Prob_CT"] <- "50% of\nContacts Traced"
data.synergy[data.synergy$Prob_CT == 0.9,"Prob_CT"] <- "90% of\nContacts Traced"

data.synergy.melt <- melt(data.synergy[,c("R_0", "R_s", "R_q", "Prob_CT")], id.vars = c("Prob_CT", "R_0"))

# View(data.synergy)

figure_synergy <- ggplot(data.synergy.melt, aes(x = R_0, y = value, color = variable)) +
  facet_grid(Prob_CT ~.) +
  geom_abline(intercept = 0, slope = 1, col = "grey") +
  geom_hline(yintercept = 1, col = "grey") +
  geom_point(alpha = 0.3) + 
  stat_smooth(method = "lm") +
  scale_color_manual(values = c("R_s" = "yellow", "R_q" = "blue"), name = "Intervention", labels = c("Active Monitoring" , "Individual Quarantine")) +
  ylab("Reproductive Number adding\nContact-Based Interventions") +
  xlab("Reproductive Number\nwith Social Distancing") +
  theme_bw() +
  theme(panel.grid = element_blank())
figure_synergy

pdf(file=paste(root, "_Plot_Synergy.pdf", sep=""), width = 6, height = 4)
plot(figure_synergy)
dev.off()


#### Save image ####
save.image(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", root, "_Synergy.RData", sep=""))
                                                                                                                                                                                                                                                                    

