#### Header ####
# Create a figure that assesses the tolerance for gamma_q.
# gamm_q = the effectiveness of quarantine at reducing the infectiousness of people who haven't developed symtpoms yet
# How low can gamma_q get before SM beats Q?

# First, set R_0 = 2.25, gamma = 1, and find a stable R_SM
# Next, measure R_Q when gamma_q = 1 and when gamma_q = 0. 
  # If R_Q when gamma_q = 0 is worse than R_SM, then explore values of gamma_q between 0 and 1 to plot the tipping point
  # If R_Q when gamma_q = 0 is nearly the same as R_Q qhen gamma_q = 1, then this variable doesn't matter much for the disease. However, we expect it should since people are spending 3-4 days in quarantine before symptoms and this could include roughly 30-50% of their infectiousness, so gamma_q should matter.

# Plot gamma_q on x axis, R on y axis, including horizontal lines denoting R_0 and R_SM


# Interventions
background_intervention = "u"

prob_CT <- 0.9

gamma <- 0.9

gamma_q <- 0

parms_epsilon = list("uniform", 0, 1, 999, "independent", "independent")
names(parms_epsilon) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_CT_delay = list("uniform", 0, 1, 999, "independent", "independent")
names(parms_CT_delay) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")


parms_d_inf = list("uniform", 1, 8, 999, "independent", "independent")
names(parms_d_inf) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

# Settings
n_pop = 500
num_generations <- 5
times <- 90
names <- c("R_0", "R_hsb", "R_s", "R_q", "Abs_Benefit","Rel_Benefit","NNQ","obs_to_iso_q","Abs_Benefit_per_Qday", "gamma_q")
data.gamma_q <- data.frame(matrix(rep(NA, length(names)*times), nrow=times))
names(data.gamma_q) <- names

options <- c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1)
gamma_q_seq <- c(rep(options, each = times/length(options)))

# sample from joint posterior distribution
sample <- sample(x = row.names(data), size = times, replace = FALSE)
params.set <- cbind(
  T_lat_offset = data[sample, "T_lat_offset"],
  d_inf = data[sample, "d_inf"],
  pi_t_triangle_center = data[sample, "pi_t_triangle_center"],
  dispersion = runif(n=times, min = 1, max = 1),
  R_0 = runif(n = times, min = 2.2, max = 2.2)) # note this is changed

for (i in 1:times){
  cat(".")
  if (i%%10 == 0){cat("|")}
  if (i%%100 == 0){cat("\n")}
  
  gamma_q <- gamma_q_seq[i]
  
  parms_T_lat$anchor_value <- params.set[i,"T_lat_offset"]
  parms_d_inf$parm2 <- params.set[i,"d_inf"]
  parms_pi_t$triangle_center <- params.set[i,"pi_t_triangle_center"]
  parms_R_0$parm1 <- params.set[i,"R_0"]
  parms_R_0$parm2 <- params.set[i,"R_0"]
  dispersion <- params.set[i, "dispersion"]
  
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
      data.gamma_q[i,"R_0"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"R"], w=In_Out$output[3:nrow(In_Out$output),"n"])
    }
    if (subseq_interventions == "s"){
      data.gamma_q[i,"R_s"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"R"], w=In_Out$output[3:nrow(In_Out$output),"n"])
    }
    if (subseq_interventions == "q"){
      data.gamma_q[i,"R_q"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"R"], w=In_Out$output[3:nrow(In_Out$output),"n"])
      data.gamma_q[i,"obs_to_iso_q"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"obs_to_iso"], w=In_Out$output[3:nrow(In_Out$output),"n"]) / 24
    }
  }
}

data.gamma_q[,"gamma_q"] <- gamma_q_seq
data.gamma_q[,"Abs_Benefit"] <- data.gamma_q[,"R_s"] - data.gamma_q[,"R_q"]
data.gamma_q[,"Rel_Benefit"] <- data.gamma_q[,"Abs_Benefit"] / data.gamma_q[,"R_s"]
data.gamma_q[,"NNQ"] <- 1 / data.gamma_q[,"Abs_Benefit"]
data.gamma_q[data.gamma_q$NNQ < 1,"NNQ"] <- 1
data.gamma_q[data.gamma_q$NNQ > 9999,"NNQ"] <- 9999
data.gamma_q[data.gamma_q$NNQ == Inf,"NNQ"] <- 9999
data.gamma_q[,"Abs_Benefit_per_Qday"] <- data.gamma_q[,"Abs_Benefit"] / data.gamma_q[,"obs_to_iso_q"]
data.gamma_q$d_inf <- params.set[,"d_inf"]
data.gamma_q$pi_t_triangle_center <- params.set[,"pi_t_triangle_center"]
data.gamma_q$R_0_input <- params.set[,"R_0"]
data.gamma_q$T_lat_offset <- params.set[,"T_lat_offset"]
data.gamma_q$dispersion <- params.set[,"dispersion"]

summary(data.gamma_q$R_s)
summary(data.gamma_q$R_q)

ggplot(data.gamma_q, aes(x = gamma_q)) +
  # geom_point(aes(y = R_s), color = "red") +
  # geom_point(aes(y = R_q), color = "black") +
  # geom_hline(yintercept = 2.25) +
  stat_smooth(aes(y = R_s), color = "red") +
  stat_smooth(aes(y = R_q), color = "black")
  
hist(data.gamma_q$obs_to_iso_q)
