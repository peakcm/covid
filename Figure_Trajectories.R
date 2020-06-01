#### Header ####
# Make a spaghetti plot with many trials starting with one individual and tracing how many are in her infection tree as a function of generations. Then apply an intervention in generation 4 or so. Color code by u, hsb, s, q

#### Load Libraries ####
library(ggplot2)

#### Load PlotTrajectories Workspace ####
desired_root <- "20200202_2019nCoV"
load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "/", desired_root, "_PlotTrajectories.RData", sep=""))

#### Load SMC Workspaces ####
desired_root <- "20200128_2019nCoV"
load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "/", desired_root, "_SMC.RData", sep=""))

#### Disease: Ebola ####
# Name the trial
date <- format(Sys.time(), "%Y%m%d")
disease <- "Ebola"
root <- paste(date, disease, sep = "_")

# Fixed Disease Parameters
parms_serial_interval <- list("gamma", 2.5, 0.2) # approximation from WHO
names(parms_serial_interval) <- c("dist","parm1","parm2")

parms_T_inc = list("gamma", 1.75, 0.182, 999, "independent", "independent") # approximation from WHO
names(parms_T_inc) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

parms_R_0 = list("uniform", 1.83, 1.83, 999, "independent", "independent") 
names(parms_R_0) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

# Variable Disease Parameters
parms_pi_t <- list("triangle", 0.50)
names(parms_pi_t) <- c("distribution","triangle_center")

parms_T_lat = list("triangle", 999, 999, 999, 0, "T_inc")
names(parms_T_lat) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_d_inf = list("uniform", 1, 8, 999, "independent", "independent")
names(parms_d_inf) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_d_symp = list("uniform", 1, 8, 999, 0, "d_inf")
names(parms_d_symp) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")
background_intervention = "u"

#### Disease: InfluenzaA ####

# Name the trial
date <- format(Sys.time(), "%Y%m%d")
disease <- "InfluenzaA"
root <- paste(date, disease, sep = "_")

# Fixed Disease Parameters
parms_serial_interval <- list("normal", 2.2, 0.8) 
names(parms_serial_interval) <- c("dist","parm1","parm2")

parms_T_inc = list("lognormal", 1.4, 1.5, 999, "independent", "independent") # Using Vink 2014, cf Fine 2003
names(parms_T_inc) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

parms_R_0 = list("uniform", 1.46, 1.46, 999, "independent", "independent") 
names(parms_R_0) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

# Variable Disease Parameters
parms_pi_t <- list("triangle", 0.50)
names(parms_pi_t) <- c("distribution","triangle_center")

parms_T_lat = list("triangle", 999, 999, 999, 0, "T_inc")
names(parms_T_lat) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_d_inf = list("uniform", 1, 8, 999, "independent", "independent")
names(parms_d_inf) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_d_symp = list("uniform", 1, 8, 999, 0, "d_inf")
names(parms_d_symp) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

#### Intervention ####

# Interventions
background_intervention <- "u"

prob_CT <- 0.9

gamma <- 0.9

gamma_q <- 0.75

parms_epsilon = list("uniform", 0, 1, 999, "independent", "independent")
names(parms_epsilon) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_CT_delay = list("uniform", 0, 1, 999, "independent", "independent")
names(parms_CT_delay) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

#### Model settings ####
n_pop = round(100/2.2)
num_generations_intervention <- 6
times = 10
names <- c("generation", "count", "trial", "intervention")
data_traj <- data.frame(matrix(rep(NA, length(names)*times*(num_generations_u + num_generations_intervention - 1)*4), nrow=times*(num_generations_u + num_generations_intervention - 1)*4))
names(data_traj) <- names
data_traj$generation <- rep(seq(1, (num_generations_u + num_generations_intervention - 1)), times = times*4)
data_traj$trial <- rep(seq(1:(times*4)), each = (num_generations_u + num_generations_intervention - 1))
data_traj$intervention <- rep(c("u", "hsb", "s", "q"), each = times * (num_generations_u + num_generations_intervention - 1))
data_traj$intervention <- factor(data_traj$intervention, levels = c("u", "hsb", "s", "q"))

#### Draw the median from SMC to represent each iteration ####
data <- outputs$data
# This is better than each iteration having a different set of attributes because the stochasticity in the model can be shown instead of just the uncertainty in the model inputs.
parms_T_lat$anchor_value <- median(data[,"T_lat_offset"])
parms_d_inf$parm1 <- median(data[,"d_inf"])
parms_d_inf$parm2 <- median(data[,"d_inf"])
parms_pi_t$distribution <- "uniform"
parms_pi_t$triangle_center<- median(data[,"pi_t_triangle_center"])


#### Run model (NEW) ####
individuals <- TRUE
data.individuals.traj <- list()

for (i in 1:times){
  cat("\ni = ", i)
  
  for (subseq_interventions in c("q")){      
    
    In_Out.post <- repeat_call_fcn(n_pop = n_pop, 
                                   parms_T_inc = parms_T_inc, 
                                   parms_T_lat = parms_T_lat, 
                                   parms_d_inf = parms_d_inf, 
                                   parms_d_symp = parms_d_symp, 
                                   parms_R_0 = parms_R_0, 
                                   parms_epsilon = parms_epsilon, 
                                   parms_pi_t = parms_pi_t,
                                   num_generations = num_generations_intervention,
                                   background_intervention = background_intervention,
                                   subseq_interventions = subseq_interventions,
                                   gamma = gamma,
                                   gamma_q = gamma_q,
                                   prob_CT = prob_CT,
                                   parms_CT_delay = parms_CT_delay,
                                   parms_serial_interval = parms_serial_interval,
                                   cap_pop = FALSE,
                                   min_infections = 1,
                                   printing = FALSE,
                                   dispersion = 1,
                                   individuals = TRUE)
    post <- In_Out.post$output[,"n"]
    
    if (individuals == TRUE){
      data.individuals.traj[[length(data.individuals.traj)+1]] <- cbind(In_Out.post$output_individuals, "simulation" = i, "intervention" = subseq_interventions)
    }
  }
}

#### Transform data ####
data.individuals.traj <- bind_rows(data.individuals.traj)
unique(data.individuals.traj$intervention)
tail(data.individuals.traj)

# Subset to only generation 2
nrow(data.individuals.traj[data.individuals.traj$generation ==2,]) # Check to see if number of infections is around 1,000
data.individuals.traj <- data.individuals.traj[data.individuals.traj$generation > 1,]

data.individuals.traj.histories <- data.frame(day = (1:(max(data.individuals.traj$T_inc)/24 + 14)), incident_infections = 0, incident_intervention = 0, under_intervention = 0, excess_1 = 0, excess_10 = 0, excess_6250 = 0)

# Calculate number of incident infections per hour
data.individuals.traj.histories[, "incident_infections"] <- apply(data.individuals.traj.histories, 1, function(x) nrow(data.individuals.traj[data.individuals.traj$T_inc >= ((x["day"]-1)*24) &
                                                                                                                                              data.individuals.traj$T_inc < ((x["day"])*24),]))

data.individuals.traj.histories[, "incident_intervention"] <- apply(data.individuals.traj.histories, 1, function(x) nrow(data.individuals.traj[data.individuals.traj$t_obs >= ((x["day"]-1)*24) &
                                                                                                                                                 data.individuals.traj$t_obs < ((x["day"])*24),]))

data.individuals.traj.histories[,"under_intervention"] <- apply(data.individuals.traj.histories, 1, function(x) nrow(data.individuals.traj[data.individuals.traj$t_obs <= as.numeric((x["day"]-1)*24) &
                                                                                                                                             data.individuals.traj$t_iso >= as.numeric((x["day"])*24),]))
for (day in 1:(nrow(data.individuals.traj.histories) - 14)){
  data.individuals.traj.histories[day:(day + 13), "excess_1"] <-  data.individuals.traj.histories[day:(day + 13), "excess_1"] + 1*data.individuals.traj.histories[day, "incident_intervention"]
  data.individuals.traj.histories[day:(day + 13), "excess_10"] <-  data.individuals.traj.histories[day:(day + 13), "excess_10"] + 99*data.individuals.traj.histories[day, "incident_intervention"]
  data.individuals.traj.histories[day:(day + 13), "excess_6250"] <-  data.individuals.traj.histories[day:(day + 13), "excess_6250"] + 6249*data.individuals.traj.histories[day, "incident_intervention"]
}

# Measure cumulative infections instead
data.individuals.traj.histories$cumulative_infections <- cumsum(data.individuals.traj.histories$incident_infections) + 1000

# 
# plot(data.individuals.traj.histories$day, data.individuals.traj.histories$incident_infections)
# plot(data.individuals.traj.histories$day, data.individuals.traj.histories$incident_intervention)
# plot(data.individuals.traj.histories$day, data.individuals.traj.histories$under_intervention)
# plot(data.individuals.traj.histories$day, data.individuals.traj.histories$excess_10)
# plot(data.individuals.traj.histories$day, data.individuals.traj.histories$excess_100)
# plot(data.individuals.traj.histories$day, data.individuals.traj.histories$excess_6250)

data.individuals.traj.histories.melt <- melt(data.individuals.traj.histories[,c("day", "incident_infections", "under_intervention", "excess_1")], id.vars = c("day"))

figure_under_intervention <- ggplot(data.individuals.traj.histories.melt, aes(x = day, y = value, color = variable)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = c(incident_infections = "red", under_intervention = "blue", excess_1 = "forestgreen"), labels = c("\nIncident Infections\n", "\nCurrently Under Quarantine\n(Truly Infected)\n", "\nCurrently Under Quarantine\n(1:1 Not Infected)\n"), name = c()) +
  # scale_x_continuous(limits = c(0, 20)) +
  xlab("Day") +
  ylab("Count") +
  ggtitle("Serial Interval = 4.8 Days") +
  theme_bw() +
  theme(legend.position = "right") 
figure_under_intervention

pdf(file=paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", root, "_PlotUnderIntervention.pdf", sep=""), width = 7, height = 4) # Dimensions are in inches
# pdf(file=paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", root, "_PlotTraj_legend.pdf", sep=""), width = 4, height = 1) # Dimensions are in inches
plot(figure_under_intervention)
dev.off()


#### Intervention (LR) ####

# Interventions
background_intervention <- "u"

prob_CT <- 0.25

gamma <- 0.1

gamma_q <- 0.1

parms_epsilon = list("uniform", 0, 4, 999, "independent", "independent")
names(parms_epsilon) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_CT_delay = list("uniform", 0, 4, 999, "independent", "independent")
names(parms_CT_delay) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

#### Model settings ####
n_pop = round(100/2.2)
num_generations_intervention <- 6
times = 10
names <- c("generation", "count", "trial", "intervention")
data_traj <- data.frame(matrix(rep(NA, length(names)*times*(num_generations_u + num_generations_intervention - 1)*4), nrow=times*(num_generations_u + num_generations_intervention - 1)*4))
names(data_traj) <- names
data_traj$generation <- rep(seq(1, (num_generations_u + num_generations_intervention - 1)), times = times*4)
data_traj$trial <- rep(seq(1:(times*4)), each = (num_generations_u + num_generations_intervention - 1))
data_traj$intervention <- rep(c("u", "hsb", "s", "q"), each = times * (num_generations_u + num_generations_intervention - 1))
data_traj$intervention <- factor(data_traj$intervention, levels = c("u", "hsb", "s", "q"))

#### Draw the median from SMC to represent each iteration ####
data <- outputs$data
# This is better than each iteration having a different set of attributes because the stochasticity in the model can be shown instead of just the uncertainty in the model inputs.
parms_T_lat$anchor_value <- median(data[,"T_lat_offset"])
parms_d_inf$parm1 <- median(data[,"d_inf"])
parms_d_inf$parm2 <- median(data[,"d_inf"])
parms_pi_t$distribution <- "uniform"
parms_pi_t$triangle_center<- median(data[,"pi_t_triangle_center"])

#### Run model (NEW) ####
individuals <- TRUE
data.individuals.traj.lr <- list()

for (i in 1:times){
  cat("\ni = ", i)
  
  for (subseq_interventions in c("q")){      
    
    In_Out.post <- repeat_call_fcn(n_pop = n_pop, 
                                   parms_T_inc = parms_T_inc, 
                                   parms_T_lat = parms_T_lat, 
                                   parms_d_inf = parms_d_inf, 
                                   parms_d_symp = parms_d_symp, 
                                   parms_R_0 = parms_R_0, 
                                   parms_epsilon = parms_epsilon, 
                                   parms_pi_t = parms_pi_t,
                                   num_generations = num_generations_intervention,
                                   background_intervention = background_intervention,
                                   subseq_interventions = subseq_interventions,
                                   gamma = gamma,
                                   gamma_q = gamma_q,
                                   prob_CT = prob_CT,
                                   parms_CT_delay = parms_CT_delay,
                                   parms_serial_interval = parms_serial_interval,
                                   cap_pop = FALSE,
                                   min_infections = 1,
                                   printing = FALSE,
                                   dispersion = 1,
                                   individuals = TRUE)
    post <- In_Out.post$output[,"n"]
    
    if (individuals == TRUE){
      data.individuals.traj.lr[[length(data.individuals.traj.lr)+1]] <- cbind(In_Out.post$output_individuals, "simulation" = i, "intervention" = subseq_interventions)
    }
  }
}

#### Transform data ####
data.individuals.traj.lr <- bind_rows(data.individuals.traj.lr)
unique(data.individuals.traj.lr$intervention)
tail(data.individuals.traj.lr)

# Subset to only generation 2
nrow(data.individuals.traj.lr[data.individuals.traj.lr$generation ==2,]) # Check to see if number of infections is around 1,000
data.individuals.traj.lr <- data.individuals.traj.lr[data.individuals.traj.lr$generation > 1,]

data.individuals.traj.histories.lr <- data.frame(day = (1:(max(data.individuals.traj.lr$T_inc)/24 + 14)), incident_infections = 0, incident_intervention = 0, under_intervention = 0, excess_1 = 0, excess_10 = 0, excess_6250 = 0)

# Calculate number of incident infections per hour
data.individuals.traj.histories.lr[, "incident_infections"] <- apply(data.individuals.traj.histories.lr, 1, function(x) nrow(data.individuals.traj.lr[data.individuals.traj.lr$T_inc >= ((x["day"]-1)*24) &
                                                                                                                                                        data.individuals.traj.lr$T_inc < ((x["day"])*24),]))

data.individuals.traj.histories.lr[, "incident_intervention"] <- apply(data.individuals.traj.histories.lr, 1, function(x) nrow(data.individuals.traj.lr[data.individuals.traj.lr$t_obs >= ((x["day"]-1)*24) &
                                                                                                                                                       data.individuals.traj.lr$t_obs < ((x["day"])*24),]))

data.individuals.traj.histories.lr[,"under_intervention"] <- apply(data.individuals.traj.histories.lr, 1, function(x) nrow(data.individuals.traj.lr[data.individuals.traj.lr$t_obs <= as.numeric((x["day"]-1)*24) &
                                                                                                                                                      data.individuals.traj.lr$t_iso >= as.numeric((x["day"])*24),]))
for (day in 1:(nrow(data.individuals.traj.histories.lr) - 14)){
  data.individuals.traj.histories.lr[day:(day + 13), "excess_1"] <-  data.individuals.traj.histories.lr[day:(day + 13), "excess_1"] + 1*data.individuals.traj.histories.lr[day, "incident_intervention"]
  data.individuals.traj.histories.lr[day:(day + 13), "excess_10"] <-  data.individuals.traj.histories.lr[day:(day + 13), "excess_10"] + 9*data.individuals.traj.histories.lr[day, "incident_intervention"]
  data.individuals.traj.histories.lr[day:(day + 13), "excess_6250"] <-  data.individuals.traj.histories.lr[day:(day + 13), "excess_6250"] + 6249*data.individuals.traj.histories.lr[day, "incident_intervention"]
}

# Measure cumulative infections instead
data.individuals.traj.histories.lr$cumulative_infections <- cumsum(data.individuals.traj.histories.lr$incident_infections) + 1000

# 
# plot(data.individuals.traj.histories$day, data.individuals.traj.histories$incident_infections)
# plot(data.individuals.traj.histories$day, data.individuals.traj.histories$incident_intervention)
# plot(data.individuals.traj.histories$day, data.individuals.traj.histories$under_intervention)
# plot(data.individuals.traj.histories$day, data.individuals.traj.histories$excess_10)
# plot(data.individuals.traj.histories$day, data.individuals.traj.histories$excess_100)
# plot(data.individuals.traj.histories$day, data.individuals.traj.histories$excess_6250)

data.individuals.traj.histories.lr.melt <- melt(data.individuals.traj.histories.lr[,c("day", "cumulative_infections", "under_intervention", "excess_1", "excess_10")], id.vars = c("day"))

figure_under_intervention_lr <- ggplot(data.individuals.traj.histories.lr.melt, aes(x = day, y = value, color = variable)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = c(cumulative_infections = "red", under_intervention = "blue", excess_1 = "forestgreen", excess_10 = "green"), labels = c("\nCumulative Infections\n", "\nCurrently Under Quarantine\n(Truly Infected)\n", "\nCurrently Under Quarantine\n(1:1 Not Infected)\n", "\nCurrently Under Quarantine\n(9:1 Not Infected)\n"), name = c()) +
  scale_x_continuous(limits = c(0, 20)) +
  scale_y_continuous(limits = c(0, 40000)) +
  geom_label(data = data.individuals.traj.histories.lr.melt[data.individuals.traj.histories.lr.melt$day == 20,], aes(x = day, y = value,  label = value), color = "black", size = 3) +
  geom_label(data = data.individuals.traj.histories.lr.melt[data.individuals.traj.histories.lr.melt$day == 1 ,], aes(x = day, y = value,  label = value), color = "black", size = 3) +
  geom_label(data = data.individuals.traj.histories.lr.melt[data.individuals.traj.histories.lr.melt$day == 1 & data.individuals.traj.histories.lr.melt$variable == "cumulative_infections",], aes(x = day, y = value,  label = value), color = "black", size = 3) +
  xlab("Day") +
  ylab("Count") +
  ggtitle("Serial Interval = 4.8 Days") +
  theme_bw() +
  theme(legend.position = "right") 
figure_under_intervention_lr

data.individuals.traj.histories.lr[20,] #starts at 1,000 cases

pdf(file=paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", root, "_PlotUnderIntervention_lr.pdf", sep=""), width = 7, height = 4) # Dimensions are in inches
# pdf(file=paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", root, "_PlotTraj_legend.pdf", sep=""), width = 4, height = 1) # Dimensions are in inches
plot(figure_under_intervention_lr)
dev.off()


#### Run model (OLD) ####
for (i in 1:times){
  cat("\ni = ", i)
  
  for (subseq_interventions in c("u", "hsb", "s","q")){      
    In_Out.pre <- repeat_call_fcn(n_pop = n_pop, 
                                  parms_T_inc = parms_T_inc, 
                                  parms_T_lat = parms_T_lat, 
                                  parms_d_inf = parms_d_inf, 
                                  parms_d_symp = parms_d_symp, 
                                  parms_R_0 = parms_R_0, 
                                  parms_epsilon = parms_epsilon, 
                                  parms_pi_t = parms_pi_t,
                                  num_generations = num_generations_u,
                                  background_intervention = background_intervention,
                                  subseq_interventions = "u",
                                  gamma = gamma,
                                  gamma_q = gamma_q,
                                  prob_CT = prob_CT,
                                  parms_CT_delay = parms_CT_delay,
                                  parms_serial_interval = parms_serial_interval,
                                  cap_pop = FALSE,
                                  min_infections = 1,
                                  printing = FALSE,
                                  dispersion = 1,
                                  individuals = TRUE)
    pre <- In_Out.pre$output[,"n"]
    
    In_Out.post <- repeat_call_fcn(n_pop = In_Out.pre$output[num_generations_u,"n"], 
                                   parms_T_inc = parms_T_inc, 
                                   parms_T_lat = parms_T_lat, 
                                   parms_d_inf = parms_d_inf, 
                                   parms_d_symp = parms_d_symp, 
                                   parms_R_0 = parms_R_0, 
                                   parms_epsilon = parms_epsilon, 
                                   parms_pi_t = parms_pi_t,
                                   num_generations = num_generations_intervention,
                                   background_intervention = background_intervention,
                                   subseq_interventions = subseq_interventions,
                                   gamma = gamma,
                                   gamma_q = gamma_q,
                                   prob_CT = prob_CT,
                                   parms_CT_delay = parms_CT_delay,
                                   parms_serial_interval = parms_serial_interval,
                                   cap_pop = FALSE,
                                   min_infections = 1,
                                   printing = FALSE,
                                   dispersion = 1,
                                   individuals = TRUE)
    post <- In_Out.post$output[,"n"]
    
    count <- c(pre, post[-1])
    if (length(count) < (num_generations_u + num_generations_intervention - 1)){
      cat(" [Error: count is too short]")
      diff <- (num_generations_u + num_generations_intervention - 1) - length(count)
      count <- c(count, rep(0, diff))
    }
    data_traj[is.na(data_traj$count) == 1 & data_traj$intervention == subseq_interventions, "count"][1:length(count)] <- count
    
  }
}

#### Plot Results (OLD) ####
plot_traj <- ggplot(data_traj, aes(x=generation, y=count, group=trial, color=intervention)) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_vline(xintercept=(num_generations_u+1), col = "grey") +
  geom_line(alpha = 0.15) +
  stat_smooth(aes(group = intervention), size = .5, method = "loess", se = FALSE, alpha = 0.5) +
  scale_x_continuous(breaks = seq(1:5), labels = c(-2, -1, 0, 1, 2)) +
  xlab("Generation") + ylab("Incident Cases") +
  scale_color_manual(name="Intervention:\n",
                     values = c("coral3", "mediumturquoise", "darkgoldenrod", "cornflowerblue"),
                       breaks=c("u", "hsb", "s", "q"),
                       labels=c("\nNone\n", "\nHealth-\nSeeking\nBehavior\n", "\nActive\nMonitoring\n", "\nQuarantine\n")) +
  theme(legend.direction = "vertical", 
        legend.position = "right",
        legend.key=element_rect(size=5, color="white"))+
  # guides(colour = FALSE) +
  # theme(legend.position="bottom",legend.direction="horizontal") +
  # scale_y_continuous(breaks = seq(0, 1400, by=200)) + ggtitle("Ebola") +
  scale_y_continuous(breaks = seq(0, 600, by=200)) + ggtitle("COVID-19") +
  theme(text = element_text(size=8))
plot_traj

pdf(file=paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", root, "_PlotTraj.pdf", sep=""), width = 4, height = 4) # Dimensions are in inches
# pdf(file=paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", root, "_PlotTraj_legend.pdf", sep=""), width = 4, height = 1) # Dimensions are in inches
plot(plot_traj)
dev.off()

#### Save Workspace ####
save.image(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", root, "_PlotTrajectories.RData", sep=""))



