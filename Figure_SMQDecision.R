#### Header ####
# Code for SMQDecision plot

#### Load Libraries ####
library(ggplot2)
library(reshape2)
library(dplyr)

#### Load Workspace ####
desired_root <- "20200202_2019nCoV" # Paste the desired root here "YYYYMMDD_DISEASE"

# If workspaces are in main folder
# load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "_HR.RData", sep=""))
# load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "_MR.RData", sep=""))
# load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "_LR.RData", sep=""))
# load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "_Plots.RData", sep=""))

# If workspaces are in their own folder, named the same as the root
load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "/", desired_root, "_HR.RData", sep=""))
load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "/", desired_root, "_MR.RData", sep=""))
load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "/", desired_root, "_LR.RData", sep=""))
load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "/", desired_root, "_Plots.RData", sep=""))

load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "_SMQDecision.RData", sep=""))

#### Create dataframe ####
R_0_relevant.min <- 1.4
R_0_relevant.max <- 3.9

data.hr$Setting <- "HR"
data.mr$Setting <- "MR"
data.lr$Setting <- "LR"
data.hr.mr.lr <- rbind(data.hr, data.mr, data.lr)

data.hr.mr.lr.melt <- melt(data.hr.mr.lr, id = c("Setting"))

Setting  = c("HR", "MR", "LR")
data.SMQ <- data.frame(Setting = Setting)
data.SMQ$R_s <- NA
data.SMQ$R_q <- NA
data.SMQ$R_s.var <- NA
data.SMQ$R_q.var<- NA

for (i in 1:nrow(data.SMQ)){
  q <- mean(data.hr.mr.lr[data.hr.mr.lr$R_0_input < R_0_relevant.max &
                       data.hr.mr.lr$R_0_input > R_0_relevant.min &
                       data.hr.mr.lr$Setting == data.SMQ[i, "Setting"], "R_q"])
  q.var <- var(data.hr.mr.lr[data.hr.mr.lr$R_0_input < R_0_relevant.max &
                            data.hr.mr.lr$R_0_input > R_0_relevant.min &
                            data.hr.mr.lr$Setting == data.SMQ[i, "Setting"], "R_q"])
  sm <- mean(data.hr.mr.lr[data.hr.mr.lr$R_0_input < R_0_relevant.max &
                            data.hr.mr.lr$R_0_input > R_0_relevant.min &
                            data.hr.mr.lr$Setting == data.SMQ[i, "Setting"], "R_s"])
  sm.var <- var(data.hr.mr.lr[data.hr.mr.lr$R_0_input < R_0_relevant.max &
                             data.hr.mr.lr$R_0_input > R_0_relevant.min &
                             data.hr.mr.lr$Setting == data.SMQ[i, "Setting"], "R_s"])
  data.SMQ[i, "R_s"] <- sm
  data.SMQ[i, "R_s.var"] <- sm.var
  data.SMQ[i, "R_q"] <- q
  data.SMQ[i, "R_q.var"] <- q.var
}
data.SMQ.melt <- melt(data.SMQ, id = c("Setting", "R_q.var", "R_s.var"))
data.SMQ.melt$variance <- NA
data.SMQ.melt[data.SMQ.melt$variable == "R_q", "variance"] <- data.SMQ.melt[data.SMQ.melt$variable == "R_q", "R_q.var"]
data.SMQ.melt[data.SMQ.melt$variable == "R_s", "variance"] <- data.SMQ.melt[data.SMQ.melt$variable == "R_s", "R_s.var"]
data.SMQ.melt <- data.SMQ.melt[,c("Setting", "variable", "value", "variance")]

data.SMQ.melt$variable <- factor(data.SMQ.melt$variable, levels = c("R_q", "R_s"), ordered = TRUE)

#### Plot 1 ####
plot1 <- ggplot(data.SMQ.melt, aes(Setting, y = value, fill = variable)) +
  theme_bw() +
  theme(legend.justification=c(0,0), legend.position=c(0.1,0.85)) +
  geom_hline(yintercept=1, color = "lightgreen", size=2) +
  geom_bar(stat = "identity", position="dodge", width = 0.8, size = 1) +
  geom_errorbar(aes(ymin = value - 1.96*sqrt(variance), ymax = value + 1.96*sqrt(variance)), width=.2, position=position_dodge(.9)) +
  scale_fill_manual(values = c("cornflowerblue", "darkgoldenrod"), labels = c("Quarantine", "Symptom Monitoring")) +
  scale_x_discrete(limits = rev(c("LR", "MR", "HR")), labels = rev(c("Low Resource", "Mid Resource", "High Resource")), name = "") +
  ylab(expression(paste(R[e]))) +
  guides(fill = guide_legend(title = element_blank()))
plot1

pdf(file=paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", root, "_PlotSMQ1.pdf", sep=""), height = 7, width = 4)
plot(plot1)
dev.off()

#### Calculate RBQD and NQD ####
data.hr.mr.lr$NQD <- 1/data.hr.mr.lr$Abs_Benefit_per_Qday
data.hr.mr.lr$Rel_Benefit_per_Qday <- data.hr.mr.lr$Rel_Benefit / data.hr.mr.lr$obs_to_iso_q

summary(data.hr.mr.lr[data.hr.mr.lr$Setting == "HR", "NQD"])
summary(data.hr.mr.lr[data.hr.mr.lr$Setting == "MR", "NQD"])
summary(data.hr.mr.lr[data.hr.mr.lr$Setting == "LR", "NQD"])

summary(data.hr.mr.lr[data.hr.mr.lr$Setting == "HR", "Rel_Benefit_per_Qday"])
summary(data.hr.mr.lr[data.hr.mr.lr$Setting == "MR", "Rel_Benefit_per_Qday"])
summary(data.hr.mr.lr[data.hr.mr.lr$Setting == "LR", "Rel_Benefit_per_Qday"])

#### Calculate number of days in quarantine or symptom monitoring before isolation ####
# D1 indicates those people who truly are infected
D1_HR_q <- data.hr$obs_to_iso_q
D1_MR_q <- data.hr.mr.lr[data.hr.mr.lr$Setting == "MR", "obs_to_iso_q"]
D1_LR_q <- data.hr.mr.lr[data.hr.mr.lr$Setting == "LR", "obs_to_iso_q"]

D1_HR_s <- data.hr[data.hr$Setting == "HR", "obs_to_iso_s"]
D1_MR_s <- data.hr.mr.lr[data.hr.mr.lr$Setting == "MR", "obs_to_iso_s"]
D1_LR_s <- data.hr.mr.lr[data.hr.mr.lr$Setting == "LR", "obs_to_iso_s"]

# D0 are those people who truly aren't infected
D0 <- 14 # days for 2019-nCoV

#### Explore number of days in quarantine and symptom monitoring until isolationNOTE THAT THIS USES THE MEAN OBS_TO_ISO FROM EACH SIMULATION. NEED TO REPEAT USING INDIVIDUALS. ####
pdf(paste(root, "_Days_under_SM_or_Q.pdf", sep = ""))
layout(cbind(c(1, 2, 3), c(4,5,6)))
hist(data.hr$obs_to_iso_s, xlim = c(0, max(data.hr$obs_to_iso_s)), main = "High Feasibility Setting", xlab = "", breaks = 10, freq = FALSE)
hist(data.mr$obs_to_iso_s, xlim = c(0, max(data.hr$obs_to_iso_s)), main = "Medium Feasibility Setting", xlab = "", breaks = 10, freq = FALSE)
hist(data.lr$obs_to_iso_s, xlim = c(0, max(data.hr$obs_to_iso_s)), main = "Low Feasibility Setting", xlab = "Days under Active Monitoring", breaks = 20, freq = FALSE)

hist(data.hr$obs_to_iso_q, xlim = c(0, max(data.hr$obs_to_iso_s)), main = "High Feasibility Setting", xlab = "", breaks = 10, freq = FALSE)
hist(data.mr$obs_to_iso_q, xlim = c(0, max(data.hr$obs_to_iso_s)), main = "Medium Feasibility Setting", xlab = "", breaks = 10, freq = FALSE)
hist(data.lr$obs_to_iso_q, xlim = c(0, max(data.hr$obs_to_iso_s)), main = "Low Feasibility Setting", xlab = "Days under Quarantine", breaks = 10, freq = FALSE)

dev.off()



#### HR: obs_to_iso using individuals rather than simulation means ####

data.hr.individuals.qdays <- bind_rows(data.hr.individuals)

# OPTIONAL - reduce the size of the dataset
# rows_q <- which(data.hr.individuals.qdays$intervention == "q")
# rows_s <- which(data.hr.individuals.qdays$intervention == "s")
# data.hr.individuals.qdays <- rbind(data.hr.individuals.qdays[sample(rows_q, 100000),], data.hr.individuals.qdays[sample(rows_s, 100000),])

# Remove individuals were "traced" at the end of their disease cycle
data.hr.individuals.qdays$keep <- NA
data.hr.individuals.qdays$keep <- apply(data.hr.individuals.qdays, 1, function(x) (as.numeric(x["t_obs"]) != max(as.numeric(x["T_inc"]) + as.numeric(x["d_symp"]), as.numeric(x["T_lat"]) + as.numeric(x["d_inf"]))))


# Remove individuals who were infected by people that weren't traced until the end of their disease cycle
data.hr.individuals.qdays$keep_2 <- NA

nrow(data.hr.individuals.qdays)

for (sim in 1:max(data.hr.individuals.qdays$simulation)){
  cat(".")
  temp <- data.hr.individuals.qdays[data.hr.individuals.qdays$simulation == sim & data.hr.individuals.qdays$intervention == "s",]
  data.hr.individuals.qdays[data.hr.individuals.qdays$simulation == sim & data.hr.individuals.qdays$intervention == "s" & data.hr.individuals.qdays$generation > 1,"keep_2"] <-
    unlist(apply(temp, 1, function(x) temp[temp$generation == (as.numeric(x["generation"]) - 1) & temp$ID == as.numeric(x["infector"]),"keep"]))
  
  temp <- data.hr.individuals.qdays[data.hr.individuals.qdays$simulation == sim & data.hr.individuals.qdays$intervention == "q",]
  data.hr.individuals.qdays[data.hr.individuals.qdays$simulation == sim & data.hr.individuals.qdays$intervention == "q" & data.hr.individuals.qdays$generation > 1,"keep_2"] <-
    unlist(apply(temp, 1, function(x) temp[temp$generation == (as.numeric(x["generation"]) - 1) & temp$ID == as.numeric(x["infector"]),"keep"]))
}

# Remove first two generations where intervention isn't settled in yet
data.hr.individuals.qdays <- data.hr.individuals.qdays[data.hr.individuals.qdays$generation > 2,] 
nrow(data.hr.individuals.qdays[data.hr.individuals.qdays$keep == FALSE,]) / nrow(data.hr.individuals.qdays) # Should be around the 1 - prob_CT
nrow(data.hr.individuals.qdays[data.hr.individuals.qdays$keep_2 == FALSE,]) / nrow(data.hr.individuals.qdays) 

data.hr.individuals.qdays <- data.hr.individuals.qdays[data.hr.individuals.qdays$keep == TRUE,]
data.hr.individuals.qdays <- data.hr.individuals.qdays[data.hr.individuals.qdays$keep_2 == TRUE,]

# Define new variables
data.hr.individuals.qdays$T_lat_offset <- data.hr.individuals.qdays$T_lat - data.hr.individuals.qdays$T_inc
data.hr.individuals.qdays$obs_to_onset <- (data.hr.individuals.qdays$T_inc - data.hr.individuals.qdays$t_obs)/24

data.hr.individuals.qdays$intervention_days <- (data.hr.individuals.qdays$t_iso - data.hr.individuals.qdays$t_obs)/24

# Remove the epsilon effect
data.hr.individuals.qdays$intervention_days_s_adj <- data.hr.individuals.qdays$intervention_days
data.hr.individuals.qdays[data.hr.individuals.qdays$intervention == "s",]$intervention_days_s_adj <- data.hr.individuals.qdays[data.hr.individuals.qdays$intervention == "s",]$intervention_days - data.hr.individuals.qdays[data.hr.individuals.qdays$intervention == "s",]$epsilon/24


figure_qdays_individuals <- ggplot(data.hr.individuals.qdays, aes(fill = intervention, x = intervention_days, alpha = intervention, stat(density))) +
  geom_histogram(position = "identity") +
  # facet_grid(generation ~ .) +
  scale_alpha_manual(values = c( q = 0.3, s = 0.7), name = "Intervention", labels = c("Quarantine", "Active Monitoring" )) +
  scale_fill_manual(values = c( q = "blue", s = "yellow"), name = "Intervention", labels = c("Quarantine", "Active Monitoring")) +
  scale_x_continuous(limits = c(0, 14), name = "Days under Intervention", breaks = seq(0, 14, 2)) +
  ylab("Density") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  ggtitle("Mean Serial Interval = 7.5 days") +
  theme(legend.position = "bottom")
figure_qdays_individuals

pdf(file=paste(root, "_Plot_Days_under_Intervention_Individuals.pdf", sep=""), width = 4, height = 4)
plot(figure_qdays_individuals)
dev.off()

summary(data.hr.individuals.qdays[data.hr.individuals.qdays$intervention == "q", "intervention_days"])
nrow(data.hr.individuals.qdays[data.hr.individuals.qdays$intervention == "q" & data.hr.individuals.qdays$intervention_days == 0,]) / nrow(data.hr.individuals.qdays[data.hr.individuals.qdays$intervention == "q" ,]) # 8.4% of traced contacts are symptomatic at time of tracing for 4.8 SI

summary(data.hr.individuals.qdays[data.hr.individuals.qdays$intervention == "s", "intervention_days"])
nrow(data.hr.individuals.qdays[data.hr.individuals.qdays$intervention == "s" & data.hr.individuals.qdays$intervention_days == 0,]) / nrow(data.hr.individuals.qdays[data.hr.individuals.qdays$intervention == "s" ,]) # 3.3% of traced contacts are symptomatic at time of tracing for 4.8 SI

ggplot(data.hr.individuals.qdays, aes(fill = intervention, x = (t_obs - t_infection)/24, alpha = intervention, stat(density))) +
  geom_histogram(position = "identity") +
  # facet_grid(Setting ~ .) +
  scale_alpha_manual(values = c( q = 0.3, s = 0.7), name = "Intervention", labels = c("Quarantine", "Active Monitoring" )) +
  scale_fill_manual(values = c( q = "blue", s = "yellow"), name = "Intervention", labels = c("Quarantine", "Active Monitoring")) +
  xlab("Days from Infection to Tracing") +
  ylab("Density") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  ggtitle("Mean Serial Interval = 4.8 days") +
  theme(legend.position = "bottom")

#### LR: obs_to_iso using individuals rather than simulation means ####

data.lr.individuals.qdays <- bind_rows(data.lr.individuals)

# OPTIONAL - reduce the size of the dataset
# rows_q <- which(data.lr.individuals.qdays$intervention == "q")
# rows_s <- which(data.lr.individuals.qdays$intervention == "s")
# data.lr.individuals.qdays <- rbind(data.lr.individuals.qdays[sample(rows_q, 100000),], data.lr.individuals.qdays[sample(rows_s, 100000),])

# Remove individuals were "traced" at the end of their disease cycle
data.lr.individuals.qdays$keep <- NA
data.lr.individuals.qdays$keep <- apply(data.lr.individuals.qdays, 1, function(x) (as.numeric(x["t_obs"]) != max(as.numeric(x["T_inc"]) + as.numeric(x["d_symp"]), as.numeric(x["T_lat"]) + as.numeric(x["d_inf"]))))


# Remove individuals who were infected by people that weren't traced until the end of their disease cycle
data.lr.individuals.qdays$keep_2 <- NA

nrow(data.lr.individuals.qdays)

for (sim in 1:max(data.lr.individuals.qdays$simulation)){
  cat(".")
  temp <- data.lr.individuals.qdays[data.lr.individuals.qdays$simulation == sim & data.lr.individuals.qdays$intervention == "s",]
  data.lr.individuals.qdays[data.lr.individuals.qdays$simulation == sim & data.lr.individuals.qdays$intervention == "s" & data.lr.individuals.qdays$generation > 1,"keep_2"] <-
    unlist(apply(temp, 1, function(x) temp[temp$generation == (as.numeric(x["generation"]) - 1) & temp$ID == as.numeric(x["infector"]),"keep"]))
  
  temp <- data.lr.individuals.qdays[data.lr.individuals.qdays$simulation == sim & data.lr.individuals.qdays$intervention == "q",]
  data.lr.individuals.qdays[data.lr.individuals.qdays$simulation == sim & data.lr.individuals.qdays$intervention == "q" & data.lr.individuals.qdays$generation > 1,"keep_2"] <-
    unlist(apply(temp, 1, function(x) temp[temp$generation == (as.numeric(x["generation"]) - 1) & temp$ID == as.numeric(x["infector"]),"keep"]))
}

# Remove first two generations where intervention isn't settled in yet
data.lr.individuals.qdays <- data.lr.individuals.qdays[data.lr.individuals.qdays$generation > 2,] 
nrow(data.lr.individuals.qdays[data.lr.individuals.qdays$keep == FALSE,]) / nrow(data.lr.individuals.qdays) # Should be around the 1 - prob_CT
nrow(data.lr.individuals.qdays[data.lr.individuals.qdays$keep_2 == FALSE,]) / nrow(data.lr.individuals.qdays) 

data.lr.individuals.qdays <- data.lr.individuals.qdays[data.lr.individuals.qdays$keep == TRUE,]
data.lr.individuals.qdays <- data.lr.individuals.qdays[data.lr.individuals.qdays$keep_2 == TRUE,]

# Define new variables
data.lr.individuals.qdays$T_lat_offset <- data.lr.individuals.qdays$T_lat - data.lr.individuals.qdays$T_inc
data.lr.individuals.qdays$obs_to_onset <- (data.lr.individuals.qdays$T_inc - data.lr.individuals.qdays$t_obs)/24

data.lr.individuals.qdays$intervention_days <- (data.lr.individuals.qdays$t_iso - data.lr.individuals.qdays$t_obs)/24

# Remove the epsilon effect
data.lr.individuals.qdays$intervention_days_s_adj <- data.lr.individuals.qdays$intervention_days
data.lr.individuals.qdays[data.lr.individuals.qdays$intervention == "s",]$intervention_days_s_adj <- data.lr.individuals.qdays[data.lr.individuals.qdays$intervention == "s",]$intervention_days - data.lr.individuals.qdays[data.lr.individuals.qdays$intervention == "s",]$epsilon/24


figure_qdays_individuals_lr <- ggplot(data.lr.individuals.qdays, aes(fill = intervention, x = intervention_days, alpha = intervention, stat(density))) +
  geom_histogram(position = "identity") +
  # facet_grid(generation ~ .) +
  scale_alpha_manual(values = c( q = 0.3, s = 0.7), name = "Intervention", labels = c("Quarantine", "Active Monitoring" )) +
  scale_fill_manual(values = c( q = "blue", s = "yellow"), name = "Intervention", labels = c("Quarantine", "Active Monitoring")) +
  scale_x_continuous(limits = c(0, 14), name = "Days under Intervention", breaks = seq(0, 14, 2)) +
  ylab("Density") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  ggtitle("Mean Serial Interval = 7.5 days") +
  theme(legend.position = "bottom")
figure_qdays_individuals_lr

pdf(file=paste(root, "_Plot_Days_under_Intervention_Individuals_lr.pdf", sep=""), width = 4, height = 4)
plot(figure_qdays_individuals_lr)
dev.off()

summary(data.lr.individuals.qdays[data.lr.individuals.qdays$intervention == "q", "intervention_days"])
nrow(data.lr.individuals.qdays[data.lr.individuals.qdays$intervention == "q" & data.lr.individuals.qdays$intervention_days == 0,]) / nrow(data.lr.individuals.qdays[data.lr.individuals.qdays$intervention == "q" ,]) # 8.4% of traced contacts are symptomatic at time of tracing for 4.8 SI

summary(data.lr.individuals.qdays[data.lr.individuals.qdays$intervention == "s", "intervention_days"])
nrow(data.lr.individuals.qdays[data.lr.individuals.qdays$intervention == "s" & data.lr.individuals.qdays$intervention_days == 0,]) / nrow(data.lr.individuals.qdays[data.lr.individuals.qdays$intervention == "s" ,]) # 3.3% of traced contacts are symptomatic at time of tracing for 4.8 SI

ggplot(data.lr.individuals.qdays, aes(fill = intervention, x = (t_obs - t_infection)/24, alpha = intervention, stat(density))) +
  geom_histogram(position = "identity") +
  # facet_grid(Setting ~ .) +
  scale_alpha_manual(values = c( q = 0.3, s = 0.7), name = "Intervention", labels = c("Quarantine", "Active Monitoring" )) +
  scale_fill_manual(values = c( q = "blue", s = "yellow"), name = "Intervention", labels = c("Quarantine", "Active Monitoring")) +
  xlab("Days from Infection to Tracing") +
  ylab("Density") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  ggtitle("Mean Serial Interval = 4.8 days") +
  theme(legend.position = "bottom")


#### Plot NQD as fraction of contacts traced truly infected decreases ####
NNQ_scenario1 <- 23.41098
Qdays_scenario1 <- 3.551424

NNQ_scenario2 <- 1.07
Qdays_scenario2 <- 3.164713

days_noninfected <- 14

prob_infected_vector <- 1/c(1e6, 1e5, 1e4, 1e3, 1e2, 10, 1)

df_NQD <- data.frame(prob_infected = rep(prob_infected_vector, times = 2), NNQ = NA, NQD = NA, scenario = NA)
df_NQD$scenario <- rep(c(1,2), each = length(prob_infected_vector))

for (i in 1:nrow(df_NQD)){
  if (df_NQD[i, "scenario"] == 1){
    NNQ <- NNQ_scenario1
    Qdays <- Qdays_scenario1
  } else {
    NNQ <- NNQ_scenario2
    Qdays <- Qdays_scenario2
  }
  
  df_NQD[i, "NNQ"] <- NNQ * 1/(df_NQD[i, "prob_infected"])
  df_NQD[i, "NQD"] <- NNQ * Qdays + ((1/df_NQD[i, "prob_infected"]) - 1) * NNQ * days_noninfected
}

ggplot(df_NQD, aes(x = prob_infected, y = NNQ, color = scenario)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  xlab()


ggplot(df_NQD, aes(x = prob_infected, y = NQD, color = scenario)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10()

#### If using two resource scenarios ####
data.hr$Setting <- "HR"
data.lr$Setting <- "LR"
data.hr.lr <- rbind(data.hr, data.lr)

data.hr.lr.melt <- melt(data.hr.lr[,c("Setting", "obs_to_iso_s", "obs_to_iso_q")], id.vars = c("Setting"))

data.hr.lr.melt[data.hr.lr.melt$Setting == "HR","Setting"] <- "High Feasiblity Setting"
data.hr.lr.melt[data.hr.lr.melt$Setting == "LR","Setting"] <- "Low Feasiblity Setting"


figure_q_days <- ggplot(data.hr.lr.melt, aes(fill = variable, alpha = variable, x = value)) +
  # geom_histogram(aes(color = "white"), position = "identity") +
  geom_histogram(position = "identity") +
  facet_grid(Setting ~ .) +
  scale_alpha_manual(values = c("obs_to_iso_s" = 0.7, "obs_to_iso_q" = 0.3), name = "Intervention", labels = c("Active Monitoring" , "Quarantine")) +
  scale_fill_manual(values = c("obs_to_iso_s" = "yellow", "obs_to_iso_q" = "blue"), name = "Intervention", labels = c("Active Monitoring" , "Quarantine")) +
  xlab("Days under Intervention") +
  ylab("Count") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  ggtitle("Mean Serial Interval = 7.5 days") +
  theme(legend.position = "bottom")
figure_q_days

pdf(file=paste(root, "_Plot_Days_under_Intervention.pdf", sep=""), width = 4, height = 4)
plot(figure_q_days)
dev.off()

#### If using two diseases, one resource scenario ####


data.scenario1$Setting <- "Scenario 1"
data.scenario2$Setting <- "Scenario 2"
data.scenario1_2 <- rbind(data.scenario1, data.scenario2)

data.scenario1_2.ReRq.melt <- melt(data.scenario1_2[,c("R_s", "R_q", "R_0_input", "Setting")], id.vars = c("R_s", "Setting", "R_0_input"))


figure_q_days <- ggplot(data.hr.lr.melt, aes(fill = variable, alpha = variable, x = value)) +
  # geom_histogram(aes(color = "white"), position = "identity") +
  geom_histogram(position = "identity") +
  facet_grid(Setting ~ .) +
  scale_alpha_manual(values = c("obs_to_iso_s" = 0.7, "obs_to_iso_q" = 0.3), name = "Intervention", labels = c("Active Monitoring" , "Quarantine")) +
  scale_fill_manual(values = c("obs_to_iso_s" = "yellow", "obs_to_iso_q" = "blue"), name = "Intervention", labels = c("Active Monitoring" , "Quarantine")) +
  xlab("Days under Intervention") +
  ylab("Count") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  ggtitle("Mean Serial Interval = 4.8 days") +
  theme(legend.position = "bottom")
figure_q_days

pdf(file=paste(root, "_Plot_Days_under_Intervention.pdf", sep=""), width = 4, height = 4)
plot(figure_q_days)
dev.off()

#### Calculate NQD50 ####
# NQDX = [ (X*D1) + ((1-X)*D0) ] / (Rs – Rq)

x <- 0.5
data.hr.mr.lr[data.hr.mr.lr$Setting == "HR", "NQD50"] <- (D1_HR_q + (1/x - 1)*(D0)) / (data.hr.mr.lr[data.hr.mr.lr$Setting == "HR", "Abs_Benefit"])
data.hr.mr.lr[data.hr.mr.lr$Setting == "MR", "NQD50"] <- (D1_MR_q + (1/x - 1)*(D0)) / (data.hr.mr.lr[data.hr.mr.lr$Setting == "MR", "Abs_Benefit"])
data.hr.mr.lr[data.hr.mr.lr$Setting == "LR", "NQD50"] <- (D1_LR_q + (1/x - 1)*(D0)) / (data.hr.mr.lr[data.hr.mr.lr$Setting == "LR", "Abs_Benefit"])

median(data.hr.mr.lr[data.hr.mr.lr$Setting == "HR", "NQD50"])
median(data.hr.mr.lr[data.hr.mr.lr$Setting == "MR", "NQD50"])
median(data.hr.mr.lr[data.hr.mr.lr$Setting == "LR", "NQD50"])

#### Calculate RBQD50 ####
# Rel_Benefit_per_QdayX = [ (Rs – Rq)/Rs ] / [ (X*D1) + ((1-X)*D0) ]

data.hr.mr.lr[data.hr.mr.lr$Setting == "HR", "RBQD50"] <- data.hr.mr.lr[data.hr.mr.lr$Setting == "HR", "Rel_Benefit"] / (D1_HR_q + (1/x - 1)*(D0))
data.hr.mr.lr[data.hr.mr.lr$Setting == "MR", "RBQD50"] <- data.hr.mr.lr[data.hr.mr.lr$Setting == "MR", "Rel_Benefit"] / (D1_MR_q + (1/x - 1)*(D0))
data.hr.mr.lr[data.hr.mr.lr$Setting == "LR", "RBQD50"] <- data.hr.mr.lr[data.hr.mr.lr$Setting == "LR", "Rel_Benefit"] / (D1_LR_q + (1/x - 1)*(D0))

median(data.hr.mr.lr[data.hr.mr.lr$Setting == "HR", "RBQD50"])
median(data.hr.mr.lr[data.hr.mr.lr$Setting == "MR", "RBQD50"])
median(data.hr.mr.lr[data.hr.mr.lr$Setting == "LR", "RBQD50"])


#### Explore Data ####
layout(cbind(c(1, 2, 3)))
hist(data.hr.mr.lr[data.hr.mr.lr$Setting == "HR", "Rel_Benefit"], xlim = c(0, max(data.hr.mr.lr$Rel_Benefit)), main = "Relative Benefit in HR setting")
hist(data.hr.mr.lr[data.hr.mr.lr$Setting == "MR", "Rel_Benefit"], xlim = c(0, max(data.hr.mr.lr$Rel_Benefit)), main = "Relative Benefit in MR setting")
hist(data.hr.mr.lr[data.hr.mr.lr$Setting == "LR", "Rel_Benefit"], xlim = c(0, max(data.hr.mr.lr$Rel_Benefit)), main = "Relative Benefit in LR setting")


#### Save Workspace ####
save.image(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", root, "_SMQDecision.RData", sep=""))

