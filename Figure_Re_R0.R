#### Header ####
# Code for the R_e vs R_0 plots

#### Load Libraries ####
library(ggplot2)
library(reshape)

#### Load HR and LR workspaces ####
desired_root <- "20200205_2019nCoV" # Paste the desired root here "YYYYMMDD_DISEASE"

# If workspaces are in main folder
# load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "_HR.RData", sep=""))
# load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "_LR.RData", sep=""))
# load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "_Plots.RData", sep=""))

# If workspaces are in their own folder, named the same as the root
load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "/", desired_root, "_HR.RData", sep=""))
load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "/", desired_root, "_LR.RData", sep=""))
load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "/", desired_root, "_Plots.RData", sep=""))

#### Load two disease scenarios ####

root_1 <- "20200205_2019nCoV"
# load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", root_1, "_HR.RData", sep=""))
data.scenario1.hr <- data.hr
data.scenario1.lr <- data.lr

root_2 <- "20200217_2019nCoV_SI4.8"
# load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", root_2, "_HR.RData", sep=""))
data.scenario2.hr <- data.hr
data.scenario2.lr <- data.lr

#### Plot for two diseases ####
data.scenario1.hr$Setting <- "HR"
data.scenario1.lr$Setting <- "LR"

data.scenario2.hr$Setting <- "HR"
data.scenario2.lr$Setting <- "LR"

data.scenario1.hr$Serial_Interval <- "Mean 7.5 Days"
data.scenario1.lr$Serial_Interval <- "Mean 7.5 Days"

data.scenario2.hr$Serial_Interval <- "Mean 4.8 Days"
data.scenario2.lr$Serial_Interval <- "Mean 4.8 Days"

data.scenario1.2 <- rbind(data.scenario1.hr, data.scenario2.hr)

data.scenario1.2.melt <- melt(data.scenario1.2[,c("R_s", "R_q", "R_0_input", "Setting", "Serial_Interval")], id.vars = c("R_s", "Setting", "R_0_input", "Serial_Interval"))

plot1 <- ggplot(data = data.scenario1.2.melt, aes(x = R_s, y = value, shape = Serial_Interval, color = R_0_input)) +
  geom_vline(xintercept=1, col="grey") + 
  geom_hline(yintercept=1, col="grey") +
  annotate("rect", xmin = 0, xmax = 1, ymin = 1, ymax = 4, alpha = .1, fill = "yellow") + 
  annotate("rect", xmin = 1, xmax = 4, ymin = 0, ymax = 1, alpha = .1, fill = "blue") +
  annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1, alpha = .1, fill = "green") +
  annotate("text", x = 2.5, y = 0.1, label = "Control with Quarantine", col = "blue") +
  annotate("text", x = 0.5, y = 3.5, label = "Control with\nActive Monitoring", col = "orange") +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, col = "grey") +
  xlim(0,4) + ylim(0,4) +
  scale_color_gradient(low="yellow", high="darkred", name = "Basic\nReproductive\nNumber") +
  scale_shape_manual(name = "Setting",
                     values = c(17, 1),
                     labels = c("Mean 7.5 Days", "Mean 4.8 Days")) +
  guides(shape = guide_legend(reverse=TRUE)) +
  xlab("Effective Reproductive Number under Active Monitoring") +
  ylab("Effective Reproductive Number under Quarantine") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle("High Feasibility Setting")
plot1

pdf(file=paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", root, "_Plot1_SIscenarios.pdf", sep=""), width = 8, height = 6)
plot(plot1)
dev.off()


#### Plot 1: two diseases, one resource setting ####

data.hr$Setting <- "HR"
data.lr$Setting <- "LR"
data.hr.lr.ReRq <- rbind(data.hr, data.lr)

data.hr.lr.ReRq.melt <- melt(data.hr.lr.ReRq[,c("R_s", "R_q", "R_0_input", "Setting")], id.vars = c("R_s", "Setting", "R_0_input"))

plot1 <- ggplot(data = data.hr.lr.ReRq.melt, aes(x = R_s, y = value, shape = Setting, color = R_0_input)) +
  geom_vline(xintercept=1, col="grey") + 
  geom_hline(yintercept=1, col="grey") +
  annotate("rect", xmin = 0, xmax = 1, ymin = 1, ymax = 4, alpha = .1, fill = "yellow") + 
  annotate("rect", xmin = 1, xmax = 4, ymin = 0, ymax = 1, alpha = .1, fill = "blue") +
  annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1, alpha = .1, fill = "green") +
  annotate("text", x = 2.5, y = 0.1, label = "Control with Quarantine", col = "blue") +
  annotate("text", x = 0.5, y = 3.5, label = "Control with\nActive Monitoring", col = "orange") +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, col = "grey") +
  xlim(0,4) + ylim(0,4) +
  scale_color_gradient(low="yellow", high="darkred", name = "Basic\nReproductive\nNumber") +
  scale_shape_manual(name = "Setting",
                     values = c(17, 1),
                     labels = c("High Feasibility", "Low Feasibility")) +
  guides(shape = guide_legend(reverse=TRUE)) +
  xlab("Effective Reproductive Number under Active Monitoring") +
  ylab("Effective Reproductive Number under Quarantine") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle("Mean Serial Interval = 4.8 days")
plot1

pdf(file=paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", root, "_Plot1.pdf", sep=""), width = 8, height = 6)
plot(plot1)
dev.off()

#### Plot 2 ####
plot2 <- ggplot(data = data.hr.lr) +
  annotate("rect", xmin = 0, xmax = 4, ymin = 0, ymax = 1, alpha = .1, fill = "green") +
  geom_hline(yintercept=1, col = "grey") +
  geom_point(aes(x=R_0_input, y=R_s, shape = Setting), col = "darkgreen", alpha = 0.7) +
  geom_point(aes(x=R_0_input, y=R_q, shape = Setting), col = "blue", alpha = 0.7) +
  stat_smooth(aes(x=R_0_input, y=R_s, shape = Setting), method = "loess", color="darkgreen", size = 1.2) +
  stat_smooth(aes(x=R_0_input, y=R_q, shape = Setting), method = "loess", color="blue", size = 1.2) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_shape_manual(name = "Setting",
                     values = c(17, 1),
                     labels = c("High Resource", "Low Resource")) +
  guides(shape = guide_legend(reverse=TRUE)) +
  xlab(expression("Basic Reproductive Number R" [0])) + ylab(expression("Effective Reproductive Number R" [e])) +
  ggtitle(paste("Disease: ", disease))
plot2

pdf(file=paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", root, "_Plot2.pdf", sep=""))
plot(plot2)
dev.off()

#### Plot 3 ####
point_size = 0.1
line_size = 1
plot3 <- ggplot(data = data.hr.lr[data.hr.lr$Setting == "HR",]) +
  geom_point(aes(x=R_0_input, y=R_hsb), col = "mediumturquoise", alpha = 0.5, size = point_size) +
  geom_point(aes(x=R_0_input, y=R_hsb), col = "mediumturquoise", alpha = 0.5, size = point_size) +
  geom_point(aes(x=R_0_input, y=R_s), col = "darkgoldenrod", alpha = 0.5, size = point_size) +
  geom_point(aes(x=R_0_input, y=R_q), col = "cornflowerblue", alpha = 0.5, size = point_size) +
  stat_smooth(aes(x=R_0_input, y=R_hsb), method = "loess", color="mediumturquoise", size = line_size, se=FALSE, alpha = 0.5) +
  stat_smooth(aes(x=R_0_input, y=R_s), method = "loess", color="darkgoldenrod", size = line_size, se=FALSE, alpha = 0.5) +
  stat_smooth(aes(x=R_0_input, y=R_q), method = "loess", color="cornflowerblue", size = line_size, se=FALSE, alpha = 0.5) +
  geom_hline(yintercept=1, col = "grey", size = 0.5) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab(expression("R" [0])) + ylab(expression("R" [e])) +
  geom_text(data = NULL, x = 1.46, y = 0, label = "|", size = 2) +
  # geom_text(data = NULL, x = 1.83, y = 0, label = "|", size = 2) +
  # ggtitle(paste("Influenza A")) +
  # ggtitle(paste("Ebola")) +
  theme(text = element_text(size=8))
plot3

pdf(file=paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", root, "_Plot3.pdf", sep=""), width = 2, height = 2)
plot(plot3)
dev.off()


#### Plot 3 FOR POSTER ####
plot3_poster <- ggplot(data = data.hr.lr[data.hr.lr$R_0 > R_0_relevant.min &
                                    data.hr.lr$R_0 < R_0_relevant.max &
                                    data.hr.lr$Setting == "HR",]) +
  annotate("rect", xmin = R_0_relevant.min, xmax = R_0_relevant.max, ymin = 0, ymax = 1, alpha = .1, fill = "green") +
  geom_hline(yintercept=1, col = "grey") +
  geom_point(aes(x=R_0, y=R_hsb), col = "mediumturquoise", alpha = 0.5, size = 0.5) +
  geom_point(aes(x=R_0, y=R_hsb), col = "mediumturquoise", alpha = 0.5, size = 0.5) +
  geom_point(aes(x=R_0, y=R_s), col = "darkgoldenrod", alpha = 0.5, size = 0.5) +
  geom_point(aes(x=R_0, y=R_q), col = "cornflowerblue", alpha = 0.5, size = 0.5) +
  stat_smooth(aes(x=R_0, y=R_hsb), method = "loess", color="mediumturquoise", size = 1, se=FALSE) +
  stat_smooth(aes(x=R_0, y=R_s), method = "loess", color="darkgoldenrod", size = 1, se=FALSE) +
  stat_smooth(aes(x=R_0, y=R_q), method = "loess", color="cornflowerblue", size = 1, se=FALSE) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab(expression("Basic Reproductive Number R" [0])) + ylab(expression("Effective Reproductive Number R" [e])) +
  geom_text(data = NULL, x = 1.46, y = 0, label = "*") +
  # geom_text(data = NULL, x = 1.83, y = 0, label = "*") +
  theme(text = element_text(size=18)) +
  ggtitle(paste("Influenza A"))
# ggtitle(paste("Ebola"))
plot3_poster

pdf(file=paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", root, "_Plot3_5by5.pdf", sep=""), width = 5, height = 5)
plot(plot3_poster)
dev.off()




#### Save Workspace ####
save.image(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", root, "_Plots.RData", sep=""))
