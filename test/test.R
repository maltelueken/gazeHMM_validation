### gazeHMM test ###

# Author: Malte Lüken
# Date: 27.04.2020


library(tidyverse)
library(depmixS4)
source("algorithm/gazeHMM.R")
source("algorithm/plotting_functions.R")


# Test on static data -----------------------------------------------------

# Eprep data set

load("data/Eprep")

dim <- c(473, 296) # screen dimensions
res <- c(1680, 1050) # screen resolution
dist <- 600 # distance to screen
fr <- 500 # frame rate

set.seed(123)

Eprep.fit.3 <- Eprep %>%
  dplyr::filter(block == sample(unique(Eprep$block), 1)) %>%
  group_by(subject) %>%
  summarise(fit = list(try(gazeHMM(x = xp, y = yp, t = (time - time[1])/1e3, 
                                              unit = "px", res = res, dim = dim, fr = fr, dist = dist, 
                                              nstates = as.integer(3), random.respstart = F)))) %>%
  pull(fit)


# Inspect parameters

lapply(Eprep.fit.3, function(x) {summary(x$model)})


# Plot classification result

for (i in 1:length(Eprep.fit.3)) {
  
  set.seed(i*25)
  
  tp <- sample(Eprep.fit.3[[i]]$samples$t[Eprep.fit.3[[i]]$samples$t > 0.5], 1)
  
  int <- c(tp - 0.5, tp) # 0.5 second interval
  
  plot_samples_time(Eprep.fit.3[[i]], from = int[1], to = int[2]) + 
    ggsave(filename = paste("Eprep_3comp_subj_", i, "_overtime.png", sep = ""), 
           width = 12, height = 6)
  
  plot_samples_xy(Eprep.fit.3[[i]], from = int[1], to = int[2]) + 
    ggsave(filename = paste("Eprep_3comp_subj_", i, "_xy.png", sep = ""), 
           width = 12, height = 6)
  
}


# Plot bivariate classification 

for (i in 1:length(Eprep.fit.3)) {
  
  ggplot(Eprep.fit.3[[i]]$samples, aes(vel, acc, color = as.factor(state))) + 
    geom_point() + scale_color_discrete(name = "state", labels = c("Fix", "Sac", "PSO")) +
    ggsave(filename = paste("Eprep_3comp_subj_", i, "_bivariate_vel_acc.png", sep = ""), 
           width = 6, height = 6)
  
  ggplot(Eprep.fit.3[[i]]$samples, aes(vel, angle, color = as.factor(state))) + 
    geom_point() + scale_color_discrete(name = "state", labels = c("Fix", "Sac", "PSO")) +
    ggsave(filename = paste("Eprep_3comp_subj_", i, "_bivariate_vel_angle.png", sep = ""), 
           width = 6, height = 6)
  
  ggplot(Eprep.fit.3[[i]]$samples, aes(acc, angle, color = as.factor(state))) + 
    geom_point() + scale_color_discrete(name = "state", labels = c("Fix", "Sac", "PSO")) +
    ggsave(filename = paste("Eprep_3comp_subj_", i, "_bivariate_acc_angle.png", sep = ""), 
           width = 6, height = 6)
  
}


# Test on dynamic data ----------------------------------------------------

# D2010 data set

load("data/D2010")

dim <- c(400, 300) # screen dimension
res <- c(1280, 960) # screen resolution
dist <- 450 # distance to screen
fr <- 250 # frame rate

set.seed(123)

D2010.fit.4 <- D2010 %>%
  dplyr::filter(subject %in% sample(unique(D2010$subject), 10, replace = F), 
                trial == sample(unique(D2010$trial), 1), 
                condition == sample(unique(D2010$condition), 1)) %>%
  group_by(subject) %>%
  summarise(fit = list(try(gazeHMM(x = x, y = y, t = (timestamp - timestamp[1])/1e6, 
                                          unit = "px", res = res, dim = dim, fr = fr, dist = dist,
                                          nstates = as.integer(4), random.respstart = F)))) %>%
  pull(fit)


# Inspect parameters

lapply(D2010.fit.4, function(x) {summary(x$model)})


# Plot classification results

for (i in 1:length(D2010.fit.4)) {
  
  set.seed(i*25)
  
  tp <- sample(D2010.fit.4[[i]]$samples$t[D2010.fit.4[[i]]$samples$t > 1], 1)
  
  int <- c(tp - 1, tp) # 1 second interval
  
  plot_samples_time(D2010.fit.4[[i]], from = int[1], to = int[2]) + 
    ggsave(filename = paste("D2010_4comp_subj_", i, "_overtime.png", sep = ""), 
           width = 12, height = 6)
  
  plot_samples_xy(D2010.fit.4[[i]], from = int[1], to = int[2]) + 
    ggsave(filename = paste("D2010_4comp_subj_", i, "_xy.png", sep = ""), 
           width = 12, height = 6)
}


# Plot bivariate classification

for (i in 1:length(D2010.fit.4)) {
  
  ggplot(D2010.fit.4[[i]]$samples, aes(vel, acc, color = as.factor(state))) + 
    geom_point() + scale_color_discrete(name = "state", labels = c("Fix", "Sac", "PSO")) +
    ggsave(filename = paste("D2010_4comp_subj_", i, "_bivariate_vel_acc.png", sep = ""), 
           width = 6, height = 6)
  
  ggplot(D2010.fit.4[[i]]$samples, aes(vel, angle, color = as.factor(state))) + 
    geom_point() + scale_color_discrete(name = "state", labels = c("Fix", "Sac", "PSO")) +
    ggsave(filename = paste("D2010_4comp_subj_", i, "_bivariate_vel_angle.png", sep = ""), 
           width = 6, height = 6)
  
  ggplot(D2010.fit.4[[i]]$samples, aes(acc, angle, color = as.factor(state))) + 
    geom_point() + scale_color_discrete(name = "state", labels = c("Fix", "Sac", "PSO")) +
    ggsave(filename = paste("D2010_4comp_subj_", i, "_bivariate_acc_angle.png", sep = ""), 
           width = 6, height = 6)
  
}
