library(tidyverse)
library(here)
library(patchwork)
theme_set(theme_classic(base_size = 12.5))

load(here("validation", "Andersson2017_raw.Rdata"))
load(here("validation", "Andersson2017_fitted.Rdata"))

merged <- bind_cols(A2017.fit[[3]][[3]][[4]]$samples, 
                coderMN = A2017[[3]][[3]]$coderMN,
                coderRA = A2017[[3]][[3]]$coderRA)

rm(A2017, A2017.fit)

df <- merged %>%
  mutate(
    gazeHMM = case_when(label == 0 ~ "Noise",
                        label == 1 ~ "Fixation",
                        label == 2 ~ "Saccade",
                        label == 3 ~ "PSO",
                        label == 4 ~ "Smooth pursuit"),
    coderMN = case_when(coderMN == 1 ~ "Fixation",
                        coderMN == 2 ~ "Saccade",
                        coderMN == 3 ~ "PSO",
                        coderMN == 4 ~ "Smooth pursuit",
                        coderMN == 5 ~ "Blink",
                        coderMN == 6 ~ "Other"),
    coderRA = case_when(coderRA == 1 ~ "Fixation",
                        coderRA == 2 ~ "Saccade",
                        coderRA == 3 ~ "PSO",
                        coderRA == 4 ~ "Smooth pursuit",
                        coderRA == 5 ~ "Blink",
                        coderRA == 6 ~ "Other")
  ) %>%
  mutate(
    gazeHMM = factor(gazeHMM, levels = c("Fixation", "Saccade", "PSO", "Smooth pursuit", "Blink", "Noise", "Other")),
    coderMN = factor(coderMN, levels = c("Fixation", "Saccade", "PSO", "Smooth pursuit", "Blink", "Noise", "Other")),
    coderRA = factor(coderRA, levels = c("Fixation", "Saccade", "PSO", "Smooth pursuit", "Blink", "Noise", "Other"))
  ) %>%
  filter(t > 4.3, t <= 4.5)


df_scarf <- df %>%
  pivot_longer(cols = c("coderMN", "coderRA", "gazeHMM"), names_to = "method", values_to = "classification") %>%
  mutate(
    method = factor(method, levels = c("coderRA", "coderMN", "gazeHMM"))
  )


scarf <- ggplot(df_scarf, aes(x=t, y=method, fill=classification)) + 
  geom_tile() + ylab(NULL) + xlab(NULL) +
  scale_fill_brewer(palette="Paired")

df_long <- df %>%
  pivot_longer(cols = c("x", "y", "vel", "acc", "angle"), names_to = "feature", values_to = "value") %>%
  mutate(
    feature = factor(
      case_when(feature == "x" ~ "x-coordinate",
                feature == "y" ~ "y-coordinate",
                feature == "vel" ~ "Velocity",
                feature == "acc" ~ "Acceleration",
                feature == "angle" ~ "Rel. angle"),
      levels = c("x-coordinate", "y-coordinate", "Velocity", "Acceleration", "Rel. angle")
    )
  )

data_plot <- ggplot(df_long, aes(x=t, y=value)) + 
  facet_wrap(vars(feature), ncol = 1, scales = "free", strip.position = "left") +
  geom_path() + xlab("Time (sec)") + ylab("") +
  theme(strip.placement = "outside",
        axis.text = element_text(size = 8))

data_plots <- list()
data_plots[["scarf"]] <- scarf

features <- levels(df_long$feature)
for(feat in features) {
  data_plots[[feat]] <- df_long %>%
    filter(feature == feat) %>%
    ggplot(aes(x=t, y=value)) +
    geom_path() + 
    xlab(if(feat != features[[length(features)]]) NULL else "Time (sec)") + 
    ylab(feat)
}

wrap_plots(data_plots, ncol = 1)
ggsave(path = "manuscript", filename = "scarf.png", width = 4.5, height = 8)
