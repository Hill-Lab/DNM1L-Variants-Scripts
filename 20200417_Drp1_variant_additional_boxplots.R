#
# Calculating and plotting Kcat,K0.5, Kcat/K0.5, and vmax for the following: 
# Drp1 WT, L230dup, G363D, G401S, and R710G GTPase results
# Experiment initially performed with recombinant L230dup
# However, later MS + SDS-PAGE confirmed protein was truncated and not usable

# Load libraries ----

library(tidyverse)
library(broom)
library(readxl)
library(minpack.lm)
library(gridExtra)
library(scales)
library(RColorBrewer)

theme_set(theme_bw() +
            theme(axis.text = element_text(size = 12, color = "black"),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank())
)

# Import and tidy data ----
pool <- read_csv("Drp1_variant_NADH_depletion_rates.csv")

# Determine and plot GTPase activity rate ----

### First, convert NADH oxidation rate into Drp1 activity (converting to min-1)
rates <- pool %>%
  na.omit() %>%
  mutate(., activity = (-a340_slope / (6220 * 0.4649 / 1e6) / 1)) %>%
  group_by(., mutant, gtp) %>%
  summarise(., avg_activity = mean(activity),
            stdev = sd(activity)) %>%
  mutate(., activity = avg_activity - avg_activity[gtp == 0]) %>%
  ungroup()

### activity is now in GTP hydrolyzed (µmol/min)
write_csv(rates, "Drp1_n3_activity_perMIN_KAM.csv")

### Then, plot activity against [GTP] (activity at 0 µM GTP)
kinetic_plot <- rates %>%
  ggplot(., aes(x = gtp, y = avg_activity, color = mutant)) +
  geom_point(shape = 15) +
  geom_errorbar(aes(ymin = avg_activity - stdev,
                    ymax = avg_activity + stdev),
                width = 15) +
  geom_line(stat = "smooth", method = "nlsLM",
            formula = y ~ (vmax*(x/(x+km))),
            method.args = list(start = c(vmax = 0.6,
                                         km = 200),
                               control = nls.control(maxiter = 100, tol = 1e-6)),
            se = FALSE,        
            fullrange = T,  
            size = 1,
            alpha = 0.6) +
  scale_x_continuous(limits = c(0, 2020),
                     breaks = c(0, 500, 1000, 1500, 2000)) +
  labs(title = "",
       x = "[GTP] (µM)",
       y = "Activity (min-1)") +
  theme(axis.text = element_text(size = 14, color = "black"),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_blank()
  )

kinetic_plot

ggsave("Drp1_n3_prior_correct_background_KAM.pdf", kinetic_plot,
       width = 10, height = 8, units = "cm", dpi = 300)

kinetic_plot2 <- rates %>%
  mutate(mutant = 
           factor(mutant,
                  levels = c("G363D", "G401S", "WT", "R710G", "L230dup"))) %>%
  ggplot(., aes(x = gtp, y = activity, color = mutant)) +
  geom_point(shape = 15) +
  geom_errorbar(aes(ymin = activity - stdev,
                    ymax = activity + stdev),
                width = 15) +
  geom_line(stat = "smooth", method = "nlsLM",
            formula = y ~ (vmax*(x/(x+km))),
            method.args = list(start = c(vmax = 0.6,
                                         km = 200),
                               control = nls.control(maxiter = 100, tol = 1e-6)),
            se = FALSE,        
            fullrange = T,  
            size = 1,
            alpha = 0.6) +
  scale_x_continuous(limits = c(0, 2020),
                     breaks = c(0, 500, 1000, 1500, 2000)) +
  scale_color_manual(values = c("WT" = "red",
                                "L230dup" = "blue",
                                "G363D" = "forestgreen",
                                "G401S" = "purple",
                                "R710G" = "orange")) +
  labs(title = "",
       x = "[GTP] (µM)",
       y = "Activity (min-1)") +
  theme(axis.text = element_text(size = 14, color = "black"),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_blank()
  )

kinetic_plot2

ggsave("Drp1_n3_activity_plot_KAM.pdf", kinetic_plot,
       width = 10, height = 8, units = "cm", dpi = 300)

# Kinetic plot - no L230dup due to protein being truncated
# Dark2 color scheme
kinetic_plot2_noL230dup_dark2 <- rates %>%
  filter(., mutant != "L230dup") %>%
  mutate(mutant = 
           factor(mutant,
                  levels = c("G363D", "G401S", "WT", "R710G"))) %>%
  ggplot(., aes(x = gtp, y = activity, color = mutant)) +
  geom_point(shape = 15) +
  geom_errorbar(aes(ymin = activity - stdev,
                    ymax = activity + stdev),
                width = 15) +
  geom_line(stat = "smooth", method = "nlsLM",
            formula = y ~ (vmax*(x/(x+km))),
            method.args = list(start = c(vmax = 0.6,
                                         km = 200),
                               control = nls.control(maxiter = 100, tol = 1e-6)),
            se = FALSE,        
            fullrange = T,  
            size = 1,
            alpha = 0.6) +
  scale_x_continuous(limits = c(0, 2020),
                     breaks = c(0, 500, 1000, 1500, 2000)) +
  scale_color_brewer(palette = "Dark2") +
  labs(title = "",
       x = "[GTP] (µM)",
       y = "Activity (min-1)") +
  theme(axis.text = element_text(size = 14, color = "black"),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_blank()
    )

kinetic_plot2_noL230dup_dark2

ggsave("Drp1_n3_activity_plot_noL230dup_dark2_KAM.pdf", kinetic_plot,
       width = 10, height = 8, units = "cm", dpi = 300)

### Calculate Vmax and K0.5 
vmax_km <- rates %>%
  group_by(mutant) %>%
  do(tidy(nlsLM(formula = activity ~ (vmax*(gtp/(gtp+km))),
                start = list(vmax = 0.6, km = 200),
                trace = TRUE,
                data = .)))
vmax_km

### vmax units = µmol/min
### km units = µM
#### Does the G401S mutation shift G401S to an unfavored ramachandran angle?
#### look at crystal structure and see if this violates it
write_csv(vmax_km, "Drp1_n3_activity_KAM.csv")

### Add residual plot and resume here
vmax_km_residual <- rates %>%
  group_by(mutant) %>%
  do(augment(nlsLM(formula = activity ~ (vmax*(gtp/(gtp+km))),
                   start = list(vmax = 0.6, km = 200),
                   trace = TRUE, 
                   data = .))) %>%
  ungroup()

  resid_plot <- vmax_km_residual %>%
  mutate(mutant = 
           factor(mutant,
                  levels = c("G363D", "G401S", "WT", "R710G", "L230dup"))) %>%
  ggplot(aes(x = gtp, y = .resid, color = mutant)) +
  geom_point(shape = 15) +
  geom_hline(yintercept = 0, linetype = 2, color = "red") +
  scale_y_continuous(breaks = c(-0.03, 0, 0.03, 0.06, 0.09)) +
  scale_color_manual(values = c("WT" = "red",
                                "L230dup" = "blue",
                                "G363D" = "forestgreen",
                                "G401S" = "purple",
                                "R710G" = "orange")) +
  labs(x = "", 
       y = "Residuals") +
  theme_bw() +
  theme(axis.text = element_text(size = 14, color = "black"),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_blank()
  )
resid_plot

### Add residual plot - no L230dup due to protein being truncated
# Dark color scheme
resid_plot_noL230dup <- vmax_km_residual %>%
  filter(., mutant != "L230dup") %>%
  mutate(mutant = 
           factor(mutant,
                  levels = c("G363D", "G401S", "WT", "R710G"))) %>%
  ggplot(aes(x = gtp, y = .resid, color = mutant)) +
  geom_point(shape = 15) +
  geom_hline(yintercept = 0, linetype = 2, color = "red") +
  scale_y_continuous(limits =c(-0.1, 0.1),
                     breaks = c(-0.1, 0, 0.1)) +
  scale_color_manual(values = c("WT" = "#7570B3",
                                "G363D" = "#1B9E77",
                                "G401S" = "#D95F02",
                                "R710G" = "#E7298A")) +
  labs(x = "", 
       y = "Residuals") +
  theme_bw() +
  theme(axis.text = element_text(size = 14, color = "black"),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position = "none",
  )
resid_plot_noL230dup

# Combine residual plot with activity plot
### generate a layout matrix to plot residuals at 1/3 size of regression plot

lay <- cbind(c(1, 2, 2))

final_plot <- grid.arrange(grobs = list(resid_plot_noL230dup, kinetic_plot2_noL230dup_dark2),
                           layout_matrix = lay,
                           top = "Drp1 GTPase Activity")

ggsave("Drp1_n3_activity_final_plot__noL230dup_KAM_v2.pdf", final_plot,
       width = 14, height = 10, units = "cm", dpi = 300)

# L230dup appears to be inactive or enzymaticaly dead - makes sense given
# MS confirmed protein to be truncated in GTPase domain
# G363D appears hyperactive
# R710G has decreased Vmax

# Determine if difference between mutants and WT activity is significant

pooled <-  pool %>%
  mutate(., activity = (-a340_slope / (6220 * 0.4649 / 1e6) / 1)) %>%
  group_by(., mutant, gtp, rep) %>%
  summarise(., avg_activity = mean(activity),
            stdev = sd(activity)) %>%
  group_by(mutant, rep) %>%
  mutate(., activity = avg_activity - avg_activity[gtp == 0]) %>%
  ungroup()

pooled

# calculate vmax and km for each replicate and then do ANOVA between mutants
vmax_km_reps <- pooled %>%
  group_by(., mutant, rep) %>%
  do(tidy(nlsLM(formula = activity ~ (vmax*(gtp/(gtp+km))),
                start = list(vmax = 0.6, km = 200),
                trace = TRUE,
                data = .))) %>%
  spread(key = term, value = estimate) %>%
  ungroup()

vmax_km_reps # results are from a global fit of n = 3 data

# calculate kcat and kcat/K0.5
kcat_reps = subset(vmax_km_reps, select = -c(std.error, statistic, p.value)) %>%
  group_by(., mutant, rep) %>%
  summarise_all(na.omit) %>%
  mutate(kcat = vmax/1000000) %>%
  mutate(kcatkm = kcat/km) %>%
  ungroup()

kcat_reps # kcat and kcat/km values calculated from 
          # vmax and km values in vmax_km_reps
 
# Plot of kcat values
kcat_reps %>%
  mutate(mutant = 
           factor(mutant,
                  levels = c("WT", "L230dup", "G363D", "G401S", "R710G"))) %>%
  ggplot(aes(x = mutant, y = kcat)) +
  geom_boxplot() +
  geom_jitter(aes(y = kcat, color = mutant), size = 1.5, shape = 15) +
  scale_color_manual(values = c("WT" = "red",
                                "L230dup" = "blue",
                                "G363D" = "forestgreen",
                                "G401S" = "purple",
                                "R710G" = "orange")) +
  labs(x = "",
       y = "kcat min-1") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        legend.position = "none")

ggsave("Drp1_kcat_boxplot_KAM.pdf",
       width = 8, height = 8, units = "cm", dpi = 300)

# Plot of kcat values - No L230dup due to truncation
kcat_reps %>%
  filter(., mutant != "L230dup") %>%
  mutate(mutant = 
           factor(mutant,
                  levels = c("WT", "G363D", "G401S", "R710G"))) %>%
  ggplot(aes(x = mutant, y = kcat)) +
  geom_boxplot() +
  geom_jitter(aes(y = kcat, color = mutant), size = 1.5, shape = 15) +
  scale_color_manual(values = c("WT" = "#7570B3",
                                "G363D" = "#1B9E77",
                                "G401S" = "#D95F02",
                                "R710G" = "#E7298A")) +
  labs(x = "",
       y = "kcat min-1") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        legend.position = "none")

ggsave("Drp1_kcat_boxplot_noL230dup_KAM.pdf",
       width = 8, height = 8, units = "cm", dpi = 300)


# Plot of kcat/K0.5 values
kcat_reps %>%
  mutate(mutant = 
           factor(mutant,
                  levels = c("WT", "L230dup", "G363D", "G401S", "R710G"))) %>%
  ggplot(aes(x = mutant, y = kcatkm)) +
  geom_boxplot() +
  geom_jitter(aes(y = kcatkm, color = mutant), size = 1.5, shape = 15) +
  scale_color_manual(values = c("WT" = "red",
                                "L230dup" = "blue",
                                "G363D" = "forestgreen",
                                "G401S" = "purple",
                                "R710G" = "orange")) +
  labs(x = "",
       y = "kcat (min-1)/K0.5 µM") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        legend.position = "none")

ggsave("Drp1_kcatkm_boxplot_KAM.pdf",
       width = 8, height = 8, units = "cm", dpi = 300)

# Plot of kcat/K0.5 values - No L230dup due to protein being truncated
kcat_reps %>%
  filter(., mutant != "L230dup") %>%
  mutate(mutant = 
           factor(mutant,
                  levels = c("WT", "G363D", "G401S", "R710G"))) %>%
  ggplot(aes(x = mutant, y = kcatkm)) +
  geom_boxplot() +
  geom_jitter(aes(y = kcatkm, color = mutant), size = 1.5, shape = 15) +
  scale_color_manual(values = c("WT" = "#7570B3",
                                "G363D" = "#1B9E77",
                                "G401S" = "#D95F02",
                                "R710G" = "#E7298A")) +
  labs(x = "",
       y = "kcat (min-1)/K0.5 µM") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        legend.position = "none")

ggsave("Drp1_kcatkm_boxplot_noL230dup_KAM.pdf",
       width = 8, height = 8, units = "cm", dpi = 300)

# ANOVA on Kcat + post-hoc Tukey
anova_kcat <- kcat_reps %>%
  do(tidy(aov(kcat ~ mutant, data = .)))

anova_kcat %>% print(width = Inf)

anova_kcat_tukey <- kcat_reps %>%
  do(tidy(TukeyHSD(aov(kcat ~ mutant, data = .))))

# ANOVA on Kcat/Km + post-hoc Tukey
anova_kcatkm <- kcat_reps %>%
  do(tidy(aov(kcatkm ~ mutant, data = .)))

anova_kcatkm %>% print(width = Inf)

anova_kcatkm_tukey <- kcat_reps %>%
  do(tidy(TukeyHSD(aov(kcatkm ~ mutant, data = .))))

# Determine if there any outliers in Vmax values via boxplot
vmax_km_reps %>%
  mutate(mutant = 
           factor(mutant,
                  levels = c("WT", "L230dup", "G363D", "G401S", "R710G"))) %>%
  ggplot(aes(x = mutant, y = vmax)) +
  geom_boxplot() +
  geom_jitter(aes(y = vmax, color = mutant), size = 1.5, shape = 15) +
  scale_color_manual(values = c("WT" = "red",
                                "L230dup" = "blue",
                                "G363D" = "forestgreen",
                                "G401S" = "purple",
                                "R710G" = "orange")) +
  labs(x = "",
       y = "Vmax (µmol/min)") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        legend.position = "none")

ggsave("Drp1_vmax_boxplot_KAM.pdf",
       width = 8, height = 8, units = "cm", dpi = 300)

# Determine if there any outliers in Vmax values via boxplot
# No L230dup due to protein being truncated
vmax_km_reps %>%
  filter(., mutant != "L230dup") %>%
  mutate(mutant = 
           factor(mutant,
                  levels = c("WT", "G363D", "G401S", "R710G"))) %>%
  ggplot(aes(x = mutant, y = vmax)) +
  geom_boxplot() +
  geom_jitter(aes(y = vmax, color = mutant), size = 1.5, shape = 15) +
  scale_color_manual(values = c("WT" = "#7570B3",
                                "G363D" = "#1B9E77",
                                "G401S" = "#D95F02",
                                "R710G" = "#E7298A")) +
  labs(x = "",
       y = "Vmax (µmol/min)") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        legend.position = "none")

ggsave("Drp1_vmax_boxplot_noL230dup_KAM.pdf",
       width = 8, height = 8, units = "cm", dpi = 300)

# ANOVA + post-hoc Tukey on Vmax
anova_vmax <- vmax_km_reps %>%
  do(tidy(aov(vmax ~ mutant, data = .)))

anova_vmax %>% print(width = Inf)

anova_vmax_tukey <- vmax_km_reps %>%
  do(tidy(TukeyHSD(aov(vmax ~ mutant, data = .))))

# significance results of Vmax values between mutants
anova_vmax_tukey %>% print(width = Inf)

write_csv(anova_vmax_tukey, "anova_tukey_VMAX_results_KAM.csv")

# Determine if there any outliers in Km values via boxplot
vmax_km_reps %>%
  mutate(mutant = 
           factor(mutant,
                  levels = c("WT", "L230dup", "G363D", "G401S", "R710G"))) %>%
  ggplot(aes(x = mutant, y = km)) +
  geom_boxplot() +
  geom_jitter(aes(y = km, color = mutant), size = 1.5, shape = 15) +
  scale_color_manual(values = c("WT" = "red",
                                "L230dup" = "blue",
                                "G363D" = "forestgreen",
                                "G401S" = "purple",
                                "R710G" = "orange")) +
  scale_y_continuous(limits = c(-5, 300),
                     breaks = c(0, 100, 200, 300)) +
  labs(x = "",
       y = "K0.5 (µM)") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        legend.position = "none")

ggsave("Drp1_km_boxplot_KAM.pdf",
       width = 8, height = 8, units = "cm", dpi = 300)

# Determine if there any outliers in Km values via boxplot
# No L230dup due to protein being truncated
vmax_km_reps %>%
  filter(., mutant != "L230dup") %>%
  mutate(mutant = 
           factor(mutant,
                  levels = c("WT", "G363D", "G401S", "R710G"))) %>%
  ggplot(aes(x = mutant, y = km)) +
  geom_boxplot() +
  geom_jitter(aes(y = km, color = mutant), size = 1.5, shape = 15) +
  scale_color_manual(values = c("WT" = "#7570B3",
                                "G363D" = "#1B9E77",
                                "G401S" = "#D95F02",
                                "R710G" = "#E7298A")) +
  scale_y_continuous(limits = c(-5, 300),
                     breaks = c(0, 100, 200, 300)) +
  labs(x = "",
       y = "K0.5 (µM)") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        legend.position = "none")

ggsave("Drp1_km_boxplot_noL230dup_KAM.pdf",
       width = 8, height = 8, units = "cm", dpi = 300)

# ANOVA + post-hoc Tukey on Km
anova_km <- vmax_km_reps %>%
  do(tidy(aov(km ~ mutant, data = .)))

anova_km %>% print(width = Inf)

anova_km_tukey <- vmax_km_reps %>%
  do(tidy(TukeyHSD(aov(km ~ mutant, data = .))))

# significance results of Km values between mutants
anova_km_tukey %>% print(width = Inf)

write_csv(anova_km_tukey, "anova_tukey_Km_results_KAM.csv")

# Perform ANOVA + post-hoc Tukey on activity
# between GTP concentrations and mutant type
activity_aov <- pooled %>%
  group_by(gtp) %>%
  do(tidy(aov(activity ~ mutant, data = .)))

activity_aov

activity_aov_tukey <- pooled %>%
  group_by(gtp) %>%
  do(tidy(TukeyHSD(aov(activity ~ mutant, data = .))))
activity_aov_tukey

activity_aov_tukey %>%
  arrange(adj.p.value) %>%
  filter(adj.p.value <= 0.05) %>%
  print(width = Inf, n = Inf)

write_csv(activity_aov_tukey, "activity_by_GTP_conc_anova_tukey_KAM.csv")