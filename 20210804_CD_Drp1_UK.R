# Analyzing circular dichroism data
# experiment: CD of WT Drp1 + UK pathological variants
# data collected on 20210804
# data corrected for empty cuvette + buffer signal
# scaled in Excel prior to loading to provide comparable baselines
# MRE calculated for each sample at each wavelength evaluated

library(tidyverse)
library(broom)
library(readxl)
library(minpack.lm)
library(ggpmisc)
library(RColorBrewer)

# Import and tidy data  --------------------------------------------------------
### update excel sheet with proper heading info separated by "_"
# this data already corrected for signal from cuvette and buffer
CD <- read_excel("20210804_corrected_scaled.xlsx", sheet = 1) %>%
  gather(., "tmp1", "signal", 2:9) %>%
  separate(., col = tmp1, into = c("prot", "type"), 
           sep = "_")

# Visualize CD spectra  - MRE
CD %>%
  filter(., type == "MRE") %>%
  filter(., wavelength > 198 & wavelength < 260) %>%
  ggplot(., aes(x = wavelength, y = signal, color = prot)) +
  geom_point(size = 0.5) +
  geom_line(size = 0.5) +
  scale_color_manual(values = c("WT" = "#7570B3",
                                "G363D" = "#1B9E77",
                                "G401S" = "#D95F02",
                                "R710G" = "#E7298A")) +
  labs(title = "CD Spectra Drp1 WT + Variants", 
       color = "Protein (0.05 mg/mL)",
       x = "wavelength (nm)",
       y = "MRE") +
  theme_bw() +
  scale_y_continuous(limits = c(-14200, 8100), breaks = c(-12000, -8000, -4000,
                                                          0, 4000, 8000)) +
  scale_x_continuous(limits = c(198, 260), breaks = c(200, 210, 220, 230, 240, 250, 260)) +
  theme(axis.text = element_text(size = 14, color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14, color = "black")
  )

ggsave("CD_MRE_corrected_scaled.pdf", 
       width = 40, height = 21, units = "cm")

