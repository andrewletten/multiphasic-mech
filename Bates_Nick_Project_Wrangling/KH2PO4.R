## Growth across different percentage KH2PO4 of KH2PO4 minimal media (only glucose held constant)

#For e.g., NaCl, find and replace all instances of KH2PO4 with NaCl, and all the plots should still work.

install.packages("RColorBrewer")
install.packages("paletteer")
install.packages("ggthemes")

### Load required packages 
library(grow96)
library(tidyverse)
library(splines)
library(patchwork)
library(wesanderson)
library(RColorBrewer)
library(paletteer)
library(ggthemes)

### Load data with grow96
data <- processODData(specPath="specs/KH2PO4", dataPath="data/KH2PO4")

### Blank OD data
blankdata <- blankODs(data, method = "fromPlate", perTimePoint = TRUE) %>% 
  filter(WellType == "DATA") %>%
  mutate(overnight = str_split_i(strains, "_", 1),
    technical = str_split_i(strains, "_", 2)) 


### Set beginning and end time cut-offs
time_limits <- tibble(KH2PO4 = unique(blankdata$KH2PO4)) %>%
  mutate(timeclip_upper = case_when(KH2PO4 == 200 ~ 10,
                                    KH2PO4 == 150 ~ 10,
                                    KH2PO4 == 100 ~ 10,
                                    KH2PO4 == 50 ~ 12.5,
                                    KH2PO4 == 25 ~ 12.5,
                                    KH2PO4 == 0 ~ 6)) %>%
  mutate(timeclip_lower = case_when(KH2PO4 == 200 ~ 1,
                                    KH2PO4 == 150 ~ 1,
                                    KH2PO4 == 100 ~ 1,
                                    KH2PO4 == 50 ~ 1,
                                    KH2PO4 == 25 ~ 1,
                                    KH2PO4 == 0 ~ 1))

## Further filtering/joining
blankdata <- blankdata %>% 
  left_join(time_limits, by = "KH2PO4") %>% 
  mutate(KH2PO4 = factor(KH2PO4, sort(as.numeric(unique(KH2PO4)))))


## Plot dynamics
gluc_gradient_ts_ms <- blankdata %>% 
  filter(Time_h > timeclip_lower) %>% # & Time_h < timeclip_upper) %>% 
  ggplot() +
  geom_line(aes(y = blankedOD, 
                x = Time_h, 
                group = interaction(Replicate,strains,KH2PO4),
                col = KH2PO4), alpha = 0.5) +
  theme_classic() +
  theme(strip.background = element_blank()) +
  labs(y = "OD (600nm)", x = "Time (h)") +
  scale_color_paletteer_d("ggthemes::excel_Parallax") #+  facet_wrap(~ as.factor(KH2PO4), scales = "free_x", ncol = 6) #+ scale_y_log10()

gluc_gradient_ts_ms

### Obtain PGR for each independent well
pgr_dat <- blankdata %>% 
  group_by(Replicate, KH2PO4, overnight, technical) %>%
  arrange(Time_min, .by_group = TRUE) %>%
  mutate(g = log(lead(blankedOD)/lag(blankedOD)) / (lead(Time_h)-lag(Time_h))) %>%
  ungroup() %>% 
  filter(blankedOD != 0) #%>%   filter(Time_h > timeclip_lower & Time_h < timeclip_upper) 

p_pgr <- ggplot(pgr_dat, aes(x = blankedOD, y = g)) +
  geom_point(alpha = 0.2, size = 0.2) +
  stat_smooth(aes(group = interaction(overnight, technical, Replicate, KH2PO4),
                  colour = KH2PO4), 
              geom = "line",
              method = "lm",
              formula = y ~ ns(x, df = 5),
              alpha = 0.4, 
              linewidth = 0.6,
              se = FALSE) +
  theme_classic()  +
  theme(strip.background = element_blank(),
        strip.text = element_text(margin = margin(0, 0, 10, 0)),
        legend.position = "none") +
  labs(x = "OD (600nm)", y = expression("per capita growth (h"^-1*")")) +
  #  scale_y_continuous(breaks = seq(-0.2, 1.0, 0.2))  +
  scale_color_paletteer_d("ggthemes::excel_Parallax") + 
  coord_cartesian(ylim = c(-0.25, 2.0)) +
  facet_wrap(~KH2PO4, scales = "free_x", ncol = 6)

p_pgr


### Obtain mean pgr in density bins
mean_by_N <- pgr_dat %>%  
  group_by(KH2PO4) %>%
  mutate(bin = cut(blankedOD, breaks = (seq((min(blankedOD)), (max(blankedOD)), length.out = 60)))) %>%
  group_by(bin, KH2PO4) %>%
  summarise(
    blankedOD_mean = mean(blankedOD),
    blankedOD_median = median(blankedOD),
    g_mean = mean(g, na.rm = TRUE),
    g_sd = sd(g, na.rm = TRUE),
    g_median = median(g, na.rm = TRUE),
    g_75 = quantile(g, 0.75, na.rm = TRUE),
    g_25 = quantile(g, 0.25, na.rm = TRUE),
    .groups = "drop"
  ) %>% 
  filter(!is.na(bin)) 

pgr_mean <- mean_by_N %>%
  ggplot(aes(y = g_median, x = blankedOD_mean, col = KH2PO4, ymin = g_25, ymax = g_75)) +
  geom_errorbar(alpha = 0.5, 
                linewidth = 0.2, 
                show.legend = FALSE) +
  geom_point(alpha = 0.5, size = 1) + 
  scale_color_paletteer_d("RColorBrewer::Dark2") + 
  theme_classic() + 
  coord_cartesian(ylim = c(-0.05, 1.2)) +
  labs(x = "OD (600nm)", 
       y = expression("per capita growth (h"^-1*")"), 
       colour = "KH2PO4 (%)") +
  theme(legend.position = c(1.02, 1.03), 
        legend.justification = c(1, 1),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.key.height = unit(0.3, "cm"))

pgr_mean


