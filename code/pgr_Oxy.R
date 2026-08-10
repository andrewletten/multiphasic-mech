## Growth across different media volumes 

### Load required packages 
library(grow96)
library(tidyverse)
library(splines)
library(patchwork)
library(wesanderson)

### Load data with grow96
data <- processODData(specPath="specs/Oxy", dataPath="data/Oxy")

### Blank OD data
blankdata <- blankODs(data, method = "fromPlate", perTimePoint = TRUE) %>% 
  filter(!is.na(glucose)) %>% 
  filter(WellType == "DATA") %>% 
  mutate(overnight = str_split_i(overnight_technical, "_", 1),
         technical = str_split_i(overnight_technical, "_", 2)) 

## Set beginning and end time cut-offs
time_limits <- tibble(volume = unique(blankdata$volume)) %>%
  mutate(timeclip_upper = case_when(volume == 180 ~ 14,
                                    volume == 140 ~ 14,
                                    volume == 100 ~ 14,
                                    volume == 60 ~ 14,
                                    volume == 20 ~ 14,
                                    volume == 0 ~ 14)) %>%
  mutate(timeclip_lower = case_when(volume == 180 ~ 1,
                                    volume == 140 ~ 1,
                                    volume == 100 ~ 1,
                                    volume == 60 ~ 1,
                                    volume == 20 ~ 1,
                                    volume == 0 ~ 1))

## Further filtering/joining
blankdata <- blankdata %>% 
  left_join(time_limits, by = "volume") %>% 
  mutate(volume = factor(volume, sort(as.numeric(unique(volume))))) %>% 
  # fudge correct for differences in volume 
  mutate(blankedOD_correct = blankedOD * 100/as.numeric(as.character(volume)))


## Plot dynamics
gluc_gradient_ts_ms <- blankdata %>% 
  filter(Time_h > timeclip_lower & Time_h < timeclip_upper) %>% 
  filter(volume != 100 & volume != 20) %>% 
  ggplot() +
  geom_line(aes(y = blankedOD_correct, # using corrected OD
                x = Time_h, 
                group = interaction(Replicate,overnight_technical,volume, Well),
                col = volume), alpha = 0.5) +
  theme_classic() +
  theme(strip.background = element_blank()) +
  labs(y = "OD (600nm)", x = "Time (h)") +
  scale_color_viridis_d() #+ facet_wrap(~ as.factor(volume), scales = "free_x", ncol = 4) #+ scale_y_log10()

gluc_gradient_ts_ms

### Obtain PGR for each independent well
pgr_dat <- blankdata %>% 
  filter(volume != 100 & volume != 20) %>% 
  group_by(Replicate, volume, overnight, technical, Well) %>%
  arrange(Time_min, .by_group = TRUE) %>%
  mutate(g = log(lead(blankedOD)/lag(blankedOD)) / (lead(Time_h)-lag(Time_h))) %>%
  mutate(g_correct = log(lead(blankedOD_correct)/lag(blankedOD_correct)) / (lead(Time_h)-lag(Time_h))) %>%
  ungroup() %>% 
  filter(blankedOD != 0) %>% 
  filter(Time_h > timeclip_lower & Time_h < timeclip_upper) 

p_pgr <- ggplot(pgr_dat, aes(x = blankedOD, y = g)) +
  geom_point(alpha = 0.2, size = 0.5) +
  stat_smooth(aes(group = interaction(overnight, technical, Replicate, volume),
                  colour = volume), 
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
  scale_color_viridis_d() + 
  coord_cartesian(ylim = c(-0.25, 1.0)) + 
  facet_wrap(~volume, scales = "free_x", ncol = 4)

p_pgr


### Obtain mean pgr in density bins
mean_by_N <- pgr_dat %>%  
  group_by(volume) %>%
  mutate(bin = cut(blankedOD, breaks = (seq((min(blankedOD)), (max(blankedOD)), length.out = 60)))) %>%
  group_by(bin, volume) %>%
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
  ggplot(aes(y = g_median, x = blankedOD_mean, col = volume, ymin = g_25, ymax = g_75)) +
  geom_errorbar(alpha = 0.5, 
                linewidth = 0.2, 
                show.legend = FALSE) +
  geom_point(alpha = 0.5, size = 1) + 
  scale_color_viridis_d() + 
  theme_classic() + 
  coord_cartesian(ylim = c(-0.05, 0.7)) +
  labs(x = "OD (600nm)", 
       y = expression("per capita growth (h"^-1*")"), 
       colour = "volume (%)") +
  theme(legend.position = c(1.02, 1.03), 
        legend.justification = c(1, 1),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.key.height = unit(0.3, "cm"))

pgr_mean





