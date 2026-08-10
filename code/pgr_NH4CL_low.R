## Growth across different percentage concentrations of NH4CL (all else held constant - glucose at 0.12%)

### Load required packages 
library(grow96)
library(tidyverse)
library(splines)
library(patchwork)
library(wesanderson)

### Load data with grow96
data <- processODData(specPath="specs/NHCA", dataPath="data/NHCA")

### Blank OD data
blankdata <- blankODs(data %>% filter(experiment == "NH4CL"), 
                      method = "fromPlate", 
                      perTimePoint = TRUE) %>% 
  filter(!is.na(glucose)) %>% 
  filter(WellType == "DATA") %>% 
  mutate(overnight = str_split_i(overnight_technical, "_", 1),
         technical = str_split_i(overnight_technical, "_", 2)) 

## Set beginning and end time cut-offs
# time_limits <- tibble(NH4CL = unique(blankdata$NH4CL)) %>%
#   mutate(timeclip_upper = case_when(NH4CL == 200 ~ 12,
#                                     NH4CL == 150 ~ 12,
#                                     NH4CL == 100 ~ 12,
#                                     NH4CL == 50 ~ 12,
#                                     NH4CL == 20 ~ 12,
#                                     NH4CL == 0 ~ 12)) %>%
#   mutate(timeclip_lower = case_when(NH4CL == 200 ~ 1,
#                                     NH4CL == 150 ~ 1,
#                                     NH4CL == 100 ~ 1,
#                                     NH4CL == 50 ~ 1,
#                                     NH4CL == 20 ~ 1,
#                                     NH4CL == 0 ~ 1))

## Further filtering/joining
blankdata <- blankdata %>% 
#  left_join(time_limits, by = "NH4CL") %>% 
  mutate(NH4CL = factor(NH4CL, sort(as.numeric(unique(NH4CL)))))


## Plot dynamics
gluc_gradient_ts_ms <- blankdata %>% 
#  filter(Time_h > timeclip_lower) %>% # & Time_h < timeclip_upper) %>% 
  ggplot() +
  geom_line(aes(y = blankedOD, 
                x = Time_h, 
                group = interaction(Replicate,overnight_technical,NH4CL),
                col = NH4CL), alpha = 0.5) +
  theme_classic() +
  theme(strip.background = element_blank()) +
  labs(y = "OD (600nm)", x = "Time (h)") +
  scale_color_viridis_d()#  facet_wrap(~ as.factor(NH4CL), scales = "free_x", ncol = 4) #+ scale_y_log10()

gluc_gradient_ts_ms

### Obtain PGR for each independent well
pgr_dat <- blankdata %>% 
  group_by(Replicate, NH4CL, overnight, technical) %>%
  arrange(Time_min, .by_group = TRUE) %>%
  mutate(g = log(lead(blankedOD)/lag(blankedOD)) / (lead(Time_h)-lag(Time_h))) %>%
  ungroup() %>% 
  filter(blankedOD != 0) # %>% filter(Time_h > timeclip_lower & Time_h < timeclip_upper) 

p_pgr <- ggplot(pgr_dat, aes(x = blankedOD, y = g)) +
  geom_point(alpha = 0.2, size = 0.5) +
  stat_smooth(aes(group = interaction(overnight, technical, Replicate, NH4CL),
                  colour = NH4CL), 
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
  scale_color_viridis_d() + 
  coord_cartesian(ylim = c(-0.25, 1.0)) + 
  facet_wrap(~NH4CL, scales = "free_x", ncol = 4)

p_pgr


### Obtain mean pgr in density bins
mean_by_N <- pgr_dat %>%  
  group_by(NH4CL) %>%
  mutate(bin = cut(blankedOD, breaks = (seq((min(blankedOD)), (max(blankedOD)), length.out = 60)))) %>%
  group_by(bin, NH4CL) %>%
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
  ggplot(aes(y = g_median, x = blankedOD_mean, col = NH4CL, ymin = g_25, ymax = g_75)) +
  geom_errorbar(alpha = 0.5, 
                linewidth = 0.2, 
                show.legend = FALSE) +
  geom_point(alpha = 0.5, size = 1) + 
  scale_color_viridis_d() + 
  theme_classic() + 
  coord_cartesian(ylim = c(-0.05, 0.7)) +
  labs(x = "OD (600nm)", 
       y = expression("per capita growth (h"^-1*")"), 
       colour = "NH4CL (%)") +
  theme(legend.position = c(1.02, 1.03), 
        legend.justification = c(1, 1),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.key.height = unit(0.3, "cm"))

pgr_mean





