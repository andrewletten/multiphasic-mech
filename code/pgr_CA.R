## Growth across different percentage concentrations of CACL (all else held constant - glucose at 0.12%)

### Load required packages 
library(grow96)
library(tidyverse)
library(splines)
library(patchwork)
library(wesanderson)

### Load data with grow96
data <- processODData(specPath="specs/NHCA", dataPath="data/NHCA")

### Blank OD data
blankdata <- blankODs(data %>% filter(experiment == "CACL"), 
                      method = "fromPlate", 
                      perTimePoint = TRUE) %>% 
  filter(!is.na(glucose)) %>% 
  filter(WellType == "DATA") %>% 
  mutate(overnight = str_split_i(overnight_technical, "_", 1),
         technical = str_split_i(overnight_technical, "_", 2)) 

# Set beginning and end time cut-offs
time_limits <- tibble(CACL = unique(blankdata$CACL)) %>%
  mutate(timeclip_upper = case_when(CACL == 200 ~ 12,
                                    CACL == 150 ~ 12,
                                    CACL == 100 ~ 12,
                                    CACL == 50 ~ 12,
                                    CACL == 20 ~ 12,
                                    CACL == 0 ~ 12)) %>%
  mutate(timeclip_lower = case_when(CACL == 200 ~ 1,
                                    CACL == 150 ~ 1,
                                    CACL == 100 ~ 1,
                                    CACL == 50 ~ 1,
                                    CACL == 20 ~ 1,
                                    CACL == 0 ~ 1))

## Further filtering/joining
blankdata <- blankdata %>% 
  left_join(time_limits, by = "CACL") %>% 
  mutate(CACL = factor(CACL, sort(as.numeric(unique(CACL)))))


## Plot dynamics
gluc_gradient_ts_ms <- blankdata %>% 
  filter(Time_h > timeclip_lower) %>% # & Time_h < timeclip_upper) %>% 
  ggplot() +
  geom_line(aes(y = blankedOD, 
                x = Time_h, 
                group = interaction(Replicate,overnight_technical,CACL),
                col = CACL), alpha = 0.5) +
  theme_classic() +
  theme(strip.background = element_blank()) +
  labs(y = "OD (600nm)", x = "Time (h)") +
  scale_color_viridis_d()#  facet_wrap(~ as.factor(CACL), scales = "free_x", ncol = 4) #+ scale_y_log10()

gluc_gradient_ts_ms

### Obtain PGR for each independent well
pgr_dat <- blankdata %>% 
  group_by(Replicate, CACL, overnight, technical) %>%
  arrange(Time_min, .by_group = TRUE) %>%
  mutate(g = log(lead(blankedOD)/lag(blankedOD)) / (lead(Time_h)-lag(Time_h))) %>%
  ungroup() %>% 
  filter(blankedOD != 0) %>% 
  filter(Time_h > timeclip_lower & Time_h < timeclip_upper) 

p_pgr <- ggplot(pgr_dat, aes(x = blankedOD, y = g)) +
  geom_point(alpha = 0.2, size = 0.5) +
  stat_smooth(aes(group = interaction(overnight, technical, Replicate, CACL),
                  colour = CACL), 
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
  facet_wrap(~CACL, scales = "free_x", ncol = 4)

p_pgr


### Obtain mean pgr in density bins
mean_by_N <- pgr_dat %>%  
  group_by(CACL) %>%
  mutate(bin = cut(blankedOD, breaks = (seq((min(blankedOD)), (max(blankedOD)), length.out = 60)))) %>%
  group_by(bin, CACL) %>%
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
  ggplot(aes(y = g_median, x = blankedOD_mean, col = CACL, ymin = g_25, ymax = g_75)) +
  geom_errorbar(alpha = 0.5, 
                linewidth = 0.2, 
                show.legend = FALSE) +
  geom_point(alpha = 0.5, size = 1) + 
  scale_color_viridis_d() + 
  theme_classic() + 
  coord_cartesian(ylim = c(-0.05, 0.8)) +
  labs(x = "OD (600nm)", 
       y = expression("per capita growth (h"^-1*")"), 
       colour = "CACL (%)") +
  theme(legend.position = c(1.02, 1.03), 
        legend.justification = c(1, 1),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.key.height = unit(0.3, "cm"))

pgr_mean





