## Growth across different percentage concentrations of glucose (all else held constant)

### Load required packages 
library(grow96)
library(tidyverse)
library(splines)
library(patchwork)
library(wesanderson)

### Load data with grow96
data <- processODData(specPath="specs/CHO", dataPath="data/CHO")

### Blank OD data
blankdata <- blankODs(data, method = "fromPlate", perTimePoint = TRUE) %>% 
  filter(!is.na(glucose)) %>% 
  filter(WellType == "DATA") %>% 
  mutate(overnight = str_split_i(overnight_technical, "_", 1),
         technical = str_split_i(overnight_technical, "_", 2)) %>% 
  # remove single outlier (OD measurements for this individual assay well far above others) 
  filter(!(overnight_technical == "O1_T1" & glucose == "0.01" & Replicate == "1")) 

### Set beginning and end time cut-offs
time_limits <- tibble(glucose = unique(data$glucose)) %>%
  mutate(timeclip_upper = case_when(glucose == "0.4" ~ 15,
                                    glucose == "0.2" ~ 10, 
                                    glucose == "0.1" ~ 7.5, 
                                    glucose == "0.05" ~ 6,
                                    glucose == "0.01" ~ 4,
                                    glucose == "0.008" ~ 3, 
                                    glucose == "0.006" ~ 3,
                                    glucose == "0.004" ~ 2.5,
                                    glucose == "0.002" ~ 2.,
                                    glucose == "0" ~ 10)) %>% 
  mutate(timeclip_lower = case_when(glucose == "0.4" ~ 1,
                                    glucose == "0.2" ~ 1, 
                                    glucose == "0.1" ~ 1, 
                                    glucose == "0.05" ~ 1,
                                    glucose == "0.01" ~ 1,
                                    glucose == "0.008" ~ 1, 
                                    glucose == "0.006" ~ 1,
                                    glucose == "0.004" ~ 1,
                                    glucose == "0.002" ~ 0.5,
                                    glucose == "0" ~ 0.5))

## Further filtering
blankdata <- blankdata %>% 
  left_join(time_limits, by = "glucose") %>% 
  # Plate 3 excluded (only contains 0.01% and below) as run on older Epoch plate reader with coarse OD resolution
  filter(Replicate != 3) %>% 
  # Lower glucose concentrations excluded as not analysed in ms
  filter(glucose == "0.01" | glucose == "0.1" | glucose == "0.2" | glucose == "0.4") 

## Plot dynamics
gluc_gradient_ts_ms <- blankdata %>% 
  filter(Time_h > timeclip_lower & Time_h < timeclip_upper) %>% 
  ggplot() +
  geom_line(aes(y = blankedOD, 
                x = Time_h, 
                group = interaction(Replicate,overnight_technical,glucose),
                col = as.factor(glucose)), alpha = 0.5) +
  theme_classic() +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_blank()) +
  labs(y = "OD (600nm)", x = "Time (h)") +
  scale_color_manual(values = c("mediumpurple4", wes_palette("Zissou1")[c(1,5,3)])) +
  facet_wrap(~ as.factor(glucose), scales = "free_x", ncol = 4) #+ scale_y_log10()

gluc_gradient_ts_ms

### Obtain PGR for each independent well
pgr_dat <- blankdata %>% 
  group_by(Replicate, glucose, overnight, technical) %>%
  arrange(Time_min, .by_group = TRUE) %>%
  mutate(g = log(lead(blankedOD)/lag(blankedOD)) / (lead(Time_h)-lag(Time_h))) %>%
  ungroup() %>% 
  filter(blankedOD != 0) %>% 
  filter(Time_h > timeclip_lower & Time_h < timeclip_upper) 

p_pgr <- ggplot(pgr_dat, aes(x = blankedOD, y = g)) +
  geom_point(alpha = 0.2, size = 0.5) +
  stat_smooth(aes(group = interaction(overnight, technical, Replicate, glucose),
                  colour = glucose), 
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
  scale_y_continuous(breaks = seq(-0.2, 1.0, 0.2))  +
  scale_color_manual(values = c("mediumpurple4", wes_palette("Zissou1")[c(1,5,3)]),
                     guide = "none") +
  coord_cartesian(ylim = c(-0.25, 1.0)) + facet_wrap(~glucose, scales = "free_x", ncol = 4)

p_pgr


### Obtain mean pgr in density bins
mean_by_N <- pgr_dat %>%  
  group_by(glucose) %>%
  mutate(bin = cut(blankedOD, breaks = (seq((min(blankedOD)), (max(blankedOD)), length.out = 60)))) %>%
  group_by(bin, glucose) %>%
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
  filter(!is.na(bin)) %>% 
  mutate(glucose = fct_relevel(glucose, "0.01", "0.1", "0.2", "0.4")) 

pgr_dat <- pgr_dat %>% 
  mutate(glucose = factor(case_when(glucose == "0.01" ~ "0.01% glucose",
                                    glucose == "0.1" ~ "0.1% glucose",
                                    glucose == "0.2" ~ "0.2% glucose",
                                    glucose == "0.4" ~ "0.4% glucose")))

pgr_mean <- mean_by_N %>%
  ggplot(aes(y = g_median, x = blankedOD_mean, col = glucose, ymin = g_25, ymax = g_75)) +
  geom_errorbar(alpha = 0.5, 
                linewidth = 0.2, 
                show.legend = FALSE) +
  geom_point(alpha = 0.5, size = 1) + 
  scale_color_manual(values = c('mediumpurple4',wes_palette("Zissou1")[c(1,5,3)])) +
  theme_classic() + 
  coord_cartesian(ylim = c(-0.05, 0.7)) +
  labs(x = "OD (600nm)", 
       y = expression("per capita growth (h"^-1*")"), 
       colour = "glucose (%)") +
  theme(legend.position = c(1.02, 1.03), 
        legend.justification = c(1, 1),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.key.height = unit(0.3, "cm"))

pgr_mean
