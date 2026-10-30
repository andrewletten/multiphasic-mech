## Growth across different percentage Na2HPO4 of Na2HPO4 minimal media (only glucose held constant)

install.packages("RColorBrewer")
install.packages("paletteer")
install.packages("ggthemes")
install.packages("scales")

### Load required packages 
library(grow96)
library(tidyverse)
library(splines)
library(patchwork)
library(wesanderson)
library(RColorBrewer)
library(paletteer)
library(ggthemes)
library(scales)

### Load data with grow96
data <- processODData(specPath="specs/Na2HPO4", dataPath="data/Na2HPO4")

### Blank OD data
blankdata <- blankODs(data, method = "fromPlate", perTimePoint = TRUE) %>% 
  filter(WellType == "DATA") %>%
  mutate(overnight = str_split_i(strains, "_", 1),
         technical = str_split_i(strains, "_", 2)) 


### Set beginning and end time cut-offs
time_limits <- tibble(Na2HPO4 = unique(blankdata$Na2HPO4)) %>%
  mutate(timeclip_upper = case_when(Na2HPO4 == 200 ~ 10,
                                    Na2HPO4 == 150 ~ 10,
                                    Na2HPO4 == 100 ~ 10,
                                    Na2HPO4 == 50 ~ 12.5,
                                    Na2HPO4 == 25 ~ 12.5,
                                    Na2HPO4 == 0 ~ 6)) %>%
  mutate(timeclip_lower = case_when(Na2HPO4 == 200 ~ 1,
                                    Na2HPO4 == 150 ~ 1,
                                    Na2HPO4 == 100 ~ 1,
                                    Na2HPO4 == 50 ~ 1,
                                    Na2HPO4 == 25 ~ 1,
                                    Na2HPO4 == 0 ~ 1))

## Further filtering/joining
blankdata <- blankdata %>% 
  left_join(time_limits, by = "Na2HPO4") %>% 
  mutate(Na2HPO4 = factor(Na2HPO4, sort(as.numeric(unique(Na2HPO4)))))


## Plot dynamics
gluc_gradient_ts_ms <- blankdata %>% 
  filter(Time_h > timeclip_lower) %>% # & Time_h < timeclip_upper) %>% 
  ggplot() +
  geom_line(aes(y = blankedOD, 
                x = Time_h, 
                group = interaction(Replicate,strains,Na2HPO4),
                col = Na2HPO4), alpha = 0.5) +
  theme_classic() +
  theme(strip.background = element_blank()) +
  labs(y = "OD (600nm)", x = "Time (h)") +
  scale_color_paletteer_d("ggthemes::excel_Parallax") #+   facet_wrap(~ as.factor(Na2HPO4), scales = "free_x", ncol = 6) #+ scale_y_log10()

gluc_gradient_ts_ms

### Obtain PGR for each independent well
pgr_dat <- blankdata %>% 
  group_by(Replicate, Na2HPO4, overnight, technical) %>%
  arrange(Time_min, .by_group = TRUE) %>%
  mutate(g = log(lead(blankedOD)/lag(blankedOD)) / (lead(Time_h)-lag(Time_h))) %>%
  ungroup() %>% 
  filter(blankedOD != 0) #%>%   filter(Time_h > timeclip_lower & Time_h < timeclip_upper) 

p_pgr <- ggplot(pgr_dat, aes(x = blankedOD, y = g)) +
  geom_point(alpha = 0.2, size = 0.2) +
  stat_smooth(aes(group = interaction(overnight, technical, Replicate, Na2HPO4),
                  colour = Na2HPO4), 
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
  scale_x_continuous(labels = label_number(accuracy = 0.01))  +
  scale_color_paletteer_d("ggthemes::excel_Parallax") + 
  coord_cartesian(ylim = c(-0.25, 2.0)) #+  facet_wrap(~Na2HPO4, scales = "free_x", ncol = 6)

p_pgr


### Obtain mean pgr in density bins
mean_by_N <- pgr_dat %>%  
  group_by(Na2HPO4) %>%
  mutate(bin = cut(blankedOD, breaks = (seq((min(blankedOD)), (max(blankedOD)), length.out = 60)))) %>%
  group_by(bin, Na2HPO4) %>%
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
  ggplot(aes(y = g_median, x = blankedOD_mean, col = Na2HPO4, ymin = g_25, ymax = g_75)) +
  geom_errorbar(alpha = 0.5, 
                linewidth = 0.2, 
                show.legend = FALSE) +
  geom_point(alpha = 0.5, size = 1) + 
  scale_color_paletteer_d("RColorBrewer::Dark2") + 
  theme_classic() + 
  coord_cartesian(ylim = c(-0.05, 1.2)) +
  labs(x = "OD (600nm)", 
       y = expression("per capita growth (h"^-1*")"), 
       colour = "Na2HPO4 (%)") +
  theme(legend.position = c(1.02, 1.03), 
        legend.justification = c(1, 1),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.key.height = unit(0.3, "cm"))

pgr_mean


