library(grow96)
library(tidyverse)
data <- processODData(specPath="specs", dataPath="data")
blankdata <- blankODs(data, method = "fromPlate", perTimePoint = TRUE) %>% 
  filter(!is.na(glucose)) %>% 
  filter(WellType == "DATA") %>% 
  mutate(overnight = str_split_i(overnight_technical, "_", 1),
         technical = str_split_i(overnight_technical, "_", 2))

ggplot(data = blankdata %>% filter(glucose == 0 & Time_h > 0), aes(y = OD, x = Time_h)) +
  geom_point(aes(col = overnight, shape = technical)) +
  facet_wrap(~ glucose)
