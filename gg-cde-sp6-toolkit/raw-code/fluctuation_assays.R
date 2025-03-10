library(ggplot2)
library(dplyr)
library(readxl)
library(tidyr)
library(viridis)
library(agricolae)
library(tibble)
library(DescTools)

#load data and line info for easier plot labelling
#make sure to edit file names to the right one
#make sure youre on the right working directory

line_info <- read_excel("emutasp6_2024_09.xlsx", sheet = "Lines")
data <- read_excel("emutasp6_2024_09.xlsx", sheet = "Data")

data$passage <- as.factor(data$passage)
data$replicate <- as.factor(data$replicate)
data$cfu_sel <- as.numeric(data$cfu_sel)

data <- data %>%
  left_join(line_info %>% select(Line, `Line name`), by = c("population" = "Line")) %>%
  rename(line = `Line name`)

data <- data %>% mutate(non_sel_total_cfu = cfu_non_sel / df_non_sel)


data <- data %>%
  mutate(
    log_Suppressor_Frequency = if_else(
      !is.na(cfu_sel) & !is.na(df) & !is.na(non_sel_total_cfu) & 
        cfu_sel > 0 & df > 0 & non_sel_total_cfu > 0,
      log10((cfu_sel / df) / non_sel_total_cfu),
      -9
    )
  )

plotting_data <- select(data, population, replicate, passage, line, log_Suppressor_Frequency)


plotting_data <- plotting_data %>%
  rename(Suppressor_Frequency = log_Suppressor_Frequency)


################################
#at this point you can apply any filters to showcase the data
################################

plotting_data <- plotting_data %>%
  mutate(jitter_condition = ifelse(Suppressor_Frequency == -9, "jitter", "no_jitter"))


supp_freq_combined <- ggplot(plotting_data, aes(x = passage, y = Suppressor_Frequency, fill = line)) +
  
  geom_point(data = subset(plotting_data, jitter_condition == "jitter"), 
             
             position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.22), 
             size = 3, shape = 21, stroke = 1.5) +  # Apply jitter to points meeting the condition
  
  geom_point(data = subset(plotting_data, jitter_condition == "no_jitter"),
             position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.22),
             size = 3, shape = 21, stroke = 1.5) +  # No jitter for other points
  scale_color_viridis_d(option = "D") +
  scale_fill_viridis_d() + 
  labs(x = "Passage",
       y = "log(Suppressor Frequency)",
       color = "Line") +
  scale_y_continuous(limits = c(-9.5, 1)) + #make sure you always modify the limites to make sure that you are plotting appropriately
  scale_x_discrete("Passage") +
  
  theme_minimal() +
  theme(
    text = element_text(size = 15),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = "top", 
    legend.title = element_blank(),
    axis.text = element_text(color = "black", size = 14),  # Increase axis text size
    axis.title = element_text(color = "black", size = 24),  # Increase axis title size
    panel.background = element_blank(), 
    axis.ticks.length = unit(0.21, "cm"),
    axis.ticks = element_line(color = "black", linewidth = 0.9),
    axis.text.x = element_text(hjust = 0.5, vjust = 0.5),
    guides(shape = "none"),
    strip.text = element_text(size = 12, ),
    legend.text = element_text(size = 12),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.75)
  ) +
  guides(color = guide_legend(override.aes = list(size = 6))) + # Adjust legend point size if necessary
  facet_wrap(~line, nrow = 2, labeller = label_wrap_gen(width = 24))
print(supp_freq_combined)


#note that for figures in paper select passages were shown




