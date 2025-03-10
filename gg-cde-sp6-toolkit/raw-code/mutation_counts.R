library(ggplot2)
library(dplyr)
library(readxl)
library(tidyr)
library(viridis)
library(agricolae)
library(tibble)
library(DescTools)

#make sure youre on the right working directory
mutation_counts <- read_excel("04_emutasp6_2025_02_5alpha_mutation_counts.xlsx")


#filter your file here to decide what and what not to plot
#mutation_counts <- mutation_counts %>% filter(population != 'XXX')

#determine if data is normally distributed 
mut_count_anova <- aov(mutations ~ population, data = mutation_counts)
summary(mut_count_anova)
shapiro.test(residuals(mut_count_anova))
#if p value of shapiro test < 0.05, it is not normal
mut_count_tukey <- hsd_result <- HSD.test(mut_count_anova, trt = "population", alpha = 0.05, group = TRUE)
letter_groups <- mut_count_tukey$groups
#add the letters to the plot
mutation_counts$group <- hsd_result$groups[as.character(mutation_counts$population), "groups"]


#if it is not, then use a non-parametric test
kruskal_result <- kruskal.test(mutations ~ population, data = mutation_counts)
kruskal_result
dunn_result <- dunnTest(mutations ~ population, data = mutation_counts, method = 'bonferroni')
write.table(dunn_result[2], "clipboard", sep = "\t", row.names = FALSE, col.names = TRUE)

#make it look nicer
mutation_counts$population <- str_wrap(mutation_counts$population, width = 7)

#you can edit to change the order of the graph
#mutation_counts$population <- factor(mutation_counts$population, levels = c("Empty vector", "No deaminase", "eMutaSP6 (ProT7)", "eMutaSP6", "eMutaT7 (reassembled)", "eMutaT7_transition"), ordered = TRUE)


jitterplot <-ggplot(data = mutation_counts, aes(x = population, y = mutations, fill = population)) +
  #you can add a boxplot if you like
  #geom_boxplot(staplewidth = 0.5, notch = FALSE, outlier.size = 4, size = 0.8) +
  geom_point(size = 3, position = position_jitterdodge(jitter.width = 0.2), 
             stroke = 1.5, shape = 21) +
  #geom_text(aes(label = group, y = max(mutations) + 0.2), vjust = -0.5, size = 8) +
  theme_minimal() +
  labs(x = "Population", y = "Mutation count") +
  scale_fill_viridis_d() +
  scale_color_viridis_d(name = "population") + 
  scale_y_continuous() + 
  theme(text = element_text(size = 16),
        plot.title = element_text(angle = 0, hjust =  0.5, vjust = 1),
        plot.title.position = "plot",
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "right", 
        axis.text.y = element_text(color = "black", size = 14),  # Increase axis text size
        axis.title = element_text(color = "black", size = 24),  # Increase axis title size
        axis.title.y = element_text(color = "black", size = 24),
        panel.background = element_blank(), 
        axis.ticks.length = unit(0.15, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.5),
        axis.text.x = element_text(hjust = 0.5, size = 14),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.75)
        
  )
print(jitterplot)


#an alternative version is this one, only if you have multi passage data
plot_by_passage <-ggplot(data = mutation_counts, aes(x = passage, y = mutations, fill = population)) +
  #you can add a boxplot if you like
  #geom_boxplot(staplewidth = 0.5, notch = FALSE, outlier.size = 4, size = 0.8) +
  geom_point(size = 3, position = position_jitterdodge(jitter.width = 1, dodge.width = 0.3), 
             
             stroke = 1.5, shape = 21) +
  #geom_text(aes(label = group, y = max(mutations) + 0.2), vjust = -0.5, size = 8) +
  theme_minimal() +
  labs(x = "Passage", y = "Mutation count") +
  scale_fill_viridis_d() +
  scale_color_viridis_d(name = "population") + 
  scale_y_continuous() + 
  scale_x_continuous(breaks = c(0,8)) +
  theme(text = element_text(size = 16),
        plot.title = element_text(angle = 0, hjust =  0.5, vjust = 1),
        plot.title.position = "plot",
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "none", #optionally change here 
        axis.text.y = element_text(color = "black", size = 14),  # Increase axis text size
        axis.title = element_text(color = "black", size = 24),  # Increase axis title size
        axis.title.y = element_text(color = "black", size = 24),
        panel.background = element_blank(), 
        axis.ticks.length = unit(0.15, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.5),
        axis.text.x = element_text(hjust = 0.5, size = 14),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.75)
        
  ) + 
  facet_wrap(~population, nrow = 2)
print(plot_by_passage)

