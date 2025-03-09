library(tidyverse)
activity <- read_csv("nrros_roi_8region_C0divC1_combined.csv")
head(activity)
#Remove excess rows
activity <- activity %>%
  filter(!row_number() %in% c(1:4))
activity
activity <- activity %>%
  mutate(Values = factor(Values))
#Add another column specifying genotypes
activity <- activity %>%
  mutate(ID = 1:38)
activity <- activity%>%
  mutate(Genotype = case_when(
    ID <= 15 ~ "HET",
    ID <= 28 & ID > 15  ~ "HOM",
    ID > 28 ~ "WT"
  ))
activity <- activity %>%
  mutate(Genotype = factor(Genotype,
                           levels = c("WT", "HET", "HOM")))
#Change the column into numeric variables
activity$Telencephalon <- as.numeric(activity$Telencephalon)
activity$Tectum <- as.numeric(activity$Tectum)
activity$Thalamus <- as.numeric(activity$Thalamus)
activity$Hypothalamus <- as.numeric(activity$Hypothalamus)
activity$Cerebellum <- as.numeric(activity$Cerebellum)
activity$Hindbrain <- as.numeric(activity$Hindbrain)
activity$Habenula <- as.numeric(activity$Habenula)
activity$`Posterior-Tuberculum` <- as.numeric(activity$`Posterior-Tuberculum`)
#Before performing one-way MANOVA, I will need to check for outliers first.
activity_het <- activity %>%
  filter(Genotype == "HET")
boxplot(activity_het$Telencephalon)
activity_hom <- activity %>%
  filter(Genotype == "HOM")
boxplot(activity_hom$Telencephalon)
boxplot(activity_hom$Tectum)
activity_wt <- activity %>%
  filter(Genotype == "WT")
boxplot(activity_wt$Telencephalon)
boxplot(activity_wt$Tectum)
boxplot(activity_wt$Thalamus)
#I think all brains should be okay.
#Now perform a one-way MANOVA
manova_model1 <- manova(cbind(Telencephalon, Tectum, Thalamus, Hypothalamus, Cerebellum, Hindbrain,Habenula, `Posterior-Tuberculum`) ~ Genotype, data = activity)
summary(manova_model1, test = "Pillai")
summary.aov(manova_model1)
#Now perform posthoc
aov1 <- aov(Telencephalon ~ Genotype, data = activity)
TukeyHSD(aov1, "Genotype")
aov2 <- aov(Tectum ~ Genotype, data = activity)
TukeyHSD(aov2, "Genotype")
aov3 <- aov(Thalamus ~ Genotype, data = activity)
TukeyHSD(aov3, "Genotype")
aov4 <- aov(Hypothalamus ~ Genotype, data = activity)
TukeyHSD(aov3, "Genotype")
aov5 <- aov(Cerebellum ~ Genotype, data = activity)
TukeyHSD(aov5, "Genotype")
aov6 <- aov(Hindbrain ~ Genotype, data = activity)
TukeyHSD(aov6, "Genotype")
aov7 <- aov(Habenula ~ Genotype, data = activity)
TukeyHSD(aov7, "Genotype")
aov8 <- aov(`Posterior-Tuberculum` ~ Genotype, data = activity)
TukeyHSD(aov8, "Genotype")
#Now, visualize the one-way MANOVA results.
#Convert the data into long format
activity_long <- activity %>%
  pivot_longer(
    cols = c("Telencephalon","Tectum","Thalamus","Hypothalamus","Cerebellum","Hindbrain","Habenula","Posterior-Tuberculum"),    
    names_to = "Variable",
    values_to = "Value"
  )
#Compute summary statistics
summary_stats <- activity_long %>%
  group_by(Genotype, Variable) %>%
  summarize(
    Mean = mean(Value, na.rm = TRUE),
    SE = sd(Value, na.rm = TRUE) / sqrt(n()),  # standard error
    .groups = "drop"
  )
summary_stats %>% print(n = Inf)
#Now plot
activity_plot <- ggplot(activity_long, aes(x = Genotype, y = Value, fill = Genotype)) +
  geom_boxplot() +
  facet_wrap(~ Variable, scales = "fixed") +
  scale_y_continuous(limits = c(0, 2.3)) +
  labs(x = "Treatment", y = "Normalized Activity (pERK/tERK)")
activity_plot

#Now let's move on to volume.
volume <- read_csv("foxp2_vol_8region_combined.csv")
str(volume)
#Remove excess rows.
volume <- volume %>%
  filter(!row_number() %in% c(1, 2, 3, 4))
volume
#Add an another column specifying genotypes
volume <- volume %>%
  mutate(ID = 1:31)
volume <- volume %>%
  mutate(Genotype = case_when(
    ID < 10 ~ "HET",
    ID < 23 & ID > 9  ~ "HOM",
    ID > 22 ~ "WT"
  ))
volume <- volume %>%
  mutate(Genotype = factor(Genotype))
#Change the column into numeric variables
volume$Telencephalon <- as.numeric(volume$Telencephalon)
volume$Tectum <- as.numeric(volume$Tectum)
volume$Thalamus <- as.numeric(volume$Thalamus)
volume$Hypothalamus <- as.numeric(volume$Hypothalamus)
volume$Cerebellum <- as.numeric(volume$Cerebellum)
volume$Hindbrain <- as.numeric(volume$Hindbrain)
volume$Habenula <- as.numeric(volume$Habenula)
volume$`Posterior-Tuberculum` <- as.numeric(volume$`Posterior-Tuberculum`)
#Now perform a one-way MANOVA
manova_model2 <- manova(cbind(Telencephalon, Tectum, Thalamus, Hypothalamus, Cerebellum, Hindbrain,Habenula, `Posterior-Tuberculum`) ~ Genotype, data = volume)
summary(manova_model2, test = "Pillai")
summary.aov(manova_model1)
#Now perform posthoc on each significant brain regions
aov2 <- aov(Telencephalon ~ Genotype, data = volume)
TukeyHSD(aov2, "Genotype")
aov3 <- aov(Tectum ~ Genotype, data = volume)
TukeyHSD(aov3, "Genotype")
aov4 <- aov(Hypothalamus ~ Genotype, data = volume)
TukeyHSD(aov4, "Genotype")
aov5 <- aov(Hindbrain ~ Genotype, data = volume)
TukeyHSD(aov5, "Genotype")
#Now, visualize the one-way MANOVA results.
#Convert the data into long format
volume_long <- volume %>%
  pivot_longer(
    cols = c("Telencephalon","Tectum","Thalamus","Hypothalamus","Cerebellum","Hindbrain","Habenula","Posterior-Tuberculum"),    
    names_to = "Variable",
    values_to = "Value"
  )
#Compute summary statistics
summary_stats_vol <- volume_long %>%
  group_by(Genotype, Variable) %>%
  summarize(
    Mean = mean(Value, na.rm = TRUE),
    SE = sd(Value, na.rm = TRUE) / sqrt(n()),  # standard error
    .groups = "drop"
  )
summary_stats %>% print(n = Inf)
#Now plot
volume_plot <- ggplot() +
  geom_col(data = summary_stats_vol, aes(x = Genotype, y = Mean, fill = Genotype),
           position = position_dodge(width = 0.9), width = 0.7) +
  geom_errorbar(data = summary_stats_vol, aes(x = Genotype, ymin = Mean - SE, ymax = Mean + SE),
                position = position_dodge(width = 0.9), width = 0.2) +
  geom_point(data = volume_long, aes(x = Genotype, y = Value),
             position = position_jitter(width = 0.2), alpha = 0.6) +
  facet_wrap(~ Variable, scales = "free_y") +
  theme_minimal() +
  labs(x = "Genotype", y = "Relative Volume") +
  theme(legend.position = "none")
volume_plot
#Add significance to the graph
max_vals_vol <- summary_stats_vol %>%
  group_by(Variable) %>%
  summarize(y_max = max(Mean+SE))
anno_df <- max_vals %>%
  mutate(
    x_start = 2,    
    x_end = 3,      
    bracket_y = y_max + 80 , 
    label_y = y_max + 80.22,   
    label = "*"              
  )
anno_df_single <- anno_df[anno_df$Variable == "Telencephalon", ]
volume_plot <- volume_plot + 
  geom_segment(data = anno_df_single,
               aes(x = x_start, xend = x_start, y = bracket_y - 0.1, yend = bracket_y),
               inherit.aes = FALSE) +
  geom_segment(data = anno_df_single,
               aes(x = x_start, xend = x_end, y = bracket_y, yend = bracket_y),
               inherit.aes = FALSE) +
  geom_segment(data = anno_df_single,
               aes(x = x_end, xend = x_end, y = bracket_y - 0.1, yend = bracket_y),
               inherit.aes = FALSE) +
  geom_text(data = anno_df_single,
            aes(x = (x_start + x_end)/2, y = label_y, label = label),
            size = 6, inherit.aes = FALSE)
volume_plot
anno_df <- max_vals %>%
  mutate(
    x_start = 2,    
    x_end = 3,      
    bracket_y = y_max + 40 , 
    label_y = y_max + 40.22,   
    label = "*"              
  )
anno_df_single <- anno_df[anno_df$Variable == "Tectum", ]
volume_plot <- volume_plot + 
  geom_segment(data = anno_df_single,
               aes(x = x_start, xend = x_start, y = bracket_y - 0.1, yend = bracket_y),
               inherit.aes = FALSE) +
  geom_segment(data = anno_df_single,
               aes(x = x_start, xend = x_end, y = bracket_y, yend = bracket_y),
               inherit.aes = FALSE) +
  geom_segment(data = anno_df_single,
               aes(x = x_end, xend = x_end, y = bracket_y - 0.1, yend = bracket_y),
               inherit.aes = FALSE) +
  geom_text(data = anno_df_single,
            aes(x = (x_start + x_end)/2, y = label_y, label = label),
            size = 6, inherit.aes = FALSE)
volume_plot
anno_df <- max_vals %>%
  mutate(
    x_start = 1,    
    x_end = 3,      
    bracket_y = y_max + 40 , 
    label_y = y_max + 40.22,   
    label = "*"              
  )
anno_df_single <- anno_df[anno_df$Variable == "Hypothalamus", ]
volume_plot <- volume_plot + 
  geom_segment(data = anno_df_single,
               aes(x = x_start, xend = x_start, y = bracket_y - 0.1, yend = bracket_y),
               inherit.aes = FALSE) +
  geom_segment(data = anno_df_single,
               aes(x = x_start, xend = x_end, y = bracket_y, yend = bracket_y),
               inherit.aes = FALSE) +
  geom_segment(data = anno_df_single,
               aes(x = x_end, xend = x_end, y = bracket_y - 0.1, yend = bracket_y),
               inherit.aes = FALSE) +
  geom_text(data = anno_df_single,
            aes(x = (x_start + x_end)/2, y = label_y, label = label),
            size = 6, inherit.aes = FALSE)
volume_plot
anno_df <- max_vals %>%
  mutate(
    x_start = 2,    
    x_end = 3,      
    bracket_y = y_max + 35 , 
    label_y = y_max + 35.22,   
    label = "**"              
  )
anno_df_single <- anno_df[anno_df$Variable == "Hypothalamus", ]
volume_plot <- volume_plot + 
  geom_segment(data = anno_df_single,
               aes(x = x_start, xend = x_start, y = bracket_y - 0.1, yend = bracket_y),
               inherit.aes = FALSE) +
  geom_segment(data = anno_df_single,
               aes(x = x_start, xend = x_end, y = bracket_y, yend = bracket_y),
               inherit.aes = FALSE) +
  geom_segment(data = anno_df_single,
               aes(x = x_end, xend = x_end, y = bracket_y - 0.1, yend = bracket_y),
               inherit.aes = FALSE) +
  geom_text(data = anno_df_single,
            aes(x = (x_start + x_end)/2, y = label_y, label = label),
            size = 6, inherit.aes = FALSE)
volume_plot
anno_df <- max_vals %>%
  mutate(
    x_start = 2,    
    x_end = 3,      
    bracket_y = y_max + 30 , 
    label_y = y_max + 30.22,   
    label = "**"              
  )
anno_df_single <- anno_df[anno_df$Variable == "Hindbrain", ]
volume_plot <- volume_plot + 
  geom_segment(data = anno_df_single,
               aes(x = x_start, xend = x_start, y = bracket_y - 0.1, yend = bracket_y),
               inherit.aes = FALSE) +
  geom_segment(data = anno_df_single,
               aes(x = x_start, xend = x_end, y = bracket_y, yend = bracket_y),
               inherit.aes = FALSE) +
  geom_segment(data = anno_df_single,
               aes(x = x_end, xend = x_end, y = bracket_y - 0.1, yend = bracket_y),
               inherit.aes = FALSE) +
  geom_text(data = anno_df_single,
            aes(x = (x_start + x_end)/2, y = label_y, label = label),
            size = 6, inherit.aes = FALSE)
volume_plot