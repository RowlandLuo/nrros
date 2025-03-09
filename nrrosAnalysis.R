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
library(tidyverse)
volume <- read_csv("nrros_vol_8region_combined.csv")
str(volume)
#Remove excess rows.
volume <- volume %>%
  filter(!row_number() %in% c(1, 2, 3, 4))
volume
#Add an another column specifying genotypes
volume <- volume %>%
  mutate(ID = 1:38)
volume <- volume %>%
  mutate(Genotype = case_when(
    ID <= 15 ~ "HET",
    ID <= 28 & ID > 15  ~ "HOM",
    ID > 28 ~ "WT"
  ))
volume <- volume %>%
  mutate(Genotype = factor(Genotype,
                           levels = c("WT", "HET", "HOM")))
#Change the column into numeric variables
volume$Telencephalon <- as.numeric(volume$Telencephalon)
volume$Tectum <- as.numeric(volume$Tectum)
volume$Thalamus <- as.numeric(volume$Thalamus)
volume$Hypothalamus <- as.numeric(volume$Hypothalamus)
volume$Cerebellum <- as.numeric(volume$Cerebellum)
volume$Hindbrain <- as.numeric(volume$Hindbrain)
volume$Habenula <- as.numeric(volume$Habenula)
volume$`Posterior-Tuberculum` <- as.numeric(volume$`Posterior-Tuberculum`)
#Outliers?
volume_het <- volume %>%
  filter(Genotype == "HET")
boxplot(volume_het$Telencephalon)
boxplot(volume_het$Tectum)
boxplot(volume_het$Thalamus)
boxplot(volume_het$Hypothalamus)
volume_hom <- volume %>%
  filter(Genotype == "HOM")
boxplot(volume_hom$Telencephalon)
boxplot(volume_hom$Tectum)
boxplot(volume_hom$Thalamus)
boxplot(volume_hom$Hypothalamus)
volume_wt <- volume %>%
  filter(Genotype == "WT")
boxplot(volume_wt$Telencephalon)
boxplot(volume_wt$Tectum)
boxplot(volume_wt$Thalamus)
#All brains should be okay.
#Now perform a one-way MANOVA
manova_model2 <- manova(cbind(Telencephalon, Tectum, Thalamus, Hypothalamus, Cerebellum, Hindbrain,Habenula, `Posterior-Tuberculum`) ~ Genotype, data = volume)
summary(manova_model2, test = "Pillai")
summary.aov(manova_model2)
#Now perform posthoc on each significant brain regions
aov1 <- aov(Telencephalon ~ Genotype, data = volume)
TukeyHSD(aov1, "Genotype")
aov2 <- aov(Tectum ~ Genotype, data = volume)
TukeyHSD(aov2, "Genotype")
aov3 <- aov(Hypothalamus ~ Genotype, data = volume)
TukeyHSD(aov3, "Genotype")
aov4 <- aov(Thalamus ~ Genotype, data = volume)
TukeyHSD(aov4, "Genotype")
aov5 <- aov(Hindbrain ~ Genotype, data = volume)
TukeyHSD(aov5, "Genotype")
aov6 <- aov(Habenula ~ Genotype, data = volume)
TukeyHSD(aov6, "Genotype")
aov7 <- aov(Cerebellum ~ Genotype, data = volume)
TukeyHSD(aov7, "Genotype")
aov8 <- aov(`Posterior-Tuberculum` ~ Genotype, data = volume)
TukeyHSD(aov8, "Genotype")
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
summary_stats_vol %>% print(n = Inf)
#Now plot
nrros_volume_plot <- ggplot(volume_long, aes(x = Genotype, y = Value, fill = Genotype)) +
  geom_boxplot() +
  facet_wrap(~ Variable, scales = "fixed") +
  scale_y_continuous(limits = c(-250, 80)) +
  labs(x = "Treatment", y = "Normalized Volume")
nrros_volume_plot