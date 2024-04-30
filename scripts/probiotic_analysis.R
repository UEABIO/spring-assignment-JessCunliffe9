# Git Token ----
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

# Clearing environment ----
rm(list=ls())

# Loading Packages ----
library(tidyverse)
library(janitor)
library(rstatix)
library(performance)
library(here)
library(tidyr)
library(dplyr)
library(kableExtra)
library(GGally)
library(emmeans)

# Loading Data ----
probiotic <- read_csv(here("data", "probiotic.csv"))

# View Data ----
# Initial look at data
glimpse(probiotic)

# time	timepoint 1 = before, 2 = after
# gender	F or M, Female or Male
# subject	subject ID
# group	Placebo or LGG (probiotic L. rhamnosus GG)
# ruminococcus_gnavus_abund	Read count abundance of pathogenic bacteria


# Making data set public friendly: wider to remove sample ID, subject duplicates and 1 and 2 rows ----
#pivotwider

  probiotic_merged_2 <- probiotic %>% 
  group_by(subject) %>% 
  pivot_wider(names_from = "time", values_from = ruminococcus_gnavus_abund)%>% 
  select(-sample) %>% 
  rename(before_abund = `1`,
         after_abund = `2`) %>% 
  pivot_longer(cols = c(before_abund, after_abund), names_to = "Time", values_to = "Measurement") %>%
  group_by(subject, group, gender, Time) %>%
  summarise(Measurement = mean(Measurement, na.rm = TRUE)) %>%
  pivot_wider(names_from = Time, values_from = Measurement) %>% 
  mutate(difference = after_abund - before_abund)

# Removing the subject column as i feel it is not needed for analysis
  probiotic_merged_2_blind <-   probiotic_merged_2[, -which(names(probiotic_merged_2) == "subject")]
  
  # Having done this i can now easily compare before and after values
  # I am keeping in gender as i would like to explore the relationship involvement 
  
# Data clean and tidy ----

# Cleaning column names 
head(probiotic)
probiotic <- janitor::clean_names(probiotic)
# Changing ruminococcus_gnavus_ambundance to rcg_abund
probiotic <-  rename(probiotic,
                       "rcg_abund"="ruminococcus_gnavus_abund")

# Looking at number of rows
nrow(probiotic)
                       
# Checking for duplicates 
probiotic %>% 
  duplicated() %>% 
  sum()

# Checking for missing values
probiotic %>% 
  is.na() %>% 
  sum()

# check for typos by looking at distinct characters
probiotic %>% 
  distinct(group)

# Removing factors not required in analysis 
probiotic <- probiotic %>% select(-sample)

# Changing values 1 to before and 2 to after and creating column to compare before and after 
probiotic_mutated <- probiotic %>%
  mutate(time = ifelse(time == 1, 'Before', 'After'))


# Hypothesis ----
#H1 = there is an reduction in abundance of the pathogenic bacterium Ruminococcus 
#gnavus in human stool samples when there is intervention with a probiotic L. rhamnosus GG (LGG).

#H0 = there is no effect on abundance of the pathogenic bacterium Ruminococcus 
#gnavus in human stool samples when there is intervention with a probiotic L. rhamnosus GG (LGG).


# Data Visualisation ----
probiotic %>% 
  ggplot()+
  geom_histogram(aes(x=rcg_abund),
                 bins=10)
# identifies that the data distribution is into normally distributed when looking at rcg_abund combing F and M 

# jitter plot of rcg_bund Vs Gender Vs Group
ggplot(data = probiotic, aes(x = group, y = rcg_abund)) +
  geom_jitter(aes(color = gender),
              width = 0.1, 
              alpha = 0.7, 
              show.legend = FALSE)
 
# Boxplot of rcg_bund Vs group Vs gender
ggplot(data = probiotic, aes(x = group, y = rcg_abund)) +
  geom_boxplot(aes(fill = gender), 
               alpha = 0.2, 
               width = 0.5, 
               outlier.shape=NA)+
  geom_jitter(aes(),
              width=0.2)+
  theme(legend.position = "none")

# Histogram with variables
probiotic %>% 
  ggplot(aes(x=rcg_abund,
             fill=gender))+
  geom_histogram(alpha=0.6,
                 bins=40)+
  facet_wrap(~group,
             ncol=1)

# Data Analysis ----
summary(probiotic)
# Is the value 913 an outlier?

#Looking at variables
probiotic %>% 
  group_by(gender) %>% 
  summarise(n = n())

probiotic %>% 
  group_by(group) %>% 
  summarise(n = n())

#Exploring means and and SD of variables 
probiotic_summary <- probiotic %>% 
  group_by(group,gender) %>% 
  summarise(mean=mean(rcg_abund),
            sd=sd(rcg_abund))

# Making a neat table using Kable of this summary 
probiotic_summary %>% 
  kbl(caption="Summary statistics of gender and group compared with RCG abudance") %>% 
  kable_styling(bootstrap_options = "striped", full_width = T, position = "left")



# Data Statistical analysis ----

# QQ- Plot for normal distribution 
# If the residuals follow a normal distribution, they should meet to produce a 
#perfect diagonal line across the plot
ggplot(probiotic, aes(sample = rcg_abund))+
  stat_qq() + 
  stat_qq_line()

probiotic %>% 
  pull(rcg_abund) %>% 
  car::qqPlot()
# Identifies data points 13 and 41 as outliers 
# No action needed until residuals are analysed 

# Creating a linear model
lsmodel0 <- lm(formula = rcg_abund ~ 1, data = probiotic)

# Table of Coefficients
summary(lsmodel0)

# Comparing means 
lsmodel1 <- lm(height ~ type, data=darwin)


