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
library(corrplot)
library(lmtest)
library(ggridges)

# Loading Data ----
probiotic <- read_csv(here("data", "probiotic.csv"))

# Manipulating data ----

probiotic_2 <- probiotic %>% 
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
probiotic_2 <-   probiotic_2[, -which(names(probiotic_2) == "subject")]

# Hypothesis ----
#H1 = there is an reduction in abundance of the pathogenic bacterium Ruminococcus 
#gnavus in human stool samples when there is intervention with a probiotic L. rhamnosus GG (LGG).

#H0 = there is no effect on abundance of the pathogenic bacterium Ruminococcus 
#gnavus in human stool samples when there is intervention with a probiotic L. rhamnosus GG (LGG).

# Data check, clean and tidy ----
# Cleaning column names 
head(probiotic_2)

# Looking at number of rows
nrow(probiotic_2)

# Checking for duplicates 
probiotic_2 %>% 
  duplicated() %>% 
  sum()

# Checking for missing values
probiotic_2 %>% 
  is.na() %>% 
  sum()

# check for typos by looking at distinct characters
probiotic_2 %>% 
  distinct(group)

# Data Visualisation probiotc_2 ----
probiotic_2 %>% 
  ggplot()+
  geom_histogram(aes(x=before_abund),
                 bins=20)

probiotic_2 %>% 
  ggplot()+
  geom_histogram(aes(x=after_abund),
                 bins=10)

probiotic_2 %>% 
  ggplot()+
  geom_histogram(aes(x=difference),
                 bins=10)


# jitter plot of rcg abundance difference Vs Gender Vs Group
ggplot(data = probiotic_2, aes(x = group, y = difference)) +
  geom_jitter(aes(color = gender),
              width = 0.1, 
              alpha = 0.7, 
              show.legend = TRUE)

ggplot(data = probiotic_2, aes(x = after_abund, y = before_abund)) +
  geom_jitter(aes(color = group),
              width = 0.1, 
              alpha = 0.7, 
              show.legend = TRUE)



# Boxplot of rcg_bund Vs group Vs gender
ggplot(data = probiotic_2, aes(x = group, y = difference)) +
  geom_boxplot(aes(fill = gender), 
               alpha = 0.2, 
               width = 0.5, 
               outlier.shape=NA)+
  geom_jitter(aes(),
              width=0.2)+
  theme(legend.position = "none")

# Histogram with variables
probiotic_2 %>% 
  ggplot(aes(x=difference,
             fill=gender))+
  geom_histogram(alpha=0.6,
                 bins=15)+
  facet_wrap(~group,
             ncol=1)

# Creating a summary of comparisons 
GGally::ggpairs(probiotic_2,
                aes(colour = group))

# Creating Corr Plot
probiotic_2 %>% 
  group_by(difference,before_abund, after_abund) %>% 
  select(where(is.numeric)) %>% 
  cor() %>% 
  corrplot()

# Data Analysis probiotic_2 ----
summary(probiotic_2)
# Is the value 913 an outlier?

#Looking at variables counts 
probiotic_2 %>% 
  group_by(gender) %>% 
  summarise(n = n())

probiotic_2 %>% 
  group_by(group) %>% 
  summarise(n = n())

#Exploring means and and SD of variables 
probiotic_2_summary <- probiotic_2 %>% 
  group_by(group,gender) %>% 
  summarise(mean_difference_in_rcg_abund=mean(difference),
            sd=sd(difference))
probiotic_2_summary

# Making a neat table using Kable of this summary 
probiotic_2_summary %>% 
  kbl(caption="Summary statistics of gender and group compared with RCG abudance difference between before and after observations") %>% 
  kable_styling(bootstrap_options = "striped", full_width = T, position = "left")

# Data Statistical analysis probiotic_2 ----

# QQ- Plot for normal distribution 
# If the residuals follow a normal distribution, they should meet to produce a 
#perfect diagonal line across the plot
ggplot(probiotic_2, aes(sample = difference))+
  stat_qq() + 
  stat_qq_line()

ggplot(probiotic_2, aes(sample = before_abund))+
  stat_qq() + 
  stat_qq_line()

ggplot(probiotic_2, aes(sample = after_abund))+
  stat_qq() + 
  stat_qq_line()

probiotic_2 %>% 
  pull(difference) %>% 
  car::qqPlot()
# Identifies data points 14 as an outlier 
# No action needed until residuals are analysed 

# Creating a linear model
model_2 <- lm(difference ~ gender + group,
            data = probiotic_2)
model_2

# Table of Coefficients
summary(model_2)

# Looking at residuals 
par(mfrow = c(2, 2))
plot(model_2)

# Breusch Pagan test for normality
lmtest::bptest(model_2)
#  the p-value (0.1338) is greater than the typical significance 
#level of 0.05, we fail to reject the null hypothesis. This suggests 
#that there is not enough evidence to conclude that heteroscedasticity 
#is present in the model at the 0.05 significance level.
# Heteroscedasticity occurs when the variability of the errors (or residuals) 
#in a regression model is not constant across all levels

# qqplot with confidence intervals
car::qqPlot(model_2)

# shapiro wilk test for homoscedasticity
shapiro.test(residuals(model_2))
# The p-value (0.0102) is less than the typical significance level of 0.05, 
#we reject the null hypothesis. This suggests that there is sufficient evidence 
#to conclude that the residuals do not follow a normal distribution.
# It measures the degree of departure from normality. For the Shapiro-Wilk test, 
#the test statistic ranges from 0 to 1, where 1 indicates perfect normality and 
#values closer to 0 indicate departure from normality. In your result, W = 0.87613.

# performance check model_2 
performance::check_model(model_2, detrend = F)
# Looking further into Outliers
performance::check_model(model_2, check="outliers")
# Looking further into cooks distance 
plot(model_2, which=c(4,4))
# Confirms data points 14 as an outlier - will be removed and can be seen in model_3 and probiotic_3 script


