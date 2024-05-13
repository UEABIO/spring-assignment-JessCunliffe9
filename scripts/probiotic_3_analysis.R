# Git Token ----
# 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

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
library(rmarkdown)
library(knitr)
library(patchwork)


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

# Removing Outlier
probiotic_3 <- probiotic_2[probiotic_2$before_abund <= 900, ]

# Hypothesis ----
#H1 = there is an reduction in abundance of the pathogenic bacterium Ruminococcus 
#gnavus in human stool samples when there is intervention with a probiotic L. rhamnosus GG (LGG).

#H0 = there is no effect on abundance of the pathogenic bacterium Ruminococcus 
#gnavus in human stool samples when there is intervention with a probiotic L. rhamnosus GG (LGG).

# Value creation ----
# Setting universal colours for graphs
colours <- c("deeppink3", "royalblue")


# Data Visualisation checks probiotc_3 ----
# Histogram to look at new data set distribution across a three variables 
probiotic_3 %>% 
  ggplot()+
  geom_histogram(aes(x=before_abund),
                 bins=20)

probiotic_3 %>% 
  ggplot()+
  geom_histogram(aes(x=after_abund),
                 bins=10)

probiotic_3 %>% 
  ggplot()+
  geom_histogram(aes(x=difference),
                 bins=10)

# Histogram with variables
probiotic_3 %>% 
  ggplot(aes(x=difference,
             fill=gender))+
  geom_histogram(alpha=0.6,
                 bins=15)+
  facet_wrap(~group,
             ncol=1)

# jitter plot of rcg abundance difference Vs Gender Vs Group
ggplot(data = probiotic_3, aes(x = group, y = difference)) +
  geom_jitter(aes(color = gender),
              width = 0.1, 
              alpha = 0.7, 
              show.legend = TRUE)

# Boxplot of rcg_bund Vs group Vs gender
ggplot(data = probiotic_3, aes(x = group, y = difference)) +
  geom_boxplot(aes(fill = gender), 
               alpha = 0.2, 
               width = 0.5, 
               outlier.shape=NA)+
  geom_jitter(aes(),
              width=0.2)+
  theme_light()

# Creating a summary of comparisons 
GGally::ggpairs(probiotic_3,
                aes(colour = group))

# Looking at density across groups and genders 
probiotic_3 %>% 
  ggplot(aes(x = difference, y = group, fill = gender))+
  geom_density_ridges(alpha = 0.5)+
  scale_fill_manual(values = colours)+
  theme_minimal()+
  theme_light()
# No male data for LGG???

# Exploring the comparison of female and male reactions to LGG and Placebo 
probiotic_3 %>% 
  ggplot(aes(x=group, y = difference, group = gender, colour = gender))+
  geom_point( alpha = 0.6)+
  geom_smooth(method = "lm",
              se = FALSE)+
  scale_colour_manual(values = colours)+
  theme_minimal()
# Observations:
# LGG males exhibit larger increase in abundance than placebo males
# LGG males echibit an increase in abundance 
# LGG females exhibit a decrease in abundance 
# Placebo females exhibit small increase in abundance 
# In females for both groups seem to have very little difference 
# In males abundance seems to increase in both groups 

# Data Analysis probiotic_3 ----

# creating a summary of probiotic_3
summary(probiotic_3)

#Looking at variables counts 
probiotic_3 %>% 
  group_by(gender) %>% 
  summarise(n = n())

probiotic_3 %>% 
  group_by(group) %>% 
  summarise(n = n())

probiotic_3 %>% 
  group_by(group, gender) %>% 
  summarise(n = n())

#Exploring means all variables as seen in figures 
probiotic_3_summary <- probiotic_3 %>% 
  group_by(group,gender) %>% 
  summarise(mean_difference_in_rcg_abund=mean(difference),
            sd=sd(difference))
probiotic_3_summary
# negative effect in females, large increase in males 

probiotic_3_summary_group <- probiotic_3 %>% 
  group_by(group) %>% 
  summarise(mean_difference_in_rcg_abund=mean(difference),
            sd=sd(difference))
probiotic_3_summary_group
# a negative difference from LGG = 37.4 to Placebo = 34.9 (-2.5)
# shows the very little difference in effect just between groups 

probiotic_3_summary_gender <- probiotic_3 %>% 
  group_by(gender) %>% 
  summarise(mean_difference_in_rcg_abund=mean(difference),
            sd=sd(difference))
probiotic_3_summary_gender
# A positive difference from F = 6.07 to M = 110 (+104)

# Making a neat table using Kable of each summary 
probiotic_3_summary %>% 
  kbl(caption="Summary statistics of gender and group compared with RCG abudance difference between before and after observations") %>% 
  kable_styling(bootstrap_options = "striped", full_width = T, position = "left")

probiotic_3_summary_group %>% 
  kbl(caption="Summary statistics of treatment group compared with RCG abudance difference between before and after observations") %>% 
  kable_styling(bootstrap_options = "striped", full_width = T, position = "left")

probiotic_3_summary_gender %>% 
  kbl(caption="Summary statistics of gender compared with RCG abudance difference between before and after observations") %>% 
  kable_styling(bootstrap_options = "striped", full_width = T, position = "left")

# Data Statistical analysis probiotic_3 ----
model_3 <- lm(difference ~ gender + group,
              data = probiotic_3)
model_3

# Linear model of coefficients
model_3_summary <- summary(model_3)
model_3_summary <- tidy(model_3_summary)
# Print the tidy summary table
kable(model_3_summary, format = "html", digits = 2) %>%
  kable_styling(bootstrap_options = "striped", full_width = FALSE)

# Looking at residuals 
par(mfrow = c(2, 2))
plot(model_3)

# Breusch Pagan test for normality
lmtest::bptest(model_3)
#  the p-value (0.1997) is greater than the typical significance 
#level of 0.05, we fail to reject the null hypothesis. This suggests 
#that there is not enough evidence to conclude that heteroscedasticity 
#is present in the model at the 0.05 significance level.
# Heteroscedasticity occurs when the variability of the errors (or residuals) 
#in a regression model is not constant across all levels

# qqplot with confidence intervals
car::qqPlot(model_3)

# shapiro wilk test for homoscedasticity
shapiro.test(residuals(model_3))
# The p-value (0.9737) is more than the typical significance level of 0.05, 
#we accept the null hypothesis. This suggests that there is sufficient evidence 
#to conclude that the residuals do follow a normal distribution.
# It measures the degree of departure from normality. For the Shapiro-Wilk test, 
#the test statistic ranges from 0 to 1, where 1 indicates perfect normality and 
#values closer to 0 indicate departure from normality. In your result, W = 0.87613.

# performance check model_3 
performance::check_model(model_3, detrend = F)
# No problem with colinearity 

# Looking further into Outliers
performance::check_model(model_3, check="outliers")

# Looking further into cooks distance 
plot(model_3, which=c(4,4))
# No problem with any outliers 

# Model_3 summaries
broom::tidy(model_3) 

broom::glance(model_3)

#graph of the estimated mean difference with an approx 95% 
GGally::ggcoef_model(model_3,
                     show_p_values=FALSE, 
                     conf.level=0.95)
# How to interpret????
# difference within males and placebo?? 
GGally::ggcoef_model(model_3,
                     show_p_values=FALSE, 
                     conf.level=0.99)
# Higher level of confidence = All null hypothesis acceptance 

# Using eemeans to compare means
probiotic_3a_means <- emmeans::emmeans(model_3, specs = ~ group)
probiotic_3a_means
# creating a visual output from this eemeans data
probiotic_3a_means %>% 
  as_tibble() %>% 
  ggplot(aes(x=group, 
             y=emmean))+
  geom_pointrange(aes(
    ymin=lower.CL, 
    ymax=upper.CL))

probiotic_3b_means <- emmeans::emmeans(model_3, specs = ~ gender)
probiotic_3b_means

# creating a visual output from this eemeans data
probiotic_3b_means %>% 
  as_tibble() %>% 
  ggplot(aes(x=gender, 
             y=emmean))+
  geom_pointrange(aes(
    ymin=lower.CL, 
    ymax=upper.CL))

# Getting the F value
drop1(model_3, test = "F")

# making a model 4 to look at data comparing just difference and group 
model_4 <- lm(difference ~ group,
              data = probiotic_3)
model_4
summary(model_4)
drop1(model_4, test = "F")

# making a model 5 to look at data comparing just difference and gender 
model_5 <- lm(difference ~ gender,
              data = probiotic_3)
model_5
summary(model_5)
drop1(model_5, test = "F")


# Data Visualisation Analysis probiotic_3 ----
# Visual of the difference in mean abundance difference between groups 
probiotic_3 %>% 
  ggplot(aes(x=group, 
             y=difference))+
  geom_boxplot(width=0.8, 
              pch=21, 
              aes(fill=gender,
                  alpha=0.5))+
  geom_jitter()+
  theme_linedraw()+
  geom_segment(aes(x=1, xend=2, y=37.43, yend=37.43-2.5), linetype="dashed")+ # The averages of means between the groups 
  stat_summary(fun=mean, geom="crossbar", width=0.2)+
  labs(x = "Treatment Type", y = "The difference in abundance of rcg between before and after observations", fill = "Gender of Subject") 

# Visual of the difference in mean abundance difference between genders 
probiotic_3 %>% 
  ggplot(aes(x=gender, 
             y=difference)) +
  geom_boxplot(width=0.6, 
               pch=21, 
               aes(fill=group),
               alpha=0.8) +  
  geom_jitter(width=0.2,alpha=0.5) +  
  theme_bw() +
  geom_segment(aes(x=1, xend=2, y=6.067, yend=6.067+103.933), linetype="dashed") +
  stat_summary(fun=mean, geom="crossbar", width=0.2)+
  labs(x = "Gender of the Subject", 
       y = "The abundance difference of rcg between before and after observations", 
       fill = "Treatment Group")


# Graph/Table creation for report ----

# Define the annotation text
annotation_text <- "A comparison line between the mean difference in Ruminococcus gnavus abundance of LGG and Placebo treatment"
# Wrap the annotation text
wrapped_text <- str_wrap(annotation_text, width = 30)  # Adjust the width as needed

plot_1 <- probiotic_3 %>% 
  ggplot(aes(x=group, y = difference, group = gender, colour = gender))+
  geom_point(alpha = 0.8, size=2)+
  theme(axis.text = element_text(size = 20))+
  geom_smooth(method = "lm", se = FALSE)+
  scale_colour_manual(values = colours)+
  theme_minimal()+
  geom_segment(aes(x=1, xend=2, y=37.43, yend=37.43-2.5), linetype="dashed",lwd = 0.5,colour="red")+
  stat_summary(geom = "point",fun = mean,size = 3)+
  scale_y_continuous(limits = c(-200, 200))+
  labs(x = "Treatment Type", 
       y = expression("The abundance difference of " * italic("R. gnavus") * " between before and after observations"), 
       color = "Gender of Subject") +
  ggtitle(expression("The Mean Abundance Differences of "* italic("Ruminococcus gnavus") *" Between Treatment Groups and Gender"))+
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")+
  annotate("text", x = 2.3, y = 36, label = wrapped_text, color = "red", size = 3)+
  theme(legend.position = "bottom", legend.direction = "horizontal", legend.text = element_text(size = 15))+
  guides(color = guide_legend(override.aes = list(shape = 15, size=8)))
plot_1
  
# LGG males exhibit larger increase in abundance than placebo males
# LGG males exhibit an increase in abundance 
# LGG females exhibit a decrease in abundance 
# Placebo females exhibit small increase in abundance 
# In females for both groups seem to have very little difference from no difference
# In males abundance seems to increase in both groups 
# The difference between the mean abundances of LGG and placebo is null 

# Boxplot comparing gender and group and difference 
plot_2<-probiotic_3 %>% 
  ggplot(aes(x=group, 
             y=difference))+
  geom_boxplot(width=0.8, 
               aes(fill=gender),alpha=0.8)+
  scale_fill_manual(values = colours)+
  geom_jitter(shape = 1, width=0.1)+
  theme_minimal()+
  scale_y_continuous(limits = c(-200, 200))+
  labs(x = "Treatment Type", 
       y = NULL,
       fill = NULL)+
  #ggtitle(expression("The abundance differences of "* italic("Ruminococcus gnavus") *" between treatment groups and gender variances"))+
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")+
  theme(legend.position = "none")
plot_2

#(1)# Creating a patchwork of the two graphs in report 
(plot_1+plot_2)+
  plot_layout(guides = "collect") 

#(2)# Table of Participant quantities
probiotic_3_participant_summary <- probiotic_3 %>% 
  group_by(group, gender) %>% 
  summarise(quantity = n())
probiotic_3_participant_summary %>% 
  kbl(caption="Table of Participant quantities within each Gender and Treatment Group") %>% 
  kable_styling(bootstrap_options = "striped", full_width = T, position = "left")

#(3)# Table the average differences of both genders and each group showing either a positive or negative value from the original observation
probiotic_3_summary <- probiotic_3 %>% 
  group_by(group,gender) %>% 
  summarise(mean_difference_in_rcg_abund = mean(difference),
            sd = sd(difference),
            se = sd(difference) / sqrt(n()))
probiotic_3_summary %>% 
  kbl(caption="Summary statistics of gender and group compared with RCG abudance difference between before and after observations") %>% 
  kable_styling(bootstrap_options = "striped", full_width = T, position = "left")
# negative effect in females, large increase in males 

#(4)# Table comparing the average differences between the two treatment groups
probiotic_3_summary_group <- probiotic_3 %>% 
  group_by(group) %>% 
  summarise(
    mean_difference_in_rcg_abund = mean(difference),
    sd = sd(difference),
    se = sd(difference) / sqrt(n()))
probiotic_3_summary_group <- probiotic_3_summary_group %>% 
  kbl(caption="Summary statistics of treatment group compared with RCG abudance difference between before and after observations") %>% 
  kable_styling(bootstrap_options = "striped", full_width = T, position = "left")
# a negative difference from LGG = 37.4 to Placebo = 34.9 (-2.5)
# shows the very little difference in effect between treatment groups 

#(5)# Table comparing just genders and the mean differences of abundances between groups 
probiotic_3_summary_gender <- probiotic_3 %>% 
  group_by(gender) %>% 
  summarise(mean_difference_in_rcg_abund = mean(difference),
            sd = sd(difference),
            se = sd(difference) / sqrt(n()))
probiotic_3_summary_gender %>% 
  kbl(caption="Summary statistics of gender compared with RCG abudance difference between before and after observations") %>% 
  kable_styling(bootstrap_options = "striped", full_width = T, position = "left")

# A positive difference from F = 6.07 to M = 110 (+104)

#(6)# Linear model of coefficients
model_3_summary <- summary(model_3)
model_3_summary <- tidy(model_3_summary)
kable(model_3_summary, format = "html", digits = 2) %>%
  kable_styling(bootstrap_options = "striped", full_width = FALSE)
#The intercept coefficient is 7.733, but it is not statistically significant (p = 0.800), 
#This means that when both gender and group are at their reference levels (e.g., female and LGG group), 
#the average value of difference is expected to be 7.733, but this value is not significantly different from zero.

# The coefficient for genderM is 103.933, and it is statistically significant (p = 0.010).
# This suggests that, holding other variables constant, males (genderM) have, on average, 
#a difference value that is 103.933 units higher than females (genderF)

# The coefficient for groupPlacebo is -2.500, but it is not statistically significant (p = 0.943).
# This implies that there is no statistically significant difference in the difference values between the Placebo group and theLGG group

