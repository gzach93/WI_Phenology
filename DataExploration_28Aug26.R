#Libraries
library(tidyverse)
library(ggplot2)
library(readxl)
library(lubridate)
library(lme4)

#Data
pheno <- read.csv('~/Dropbox/Documents/Projects/WI_Phenology/DNR_Phenology.csv')

#Any Empty Call index observation becomes a zero 
pheno[pheno$Wood %in% c(''),'Wood'] <- 0
pheno[pheno$Boreal_Chorus %in% c(''),'Boreal_Chorus'] <- 0
pheno[pheno$Peeper %in% c(''),'Peeper'] <- 0
pheno[pheno$Leopard %in% c(''),'Leopard'] <- 0
pheno[pheno$Pickerel %in% c(''),'Pickerel'] <- 0
pheno[pheno$Am_Toad %in% c(''),'Am_Toad'] <- 0
pheno[pheno$Gray %in% c(''),'Gray'] <- 0
pheno[pheno$Copes_Gray %in% c(''),'Copes_Gray'] <- 0
pheno[pheno$B_Cricket %in% c(''),'B_Cricket'] <- 0
pheno[pheno$Mink %in% c(''),'Mink'] <- 0
pheno[pheno$Green %in% c(''),'Green'] <- 0
pheno[pheno$Bullfrog %in% c(''),'Bullfrog'] <- 0

#Format Date
pheno$Date <- as.Date(pheno$Date, format = "%m/%d/%y")
pheno$julian_day <- as.numeric(format(pheno$Date, "%j"))
pheno$year <- year(pheno$Date)
pheno$julian_day2 <- pheno$julian_day^2

#Change X to NA -- these are presence only and maybe useful for something else
pheno[pheno$Boreal_Chorus %in% c('X'),'Boreal_Chorus'] <- NA #13
pheno[pheno$Am_Toad %in% c('X'),'Am_Toad'] <- NA #27

#Change ? to NA -- no call index guessed
pheno[pheno$Boreal_Chorus %in% c('?'),'Boreal_Chorus'] <- NA #1
pheno[pheno$Am_Toad %in% c('?', '*'),'Am_Toad'] <- NA #5

#Assume numbers with questions are right
pheno[pheno$Boreal_Chorus %in% c('2?'),'Boreal_Chorus'] <- 2 #1
pheno[pheno$Boreal_Chorus %in% c('1?'),'Boreal_Chorus'] <- 1 #1

pheno[pheno$Am_Toad %in% c('2?'),'Am_Toad'] <- 2 #1

#Change call index to numeric
pheno$Boreal_Chorus <- as.numeric(pheno$Boreal_Chorus)
pheno$Am_Toad <- as.numeric(pheno$Am_Toad)

#How many years? and Sites per year?
sort(unique(pheno$year)) #44 years of data
length(unique(pheno$Site_No)) #197 Sites

sites_per_year <- pheno %>%
                  group_by(year) %>%
                  summarize(unique_sites = n_distinct(Site_No))

ggplot() + 
  geom_line(data = sites_per_year, aes(x = year , y = unique_sites)) +
  theme_classic() +
  xlab('Year') + ylab('Number of Unique Sites')


first_call_date <- pheno %>%
  group_by(Site_No, year) %>%
  arrange(Date) %>%          
  filter(Am_Toad > 0) %>%     
  dplyr::slice(1)                  

ggplot(data = first_call_date, aes(x = year, y = julian_day)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  theme_classic()

summary(lm(julian_day ~ year, data = first_call_date))





#############################################################################
### Simple Visuals of Chorus Frogs ###
#Base Data Plot
ggplot() +
  geom_jitter(data = pheno, 
              aes(x = julian_day, y = Boreal_Chorus, col = year), height = 0.1) +
  theme_classic() + ylab('Chorus Frog Call Index') + xlab('Julian Day')

### Simple Visuals of Toads ###
#Base Data Plot
ggplot() +
  geom_jitter(data = pheno, 
              aes(x = julian_day, y = Am_Toad, col = year), height = 0.1) +
  theme_classic() + ylab('American Toad Call Index') + xlab('Julian Day')

toad <- pheno %>%
  group_by(julian_day, year) %>%
  summarise(toad = mean(Am_Toad, na.rm = TRUE))

ggplot() +
  geom_jitter(data = toad, 
              aes(x = julian_day, y = toad, col = year), height = 0.1) +
  theme_classic() + ylab('American Toad Call Index') + xlab('Julian Day')

pheno$year_c <- pheno$year - mean(pheno$year, na.rm = TRUE)
pheno$julian_c <- pheno$julian_day - mean(pheno$julian_day, na.rm = TRUE)
pheno$Site_No <- as.factor(pheno$Site_No)

call.mod <- glmer(Am_Toad ~ poly(julian_c,2)*year_c + (1 | Site_No), 
                data = pheno[!(pheno$julian_c %in% c(NA)),], 
                family = 'poisson')
summary(call.mod)

new.data <- data.frame('julian_c' = rep(seq(-83, 124, by = 1), times = 46),
                       'year_c' = rep(unique(pheno$year_c), each = 208))

new.data <- cbind(new.data, data.frame(predict(call.mod, newdata = new.data, type = 'response', se.fit = TRUE, re.form = NA)))

ggplot() +
  geom_jitter(data = pheno, aes(x = julian_c, y = Am_Toad, col = year_c), height = 0.1) + #jitter overlapping points
  geom_line(data = new.data, aes(x = julian_c, y = fit, col = year_c, group = year_c)) +
  #geom_ribbon(data = new.data, aes(x = julian_c, ymin = fit - 2*se.fit, ymax = fit + 2*se.fit, fill = year_c, group = year_c),
   #           alpha = .2) +
  theme_classic() +
  ylab('Toad Call Intensity') + xlab('Julian Day') +
  coord_cartesian(ylim = c(0,3))

annual_summary <- new.data %>%
  filter(fit > 0.001) %>% #Change based on threshold -- what is above 0?
  group_by(year_c) %>%
  summarize(firstcallday = first(julian_c), # Replace 'Time' with your x-axis variable
            maxcallday = julian_c[which.max(fit)])


######################################################
library(ordinal)
pheno_test <- pheno[,c('Site_No', 'Am_Toad', 'julian_day', 'julian_day2', 'year')]
pheno_test <- na.omit(pheno_test)
pheno_test$year_c <- pheno_test$year - mean(pheno_test$year)

pheno$Am_Toad <- factor(pheno$Am_Toad, levels = 0:3, ordered = TRUE)
call.mod <- clm(Am_Toad ~  poly(julian_day,2) * year_c, 
                 data = pheno_test)
summary(call.mod)
#ctable <- coef(summary(call.mod))
#p_vals <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
nominal_test(call.mod)

newdat <- expand.grid(julian_day = seq(min(pheno_test$julian_day), max(pheno_test$julian_day), by = 1),
                      year_c = unique(pheno_test$year_c))
preds <- predict(call.mod, newdat, type = "probs")
newdat <- cbind(newdat, preds)
newdat_long <- newdat %>%
  pivot_longer(cols = as.character(0:3), 
               names_to = "Am_Toad_level", 
               values_to = "probability")

ggplot(newdat_long, aes(x = julian_day, y = probability, color = Am_Toad_level)) +
  geom_line(linewidth = 1) +
  facet_wrap(~ year_c) +
  labs(x = "Julian Day", y = "Predicted Probability", 
       color = "Calling\nIntensity",
       title = "Predicted Am. Toad Calling Intensity by Day and Year") +
  theme_minimal()

obs_summary <- pheno_test %>%
  mutate(day_bin = round(julian_day / 7) * 7) %>%  # weekly bins
  count(year, day_bin, Am_Toad) %>%
  group_by(year, day_bin) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

ggplot() +
  geom_line(data = newdat_long, 
            aes(x = julian_day, y = probability, color = Am_Toad_level),
            linewidth = 1) +
  geom_point(data = obs_summary,
             aes(x = day_bin, y = prop, color = Am_Toad), 
             alpha = 0.6, size = 2) +
  facet_wrap(~ year) +
  labs(x = "Julian Day", y = "Probability", 
       color = "Calling\nIntensity",
       title = "Model Predictions (lines) vs. Observed Proportions (points)") +
  theme_minimal()


library(VGAM)
call.mod2 <- vglm(Am_Toad ~ poly(julian_day, 2) * year_c,
                  family = cumulative(parallel = FALSE ~ poly(julian_day, 2)),
                  data = pheno_test)
summary(call.mod2)

newdat <- expand.grid(
  julian_day = seq(min(pheno_test$julian_day), max(pheno_test$julian_day), by = 1),
  year_c = sort(unique(pheno_test$year_c))
)
preds <- predict(call.mod2, newdata = newdat, type = 'response')

prob_matrix <- as.data.frame(preds)
colnames(prob_matrix) <- levels(pheno_test$Am_Toad)
newdat <- cbind(newdat, prob_matrix)

newdat_long <- newdat %>%
  pivot_longer(cols = as.character(0:3),
               names_to = "Am_Toad_level",
               values_to = "probability")

ggplot(newdat_long, aes(x = julian_day, y = probability, color = Am_Toad_level)) +
  geom_line(linewidth = 1) +
  facet_wrap(~ year_c) +
  labs(x = "Julian Day", y = "Predicted Probability",
       color = "Calling\nIntensity",
       title = "Predicted Am. Toad Calling Intensity (Partial Proportional Odds)") +
  theme_minimal()

obs_summary <- pheno_test %>%
  mutate(day_bin = round(julian_day / 7) * 7) %>%
  count(year, day_bin, Am_Toad) %>%
  group_by(year, day_bin) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  rename(Am_Toad_level = Am_Toad)

ggplot() +
  geom_line(data = newdat_long,
            aes(x = julian_day, y = probability, color = Am_Toad_level),
            linewidth = 1) +
  geom_point(data = obs_summary,
             aes(x = day_bin, y = prop, color = Am_Toad_level),
             alpha = 0.6, size = 2) +
  facet_wrap(~ year) +
  labs(x = "Julian Day", y = "Probability",
       color = "Calling\nIntensity",
       title = "Partial PO Model Predictions (lines) vs. Observed (points)") +
  theme_minimal()


pheno_test %>%
  filter(julian_day > 125, julian_day < 175) %>%  # adjust to your known peak window
  count(Am_Toad) %>%
  mutate(prop = n / sum(n))


#####################################################


library(mgcv)
pheno_test$Am_Toad <- as.integer(pheno_test$Am_Toad)

call.mod_gam <- gam(Am_Toad ~ s(julian_day, by = year_c) + s(year_c) + s(Site_No, bs = "re"),
                    family = ocat(R = 4), data = pheno_test)
summary(call.mod_gam)

newdat <- expand.grid(
  julian_day = seq(min(pheno_test$julian_day), max(pheno_test$julian_day), by = 1),
  year_c = unique(pheno_test$year_c))

newdat$Site_No <- factor(levels(pheno_test$Site_No)[1], levels = levels(pheno_test$Site_No))


preds <- predict(call.mod_gam, newdata = newdat, type = "response", 
                 exclude = "s(Site_No)")

prob_matrix <- as.data.frame(preds)
colnames(prob_matrix) <- levels(pheno_test$Am_Toad)  # "0","1","2","3"
newdat <- cbind(newdat, prob_matrix)

newdat_long <- newdat %>%
  pivot_longer(cols = as.character(1:4),
               names_to = "Am_Toad_level",
               values_to = "probability")

newdat_long$Am_Toad_level <- as.numeric(newdat_long$Am_Toad_level) - 1


ggplot(newdat_long, aes(x = julian_day, y = probability, color = Am_Toad_level)) +
  geom_line(linewidth = 1) +
  facet_wrap(~ year_c) +
  labs(x = "Julian Day", y = "Predicted Probability",
       color = "Calling\nIntensity",
       title = "Predicted Am. Toad Calling Intensity (Ordinal GAM)") +
  theme_minimal()



clmm(calls ~ poly(Day, 2) * Year + (1|Site), data = your_data)