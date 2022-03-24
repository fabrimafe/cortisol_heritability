## Heritability of chimp cortisol - repeatability #####

###
# 1 Import data, packages, and functions ####

d = read.csv("new_cort.csv", 
                  header = T, colClasses = "character")

library(tidyverse)
library(lme4)
library(brms)
library(car)
library(coda)
library(loo)
library(arm)
library(RColorBrewer)

source("C:/Users/Megaport/OneDrive/Project_data/heritability_slopes/diagnostic_fcns.r")
source("C:/Users/Megaport/OneDrive/Project_data/heritability_slopes/helpers.r")
source("C:/Users/Megaport/OneDrive/Project_data/heritability_slopes/drop1_para.r")

d = filter(d, !community == "Middle")

sort(unique(d$demographic))

str(d)

d$X=NULL

d = d %>%
  mutate(log.cort = as.numeric(as.character(log.cort)),
         time.h = as.numeric(as.character(time.h)),
         z.time = as.numeric(as.character(z.time )),
         sex_ratio = as.numeric(as.character(sex_ratio)),
         z.sex.ratio = as.numeric(as.character(z.sex.ratio)),
         z.demo.code = as.numeric(as.character(z.demo.code)),
         age_at_sample = as.numeric(as.character(age_at_sample)),
         z.age_at_sample = as.numeric(as.character(z.age_at_sample)),
         community_size = as.numeric(as.character(community_size)),
         z.comm_size = as.numeric(as.character(z.comm_size)),
         seasonDate = as.numeric(as.character(seasonDate)),
         Year = as.numeric(as.character(Year)),
         Month = as.numeric(as.character(Month)),
         Day = as.numeric(as.character(Day)),
         z.database.code = as.numeric(as.character(z.database.code)),
         z.group.code = as.numeric(as.character(z.group.code)),
         z.site.code = as.numeric(as.character(z.site.code)))


adult_males = filter(d, demographic == "adult_male")
adult_females = filter(d, demographic == "cycling" | demographic == "lactating")
juveniles = filter(d, demographic == "juvenile_female" | demographic == "juvenile_male")
all_inds = d

nrow(adult_males)+nrow(adult_females)+nrow(juveniles) # all individuals covered

# for adult males, we don't want individuals without rank

adult_males$rank = as.numeric(as.character(adult_males$rank))
xx = filter(adult_males, is.na(rank)) # 84 samples from individuals before we had 5 pant grunts
adult_males = filter(adult_males, !is.na(rank))

# for adult females, we don't want individuals without reproductive state

sort(unique(adult_females$repro.stat))

xx = filter(adult_females, is.na(repro.stat))
xx = filter(adult_females, repro.stat == "")
adult_females = filter(adult_females, !repro.stat == "")
adult_females = filter(adult_females, !is.na(repro.stat))

# 2 All individuals  ####

test_data = all_inds

# vifs
xx=lm(log.cort ~ z.time + z.sex.ratio + z.demo.code + z.age_at_sample +z.lcms.code + 
        z.group.code + z.site.code + z.comm_size +
        sin(seasonDate) + cos(seasonDate), data=test_data)
vif(xx)

# perfectly collinear in there so remove site
xx=lm(log.cort ~ z.time + z.sex.ratio + z.demo.code + z.age_at_sample +z.lcms.code + 
        z.group.code + z.comm_size +
        sin(seasonDate) + cos(seasonDate), data=test_data)
vif(xx)#146.7919

# without community
xx=lm(log.cort ~ z.time + z.sex.ratio + z.demo.code + z.age_at_sample +z.lcms.code + 
        z.site.code + z.comm_size + sin(seasonDate) + cos(seasonDate), data=test_data)
vif(xx)

xx=lm(log.cort ~ z.time + z.sex.ratio + z.demo.code + z.age_at_sample + z.lcms.code + 
        z.comm_size + sin(seasonDate) + cos(seasonDate), data=test_data)
max(vif(xx)) #3.11

# having just community size is best.
names(test_data)

hist(test_data$log.cort)
str(test_data)

xx.fe.re = fe.re.tab(fe.model=
                       "log.cort~z.time+z.sex.ratio+z.demo.code+z.age_at_sample+z.comm_size+
                     z.lcms.code+seasonDate", 
                     re="(1|name_code)+(1|ID_year)+(1|group_year)+(1|database.code)", data=test_data)

xx.fe.re$summary

nrow(test_data) #6123
nrow(xx.fe.re$data) #6123 so all ok

test_data = xx.fe.re$data

# 3 All inds intercept ####

all_test_data = test_data

mprior = get_prior(log.cort ~ 
                     I(z.time^2)*z.demo.code + z.time*z.demo.code + 
                     I(z.time^2)*z.age_at_sample + z.time*z.age_at_sample + 
                     I(z.time^2)*sin(seasonDate) + z.time*sin(seasonDate) + 
                     I(z.time^2)*cos(seasonDate) + z.time*cos(seasonDate) + 
                     I(z.time^2)*z.sex.ratio + z.time*z.sex.ratio +
                     I(z.time^2)*z.comm_size + z.time*z.comm_size + 
                     I(z.time^2)*z.lcms.code + z.time*z.lcms.code + 
                     (1|database.code) + (1|group_year), 
                   data = all_test_data, family = gaussian(link="identity"))

mprior$prior[2:24] <- "normal(0,1)"

make_stancode(log.cort ~ 
                I(z.time^2)*z.demo.code + z.time*z.demo.code + 
                I(z.time^2)*z.age_at_sample + z.time*z.age_at_sample + 
                I(z.time^2)*sin(seasonDate) + z.time*sin(seasonDate) + 
                I(z.time^2)*cos(seasonDate) + z.time*cos(seasonDate) + 
                I(z.time^2)*z.sex.ratio + z.time*z.sex.ratio +
                I(z.time^2)*z.comm_size + z.time*z.comm_size + 
                I(z.time^2)*z.lcms.code + z.time*z.lcms.code + 
                (1|database.code) + (1|group_year), 
              data = all_test_data, family = gaussian(link="identity"), prior = mprior)

all.int.null = brm(log.cort ~ 
                     I(z.time^2)*z.demo.code + z.time*z.demo.code + 
                     I(z.time^2)*z.age_at_sample + z.time*z.age_at_sample + 
                     I(z.time^2)*sin(seasonDate) + z.time*sin(seasonDate) + 
                     I(z.time^2)*cos(seasonDate) + z.time*cos(seasonDate) + 
                     I(z.time^2)*z.sex.ratio + z.time*z.sex.ratio +
                     I(z.time^2)*z.comm_size + z.time*z.comm_size + 
                     I(z.time^2)*z.lcms.code + z.time*z.lcms.code + 
                     (1|database.code) + (1|group_year), 
                   data = all_test_data, family = gaussian(link="identity"), prior = mprior,
                   chains = 4, cores = 4, warmup = 2000, iter = 4000, thin=2,
                   control=list(adapt_delta=0.99, max_treedepth = 20))

mprior2 = get_prior(log.cort ~ 
                      I(z.time^2)*z.demo.code + z.time*z.demo.code + 
                      I(z.time^2)*z.age_at_sample + z.time*z.age_at_sample + 
                      I(z.time^2)*sin(seasonDate) + z.time*sin(seasonDate) + 
                      I(z.time^2)*cos(seasonDate) + z.time*cos(seasonDate) + 
                      I(z.time^2)*z.sex.ratio + z.time*z.sex.ratio +
                      I(z.time^2)*z.comm_size + z.time*z.comm_size + 
                      I(z.time^2)*z.lcms.code + z.time*z.lcms.code + 
                      (1|database.code) + (1|group_year) + (1|name_code), 
                    data = all_test_data, family = gaussian(link="identity"))

mprior2$prior[2:24] <- "normal(0,1)"

make_stancode(log.cort ~ 
                I(z.time^2)*z.demo.code + z.time*z.demo.code + 
                I(z.time^2)*z.age_at_sample + z.time*z.age_at_sample + 
                I(z.time^2)*sin(seasonDate) + z.time*sin(seasonDate) + 
                I(z.time^2)*cos(seasonDate) + z.time*cos(seasonDate) + 
                I(z.time^2)*z.sex.ratio + z.time*z.sex.ratio +
                I(z.time^2)*z.comm_size + z.time*z.comm_size + 
                I(z.time^2)*z.lcms.code + z.time*z.lcms.code + 
                (1|database.code) + (1|group_year) + (1|name_code), 
              data = all_test_data, family = gaussian(link="identity"), prior = mprior2)

all.int.mod = brm(log.cort ~ 
                    I(z.time^2)*z.demo.code + z.time*z.demo.code + 
                    I(z.time^2)*z.age_at_sample + z.time*z.age_at_sample + 
                    I(z.time^2)*sin(seasonDate) + z.time*sin(seasonDate) + 
                    I(z.time^2)*cos(seasonDate) + z.time*cos(seasonDate) + 
                    I(z.time^2)*z.sex.ratio + z.time*z.sex.ratio +
                    I(z.time^2)*z.comm_size + z.time*z.comm_size + 
                    I(z.time^2)*z.lcms.code + z.time*z.lcms.code + 
                    (1|database.code) + (1|group_year) + (1|name_code), 
                  data = all_test_data, family = gaussian(link="identity"), prior = mprior2,
                  chains = 4, cores = 4, warmup = 2000, iter = 4000, thin=2,
                  control=list(adapt_delta=0.99,max_treedepth = 20))

# 4 All inds slopes ####

mprior3 = get_prior(log.cort ~ 
                      I(z.time^2)*z.demo.code + z.time*z.demo.code + 
                      I(z.time^2)*z.age_at_sample + z.time*z.age_at_sample + 
                      I(z.time^2)*sin(seasonDate) + z.time*sin(seasonDate) + 
                      I(z.time^2)*cos(seasonDate) + z.time*cos(seasonDate) + 
                      I(z.time^2)*z.sex.ratio + z.time*z.sex.ratio +
                      I(z.time^2)*z.comm_size + z.time*z.comm_size + 
                      I(z.time^2)*z.lcms.code + z.time*z.lcms.code + 
                      (1|name_code) + 
                      (1|ID_year) + 
                      (1|database.code) + (1|group_year), 
                    data = all_test_data, family = gaussian(link="identity"))

mprior3$prior[2:24] <- "normal(0,1)"

make_stancode(log.cort ~ 
                I(z.time^2)*z.demo.code + z.time*z.demo.code + 
                I(z.time^2)*z.age_at_sample + z.time*z.age_at_sample + 
                I(z.time^2)*sin(seasonDate) + z.time*sin(seasonDate) + 
                I(z.time^2)*cos(seasonDate) + z.time*cos(seasonDate) + 
                I(z.time^2)*z.sex.ratio + z.time*z.sex.ratio +
                I(z.time^2)*z.comm_size + z.time*z.comm_size + 
                I(z.time^2)*z.lcms.code + z.time*z.lcms.code + 
                (1|name_code) + 
                (1|ID_year) + 
                (1|database.code) + (1|group_year), 
              data = all_test_data, family = gaussian(link="identity"), prior = mprior3)

all.slo.null = brm(log.cort ~ 
                     I(z.time^2)*z.demo.code + z.time*z.demo.code + 
                     I(z.time^2)*z.age_at_sample + z.time*z.age_at_sample + 
                     I(z.time^2)*sin(seasonDate) + z.time*sin(seasonDate) + 
                     I(z.time^2)*cos(seasonDate) + z.time*cos(seasonDate) + 
                     I(z.time^2)*z.sex.ratio + z.time*z.sex.ratio +
                     I(z.time^2)*z.comm_size + z.time*z.comm_size + 
                     I(z.time^2)*z.lcms.code + z.time*z.lcms.code + 
                     (1|name_code) + 
                     (1|ID_year) + 
                     (1|database.code) + (1|group_year), 
                   data = all_test_data, family = gaussian(link="identity"), prior = mprior3, 
                   chains = 4, cores = 4, warmup = 2000, iter = 4000, thin=2,
                   control=list(adapt_delta=0.99,max_treedepth = 20))

mprior4 = get_prior(log.cort ~ 
                      I(z.time^2)*z.demo.code + z.time*z.demo.code + 
                      I(z.time^2)*z.age_at_sample + z.time*z.age_at_sample + 
                      I(z.time^2)*sin(seasonDate) + z.time*sin(seasonDate) + 
                      I(z.time^2)*cos(seasonDate) + z.time*cos(seasonDate) + 
                      I(z.time^2)*z.sex.ratio + z.time*z.sex.ratio +
                      I(z.time^2)*z.comm_size + z.time*z.comm_size + 
                      I(z.time^2)*z.lcms.code + z.time*z.lcms.code + 
                      (I(z.time^2) + z.time|name_code) + 
                      (I(z.time^2) + z.time|ID_year) + 
                      (1|database.code) + (1|group_year), 
                    data = all_test_data, family = gaussian(link="identity"))

mprior4$prior[2:24] <- "normal(0,1)"

make_stancode(log.cort ~ 
                I(z.time^2)*z.demo.code + z.time*z.demo.code + 
                I(z.time^2)*z.age_at_sample + z.time*z.age_at_sample + 
                I(z.time^2)*sin(seasonDate) + z.time*sin(seasonDate) + 
                I(z.time^2)*cos(seasonDate) + z.time*cos(seasonDate) + 
                I(z.time^2)*z.sex.ratio + z.time*z.sex.ratio +
                I(z.time^2)*z.comm_size + z.time*z.comm_size + 
                I(z.time^2)*z.lcms.code + z.time*z.lcms.code + 
                (I(z.time^2) + z.time|name_code) + 
                (I(z.time^2) + z.time|ID_year) + 
                (1|database.code) + (1|group_year), 
              data = all_test_data, family = gaussian(link="identity"), prior = mprior3)

all.slo.mod = brm(log.cort ~ 
                    I(z.time^2)*z.demo.code + z.time*z.demo.code + 
                    I(z.time^2)*z.age_at_sample + z.time*z.age_at_sample + 
                    I(z.time^2)*sin(seasonDate) + z.time*sin(seasonDate) + 
                    I(z.time^2)*cos(seasonDate) + z.time*cos(seasonDate) + 
                    I(z.time^2)*z.sex.ratio + z.time*z.sex.ratio +
                    I(z.time^2)*z.comm_size + z.time*z.comm_size + 
                    I(z.time^2)*z.lcms.code + z.time*z.lcms.code + 
                    (I(z.time^2) + z.time|name_code) + 
                    (I(z.time^2) + z.time|ID_year) + 
                    (1|database.code) + (1|group_year), 
                  data = all_test_data, family = gaussian(link="identity"), prior = mprior3, 
                  chains = 4, cores = 4, warmup = 2000, iter = 4000, thin=2,
                  control=list(adapt_delta=0.99,max_treedepth = 20))

# 5 All inds mod comps ####

fit1 = loo(all.int.null, cores = 4)
fit2 = loo(all.int.mod, cores = 4)
fit3 = loo(all.slo.null, cores = 4)
fit4 = loo(all.slo.mod, cores = 4)

comp = loo_compare(fit1, fit2, fit3, fit4)
print(comp, digits = 3)

'            elpd_diff se_diff 
all.slo.null    0.000     0.000
all.slo.mod    -1.794     4.609
all.int.mod   -82.951    13.413
all.int.null -376.514    28.706'

# support for intercept only reaction norm model

# strong support for the slo null, aka the intercept model

# 6 All inds repeatability ####

summary(all.slo.mod)
plot(all.slo.mod)
pp_check(all.slo.mod)

colnames(posterior_samples(all.slo.mod))[25:43]

var.grpyr <- posterior_samples(all.slo.mod)$"sd_group_year__Intercept"^2
var.db <- posterior_samples(all.slo.mod)$"sd_database.code__Intercept"^2

var.ind.int <- posterior_samples(all.slo.mod)$"sd_name_code__Intercept"^2
var.ind.yr.int <- posterior_samples(all.slo.mod)$"sd_ID_year__Intercept"^2

var.ind.yr.slp2 <- posterior_samples(all.slo.mod)$"sd_ID_year__Iz.timeE2"^2
var.ind.yr.slp <- posterior_samples(all.slo.mod)$"sd_ID_year__z.time"^2

var.ind.slp2 <- posterior_samples(all.slo.mod)$"sd_name_code__Iz.timeE2"^2
var.ind.slp <- posterior_samples(all.slo.mod)$"sd_name_code__z.time"^2

var.res <- posterior_samples(all.slo.mod)$"sigma"^2

VTotal = var.ind.int + var.ind.yr.int + var.grpyr + var.res + var.db
short_termR = (var.ind.int+var.ind.yr.int)/VTotal

mean(short_termR);HPDinterval(as.mcmc(short_termR),0.95)

long_termR = var.ind.int/VTotal

mean(long_termR);HPDinterval(as.mcmc(long_termR),0.95)

InterceptR = var.ind.int/(var.ind.int+var.ind.yr.int)
mean(InterceptR);HPDinterval(as.mcmc(InterceptR),0.95)

LinearR = var.ind.slp/(var.ind.slp+var.ind.yr.slp)
mean(LinearR);
round(HPDinterval(as.mcmc(LinearR),0.95),3)

QuadR = var.ind.slp2/(var.ind.slp2+var.ind.yr.slp2)
mean(QuadR)
round(HPDinterval(as.mcmc(QuadR),0.95),3)

# 7 Male only analysis ####

test_data = adult_males

# for males only we need to exclude individuals without rank info

test_data$rank = as.numeric(as.character(test_data$rank))

test_data = filter(test_data, !is.na(rank)) # all good

names(test_data)

# need to restandardize within this particular dataset

test_data$time.h=(as.numeric(as.character(test_data$time.h)))
range(test_data$time.h)
test_data$z.time = as.vector((test_data$time.h - mean(test_data$time.h, na.rm=TRUE))/(2*sd(test_data$time.h, na.rm=TRUE)))
range(test_data$z.time)

test_data$z.age_at_sample=as.vector((test_data$age_at_sample - mean(test_data$age_at_sample))/(2*sd(test_data$age_at_sample)))
range(test_data$z.age_at_sample)

test_data$z.sex.ratio = as.vector((test_data$sex_ratio - mean(test_data$sex_ratio, na.rm=TRUE))
                                  /(2*sd(test_data$sex_ratio, na.rm=TRUE)))
range(test_data$z.sex.ratio)

test_data$z.comm_size = as.vector((test_data$community_size - mean(test_data$community_size, na.rm=TRUE))
                                  /(2*sd(test_data$community_size, na.rm=TRUE)))
range(test_data$z.comm_size)

test_data$z.rank = as.vector((test_data$rank - mean(test_data$rank, na.rm=TRUE))
                             /(2*sd(test_data$rank, na.rm=TRUE)))
range(test_data$z.rank)

# check data

names(test_data)
xx.fe.re = fe.re.tab(fe.model=
                       "log.cort~z.time+z.sex.ratio+z.rank+z.age_at_sample+z.comm_size+
                     z.lcms.code+seasonDate", 
                     re="(1|name_code)+(1|ID_year)+(1|group_year)+(1|database.code)", data=test_data)

xx.fe.re$summary

nrow(test_data) #3196
nrow(xx.fe.re$data) #3196

test_data = xx.fe.re$data

# 8 Male intercept ####

male_test_data = test_data

mprior = get_prior(log.cort ~ 
                     I(z.time^2)*z.rank + z.time*z.rank + 
                     I(z.time^2)*z.age_at_sample + z.time*z.age_at_sample + 
                     I(z.time^2)*sin(seasonDate) + z.time*sin(seasonDate) + 
                     I(z.time^2)*cos(seasonDate) + z.time*cos(seasonDate) + 
                     I(z.time^2)*z.sex.ratio + z.time*z.sex.ratio +
                     I(z.time^2)*z.comm_size + z.time*z.comm_size + 
                     I(z.time^2)*z.lcms.code + z.time*z.lcms.code + 
                     (1|database.code) + (1|group_year), 
                   data = male_test_data, family = gaussian(link="identity"))

mprior$prior[2:24] <- "normal(0,1)"

make_stancode(log.cort ~ 
                I(z.time^2)*z.rank + z.time*z.rank + 
                I(z.time^2)*z.age_at_sample + z.time*z.age_at_sample + 
                I(z.time^2)*sin(seasonDate) + z.time*sin(seasonDate) + 
                I(z.time^2)*cos(seasonDate) + z.time*cos(seasonDate) + 
                I(z.time^2)*z.sex.ratio + z.time*z.sex.ratio +
                I(z.time^2)*z.comm_size + z.time*z.comm_size + 
                I(z.time^2)*z.lcms.code + z.time*z.lcms.code + 
                (1|database.code) + (1|group_year), 
              data = male_test_data, family = gaussian(link="identity"), prior = mprior)

male.int.null = brm(log.cort ~ 
                      I(z.time^2)*z.rank + z.time*z.rank +  
                      I(z.time^2)*z.age_at_sample + z.time*z.age_at_sample + 
                      I(z.time^2)*sin(seasonDate) + z.time*sin(seasonDate) + 
                      I(z.time^2)*cos(seasonDate) + z.time*cos(seasonDate) + 
                      I(z.time^2)*z.sex.ratio + z.time*z.sex.ratio +
                      I(z.time^2)*z.comm_size + z.time*z.comm_size + 
                      I(z.time^2)*z.lcms.code + z.time*z.lcms.code + 
                      (1|database.code) + (1|group_year), 
                    data = male_test_data, family = gaussian(link="identity"), prior = mprior,
                    chains = 4, cores = 4, warmup = 2000, iter = 4000, thin=2,
                    control=list(adapt_delta=0.99, max_treedepth = 20))

mprior2 = get_prior(log.cort ~ 
                      I(z.time^2)*z.rank + z.time*z.rank + 
                      I(z.time^2)*z.age_at_sample + z.time*z.age_at_sample + 
                      I(z.time^2)*sin(seasonDate) + z.time*sin(seasonDate) + 
                      I(z.time^2)*cos(seasonDate) + z.time*cos(seasonDate) + 
                      I(z.time^2)*z.sex.ratio + z.time*z.sex.ratio +
                      I(z.time^2)*z.comm_size + z.time*z.comm_size + 
                      I(z.time^2)*z.lcms.code + z.time*z.lcms.code + 
                      (1|database.code) + (1|group_year) + (1|name_code), 
                    data = male_test_data, family = gaussian(link="identity"))

mprior2$prior[2:24] <- "normal(0,1)"

make_stancode(log.cort ~ 
                I(z.time^2)*z.rank + z.time*z.rank + 
                I(z.time^2)*z.age_at_sample + z.time*z.age_at_sample + 
                I(z.time^2)*sin(seasonDate) + z.time*sin(seasonDate) + 
                I(z.time^2)*cos(seasonDate) + z.time*cos(seasonDate) + 
                I(z.time^2)*z.sex.ratio + z.time*z.sex.ratio +
                I(z.time^2)*z.comm_size + z.time*z.comm_size + 
                I(z.time^2)*z.lcms.code + z.time*z.lcms.code + 
                (1|database.code) + (1|group_year) + (1|name_code), 
              data = male_test_data, family = gaussian(link="identity"), prior = mprior2)

male.int.mod = brm(log.cort ~ 
                     I(z.time^2)*z.rank + z.time*z.rank + 
                     I(z.time^2)*z.age_at_sample + z.time*z.age_at_sample + 
                     I(z.time^2)*sin(seasonDate) + z.time*sin(seasonDate) + 
                     I(z.time^2)*cos(seasonDate) + z.time*cos(seasonDate) + 
                     I(z.time^2)*z.sex.ratio + z.time*z.sex.ratio +
                     I(z.time^2)*z.comm_size + z.time*z.comm_size + 
                     I(z.time^2)*z.lcms.code + z.time*z.lcms.code + 
                     (1|database.code) + (1|group_year) + (1|name_code), 
                   data = male_test_data, family = gaussian(link="identity"), prior = mprior2,
                   chains = 4, cores = 4, warmup = 2000, iter = 4000, thin=2,
                   control=list(adapt_delta=0.99,max_treedepth = 20))

# 9 Male slopes ####

mprior3 = get_prior(log.cort ~ 
                      I(z.time^2)*z.rank + z.time*z.rank + 
                      I(z.time^2)*z.age_at_sample + z.time*z.age_at_sample + 
                      I(z.time^2)*sin(seasonDate) + z.time*sin(seasonDate) + 
                      I(z.time^2)*cos(seasonDate) + z.time*cos(seasonDate) + 
                      I(z.time^2)*z.sex.ratio + z.time*z.sex.ratio +
                      I(z.time^2)*z.comm_size + z.time*z.comm_size + 
                      I(z.time^2)*z.lcms.code + z.time*z.lcms.code + 
                      (1|name_code) + 
                      (1|ID_year) + 
                      (1|database.code) + (1|group_year), 
                    data = male_test_data, family = gaussian(link="identity"))

mprior3$prior[2:24] <- "normal(0,1)"

make_stancode(log.cort ~ 
                I(z.time^2)*z.rank + z.time*z.rank + 
                I(z.time^2)*z.age_at_sample + z.time*z.age_at_sample + 
                I(z.time^2)*sin(seasonDate) + z.time*sin(seasonDate) + 
                I(z.time^2)*cos(seasonDate) + z.time*cos(seasonDate) + 
                I(z.time^2)*z.sex.ratio + z.time*z.sex.ratio +
                I(z.time^2)*z.comm_size + z.time*z.comm_size + 
                I(z.time^2)*z.lcms.code + z.time*z.lcms.code + 
                (1|name_code) + 
                (1|ID_year) + 
                (1|database.code) + (1|group_year), 
              data = male_test_data, family = gaussian(link="identity"), prior = mprior3)

male.slo.null = brm(log.cort ~ 
                      I(z.time^2)*z.rank + z.time*z.rank + 
                      I(z.time^2)*z.age_at_sample + z.time*z.age_at_sample + 
                      I(z.time^2)*sin(seasonDate) + z.time*sin(seasonDate) + 
                      I(z.time^2)*cos(seasonDate) + z.time*cos(seasonDate) + 
                      I(z.time^2)*z.sex.ratio + z.time*z.sex.ratio +
                      I(z.time^2)*z.comm_size + z.time*z.comm_size + 
                      I(z.time^2)*z.lcms.code + z.time*z.lcms.code + 
                      (1|name_code) + 
                      (1|ID_year) + 
                      (1|database.code) + (1|group_year), 
                    data = male_test_data, family = gaussian(link="identity"), prior = mprior3, 
                    chains = 4, cores = 4, warmup = 2000, iter = 4000, thin=2,
                    control=list(adapt_delta=0.99,max_treedepth = 20))

mprior4 = get_prior(log.cort ~ 
                      I(z.time^2)*z.rank + z.time*z.rank + 
                      I(z.time^2)*z.age_at_sample + z.time*z.age_at_sample + 
                      I(z.time^2)*sin(seasonDate) + z.time*sin(seasonDate) + 
                      I(z.time^2)*cos(seasonDate) + z.time*cos(seasonDate) + 
                      I(z.time^2)*z.sex.ratio + z.time*z.sex.ratio +
                      I(z.time^2)*z.comm_size + z.time*z.comm_size + 
                      I(z.time^2)*z.lcms.code + z.time*z.lcms.code + 
                      (I(z.time^2) + z.time|name_code) + 
                      (I(z.time^2) + z.time|ID_year) + 
                      (1|database.code) + (1|group_year), 
                    data = male_test_data, family = gaussian(link="identity"))

mprior4$prior[2:24] <- "normal(0,1)"

make_stancode(log.cort ~ 
                I(z.time^2)*z.rank + z.time*z.rank + 
                I(z.time^2)*z.age_at_sample + z.time*z.age_at_sample + 
                I(z.time^2)*sin(seasonDate) + z.time*sin(seasonDate) + 
                I(z.time^2)*cos(seasonDate) + z.time*cos(seasonDate) + 
                I(z.time^2)*z.sex.ratio + z.time*z.sex.ratio +
                I(z.time^2)*z.comm_size + z.time*z.comm_size + 
                I(z.time^2)*z.lcms.code + z.time*z.lcms.code + 
                (I(z.time^2) + z.time|name_code) + 
                (I(z.time^2) + z.time|ID_year) + 
                (1|database.code) + (1|group_year), 
              data = male_test_data, family = gaussian(link="identity"), prior = mprior3)

male.slo.mod = brm(log.cort ~ 
                     I(z.time^2)*z.rank + z.time*z.rank + 
                     I(z.time^2)*z.age_at_sample + z.time*z.age_at_sample + 
                     I(z.time^2)*sin(seasonDate) + z.time*sin(seasonDate) + 
                     I(z.time^2)*cos(seasonDate) + z.time*cos(seasonDate) + 
                     I(z.time^2)*z.sex.ratio + z.time*z.sex.ratio +
                     I(z.time^2)*z.comm_size + z.time*z.comm_size + 
                     I(z.time^2)*z.lcms.code + z.time*z.lcms.code + 
                     (I(z.time^2) + z.time|name_code) + 
                     (I(z.time^2) + z.time|ID_year) + 
                     (1|database.code) + (1|group_year), 
                   data = male_test_data, family = gaussian(link="identity"), prior = mprior3, 
                   chains = 4, cores = 4, warmup = 2000, iter = 4000, thin=2,
                   control=list(adapt_delta=0.99,max_treedepth = 20))

# 10 Male mod comps ####

fit1 = loo(male.int.null, cores = 4)
fit2 = loo(male.int.mod, cores = 4)
fit3 = loo(male.slo.null, cores = 4)
fit4 = loo(male.slo.mod, cores = 4)

comp = loo_compare(fit1, fit2, fit3, fit4)
print(comp, digits = 3)

'              elpd_diff se_diff
male.slo.null   0.000     0.000
male.slo.mod   -1.231     3.568
male.int.mod  -31.986     8.144
male.int.null -93.334    14.090'

# 11 Males repeatability ####

summary(male.slo.mod)
plot(male.slo.mod)
pp_check(male.slo.mod)


var.grpyr <- posterior_samples(male.slo.mod)$"sd_group_year__Intercept"^2
var.db <- posterior_samples(male.slo.mod)$"sd_database.code__Intercept"^2

var.ind.int <- posterior_samples(male.slo.mod)$"sd_name_code__Intercept"^2
var.ind.yr.int <- posterior_samples(male.slo.mod)$"sd_ID_year__Intercept"^2

var.ind.yr.slp2 <- posterior_samples(male.slo.mod)$"sd_ID_year__Iz.timeE2"^2
var.ind.yr.slp <- posterior_samples(male.slo.mod)$"sd_ID_year__z.time"^2

var.ind.slp2 <- posterior_samples(male.slo.mod)$"sd_name_code__Iz.timeE2"^2
var.ind.slp <- posterior_samples(male.slo.mod)$"sd_name_code__z.time"^2

var.res <- posterior_samples(male.slo.mod)$"sigma"^2

VTotal = var.ind.int + var.ind.yr.int + var.grpyr + var.res + var.db
short_termR = (var.ind.int+var.ind.yr.int)/VTotal

mean(short_termR);HPDinterval(as.mcmc(short_termR),0.95)

long_termR = var.ind.int/VTotal

mean(long_termR);HPDinterval(as.mcmc(long_termR),0.95)

InterceptR = var.ind.int/(var.ind.int+var.ind.yr.int)
mean(InterceptR);HPDinterval(as.mcmc(InterceptR),0.95)

LinearR = var.ind.slp/(var.ind.slp+var.ind.yr.slp)
mean(LinearR);round(HPDinterval(as.mcmc(LinearR),0.95),3)

QuadR = var.ind.slp2/(var.ind.slp2+var.ind.yr.slp2)
mean(QuadR);round(HPDinterval(as.mcmc(QuadR),0.95),3)

# 12 Female only analysis ####

test_data = adult_females

# z transformations

test_data$time.h=(as.numeric(as.character(test_data$time.h)))
range(test_data$time.h)
test_data$z.time = as.vector((test_data$time.h - mean(test_data$time.h, na.rm=TRUE))/(2*sd(test_data$time.h, na.rm=TRUE)))
range(test_data$z.time)

test_data$z.age_at_sample=as.vector((test_data$age_at_sample - mean(test_data$age_at_sample))/(2*sd(test_data$age_at_sample)))
range(test_data$z.age_at_sample)

test_data$z.sex.ratio = as.vector((test_data$sex_ratio - mean(test_data$sex_ratio, na.rm=TRUE))
                                  /(2*sd(test_data$sex_ratio, na.rm=TRUE)))
range(test_data$z.sex.ratio)

test_data$z.comm_size = as.vector((test_data$community_size - mean(test_data$community_size, na.rm=TRUE))
                                  /(2*sd(test_data$community_size, na.rm=TRUE)))
range(test_data$z.comm_size)

sort(unique(test_data$demographic))
repro.code=as.numeric(as.character(factor(test_data$demographic, levels=c("cycling","lactating"), 
                                          labels=c("0", "1")))) 
test_data = data.frame(test_data, repro.code)
rm(repro.code)

test_data$z.repro.code = as.vector((test_data$repro.code - mean(test_data$repro.code, na.rm=TRUE))
                                   /(2*sd(test_data$repro.code, na.rm=TRUE)))
test_data$z.repro.code = as.vector(as.factor((test_data$z.repro.code)))

names(test_data)
xx.fe.re = fe.re.tab(fe.model=
                       "log.cort~z.time+z.sex.ratio+z.repro.code+z.age_at_sample+z.comm_size+
                     z.lcms.code+seasonDate", 
                     re="(1|name_code)+(1|ID_year)+(1|group_year)+(1|database.code)", data=test_data)

xx.fe.re$summary

nrow(test_data) #1742
nrow(xx.fe.re$data) #1742

test_data = xx.fe.re$data

# 13 Female intercept ####

fem_test_data = test_data

mprior = get_prior(log.cort ~ 
                     I(z.time^2)*z.repro.code + z.time*z.repro.code + 
                     I(z.time^2)*z.age_at_sample + z.time*z.age_at_sample + 
                     I(z.time^2)*sin(seasonDate) + z.time*sin(seasonDate) + 
                     I(z.time^2)*cos(seasonDate) + z.time*cos(seasonDate) + 
                     I(z.time^2)*z.sex.ratio + z.time*z.sex.ratio +
                     I(z.time^2)*z.comm_size + z.time*z.comm_size + 
                     I(z.time^2)*z.lcms.code + z.time*z.lcms.code + 
                     (1|database.code) + (1|group_year), 
                   data = fem_test_data, family = gaussian(link="identity"))

mprior$prior[2:24] <- "normal(0,1)"

make_stancode(log.cort ~ 
                I(z.time^2)*z.repro.code + z.time*z.repro.code + 
                I(z.time^2)*z.age_at_sample + z.time*z.age_at_sample + 
                I(z.time^2)*sin(seasonDate) + z.time*sin(seasonDate) + 
                I(z.time^2)*cos(seasonDate) + z.time*cos(seasonDate) + 
                I(z.time^2)*z.sex.ratio + z.time*z.sex.ratio +
                I(z.time^2)*z.comm_size + z.time*z.comm_size + 
                I(z.time^2)*z.lcms.code + z.time*z.lcms.code + 
                (1|database.code) + (1|group_year), 
              data = fem_test_data, family = gaussian(link="identity"), prior = mprior)

fem.int.null = brm(log.cort ~ 
                     I(z.time^2)*z.repro.code + z.time*z.repro.code + 
                     I(z.time^2)*z.age_at_sample + z.time*z.age_at_sample + 
                     I(z.time^2)*sin(seasonDate) + z.time*sin(seasonDate) + 
                     I(z.time^2)*cos(seasonDate) + z.time*cos(seasonDate) + 
                     I(z.time^2)*z.sex.ratio + z.time*z.sex.ratio +
                     I(z.time^2)*z.comm_size + z.time*z.comm_size + 
                     I(z.time^2)*z.lcms.code + z.time*z.lcms.code + 
                     (1|database.code) + (1|group_year), 
                   data = fem_test_data, family = gaussian(link="identity"), prior = mprior,
                   chains = 4, cores = 4, warmup = 2000, iter = 4000, thin=2,
                   control=list(adapt_delta=0.99, max_treedepth = 20))

mprior2 = get_prior(log.cort ~ 
                      I(z.time^2)*z.repro.code + z.time*z.repro.code + 
                      I(z.time^2)*z.age_at_sample + z.time*z.age_at_sample + 
                      I(z.time^2)*sin(seasonDate) + z.time*sin(seasonDate) + 
                      I(z.time^2)*cos(seasonDate) + z.time*cos(seasonDate) + 
                      I(z.time^2)*z.sex.ratio + z.time*z.sex.ratio +
                      I(z.time^2)*z.comm_size + z.time*z.comm_size + 
                      I(z.time^2)*z.lcms.code + z.time*z.lcms.code + 
                      (1|database.code) + (1|group_year) + (1|name_code), 
                    data = fem_test_data, family = gaussian(link="identity"))

mprior2$prior[2:24] <- "normal(0,1)"

make_stancode(log.cort ~ 
                I(z.time^2)*z.repro.code + z.time*z.repro.code + 
                I(z.time^2)*z.age_at_sample + z.time*z.age_at_sample + 
                I(z.time^2)*sin(seasonDate) + z.time*sin(seasonDate) + 
                I(z.time^2)*cos(seasonDate) + z.time*cos(seasonDate) + 
                I(z.time^2)*z.sex.ratio + z.time*z.sex.ratio +
                I(z.time^2)*z.comm_size + z.time*z.comm_size + 
                I(z.time^2)*z.lcms.code + z.time*z.lcms.code + 
                (1|database.code) + (1|group_year) + (1|name_code), 
              data = fem_test_data, family = gaussian(link="identity"), prior = mprior2)

fem.int.mod = brm(log.cort ~ 
                    I(z.time^2)*z.repro.code + z.time*z.repro.code + 
                    I(z.time^2)*z.age_at_sample + z.time*z.age_at_sample + 
                    I(z.time^2)*sin(seasonDate) + z.time*sin(seasonDate) + 
                    I(z.time^2)*cos(seasonDate) + z.time*cos(seasonDate) + 
                    I(z.time^2)*z.sex.ratio + z.time*z.sex.ratio +
                    I(z.time^2)*z.comm_size + z.time*z.comm_size + 
                    I(z.time^2)*z.lcms.code + z.time*z.lcms.code + 
                    (1|database.code) + (1|group_year) + (1|name_code), 
                  data = fem_test_data, family = gaussian(link="identity"), prior = mprior2,
                  chains = 4, cores = 4, warmup = 2000, iter = 4000, thin=2,
                  control=list(adapt_delta=0.99,max_treedepth = 20))

# 14 Female slopes ####

mprior3 = get_prior(log.cort ~ 
                      I(z.time^2)*z.repro.code + z.time*z.repro.code + 
                      I(z.time^2)*z.age_at_sample + z.time*z.age_at_sample + 
                      I(z.time^2)*sin(seasonDate) + z.time*sin(seasonDate) + 
                      I(z.time^2)*cos(seasonDate) + z.time*cos(seasonDate) + 
                      I(z.time^2)*z.sex.ratio + z.time*z.sex.ratio +
                      I(z.time^2)*z.comm_size + z.time*z.comm_size + 
                      I(z.time^2)*z.lcms.code + z.time*z.lcms.code + 
                      (1|name_code) + 
                      (1|ID_year) + 
                      (1|database.code) + (1|group_year), 
                    data = fem_test_data, family = gaussian(link="identity"))

mprior3$prior[2:24] <- "normal(0,1)"

make_stancode(log.cort ~ 
                I(z.time^2)*z.repro.code + z.time*z.repro.code + 
                I(z.time^2)*z.age_at_sample + z.time*z.age_at_sample + 
                I(z.time^2)*sin(seasonDate) + z.time*sin(seasonDate) + 
                I(z.time^2)*cos(seasonDate) + z.time*cos(seasonDate) + 
                I(z.time^2)*z.sex.ratio + z.time*z.sex.ratio +
                I(z.time^2)*z.comm_size + z.time*z.comm_size + 
                I(z.time^2)*z.lcms.code + z.time*z.lcms.code + 
                (1|name_code) + 
                (1|ID_year) + 
                (1|database.code) + (1|group_year), 
              data = fem_test_data, family = gaussian(link="identity"), prior = mprior3)

fem.slo.null = brm(log.cort ~ 
                     I(z.time^2)*z.repro.code + z.time*z.repro.code + 
                     I(z.time^2)*z.age_at_sample + z.time*z.age_at_sample + 
                     I(z.time^2)*sin(seasonDate) + z.time*sin(seasonDate) + 
                     I(z.time^2)*cos(seasonDate) + z.time*cos(seasonDate) + 
                     I(z.time^2)*z.sex.ratio + z.time*z.sex.ratio +
                     I(z.time^2)*z.comm_size + z.time*z.comm_size + 
                     I(z.time^2)*z.lcms.code + z.time*z.lcms.code + 
                     (1|name_code) + 
                     (1|ID_year) + 
                     (1|database.code) + (1|group_year), 
                   data = fem_test_data, family = gaussian(link="identity"), prior = mprior3, 
                   chains = 4, cores = 4, warmup = 2000, iter = 4000, thin=2,
                   control=list(adapt_delta=0.99,max_treedepth = 20))

mprior4 = get_prior(log.cort ~ 
                      I(z.time^2)*z.repro.code + z.time*z.repro.code + 
                      I(z.time^2)*z.age_at_sample + z.time*z.age_at_sample + 
                      I(z.time^2)*sin(seasonDate) + z.time*sin(seasonDate) + 
                      I(z.time^2)*cos(seasonDate) + z.time*cos(seasonDate) + 
                      I(z.time^2)*z.sex.ratio + z.time*z.sex.ratio +
                      I(z.time^2)*z.comm_size + z.time*z.comm_size + 
                      I(z.time^2)*z.lcms.code + z.time*z.lcms.code + 
                      (I(z.time^2) + z.time|name_code) + 
                      (I(z.time^2) + z.time|ID_year) + 
                      (1|database.code) + (1|group_year), 
                    data = fem_test_data, family = gaussian(link="identity"))

mprior4$prior[2:24] <- "normal(0,1)"

make_stancode(log.cort ~ 
                I(z.time^2)*z.repro.code + z.time*z.repro.code + 
                I(z.time^2)*z.age_at_sample + z.time*z.age_at_sample + 
                I(z.time^2)*sin(seasonDate) + z.time*sin(seasonDate) + 
                I(z.time^2)*cos(seasonDate) + z.time*cos(seasonDate) + 
                I(z.time^2)*z.sex.ratio + z.time*z.sex.ratio +
                I(z.time^2)*z.comm_size + z.time*z.comm_size + 
                I(z.time^2)*z.lcms.code + z.time*z.lcms.code + 
                (I(z.time^2) + z.time|name_code) + 
                (I(z.time^2) + z.time|ID_year) + 
                (1|database.code) + (1|group_year), 
              data = fem_test_data, family = gaussian(link="identity"), prior = mprior3)

fem.slo.mod = brm(log.cort ~ 
                    I(z.time^2)*z.repro.code + z.time*z.repro.code + 
                    I(z.time^2)*z.age_at_sample + z.time*z.age_at_sample + 
                    I(z.time^2)*sin(seasonDate) + z.time*sin(seasonDate) + 
                    I(z.time^2)*cos(seasonDate) + z.time*cos(seasonDate) + 
                    I(z.time^2)*z.sex.ratio + z.time*z.sex.ratio +
                    I(z.time^2)*z.comm_size + z.time*z.comm_size + 
                    I(z.time^2)*z.lcms.code + z.time*z.lcms.code + 
                    (I(z.time^2) + z.time|name_code) + 
                    (I(z.time^2) + z.time|ID_year) + 
                    (1|database.code) + (1|group_year), 
                  data = fem_test_data, family = gaussian(link="identity"), prior = mprior3, 
                  chains = 4, cores = 4, warmup = 2000, iter = 4000, thin=2,
                  control=list(adapt_delta=0.99,max_treedepth = 20))

# 15 Female mod comps ####

fit1 = loo(fem.int.null, cores = 4)
fit2 = loo(fem.int.mod, cores = 4)
fit3 = loo(fem.slo.null, cores = 4)
fit4 = loo(fem.slo.mod, cores = 4)

comp = loo_compare(fit1, fit2, fit3, fit4)
print(comp, digits = 3)

'             elpd_diff se_diff
fem.slo.null   0.000     0.000
fem.int.mod   -1.686     2.329
fem.slo.mod   -1.830     2.400
fem.int.null -60.191    12.037'

# 16 Females repeatability ####

summary(fem.slo.mod)
plot(fem.slo.mod)
pp_check(fem.slo.mod)

var.grpyr <- posterior_samples(fem.slo.mod)$"sd_group_year__Intercept"^2
var.db <- posterior_samples(fem.slo.mod)$"sd_database.code__Intercept"^2

var.ind.int <- posterior_samples(fem.slo.mod)$"sd_name_code__Intercept"^2
var.ind.yr.int <- posterior_samples(fem.slo.mod)$"sd_ID_year__Intercept"^2

var.ind.yr.slp2 <- posterior_samples(fem.slo.mod)$"sd_ID_year__Iz.timeE2"^2
var.ind.yr.slp <- posterior_samples(fem.slo.mod)$"sd_ID_year__z.time"^2

var.ind.slp2 <- posterior_samples(fem.slo.mod)$"sd_name_code__Iz.timeE2"^2
var.ind.slp <- posterior_samples(fem.slo.mod)$"sd_name_code__z.time"^2

var.res <- posterior_samples(fem.slo.mod)$"sigma"^2

VTotal = var.ind.int + var.ind.yr.int + var.grpyr + var.res + var.db
short_termR = (var.ind.int+var.ind.yr.int)/VTotal

mean(short_termR);HPDinterval(as.mcmc(short_termR),0.95)

long_termR = var.ind.int/VTotal

mean(long_termR);HPDinterval(as.mcmc(long_termR),0.95)

InterceptR = var.ind.int/(var.ind.int+var.ind.yr.int)
mean(InterceptR);HPDinterval(as.mcmc(InterceptR),0.95)

LinearR = var.ind.slp/(var.ind.slp+var.ind.yr.slp)
mean(LinearR);round(HPDinterval(as.mcmc(LinearR),0.95),3)

QuadR = var.ind.slp2/(var.ind.slp2+var.ind.yr.slp2)
mean(QuadR);round(HPDinterval(as.mcmc(QuadR),0.95),3)

# 17 Juveniles only analysis ####

test_data = juveniles

# z transformations

test_data$time.h=(as.numeric(as.character(test_data$time.h)))
range(test_data$time.h)
test_data$z.time = as.vector((test_data$time.h - mean(test_data$time.h, na.rm=TRUE))/(2*sd(test_data$time.h, na.rm=TRUE)))
range(test_data$z.time)

test_data$z.age_at_sample=as.vector((test_data$age_at_sample - mean(test_data$age_at_sample))/(2*sd(test_data$age_at_sample)))
range(test_data$z.age_at_sample)

test_data$z.sex.ratio = as.vector((test_data$sex_ratio - mean(test_data$sex_ratio, na.rm=TRUE))
                                  /(2*sd(test_data$sex_ratio, na.rm=TRUE)))
range(test_data$z.sex.ratio)

test_data$z.comm_size = as.vector((test_data$community_size - mean(test_data$community_size, na.rm=TRUE))
                                  /(2*sd(test_data$community_size, na.rm=TRUE)))
range(test_data$z.comm_size)

levels(as.factor(test_data$demographic))
sex_code=as.numeric(as.character(factor(test_data$demographic, levels=c("juvenile_female","juvenile_male"), 
                                        labels=c("0", "1")))) 
length(sex_code[is.na(sex_code)])
test_data=data.frame(test_data, sex_code)
rm(sex_code)
test_data$z.sex.code = as.vector((test_data$sex_code - mean(test_data$sex_code, na.rm=TRUE))
                                 /(2*sd(test_data$sex_code, na.rm=TRUE)))
test_data$z.sex.code = as.vector(as.factor((test_data$z.sex.code)))
range(test_data$z.sex.code)

names(test_data)
xx.fe.re = fe.re.tab(fe.model=
                       "log.cort~z.time+z.sex.ratio+z.sex.code+z.age_at_sample+z.comm_size+
                     z.lcms.code+seasonDate", 
                     re="(1|name_code)+(1|ID_year)+(1|group_year)+(1|database.code)", data=test_data)

xx.fe.re$summary

nrow(test_data) #1138
nrow(xx.fe.re$data) #1138

test_data = xx.fe.re$data

# 18 Juvenile intercept ####

juv_test_data = test_data

mprior = get_prior(log.cort ~ 
                     I(z.time^2)*z.sex.code + z.time*z.sex.code + 
                     I(z.time^2)*z.age_at_sample + z.time*z.age_at_sample + 
                     I(z.time^2)*sin(seasonDate) + z.time*sin(seasonDate) + 
                     I(z.time^2)*cos(seasonDate) + z.time*cos(seasonDate) + 
                     I(z.time^2)*z.sex.ratio + z.time*z.sex.ratio +
                     I(z.time^2)*z.comm_size + z.time*z.comm_size + 
                     I(z.time^2)*z.lcms.code + z.time*z.lcms.code + 
                     (1|database.code) + (1|group_year), 
                   data = juv_test_data, family = gaussian(link="identity"))

mprior$prior[2:24] <- "normal(0,1)"

make_stancode(log.cort ~ 
                I(z.time^2)*z.sex.code + z.time*z.sex.code +
                I(z.time^2)*z.age_at_sample + z.time*z.age_at_sample + 
                I(z.time^2)*sin(seasonDate) + z.time*sin(seasonDate) + 
                I(z.time^2)*cos(seasonDate) + z.time*cos(seasonDate) + 
                I(z.time^2)*z.sex.ratio + z.time*z.sex.ratio +
                I(z.time^2)*z.comm_size + z.time*z.comm_size + 
                I(z.time^2)*z.lcms.code + z.time*z.lcms.code + 
                (1|database.code) + (1|group_year), 
              data = juv_test_data, family = gaussian(link="identity"), prior = mprior)

juv.int.null = brm(log.cort ~ 
                     I(z.time^2)*z.sex.code + z.time*z.sex.code +
                     I(z.time^2)*z.age_at_sample + z.time*z.age_at_sample + 
                     I(z.time^2)*sin(seasonDate) + z.time*sin(seasonDate) + 
                     I(z.time^2)*cos(seasonDate) + z.time*cos(seasonDate) + 
                     I(z.time^2)*z.sex.ratio + z.time*z.sex.ratio +
                     I(z.time^2)*z.comm_size + z.time*z.comm_size + 
                     I(z.time^2)*z.lcms.code + z.time*z.lcms.code + 
                     (1|database.code) + (1|group_year), 
                   data = juv_test_data, family = gaussian(link="identity"), prior = mprior,
                   chains = 4, cores = 4, warmup = 2000, iter = 4000, thin=2,
                   control=list(adapt_delta=0.99, max_treedepth = 20))

mprior2 = get_prior(log.cort ~ 
                      I(z.time^2)*z.sex.code + z.time*z.sex.code +
                      I(z.time^2)*z.age_at_sample + z.time*z.age_at_sample + 
                      I(z.time^2)*sin(seasonDate) + z.time*sin(seasonDate) + 
                      I(z.time^2)*cos(seasonDate) + z.time*cos(seasonDate) + 
                      I(z.time^2)*z.sex.ratio + z.time*z.sex.ratio +
                      I(z.time^2)*z.comm_size + z.time*z.comm_size + 
                      I(z.time^2)*z.lcms.code + z.time*z.lcms.code + 
                      (1|database.code) + (1|group_year) + (1|name_code), 
                    data = juv_test_data, family = gaussian(link="identity"))

mprior2$prior[2:24] <- "normal(0,1)"

make_stancode(log.cort ~ 
                I(z.time^2)*z.sex.code + z.time*z.sex.code +
                I(z.time^2)*z.age_at_sample + z.time*z.age_at_sample + 
                I(z.time^2)*sin(seasonDate) + z.time*sin(seasonDate) + 
                I(z.time^2)*cos(seasonDate) + z.time*cos(seasonDate) + 
                I(z.time^2)*z.sex.ratio + z.time*z.sex.ratio +
                I(z.time^2)*z.comm_size + z.time*z.comm_size + 
                I(z.time^2)*z.lcms.code + z.time*z.lcms.code + 
                (1|database.code) + (1|group_year) + (1|name_code), 
              data = juv_test_data, family = gaussian(link="identity"), prior = mprior2)

juv.int.mod = brm(log.cort ~ 
                    I(z.time^2)*z.sex.code + z.time*z.sex.code +
                    I(z.time^2)*z.age_at_sample + z.time*z.age_at_sample + 
                    I(z.time^2)*sin(seasonDate) + z.time*sin(seasonDate) + 
                    I(z.time^2)*cos(seasonDate) + z.time*cos(seasonDate) + 
                    I(z.time^2)*z.sex.ratio + z.time*z.sex.ratio +
                    I(z.time^2)*z.comm_size + z.time*z.comm_size + 
                    I(z.time^2)*z.lcms.code + z.time*z.lcms.code + 
                    (1|database.code) + (1|group_year) + (1|name_code), 
                  data = juv_test_data, family = gaussian(link="identity"), prior = mprior2,
                  chains = 4, cores = 4, warmup = 2000, iter = 4000, thin=2,
                  control=list(adapt_delta=0.99,max_treedepth = 20))

# 19 Juvenile slopes ####

mprior3 = get_prior(log.cort ~ 
                      I(z.time^2)*z.sex.code + z.time*z.sex.code +
                      I(z.time^2)*z.age_at_sample + z.time*z.age_at_sample + 
                      I(z.time^2)*sin(seasonDate) + z.time*sin(seasonDate) + 
                      I(z.time^2)*cos(seasonDate) + z.time*cos(seasonDate) + 
                      I(z.time^2)*z.sex.ratio + z.time*z.sex.ratio +
                      I(z.time^2)*z.comm_size + z.time*z.comm_size + 
                      I(z.time^2)*z.lcms.code + z.time*z.lcms.code + 
                      (1|name_code) + 
                      (1|ID_year) + 
                      (1|database.code) + (1|group_year), 
                    data = juv_test_data, family = gaussian(link="identity"))

mprior3$prior[2:24] <- "normal(0,1)"

make_stancode(log.cort ~ 
                I(z.time^2)*z.sex.code + z.time*z.sex.code +
                I(z.time^2)*z.age_at_sample + z.time*z.age_at_sample + 
                I(z.time^2)*sin(seasonDate) + z.time*sin(seasonDate) + 
                I(z.time^2)*cos(seasonDate) + z.time*cos(seasonDate) + 
                I(z.time^2)*z.sex.ratio + z.time*z.sex.ratio +
                I(z.time^2)*z.comm_size + z.time*z.comm_size + 
                I(z.time^2)*z.lcms.code + z.time*z.lcms.code + 
                (1|name_code) + 
                (1|ID_year) + 
                (1|database.code) + (1|group_year), 
              data = juv_test_data, family = gaussian(link="identity"), prior = mprior3)

juv.slo.null = brm(log.cort ~ 
                     I(z.time^2)*z.sex.code + z.time*z.sex.code +
                     I(z.time^2)*z.age_at_sample + z.time*z.age_at_sample + 
                     I(z.time^2)*sin(seasonDate) + z.time*sin(seasonDate) + 
                     I(z.time^2)*cos(seasonDate) + z.time*cos(seasonDate) + 
                     I(z.time^2)*z.sex.ratio + z.time*z.sex.ratio +
                     I(z.time^2)*z.comm_size + z.time*z.comm_size + 
                     I(z.time^2)*z.lcms.code + z.time*z.lcms.code + 
                     (1|name_code) + 
                     (1|ID_year) + 
                     (1|database.code) + (1|group_year), 
                   data = juv_test_data, family = gaussian(link="identity"), prior = mprior3, 
                   chains = 4, cores = 4, warmup = 2000, iter = 4000, thin=2,
                   control=list(adapt_delta=0.99,max_treedepth = 20))

mprior4 = get_prior(log.cort ~ 
                      I(z.time^2)*z.sex.code + z.time*z.sex.code +
                      I(z.time^2)*z.age_at_sample + z.time*z.age_at_sample + 
                      I(z.time^2)*sin(seasonDate) + z.time*sin(seasonDate) + 
                      I(z.time^2)*cos(seasonDate) + z.time*cos(seasonDate) + 
                      I(z.time^2)*z.sex.ratio + z.time*z.sex.ratio +
                      I(z.time^2)*z.comm_size + z.time*z.comm_size + 
                      I(z.time^2)*z.lcms.code + z.time*z.lcms.code + 
                      (I(z.time^2) + z.time|name_code) + 
                      (I(z.time^2) + z.time|ID_year) + 
                      (1|database.code) + (1|group_year), 
                    data = juv_test_data, family = gaussian(link="identity"))

mprior4$prior[2:24] <- "normal(0,1)"

make_stancode(log.cort ~ 
                I(z.time^2)*z.sex.code + z.time*z.sex.code +
                I(z.time^2)*z.age_at_sample + z.time*z.age_at_sample + 
                I(z.time^2)*sin(seasonDate) + z.time*sin(seasonDate) + 
                I(z.time^2)*cos(seasonDate) + z.time*cos(seasonDate) + 
                I(z.time^2)*z.sex.ratio + z.time*z.sex.ratio +
                I(z.time^2)*z.comm_size + z.time*z.comm_size + 
                I(z.time^2)*z.lcms.code + z.time*z.lcms.code + 
                (I(z.time^2) + z.time|name_code) + 
                (I(z.time^2) + z.time|ID_year) + 
                (1|database.code) + (1|group_year), 
              data = juv_test_data, family = gaussian(link="identity"), prior = mprior3)

juv.slo.mod = brm(log.cort ~ 
                    I(z.time^2)*z.sex.code + z.time*z.sex.code +
                    I(z.time^2)*z.age_at_sample + z.time*z.age_at_sample + 
                    I(z.time^2)*sin(seasonDate) + z.time*sin(seasonDate) + 
                    I(z.time^2)*cos(seasonDate) + z.time*cos(seasonDate) + 
                    I(z.time^2)*z.sex.ratio + z.time*z.sex.ratio +
                    I(z.time^2)*z.comm_size + z.time*z.comm_size + 
                    I(z.time^2)*z.lcms.code + z.time*z.lcms.code + 
                    (I(z.time^2) + z.time|name_code) + 
                    (I(z.time^2) + z.time|ID_year) + 
                    (1|database.code) + (1|group_year), 
                  data = juv_test_data, family = gaussian(link="identity"), prior = mprior3, 
                  chains = 4, cores = 4, warmup = 2000, iter = 4000, thin=2,
                  control=list(adapt_delta=0.99,max_treedepth = 20))

# 20 Juvenile mod comps ####

fit1 = loo(juv.int.null, cores = 4)
fit2 = loo(juv.int.mod, cores = 4)
fit3 = loo(juv.slo.null, cores = 4)
fit4 = loo(juv.slo.mod, cores = 4)

comp = loo_compare(fit1, fit2, fit3, fit4)
print(comp, digits = 3)

'             elpd_diff se_diff
juv.slo.mod    0.000     0.000
juv.slo.null -14.521     8.006
juv.int.mod  -23.768     8.980
juv.int.null -61.897    13.775'


# 21 Juv repeatability ####

summary(juv.slo.mod)
plot(juv.slo.mod)
pp_check(juv.slo.mod)

var.grpyr <- posterior_samples(juv.slo.mod)$"sd_group_year__Intercept"^2
var.db <- posterior_samples(juv.slo.mod)$"sd_database.code__Intercept"^2

var.ind.int <- posterior_samples(juv.slo.mod)$"sd_name_code__Intercept"^2
var.ind.yr.int <- posterior_samples(juv.slo.mod)$"sd_ID_year__Intercept"^2

var.ind.yr.slp2 <- posterior_samples(juv.slo.mod)$"sd_ID_year__Iz.timeE2"^2
var.ind.yr.slp <- posterior_samples(juv.slo.mod)$"sd_ID_year__z.time"^2

var.ind.slp2 <- posterior_samples(juv.slo.mod)$"sd_name_code__Iz.timeE2"^2
var.ind.slp <- posterior_samples(juv.slo.mod)$"sd_name_code__z.time"^2

var.res <- posterior_samples(juv.slo.mod)$"sigma"^2

VTotal = var.ind.int + var.ind.yr.int + var.grpyr + var.res + var.db
short_termR = (var.ind.int+var.ind.yr.int)/VTotal

mean(short_termR);HPDinterval(as.mcmc(short_termR),0.95)

long_termR = var.ind.int/VTotal

mean(long_termR);HPDinterval(as.mcmc(long_termR),0.95)

InterceptR = var.ind.int/(var.ind.int+var.ind.yr.int)
mean(InterceptR);HPDinterval(as.mcmc(InterceptR),0.95)

LinearR = var.ind.slp/(var.ind.slp+var.ind.yr.slp)
mean(LinearR);round(HPDinterval(as.mcmc(LinearR),0.95),3)

QuadR = var.ind.slp2/(var.ind.slp2+var.ind.yr.slp2)
mean(QuadR);round(HPDinterval(as.mcmc(QuadR),0.95),3)

# XX Backup workspace ####

save.image("C:/Users/Megaport/OneDrive/Project_data/heritability_slopes/repeatability_analysis.RData")
#load("C:/Users/Megaport/OneDrive/Project_data/heritability_slopes/repeatability_analysis.RData")

