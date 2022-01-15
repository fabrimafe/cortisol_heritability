library(dplyr)
library(stringr)
library(tidyverse)
library(lme4)

library(MCMCglmm)
library(GeneticsPed)
library(brms)
library(data.table)
library(lme4)
library(MasterBayes)
library(invgamma)
library(pedigreemm)
library(AGHmatrix)

###############################
# data preparation from Paddy
###############################
{

d = read.csv("rank_cat_data.csv", header = T, colClasses = "character")

# 2 Separate sheets for repeatability ####

sort(unique(d$demographic))

adult_males = filter(d, demographic == "adult_male")
adult_females = filter(d, demographic == "cycling" | demographic == "lactating")
juveniles = filter(d, demographic == "juvenile_female" | demographic == "juvenile_male")

nrow(adult_males)+nrow(adult_females)+nrow(juveniles) # all individuals covered

# 3 Adult males sample check ####

# for adult males, we don't want individuals without rank #I remove this, assuming that Paddy already filtered them out -> I should ask

#adult_males$rank = as.numeric(as.character(adult_males$rank))
#xx = filter(adult_males, is.na(rank)) # individuals before we had 5 pant grunts
#adult_males = filter(adult_males, !is.na(rank))

# going to double check all individuals have at least 1 year with 3 samples when
# split up like this

list.ad.m.samp=unique(adult_males$name_code) #49 to start with

#create a matrix with name 

tab.ad.male=data.frame(matrix(0,nrow=length(list.ad.m.samp),ncol=1))

colnames(tab.ad.male)=c("name_code")

tab.ad.male$name_code=list.ad.m.samp

# Calculate the sample size for each individual 

tab.ad.male$n.samples=NA
tab.ad.male$n.years=NA

##create an individual sample list

samp.data.adu.m=adult_males
names(samp.data.adu.m)

for (i in 1:nrow(tab.ad.male)){
  
  xx=which(samp.data.adu.m$name_code==tab.ad.male$name_code[i])
  
  tab.ad.male$n.samples[i]=length(xx)
  tab.ad.male$n.years[i]=length(unique(adult_males$Year[xx]))
  
  
}

##check ratio samples per year (average)

tab.ad.male$av.sample.p.year=tab.ad.male$n.samples/tab.ad.male$n.years

###check how many individuals have 3 or more average samples per year 
yy = nrow(tab.ad.male)
zz = length(which(tab.ad.male$av.sample.p.year>=3))#47

### then with more than 3 per year for at least one years
zz = length(which(tab.ad.male$av.sample.p.year>=3 & tab.ad.male$n.years>=1)) #47

# recalculate the time in hours
# remove samples with no time
xx=which(samp.data.adu.m$time_corr%in%c("","Check tube","1", "unk"))

#samp.data.adu.m=samp.data.adu.m[-xx,] not needed as xx is empty
samp.data.adu.m=droplevels(subset(samp.data.adu.m, !coll_hour == "unk"))

range(as.numeric(as.character(samp.data.adu.m$coll_hour)))  ## ok
range(as.numeric(as.character(samp.data.adu.m$coll_min))) ##ok
samp.data.adu.m$time.h=as.numeric(as.character(samp.data.adu.m$coll_hour)) + as.numeric(as.character(samp.data.adu.m$coll_min))/60
range(samp.data.adu.m$time.h)

##remove the smaples for which we don't have a time

samp.data.adu.m=subset(samp.data.adu.m,is.na(samp.data.adu.m$time.h)==
                         F&samp.data.adu.m$name_code%in%tab.ad.male$name_code) #3214 lines

tab.samp.size=aggregate(samp.data.adu.m$name_code,by=list(samp.data.adu.m$Year,samp.data.adu.m$name_code),
                        length)

colnames(tab.samp.size)=c("year","name_code","n.samp")

#get the minimum and maximum time a sample was collected

xx1=aggregate(samp.data.adu.m$time.h,by=list(samp.data.adu.m$Year,samp.data.adu.m$name_code),min)
xx2=aggregate(samp.data.adu.m$time.h,by=list(samp.data.adu.m$Year,samp.data.adu.m$name_code),max)


tab.samp.size$min.time.h=xx1$x
tab.samp.size$max.time.h=xx2$x

###calculate years for which we have at least 6 hours span between the earliest and latests sample

tab.samp.size$time.spread=tab.samp.size$max.time.h-tab.samp.size$min.time.h  ##range from 0 to 12.15

tab.samp.size.6h=subset(tab.samp.size,tab.samp.size$time.spread>6&tab.samp.size$n.samp>2)

range(tab.samp.size.6h$n.samp) #3 to 110

#total number of inds
length(unique(tab.samp.size.6h$name_code)) ##48

####create a subset of the dataset with only the right males and the right years

samp.data.adu.m$name_year=NA

for (i in 1:nrow(samp.data.adu.m)){
  samp.data.adu.m$name_year[i]=paste(samp.data.adu.m$name_code[i],samp.data.adu.m$Year[i],sep="_")
}

tab.samp.size.6h$name_year=NA

for (i in 1:nrow(tab.samp.size.6h)){
  tab.samp.size.6h$name_year[i]=paste(tab.samp.size.6h$name_code[i],tab.samp.size.6h$year[i],sep="_")
}

# now we need to add our threshold for number of years

tab.samp.size.6h$no.of.years = NA

tab.samp.size.6h <- within(tab.samp.size.6h, years <- ave(name_code, list(name_code), FUN=length))
# subset those without enough years
tab.samp.size.6h = droplevels(subset(tab.samp.size.6h, years >=1))

samp.data.adu.m2=subset(samp.data.adu.m, samp.data.adu.m$name_year%in%unique(tab.samp.size.6h$name_year))

nrow(samp.data.adu.m2) #3192
xx = levels(as.factor(samp.data.adu.m2$name_code)) #48
adult_males = samp.data.adu.m2

# 3 Adult females sample check ####

# for adult females, we don't want individuals without reproductive state

sort(unique(adult_females$repro.stat))

adult_females = filter(adult_females, !is.na(repro.stat))
adult_females = filter(adult_females, !repro.stat == "")

# going to double check all individuals have at least 1 year with 3 samples when
# split up like this

list.ad.m.samp=unique(adult_females$name_code) #70 to start with

#create a matrix with name 

tab.ad.male=data.frame(matrix(0,nrow=length(list.ad.m.samp),ncol=1))

colnames(tab.ad.male)=c("name_code")

tab.ad.male$name_code=list.ad.m.samp

# Calculate the sample size for each individual 

tab.ad.male$n.samples=NA
tab.ad.male$n.years=NA

##create an individual sample list

samp.data.adu.m=adult_females
names(samp.data.adu.m)

for (i in 1:nrow(tab.ad.male)){
  
  xx=which(samp.data.adu.m$name_code==tab.ad.male$name_code[i])
  
  tab.ad.male$n.samples[i]=length(xx)
  tab.ad.male$n.years[i]=length(unique(adult_males$Year[xx]))
  
  
}

##check ratio samples per year (average)

tab.ad.male$av.sample.p.year=tab.ad.male$n.samples/tab.ad.male$n.years

###check how many individuals have 3 or more average samples per year 
yy = nrow(tab.ad.male)
zz = length(which(tab.ad.male$av.sample.p.year>=3))#60

### then with more than 3 per year for at least one years
zz = length(which(tab.ad.male$av.sample.p.year>=3 & tab.ad.male$n.years>=1)) #60

# recalculate the time in hours
# remove samples with no time
xx=which(samp.data.adu.m$time_corr%in%c("","Check tube","1", "unk"))

#samp.data.adu.m=samp.data.adu.m[-xx,] not needed as xx is empty
samp.data.adu.m=droplevels(subset(samp.data.adu.m, !coll_hour == "unk"))

range(as.numeric(as.character(samp.data.adu.m$coll_hour)))  ## ok
range(as.numeric(as.character(samp.data.adu.m$coll_min))) ##ok
samp.data.adu.m$time.h=as.numeric(as.character(samp.data.adu.m$coll_hour)) + as.numeric(as.character(samp.data.adu.m$coll_min))/60
range(samp.data.adu.m$time.h)

##remove the smaples for which we don't have a time

samp.data.adu.m=subset(samp.data.adu.m,is.na(samp.data.adu.m$time.h)==
                         F&samp.data.adu.m$name_code%in%tab.ad.male$name_code) #1758 lines

tab.samp.size=aggregate(samp.data.adu.m$name_code,by=list(samp.data.adu.m$Year,samp.data.adu.m$name_code),
                        length)

colnames(tab.samp.size)=c("year","name_code","n.samp")

#get the minimum and maximum time a sample was collected

xx1=aggregate(samp.data.adu.m$time.h,by=list(samp.data.adu.m$Year,samp.data.adu.m$name_code),min)
xx2=aggregate(samp.data.adu.m$time.h,by=list(samp.data.adu.m$Year,samp.data.adu.m$name_code),max)


tab.samp.size$min.time.h=xx1$x
tab.samp.size$max.time.h=xx2$x

###calculate years for which we have at least 6 hours span between the earliest and latests sample

tab.samp.size$time.spread=tab.samp.size$max.time.h-tab.samp.size$min.time.h  ##range from 0 to 12.15

tab.samp.size.6h=subset(tab.samp.size,tab.samp.size$time.spread>6&tab.samp.size$n.samp>2)

range(tab.samp.size.6h$n.samp) #3 to 34

#total number of inds
length(unique(tab.samp.size.6h$name_code)) ##69

####create a subset of the dataset with only the right males and the right years

samp.data.adu.m$name_year=NA

for (i in 1:nrow(samp.data.adu.m)){
  samp.data.adu.m$name_year[i]=paste(samp.data.adu.m$name_code[i],samp.data.adu.m$Year[i],sep="_")
}

tab.samp.size.6h$name_year=NA

for (i in 1:nrow(tab.samp.size.6h)){
  tab.samp.size.6h$name_year[i]=paste(tab.samp.size.6h$name_code[i],tab.samp.size.6h$year[i],sep="_")
}

# now we need to add our threshold for number of years

tab.samp.size.6h$no.of.years = NA

tab.samp.size.6h <- within(tab.samp.size.6h, years <- ave(name_code, list(name_code), FUN=length))
# subset those without enough years
tab.samp.size.6h = droplevels(subset(tab.samp.size.6h, years >=1))

samp.data.adu.m2=subset(samp.data.adu.m, samp.data.adu.m$name_year%in%unique(tab.samp.size.6h$name_year))

nrow(samp.data.adu.m2) #1756
xx = levels(as.factor(samp.data.adu.m2$name_code)) #69
adult_females = samp.data.adu.m2

# 5 Juvenile sample check ####

list.ad.m.samp=unique(juveniles$name_code) #174 to start with

#create a matrix with name 

tab.ad.male=data.frame(matrix(0,nrow=length(list.ad.m.samp),ncol=1))

colnames(tab.ad.male)=c("name_code")

tab.ad.male$name_code=list.ad.m.samp

# Calculate the sample size for each individual 

tab.ad.male$n.samples=NA
tab.ad.male$n.years=NA

##create an individual sample list

samp.data.adu.m=juveniles
names(samp.data.adu.m)

for (i in 1:nrow(tab.ad.male)){
  
  xx=which(samp.data.adu.m$name_code==tab.ad.male$name_code[i])
  
  tab.ad.male$n.samples[i]=length(xx)
  tab.ad.male$n.years[i]=length(unique(juveniles$Year[xx]))
  
  
}

##check ratio samples per year (average)

tab.ad.male$av.sample.p.year=tab.ad.male$n.samples/tab.ad.male$n.years

###check how many individuals have 3 or more average samples per year 
yy = nrow(tab.ad.male)
zz = length(which(tab.ad.male$av.sample.p.year>=3))#67

### then with more than 3 per year for at least two years
zz = length(which(tab.ad.male$av.sample.p.year>=3 & tab.ad.male$n.years>=1)) #67

# recalculate the time in hours
# remove samples with no time
xx=which(samp.data.adu.m$time_corr%in%c("","Check tube","1", "unk"))

#samp.data.adu.m=samp.data.adu.m[-xx,] not needed as xx is empty
samp.data.adu.m=droplevels(subset(samp.data.adu.m, !coll_hour == "unk"))

range(as.numeric(as.character(samp.data.adu.m$coll_hour)))  ## ok
range(as.numeric(as.character(samp.data.adu.m$coll_min))) ##ok
samp.data.adu.m$time.h=as.numeric(as.character(samp.data.adu.m$coll_hour)) + as.numeric(as.character(samp.data.adu.m$coll_min))/60
range(samp.data.adu.m$time.h)

##remove the smaples for which we don't have a time

samp.data.adu.m=subset(samp.data.adu.m,is.na(samp.data.adu.m$time.h)==F&samp.data.adu.m$name_code%in%tab.ad.male$name_code) #2678 lines

tab.samp.size=aggregate(samp.data.adu.m$name_code,by=list(samp.data.adu.m$Year,samp.data.adu.m$name_code),length)

colnames(tab.samp.size)=c("year","name_code","n.samp")

#get the minimum and maximum time a sample was collected

xx1=aggregate(samp.data.adu.m$time.h,by=list(samp.data.adu.m$Year,samp.data.adu.m$name_code),min)
xx2=aggregate(samp.data.adu.m$time.h,by=list(samp.data.adu.m$Year,samp.data.adu.m$name_code),max)


tab.samp.size$min.time.h=xx1$x
tab.samp.size$max.time.h=xx2$x

###calculate years for which we have at least 6 hours span between the earliest and latests sample

tab.samp.size$time.spread=tab.samp.size$max.time.h-tab.samp.size$min.time.h  ##range from 0 to 12.15

tab.samp.size.6h=subset(tab.samp.size,tab.samp.size$time.spread>6&tab.samp.size$n.samp>2)

range(tab.samp.size.6h$n.samp) #3 to 45

#total number of inds
length(unique(tab.samp.size.6h$name_code)) ##65

####create a subset of the dataset with only the right males and the right years

samp.data.adu.m$name_year=NA

for (i in 1:nrow(samp.data.adu.m)){
  samp.data.adu.m$name_year[i]=paste(samp.data.adu.m$name_code[i],samp.data.adu.m$Year[i],sep="_")
}

tab.samp.size.6h$name_year=NA

for (i in 1:nrow(tab.samp.size.6h)){
  tab.samp.size.6h$name_year[i]=paste(tab.samp.size.6h$name_code[i],tab.samp.size.6h$year[i],sep="_")
}

# now we need to add our threshold for number of years

tab.samp.size.6h$no.of.years = NA

tab.samp.size.6h <- within(tab.samp.size.6h, years <- ave(name_code, list(name_code), FUN=length))
# subset those without enough years
tab.samp.size.6h = droplevels(subset(tab.samp.size.6h, years >=2))

samp.data.adu.m2=subset(samp.data.adu.m, samp.data.adu.m$name_year%in%unique(tab.samp.size.6h$name_year))

nrow(samp.data.adu.m2) #775
xx = levels(as.factor(samp.data.adu.m2$name_code)) #32
juveniles = samp.data.adu.m2

all_inds = d

#FM:A few warnings of objects never created. It should not be a problem
#rm(samp.data.adu.m, samp.data.adu.m2, tab.ad.male, tab.samp.size, tab.samp.size.6h, xx1, xx2,
#   f_sampled, fathers, fathers_with_multiple_siblings_sampled, i, inds_sampled, list.ad.m.samp,
#   m_sampled, mothers, mothers_with_multiple_siblings_sampled, ny, sampled, x, xx, yy, zz)

   res.Int<-function(x){
  smod<-sim(x, 1000)
  R<-data.frame(smod@fixef)
  R$V.fem=apply(smod@ranef$name.code,1,var)  
  R$V.fem.year<-apply(smod@ranef$fem_year,1,var)
  R$V.year<-apply(smod@ranef$Year,1,var)
  R$V.res<-smod@sigma^2
  V.total<-R$V.fem + R$V.fem.year + R$V.year + R$V.res 
  R$long.term.R<-R$V.fem/(V.total)
  R$short.term.R<-(R$V.fem + R$V.fem.year)/(V.total)
  V.model1=data.frame(mean=round(apply(R, 2, mean),2), lCI=round(apply(R, 2, quantile, 0.025),2), uCI=round(apply(R, 2, quantile, 0.975),2))
  V.model1
}

res.Slope<-function(x){
  smod<-sim(x, 1000)
  R<-data.frame(smod@fixef)
  R$V.fem.Int=apply(smod@ranef$name.code[,,1],1,var)
  R$V.fem.Slope=apply(smod@ranef$name.code[,,2],1,var)
  R$V.fem.Slope2=apply(smod@ranef$name.code[,,3],1,var)
  R$V.fem.year.Int<-apply(smod@ranef$fem_year[,,1],1,var)
  R$V.fem.year.Slope<-apply(smod@ranef$fem_year[,,2],1,var)
  R$V.fem.year.Slope2<-apply(smod@ranef$fem_year[,,3],1,var)
  R$V.year<-apply(smod@ranef$Year,1,var)
  R$V.res<-smod@sigma^2
  V.total<-R$V.fem.Int + R$V.fem.year.Int + R$V.year + R$V.res
  R$long.term.R<-R$V.fem.Int/(V.total)
  R$short.term.R<-(R$V.fem.Int + R$V.fem.year.Int)/(V.total)
  R$R.Intercept<-R$V.fem.Int/(R$V.fem.Int + R$V.fem.year.Int)  ##Reaction norm repeatability for intercepts
  R$R.Slope<-R$V.fem.Slope/(R$V.fem.Slope + R$V.fem.year.Slope) ##Reaction norm repeatability for slopes
  #here additional calculation for slope over time squared:
  R$R.Slope2<-R$V.fem.Slope2/(R$V.fem.Slope2 + R$V.fem.year.Slope2)
  V.model1=data.frame(mean=round(apply(R, 2, mean),2), lCI=round(apply(R, 2, quantile, 0.025),2), uCI=round(apply(R, 2, quantile, 0.975),2))
  V.model1
}

res.Slope3<-function(x){
  smod<-sim(x, 1000)
  R<-data.frame(smod@fixef)
  R$V.fem.Int=apply(smod@ranef$name.code[,,1],1,var)
  R$V.fem.Slope=apply(smod@ranef$name.code[,,2],1,var)
  R$V.fem.year.Int<-apply(smod@ranef$fem_year[,,1],1,var)
  R$V.fem.year.Slope<-apply(smod@ranef$fem_year[,,2],1,var)
  R$V.year<-apply(smod@ranef$Year,1,var)
  R$V.res<-smod@sigma^2
  V.total<-R$V.fem.Int + R$V.fem.year.Int + R$V.year + R$V.res
  R$long.term.R<-R$V.fem.Int/(V.total)
  R$short.term.R<-(R$V.fem.Int + R$V.fem.year.Int)/(V.total)
  R$R.Intercept<-R$V.fem.Int/(R$V.fem.Int + R$V.fem.year.Int)  ##Reaction norm repeatability for intercepts
  R$R.Slope<-R$V.fem.Slope/(R$V.fem.Slope + R$V.fem.year.Slope) ##Reaction norm repeatability for slopes
  #here additional calculation for slope over time squared:
  V.model1=data.frame(mean=round(apply(R, 2, mean),2), lCI=round(apply(R, 2, quantile, 0.025),2), uCI=round(apply(R, 2, quantile, 0.975),2))
  V.model1
}

test_data = all_inds
names(test_data)
test_data$community_size = as.numeric(as.character(test_data$group_size))

# Time
test_data$time.h = as.numeric(as.character(test_data$time.h))
test_data = droplevels(subset(test_data, !is.na(time.h)))
range(test_data$time.h)

# Date

test_data$rrdate = as.numeric(as.Date(as.character(test_data$rdate)))
range(test_data$rrdate)
test_data$seasonDate=2*pi*as.numeric(test_data$rrdate)/365.25
range(test_data$seasonDate)

# control variables

test_data$sex_ratio = as.numeric(as.character(test_data$male_to_fems))
range(test_data$sex_ratio)
test_data$age_at_sample = as.numeric(as.character(test_data$age_at_sample))
range(test_data$age_at_sample)
test_data$community_size = as.numeric(as.character(test_data$community_size))
range(test_data$community_size)
test_data$Year = as.numeric(as.character(test_data$Year))
range(test_data$Year)
test_data$Month = as.numeric(as.character(test_data$Month))
range(test_data$Month)

test_data$project_person = as.factor(as.character(test_data$project_person))

test_data$community = as.factor(as.character(test_data$community))
test_data$site = as.factor(as.character(test_data$site))

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

# dummy the factors #FM: never clear why he does this. They are factors, why ever z rescaling them?! I should ask him
names(d)
sort(unique(test_data$project_person))
database.code=as.numeric(as.character(factor(test_data$project_person, levels=c("Anna","Cedric","Christina","Coco","ERC",
                                                                                "Erin","Liran","Paddy","Pawel","Prince", "Resi",
                                                                                "RomanCathy","Tobias","Verena","Virgile","Zinta"), 
                                             labels=c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13","14","15")))) 
test_data = data.frame(test_data, database.code)
rm(database.code)

test_data$z.database.code = as.vector((test_data$database.code - mean(test_data$database.code, na.rm=TRUE))
                                      /(2*sd(test_data$database.code, na.rm=TRUE)))
test_data$z.database.code = as.vector(as.factor((test_data$z.database.code)))


test_data<-test_data[test_data$community!="Middle",]
#Exclude Middle
levels(as.factor(test_data$community))
group_code=as.numeric(as.character(factor(test_data$community, levels=c("East","North","Sonso","South","Waibira"), 
                                          labels=c("0", "1", "2", "3", "4")))) 
test_data=data.frame(test_data, group_code)
rm(group_code)
test_data$z.group.code = as.vector((test_data$group_code - mean(test_data$group_code, na.rm=TRUE))
                                   /(2*sd(test_data$group_code, na.rm=TRUE)))
test_data$z.group.code = as.vector(as.factor((test_data$z.group.code)))

levels(as.factor(test_data$LCMS))
lcms_code=as.numeric(as.character(factor(test_data$LCMS, levels=c("new","old"), 
                                         labels=c("0", "1")))) 
test_data=data.frame(test_data, lcms_code)
rm(lcms_code)
test_data$z.lcms.code = as.vector((test_data$lcms_code - mean(test_data$lcms_code, na.rm=TRUE))
                                  /(2*sd(test_data$lcms_code, na.rm=TRUE)))
test_data$z.lcms.code = as.vector(as.factor((test_data$z.lcms.code)))

site.code=as.numeric(as.character(factor(test_data$site, levels=c("Budongo", "Tai"), labels=c("0", "1")))) 

test_data=data.frame(test_data, site.code)
rm(site.code)

test_data$z.site.code = as.vector((test_data$site.code - mean(test_data$site.code, na.rm=TRUE))/
                                    (2*sd(test_data$site.code, na.rm=TRUE)))
test_data$z.site.code = as.vector(as.factor((test_data$z.site.code)))

names(test_data)

test_data$group_year=as.character(apply(test_data[, c("community", "Year")], 1, paste, collapse="-"))
test_data$group_year = as.factor(test_data$group_year)

test_data$ID_year=as.character(apply(test_data[, c("name_code", "Year")], 1, paste, collapse="-"))
test_data$ID_year = as.factor(test_data$ID_year)

sort(unique(test_data$demographic))
demo.code=as.numeric(as.character(factor(test_data$demographic, levels=c("adult_male","cycling","lactating","juvenile_female","juvenile_male"), 
                                             labels=c("0", "1", "2", "3", "4")))) 
test_data = data.frame(test_data, demo.code)
rm(demo.code)

test_data$z.demo.code = as.vector((test_data$demo.code - mean(test_data$demo.code, na.rm=TRUE))
                                      /(2*sd(test_data$demo.code, na.rm=TRUE)))
test_data$z.demo.code = as.vector(as.factor((test_data$z.demo.code)))

#response variable:

test_data$Cortisol_ng_ml_SG = as.numeric(as.character(test_data$Cortisol_ng_ml_SG))
range(test_data$Cortisol_ng_ml_SG)

test_data$log.cort = log(test_data$Cortisol_ng_ml_SG)
hist(test_data$log.cort)

library(car)
# vifs
xx=lm(log.cort ~ z.time + z.sex.ratio + z.demo.code + z.age_at_sample + z.group.code + z.lcms.code + z.site.code + z.comm_size +
        sin(seasonDate) + cos(seasonDate), data=test_data)
#vif(xx)
if (FALSE)
{
# perfectly collinear in there so remove site
xx=lm(log.cort ~ z.time + z.sex.ratio + z.demo.code + z.age_at_sample + z.group.code + z.lcms.code + z.comm_size +
        sin(seasonDate) + cos(seasonDate), data=test_data)
#max(vif(xx)) #172.21

xx=lm(log.cort ~ z.time + z.sex.ratio + z.demo.code + z.age_at_sample + z.group.code + z.lcms.code + 
        sin(seasonDate) + cos(seasonDate), data=test_data)
#max(vif(xx)) #5.44

xx=lm(log.cort ~ z.time + z.sex.ratio + z.demo.code + z.age_at_sample + z.comm_size + z.lcms.code + 
        sin(seasonDate) + cos(seasonDate), data=test_data)
#max(vif(xx)) #3.04

xx=lm(log.cort ~ z.time + z.sex.ratio + z.demo.code + z.age_at_sample + z.comm_size + z.lcms.code + z.site.code +
        sin(seasonDate) + cos(seasonDate), data=test_data)
#max(vif(xx))#12.89
}
}
###############################
#define functions
###############################
{

#mydata<-mydata[mydata$community=="South" | mydata$community=="North",]
data2pedigree<-function(mydata) {
pedigree<-data.frame(animal=mydata$animal,father=mydata$father_name,mother=mydata$mother_name)
pedigree$animal<-as.character(pedigree$animal)
pedigree$father<-as.character(pedigree$father)
pedigree$mother<-as.character(pedigree$mother)
pedigree[pedigree==""]<-NA
pedigree0<-pedigree[is.na(pedigree$father) & is.na(pedigree$mother),]
#properly formatting pedigree
pedigree0<-rbind(pedigree0,pedigree[is.na(pedigree$father) & !is.na(pedigree$mother),])
pedigree0<-rbind(pedigree0,pedigree[!is.na(pedigree$father) & is.na(pedigree$mother),])
pedigree<-rbind(pedigree0,pedigree[!is.na(pedigree$father) & !is.na(pedigree$mother),])
pedigree<-unique(pedigree)
pedigree<-pedigree[pedigree$father%in%mydata$animal || is.na(pedigree$father),]
pedigree<-pedigree[pedigree$mother%in%mydata$animal || is.na(pedigree$mother),]
#pedigree0<-cbind(pedigree$father[!pedigree$father %in% pedigree$animal & !is.na(pedigree$father)],NA,NA)
if (length(pedigree$father[!pedigree$father %in% pedigree$animal & !is.na(pedigree$father)])!=0){
pedigree0<-cbind(pedigree$father[!pedigree$father %in% pedigree$animal & !is.na(pedigree$father)],NA,NA)
}
colnames(pedigree0)<-c("animal","father","mother")
pedigree<-rbind(pedigree0,pedigree)
pedigree$animal<-as.character(pedigree$animal)
pedigree$father<-as.character(pedigree$father)
pedigree$mother<-as.character(pedigree$mother)
pedigree0<-cbind(pedigree$mother[!pedigree$mother %in% pedigree$animal & !is.na(pedigree$mother)],NA,NA)
colnames(pedigree0)<-c("animal","father","mother")
pedigree<-rbind(pedigree0,pedigree)
pedigree<-unique(pedigree)
pedigree<-orderPed(pedigree)
pedigree
}
data2randompedigree_allindividuals<-function(mydata, nestingcategory) 
{
#this version shuffle within community both parents and offspring. Therefore it won't work if migration between communities..and in this dataset there are some!
count_categories<-1
for (ic in unique(mydata[[nestingcategory]]) )
    {
    mydata_temp<-mydata[mydata[[nestingcategory]]==ic,]
    result<-"error";class(result)<-"try-error"
    count<-1
    while (is(result)[1]=="try-error")
        {
        pedigree<-data.frame(animal=mydata_temp$animal,father=mydata_temp$father_name,mother=mydata_temp$mother_name)
        pedigree$animal<-as.character(pedigree$animal)
        pedigree$father<-as.character(pedigree$father)
        pedigree$mother<-as.character(pedigree$mother)
        pedigree[pedigree==""]<-NA
        pedigree0<-pedigree[is.na(pedigree$father) & is.na(pedigree$mother),]
        #properly formatting pedigree
        pedigree0<-rbind(pedigree0,pedigree[is.na(pedigree$father) & !is.na(pedigree$mother),])
        pedigree0<-rbind(pedigree0,pedigree[!is.na(pedigree$father) & is.na(pedigree$mother),])
        pedigree<-rbind(pedigree0,pedigree[!is.na(pedigree$father) & !is.na(pedigree$mother),])
        pedigree<-unique(pedigree)
        pedigree<-pedigree[pedigree$father%in%mydata$animal || is.na(pedigree$father),]
        pedigree<-pedigree[pedigree$mother%in%mydata$animal || is.na(pedigree$mother),]
        if (length(pedigree$father[!pedigree$father %in% pedigree$animal & !is.na(pedigree$father)])!=0){
        pedigree0<-cbind(pedigree$father[!pedigree$father %in% pedigree$animal & !is.na(pedigree$father)],NA,NA)
        }
        colnames(pedigree0)<-c("animal","father","mother")
        pedigree<-rbind(pedigree0,pedigree)
        pedigree$animal<-as.character(pedigree$animal)
        pedigree$father<-as.character(pedigree$father)
        pedigree$mother<-as.character(pedigree$mother)
        if (length(pedigree$mother[!pedigree$mother %in% pedigree$animal & !is.na(pedigree$mother)])!=0){
        pedigree0<-cbind(pedigree$mother[!pedigree$mother %in% pedigree$animal & !is.na(pedigree$mother)],NA,NA)
        }
        colnames(pedigree0)<-c("animal","father","mother")
        pedigree<-rbind(pedigree0,pedigree)
        pedigree<-unique(pedigree)
        result<-"error"
        result<-try( { pedigree$father<-sample(pedigree$father); pedigree$mother<-sample(pedigree$mother); pedigree<-orderPed(pedigree)} )
        print(paste("try ",count))
        count<-count+1
        }
    if (count_categories==1) { res_pedigree<-pedigree } else { res_pedigree<-rbind(res_pedigree,pedigree) }
    count_categories<-count_categories+1
    }
    res_pedigree<-orderPed(res_pedigree)
}
data2randompedigree<-function(mydata, nestingcategory,pedigree)
{
#this version only shuffles sampled individuals, retaining only one community. More approximate but it works also with migrations.
pedigree_t<-pedigree
pedigree_t[[nestingcategory]]<-NA
for ( i in 1:length(pedigree[,1]) ){
pedigree_t[[nestingcategory]][i]<-mydata[[nestingcategory]][mydata$animal==as.character(pedigree[i,1])][1]
}
mylevels<-unique(pedigree_t[[nestingcategory]])
mylevels[is.na(mylevels)]<-"NA"
pedigree_t$mother<-as.character(pedigree_t$mother)
pedigree_t$father<-as.character(pedigree_t$father)
pedigree_t[is.na(pedigree_t)]<-"NA"
#library(tidyr); #replace_na(pedigree_tt,"NA") #pedigree_tt$mother %>% replace_na("NA")
result<-"error";class(result)<-"try-error"
count<-1
while (is(result)[1]=="try-error" || is(result2)[1]=="try-error")
    {
    pedigree_tt<-pedigree_t
    pedigree_tt$father<-unname(unlist(sapply(mylevels, function(x) sample(as.character(pedigree_tt$father)[pedigree_tt[[nestingcategory]]==x]))))
    pedigree_tt$mother<-unname(unlist(sapply(mylevels, function(x) sample(as.character(pedigree_tt$mother)[pedigree_tt[[nestingcategory]]==x]))))
    pedigree_tt$father<-as.factor(pedigree_tt$father)
    pedigree_tt$mother<-as.factor(pedigree_tt$mother)
    pedigree_tt$animal<-as.factor(pedigree_tt$animal)
    pedigree_tt[pedigree_tt=="NA"]<-NA
    result<-try( { pedigree_tt<-orderPed(pedigree_tt) } )
    result2<-try( { mypedigree.random<-pedigree(sire=as.character(pedigree_tt$father),dam=as.character(pedigree_tt$mother),label=as.character(pedigree_tt$animal)) } )
    print(paste("try ",count))
    count<-count+1
    }
pedigree_tt<-pedigree_tt[,1:3]
pedigree_tt
}
data2randomdata<-function(mydata, nestingcategory)
{
#this version shuffles sampled names of individuals samples. Now overall relationships should be identical. It should be the final version for the paper.
mylevels<-unique(mydata[[nestingcategory]])
mydata_t<-mydata    
newanimals<-c()
    for (x in 1:length(mylevels)){
    animal_original<-as.character(mydata_t$animal[mydata_t[[nestingcategory]]==mylevels[x]])
    animal_original_factor<-as.factor(animal_original)
    animal_original_factor_new<-animal_original_factor
    levels(animal_original_factor_new)<-sample(levels(animal_original_factor))
    newanimals<-c(newanimals,as.character(animal_original_factor_new))
    }
mydata_t$animal<-newanimals
mydata_t
}
#data2randomdata(mydata,"community")

}

mydata<-test_data
mydata$name.code<-mydata$name_code
mydata$animal<-mydata$name_code
#mothers make no sense as Paddy did (missing as ""). Not to miss data (MCMCglmm skips data if missing values for predictors), assumption is that if we do not know if siblings,
#then unrelated
mydata$z.demo.code<-as.factor(mydata$z.demo.code)
mydata$z.database.code<-as.factor(mydata$z.database.code)
mydata$z.lcms.code<-as.factor(mydata$z.lcms.code)
mydata$group.year<-paste(mydata$community,mydata$Year,sep=".")
mydata$Year<-as.factor(mydata$Year)
mydata$database.code<-as.factor(mydata$database.code)
pedigree<-data2pedigree(mydata)
mydata$mother_name[mydata$mother_name==""]<-paste0(as.character(mydata$name.code[mydata$mother_name==""]),"_mother")
#mydata$father_name[mydata$father_name==""]<-paste0(as.character(mydata$name.code[mydata$father_name==""]),"_father")
mydata$father_name[mydata$father_name==""]<-paste0(as.character(mydata$name.code[mydata$father_name==""]),"_father")
mydata$quad<-mydata$z.time^2

pedigree_Amat<-pedigree
temp<-as.character(pedigree_Amat[,2])
temp[is.na(temp)]<-0
pedigree_Amat[,2]<-temp
temp<-as.character(pedigree_Amat[,3])
temp[is.na(temp)]<-0
pedigree_Amat[,3]<-temp
temp<-as.character(pedigree_Amat[,1])
temp[is.na(temp)]<-0
pedigree_Amat[,1]<-temp
A  <- Amatrix(pedigree_Amat)
