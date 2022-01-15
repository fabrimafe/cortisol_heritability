source("model_brms_header2.R")
library(tidyverse)

res<-list()
for ( izzy in 1:100 )
	{
	mydata$group.year<-paste(mydata$community,mydata$Year,sep=".")
	mydata.random<-mydata
	A.random<-A
	old_names<-c()
	new_names<-c()
	community_names<-c()
	for ( ic in as.character(unique(mydata$community)))
		{
		old_names<-c(old_names,unique(mydata$name.code[mydata$community==ic]))
		new_names<-c(new_names,sample(unique(mydata$name.code[mydata$community==ic])))
		}
	#when an individual is in multiple communities,to avoid messing with the pedigree structure, I shuffle it only for one community 
        duplicatednames<-old_names[duplicated(old_names)]
        for (inames in duplicatednames){
                order_discarded_indexname<-sample(1:length(which(old_names==inames)))
		discarded_indexname<-which(old_names==inames)[order_discarded_indexname]
                old_names<-old_names[-discarded_indexname]
                newpositionfornametoretain<-which(new_names==inames)[order_discarded_indexname]
		new_names[newpositionfornametoretain]<-new_names[discarded_indexname]
		new_names<-new_names[-discarded_indexname]
		}
	A_random_labels<-rownames(A)
	for (irr in 1:length(old_names))
		{
		A_random_labels[rownames(A)==old_names[irr]]<-new_names[irr]
		}
	#first individuals not shuffled because only parents without data
	#A_random_labels==rownames(A)
	rownames(A.random)<-A_random_labels
	colnames(A.random)<-A_random_labels
	#A_random_labels[A_random_labels==rownames(A)]

#v7 (old v5 but improved)
	table_mot<- mydata %>% group_by(mother_name) %>% dplyr::select(name.code,mother_name) %>% unique %>% dplyr::select(mother_name) %>% table
	siblings<-mydata %>% group_by(mother_name) %>% dplyr::select(name.code,mother_name) %>% unique %>% ungroup %>% filter( ! mother_name %in% names(table_mot)[table_mot==1]) %>% arrange(by=mother_name)
	new_mother_names<-mydata$name.code
	for (isib in 2:nrow(siblings))
		{
		if (siblings$mother_name[isib]==siblings$mother_name[isib-1])
			{
			new_mother_names[new_mother_names==siblings$name.code[isib]]<-siblings$name.code[isib-1]
    			}
		}
	mydata.random$mother_name<-new_mother_names


model_brms.3.16_random <- brm(log.cort ~ 1 + (1+z.time+quad)*z.age_at_sample +
(1+z.time+quad)*sin(seasonDate) + (1+z.time+quad)*cos(seasonDate) +
(1+z.time+quad)*z.sex.ratio+(1+z.time+quad)*z.comm_size+
(1+z.time+quad)*z.demo.code + (1+z.time+quad)*z.lcms.code +
(1+z.time+quad|group.year)+(1+z.time+quad|database.code)+
(1+z.time+quad|ID_year)+(1+z.time+quad|name.code)+(1+z.time+quad|animal)+(1+z.time+quad|mother_name),
                prior = c(prior(normal(0, 1), class = Intercept),
                prior(normal(0, 1), class = b)),
                  data = mydata.random,
                  family = gaussian(),
                  cov_ranef = list(animal = A.random),
                  chains = 4,
                  cores = 4,
		  warmup=4000,
                  iter = 8000, control=list(adapt_delta=0.99,max_treedepth = 20))
res[[izzy]]<-model_brms.3.16_random
model_brms.3.16_m53_random.v9.a<-res
save(model_brms.3.16_m53_random.v9.a,file="model_brms.3.16_m53.random.v9.a.RData")
}

