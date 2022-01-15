source("model_brms_header2.R")
mydata_rank<-read.csv("rank_cat_data.csv"); mydata_rank<-mydata_rank %>% dplyr::select(sample.id,rank_cat); mydata<-left_join(mydata_rank,mydata) 

#+fixed rank
model_brms.3.16_m61 <- brm(log.cort ~ 1 + (1+z.time+quad)*z.age_at_sample +(1+z.time+quad)*z.lcms.code+
(1+z.time+quad)*sin(seasonDate) + (1+z.time+quad)*cos(seasonDate) +
(1+z.time+quad)*z.sex.ratio+(1+z.time+quad)*z.comm_size+
(1+z.time+quad)*z.demo.code + (1+z.time+quad)*rank_cat+
(1+z.time+quad|group.year)+(1+z.time+quad|database.code)+
(1+z.time+quad|ID_year)+(1+z.time+quad|name.code)+(1+z.time+quad|animal)+(1+z.time+quad|mother_name),
                prior = c(prior(normal(0, 1), class = Intercept),
                prior(normal(0, 1), class = b)),
                  data = mydata,
                  family = gaussian(),
                  cov_ranef = list(animal = A),
                  chains = 4,
                  cores = 4,
                  warmup = 4000,
                  iter = 8000, control=list(adapt_delta=0.99,max_treedepth = 20),seed=1)
save(model_brms.3.16_m61,file="model_brms.3.16_m61.RData")

