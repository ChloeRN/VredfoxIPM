library(metafor)
library(tidyverse)

set.seed(2023) # simmulation for imputing SDs 

## Read in the data

RedFox_LiteratureData <- read_csv("Data/RedFox_LiteratureData.csv")

## subset by include_MetaAnalysis

RedFox_MA<-RedFox_LiteratureData |> 
  filter(include_MetaAnalysis==1)

## We need SD but we can take the SD across all studies and use a single value for each age class
## We can improve this though!

RedFox_MA_1SD<-RedFox_MA |>
  mutate(Age0_mean=mean(Age0, na.rm=TRUE)) |>
  mutate(Age1_mean=mean(Age1, na.rm=TRUE)) |>
  mutate(Age2_mean=mean(Age2, na.rm=TRUE)) |>
  mutate(Age3_mean=mean(Age3, na.rm=TRUE)) |>
  mutate(Age4_mean=mean(Age4, na.rm=TRUE))

## imputation of SDs?

RedFox_MAsd<-RedFox_MA |> 
  rowwise() |> 
  mutate(Age0_sd=sd(TruncatedNormal::rtnorm(maxN_all,Age0, lb=0, ub=1)),
         Age1_sd=sd(TruncatedNormal::rtnorm(maxN_all,Age1, lb=0, ub=1)),
         Age2_sd=sd(TruncatedNormal::rtnorm(maxN_all,Age2, lb=0, ub=1)),
         Age3_sd=sd(TruncatedNormal::rtnorm(maxN_all,Age3, lb=0, ub=1))) 

RedFox_MAsd$Age4_sd<-NULL

for (i in 1:dim(RedFox_MAsd)[1]){
  if(is.na(RedFox_MAsd$Age4[i])){
    RedFox_MAsd$Age4_sd[i]=NA
  } else{
    RedFox_MAsd$Age4_sd[i]=sd(TruncatedNormal::rtnorm(RedFox_MAsd$maxN_all[i],RedFox_MAsd$Age4[i], lb=0, ub=1))
  }
    
}

## Weighted metaAnalysis
Age0_MA<-rma(yi=Age0,vi=Age0_mean, weights=maxN_all, data=RedFox_MA_1SD)
Age1_MA<-rma(yi=Age1,vi=Age1_mean, weights=maxN_all, data=RedFox_MA_1SD)
Age2_MA<-rma(yi=Age2,vi=Age2_mean, weights=maxN_all, data=RedFox_MA_1SD)
Age3_MA<-rma(yi=Age3,vi=Age3_mean, weights=maxN_all, data=RedFox_MA_1SD)
Age4_MA<-rma(yi=Age4,vi=Age4_mean, weights=maxN_all, data=RedFox_MA_1SD)

forest(Age0_MA)
forest(Age1_MA)
forest(Age2_MA)
forest(Age3_MA)
forest(Age4_MA)


Age0_MA_sd<-rma(yi=Age0,vi=Age0_sd, weights=maxN_all, data=RedFox_MAsd)
Age1_MA_sd<-rma(yi=Age1,vi=Age1_sd, weights=maxN_all, data=RedFox_MAsd)
Age2_MA_sd<-rma(yi=Age2,vi=Age2_sd, weights=maxN_all, data=RedFox_MAsd)
Age3_MA_sd<-rma(yi=Age3,vi=Age3_sd, weights=maxN_all, data=RedFox_MAsd)
Age4_MA_sd<-rma(yi=Age4,vi=Age4_sd, weights=maxN_all, data=RedFox_MAsd)

forest(Age0_MA_sd)
forest(Age1_MA_sd)
forest(Age2_MA_sd)
forest(Age3_MA_sd)
forest(Age4_MA_sd)

##### 
## make a data table

out_dat<-tibble(Estimate=c(Age0_MA$beta,
                               Age1_MA$beta,
                               Age2_MA$beta,
                               Age3_MA$beta,
                               Age4_MA$beta,
                               Age0_MA_sd$beta,
                               Age1_MA_sd$beta,
                               Age2_MA_sd$beta,
                               Age3_MA_sd$beta,
                               Age4_MA_sd$beta),
                    se=c(Age0_MA$se,
                         Age1_MA$se,
                         Age2_MA$se,
                         Age3_MA$se,
                         Age4_MA$se,
                         Age0_MA_sd$se,
                         Age1_MA_sd$se,
                         Age2_MA_sd$se,
                         Age3_MA_sd$se,
                         Age4_MA_sd$se),
                    ci.lb=c(
                      Age0_MA$ci.lb,
                       Age1_MA$ci.lb,
                       Age2_MA$ci.lb,
                       Age3_MA$ci.lb,
                       Age4_MA$ci.lb,
                       Age0_MA_sd$ci.lb,
                       Age1_MA_sd$ci.lb,
                       Age2_MA_sd$ci.lb,
                       Age3_MA_sd$ci.lb,
                       Age4_MA_sd$ci.lb),
                    ci.ub=c(
                      Age0_MA$ci.ub,
                      Age1_MA$ci.ub,
                      Age2_MA$ci.ub,
                      Age3_MA$ci.ub,
                      Age4_MA$ci.ub,
                      Age0_MA_sd$ci.ub,
                      Age1_MA_sd$ci.ub,
                      Age2_MA_sd$ci.ub,
                      Age3_MA_sd$ci.ub,
                      Age4_MA_sd$ci.ub),
                    SD_imputed=c(rep("SD1", 5), rep("SDimputed", 5)),
                    Age_class=c(rep(c("Age0","Age1","Age2", "Age3", "Age4"),2)))


write.csv(out_dat, "Data/weighted_MA.csv") 




