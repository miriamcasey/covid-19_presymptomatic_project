#### taking uncertainty into account: Simulation A
########### I have collected simulation outputs in particular folders to make it easier to summarise them later
############ If rerunning/testing code, the same will need to be done - that is replace my folders with your own
######### one file for serial interval or generation time outputs setwd("~/Documents/covid_19/r_files/code/sim_b_output/si")
######## one file for transmission time relative to symptom onset outputs setwd("~/Documents/covid_19/r_files/code/sim_b_output/trs")
##### one file for bootstrapping output for transmission time relative to symptom output: setwd("~/Documents/covid_19/r_files/code/sim_b_output/strap_trs")
###### summary csv file for serial interval/generation time bootstrap summary stats: "si_df.csv"
###### summary csv file for transmission time relative to symptom onset bootstrap summary stats: "trs_df.csv"

#setwd("~/Documents/covid_19/r_files/code")
rm(list=ls())
library(tidyverse)
library(janitor)
library(epitrix)
library(fitdistrplus)
library(rriskDistributions)
library(reshape2)
library(mixdist)


#### load data
d <- read_csv("final_data_19_may.csv")
d[is.na(d$study_r),]
d[is.na(d$ci_category),]

summary(d)
d1 <- d[!duplicated(d$study_r),]
(summary(as.factor(d1$ci_category)))
d1$distribution[d1$ci_category=="missing_sd_cis"]
d1$study_r[d1$ci_category=="missing_sd_cis"]
length(unique(d1$study_r))
#### remove liao pending clarification on distribution 
d <- d[d$study_r!="liao",]
##### remove cheng as cheng does not allow for pre-symptomatic transmission
d <- d[d$study_r!="cheng",]
d$central_measure

n_primary_runs <-  10  ####1000
n_primary_samples <- 10 ####1000
n_bootstrap_runs <- 100  #####10000
n_bootstrap_samples <- 10 ##### 1000

### load data and estimate transmission time relative to symptoms

#### load incubation period from McAloon et al metanalysis:
## https://www.medrxiv.org/content/10.1101/2020.04.24.20073957v1

d[d$study_r=="mca",]

###### incubation period from mcaloon et al
sims <- data.frame (IP=numeric())
for (y in 1:n_primary_runs){
  mu = rnorm(1, 1.63, (1.75-1.51)/(2*1.96))
  sigma = rnorm(1, 0.5, (0.55 - 0.45)/(2*1.96))
  A=data.frame(IP=rlnorm(n_primary_samples, mu, sigma))
  sims=bind_rows(sims, A)
}
summary(sims)
write.csv(sims, "IP_simB.csv")

#### make loop for 8 gamma distributed serial intervals / generation times with ci's for mean and standard deviation
summary(as.factor(d$ci_category))
summary(as.factor(d1$ci_category))
df_gamma_mean_sd <- d[d$ci_category=="mean_sd_gamma" & !is.na(d$ci_category),]
names_df_gamma_mean_sd <- unique(df_gamma_mean_sd$study_r)
names_df_gamma_mean_sd 
collection_si <- data.frame()
collection_trs <- data.frame()


for (i in 1:length(names_df_gamma_mean_sd))
{
  df1 <- d[d$study_r==names_df_gamma_mean_sd[i],]
  mean_data <- df1$central_measure[df1$dist_par_r=="mean"]
  mean_data
  mean_data_lci <- df1$lower_ci[df1$dist_par_r=="mean"]  ##### take 95% confidence intervals from source publication for serial interval or generation time
  mean_data_uci <- df1$upper_ci[df1$dist_par_r=="mean"]
  sd_data <- df1$central_measure[df1$dist_par_r=="sd"]
  sd_data_lci <- df1$lower_ci[df1$dist_par_r=="sd"]
  sd_data_uci <- df1$upper_ci[df1$dist_par_r=="sd"]
  
  sims1 <- data.frame ()
  for (j in 1:n_primary_runs){
    meanx <- rnorm(1, mean_data, (mean_data_uci-mean_data_lci)/(2*1.96))
    
    #######  gamma distribution for standard deviation
    #####gamma shape=mean^2/sd^2  gamma rate=mean/sd^2
    sd_sd <- (sd_data_uci-sd_data_lci)/(2*1.96)
    sd_mean <- sd_data
    sd_shape <- sd_mean^2/(sd_sd^2)
    sd_rate <- sd_mean/(sd_sd^2)
    sdx <- rgamma(1, sd_shape, sd_rate)
    #####gamma shape=mean^2/sd^2  gamma rate=mean/sd^2
    shapex <- meanx^2/(sdx^2)
    ratex <- meanx/(sdx^2)
    B=data.frame(SI=rgamma(n_primary_samples, shapex, ratex))
    summary(B)
    sims1 <- rbind(sims1, B)
  }
  summary(sims1)
  sims1x <- as.data.frame(sims1)
  sims1x$study <- names_df_gamma_mean_sd[i]
  name_csv <- paste("SI", names_df_gamma_mean_sd[i], i, ".csv", sep="_")
  name_csv
#  setwd("~/Documents/covid_19/r_files/code/sim_b_output/si")
  write.csv(sims1x, name_csv)
  strap_si <- data.frame()
  for(k in 1:n_bootstrap_runs){
    a <- sample_n(sims1, n_bootstrap_samples, replace=TRUE)
    summary(a)
    Ax = data.frame(mean=mean(a$SI, na.rm=T), 
                    sd=sd(a$SI, na.rm=T),
                    X2.5=quantile(a$SI, probs=0.025, na.rm=T),
                    X25=quantile(a$SI, probs=0.25, na.rm=T),
                    X50=quantile(a$SI, probs=0.5, na.rm=T),
                    X75=quantile(a$SI, probs=0.75, na.rm=T),
                    X97.5=quantile(a$SI, probs=0.975, na.rm=T))
    strap_si=bind_rows(strap_si, Ax)
  }
  
  melt_strap_si <- melt(strap_si)
  summary_table_si <- melt_strap_si  %>%  #### summary at estimate level - that is based on each serial interval or generation time estimate
    group_by(variable) %>%
    summarise(n_samples=n(),
              meanx=mean(value),
              sdx=sd(value),
              lci=meanx-1.96*sdx,
              uci=meanx+1.96*sdx)
  
  names_df_gamma_mean_sd[i]
  summary_table_si$study <- names_df_gamma_mean_sd[i]
  collection_si <- rbind(collection_si, summary_table_si)
  ######## transmission relative to symptoms plus bootstrap
  TRS <- sims1$SI-sims$IP
  sims2 <- data.frame(TRS=TRS)
  sims2x <- as.data.frame(sims2)
  sims2x$study <- names_df_gamma_mean_sd[i]
  name_csv <- paste("TRS", names_df_gamma_mean_sd[i], i, ".csv", sep="_")
  name_csv
 # setwd("~/Documents/covid_19/r_files/code/sim_b_output/trs")
  write.csv(sims2x, name_csv)
  sum(TRS<0, na.rm=T)
  length(a$TRS)
  length(a$TRS[!is.na(a$TRS)])
  strap_trs <- data.frame()
  
  for(l in 1:n_bootstrap_runs){
    a <- sample_n(sims2, n_bootstrap_samples, replace=TRUE)
    summary(a)
    Ax = data.frame(meanx=mean(a$TRS, na.rm=T), 
                    sdx=sd(a$TRS, na.rm=T),
                    X2.5=quantile(a$TRS, probs=0.025, na.rm=T),
                    X25=quantile(a$TRS, probs=0.25, na.rm=T),
                    X50=quantile(a$TRS, probs=0.5, na.rm=T),
                    X75=quantile(a$TRS, probs=0.75, na.rm=T),
                    X97.5=quantile(a$TRS, probs=0.975, na.rm=T),
                    n_samples=length(a$TRS[!is.na(a$TRS)]),
                    n_trs_neg=sum(a$TRS<0, na.rm=T))
    Ax$prop_trs_neg <- Ax$n_trs_neg/Ax$n_samples
    strap_trs=bind_rows(strap_trs, Ax)
  }
  ###### save strap trs
  name_csv <- paste("strap_TRS", names_df_gamma_mean_sd[i], i, ".csv", sep="_")
  name_csv
 # setwd("~/Documents/covid_19/r_files/code/sim_b_output/strap_trs")
  write.csv(strap_trs, name_csv)
  
  
  melt_strap_trs <- melt(strap_trs)
  summary_table_trs <- melt_strap_trs  %>%  
    group_by(variable) %>%
    summarise(n_samples=n(),
              meanx=mean(value),
              sdx=sd(value),
              lci=meanx-1.96*sdx,
              uci=meanx+1.96*sdx)
  
  summary_table_trs$study <- names_df_gamma_mean_sd[i]
  
  collection_trs <- rbind(collection_trs, summary_table_trs)
  
}
summary(sims1x)
summary(sims2x)
####### mean_sd_normal
#### now make loop for all 3 normally distributed serial intervals with ci's for mean and standard deviation
summary(as.factor(d$ci_category))

df_norm_mean_sd <- d[d$ci_category=="mean_sd_normal",]
names_df_norm_mean_sd <- unique(df_norm_mean_sd$study_r)
names_df_norm_mean_sd
collection_si1 <- data.frame()
collection_trs1 <- data.frame()

for (i in 1:length(names_df_norm_mean_sd))
{
  df1 <- d[d$study_r==names_df_norm_mean_sd[i],]
  mean_data <- df1$central_measure[df1$dist_par_r=="n_m"]
  mean_data_lci <- df1$lower_ci[df1$dist_par_r=="n_m"]
  mean_data_uci <- df1$upper_ci[df1$dist_par_r=="n_m"]
  sd_data <- df1$central_measure[df1$dist_par_r=="n_sd"]
  sd_data_lci <- df1$lower_ci[df1$dist_par_r=="n_sd"]
  sd_data_uci <- df1$upper_ci[df1$dist_par_r=="n_sd"]
  
  sims1 <- data.frame ()
  for (j in 1:n_primary_runs){
    meanx <- rnorm(1, mean_data, (mean_data_uci-mean_data_lci)/2*1.96)
    ####### try gamma distribution for standard deviation
    #####gamma shape=mean^2/sd^2  gamma rate=mean/sd^2
    sd_sd <- (sd_data_uci-sd_data_lci)/(2*1.96)
    sd_mean <- sd_data
    sd_shape <- sd_mean^2/(sd_sd^2)
    sd_rate <- sd_mean/(sd_sd^2)
    sdx <- rgamma(1, sd_shape, sd_rate)
    B=data.frame(SI=rnorm(n_primary_samples, meanx, sdx))
    summary(B)
    sims1 <- rbind(sims1, B)
  }
  summary(sims1)
  sims1x <- as.data.frame(sims1)
  sims1x$study <- names_df_norm_mean_sd[i]
  name_csv <- paste("SI", names_df_norm_mean_sd[i], i, ".csv", sep="_")
  name_csv
#  setwd("~/Documents/covid_19/r_files/code/sim_b_output/si")
  write.csv(sims1x, name_csv)
  
  
  
  strap_si <- data.frame()
  for(k in 1:n_bootstrap_runs){
    a <- sample_n(sims1, n_bootstrap_samples, replace=TRUE)
    summary(a)
    Ax = data.frame(mean=mean(a$SI, na.rm=T), 
                    sd=sd(a$SI, na.rm=T),
                    X2.5=quantile(a$SI, probs=0.025, na.rm=T),
                    X25=quantile(a$SI, probs=0.25, na.rm=T),
                    X50=quantile(a$SI, probs=0.5, na.rm=T),
                    X75=quantile(a$SI, probs=0.75, na.rm=T),
                    X97.5=quantile(a$SI, probs=0.975, na.rm=T))
    strap_si=bind_rows(strap_si, Ax)
  }
  
  melt_strap_si <- melt(strap_si)
  summary_table_si <- melt_strap_si  %>%  #### summary at estimate level - that is based on each serial interval or generation time estimate
    group_by(variable) %>%
    summarise(n_samples=n(),
              meanx=mean(value),
              sdx=sd(value),
              lci=meanx-1.96*sdx,
              uci=meanx+1.96*sdx)
  
  names_df_norm_mean_sd[i]
  summary_table_si$study <- names_df_norm_mean_sd[i]
  collection_si1 <- rbind(collection_si1, summary_table_si)
  ######## transmission relative to symptoms plus bootstrap
  TRS <- sims1$SI-sims$IP
  sims2 <- data.frame(TRS=TRS)
  sims2x <- as.data.frame(sims2)
  sims2x$study <- names_df_norm_mean_sd[i]
  name_csv <- paste("TRS", names_df_norm_mean_sd[i], i, ".csv", sep="_")
  name_csv
 # setwd("~/Documents/covid_19/r_files/code/sim_b_output/trs")
  write.csv(sims2x, name_csv)
  sum(TRS<0, na.rm=T)
  length(a$TRS)
  length(a$TRS[!is.na(a$TRS)])
  strap_trs <- data.frame()
  
  for(l in 1:n_bootstrap_runs){
    a <- sample_n(sims2, n_bootstrap_samples, replace=TRUE)
    summary(a)
    Ax = data.frame(meanx=mean(a$TRS, na.rm=T), 
                    sdx=sd(a$TRS, na.rm=T),
                    X2.5=quantile(a$TRS, probs=0.025, na.rm=T),
                    X25=quantile(a$TRS, probs=0.25, na.rm=T),
                    X50=quantile(a$TRS, probs=0.5, na.rm=T),
                    X75=quantile(a$TRS, probs=0.75, na.rm=T),
                    X97.5=quantile(a$TRS, probs=0.975, na.rm=T),
                    n_samples=length(a$TRS[!is.na(a$TRS)]),
                    n_trs_neg=sum(a$TRS<0, na.rm=T))
    Ax$prop_trs_neg <- Ax$n_trs_neg/Ax$n_samples
    strap_trs=bind_rows(strap_trs, Ax)
  }
  ###### save strap trs
  name_csv <- paste("strap_TRS", names_df_norm_mean_sd[i], i, ".csv", sep="_")
  name_csv
 # setwd("~/Documents/covid_19/r_files/code/sim_b_output/strap_trs")
  write.csv(strap_trs, name_csv)
  
  
  
  melt_strap_trs <- melt(strap_trs)
  summary_table_trs <- melt_strap_trs  %>%  
    group_by(variable) %>%
    summarise(n_samples=n(),
              meanx=mean(value),
              sdx=sd(value),
              lci=meanx-1.96*sdx,
              uci=meanx+1.96*sdx)
  
  summary_table_trs$study <- names_df_norm_mean_sd[i]
  
  collection_trs1 <- rbind(collection_trs1, summary_table_trs)
  
}

###### lognormal distribution with mean and standard deviation
#sigma <- sqrt(log(v/(mean^2)+1))
#mu <- log(mean) - (sigma)^2/2

####### mean_sd_lognormal - Nishiura 28

summary(as.factor(d$ci_category))

df_lnorm_mean_sd <- d[d$ci_category=="mean_sd_lognormal",]
names_df_lnorm_mean_sd <- unique(df_lnorm_mean_sd$study_r)
names_df_lnorm_mean_sd[i]
collection_si2 <- data.frame()
collection_trs2 <- data.frame()

for (i in 1:length(names_df_lnorm_mean_sd))
{
  df1 <- d[d$study_r==names_df_lnorm_mean_sd[i],]
  mean_data <- df1$central_measure[df1$dist_par_r=="mean"]
  mean_data_lci <- df1$lower_ci[df1$dist_par_r=="mean"]
  mean_data_uci <- df1$upper_ci[df1$dist_par_r=="mean"]
  sd_data <- df1$central_measure[df1$dist_par_r=="sd"]
  sd_data_lci <- df1$lower_ci[df1$dist_par_r=="sd"]
  sd_data_uci <- df1$upper_ci[df1$dist_par_r=="sd"]
  
  sims1 <- data.frame ()
  for (j in 1:n_primary_runs){
    meanx <- rnorm(1, mean_data, (mean_data_uci-mean_data_lci)/2*1.96)
    ####### try gamma distribution for standard deviation
    #####gamma shape=mean^2/sd^2  gamma rate=mean/sd^2
    sd_sd <- (sd_data_uci-sd_data_lci)/(2*1.96)
    sd_mean <- sd_data
    sd_shape <- sd_mean^2/(sd_sd^2)
    sd_rate <- sd_mean/(sd_sd^2)
    sdx <- rgamma(1, sd_shape, sd_rate)
    v <- sdx^2
    sigmax <- sqrt(log(v/(meanx^2)+1))
    mux <- log(meanx) - (sigmax)^2/2
    B=data.frame(SI=rlnorm(n_primary_samples, mux, sigmax))
    sims1 <- rbind(sims1, B)
  }
  summary(sims1)
  sims1x <- as.data.frame(sims1)
  sims1x$study <- names_df_lnorm_mean_sd[i]
  name_csv <- paste("SI", names_df_lnorm_mean_sd[i], i, ".csv", sep="_")
  name_csv
  #setwd("~/Documents/covid_19/r_files/code/sim_b_output/si")
  write.csv(sims1x, name_csv)
  
  strap_si <- data.frame()
  for(k in 1:n_bootstrap_runs){
    a <- sample_n(sims1, n_bootstrap_samples, replace=TRUE)
    summary(a)
    Ax = data.frame(mean=mean(a$SI, na.rm=T), 
                    sd=sd(a$SI, na.rm=T),
                    X2.5=quantile(a$SI, probs=0.025, na.rm=T),
                    X25=quantile(a$SI, probs=0.25, na.rm=T),
                    X50=quantile(a$SI, probs=0.5, na.rm=T),
                    X75=quantile(a$SI, probs=0.75, na.rm=T),
                    X97.5=quantile(a$SI, probs=0.975, na.rm=T))
    strap_si=bind_rows(strap_si, Ax)
  }
  
  melt_strap_si <- melt(strap_si)
  summary_table_si <- melt_strap_si  %>%  #### summary at estimate level - that is based on each serial interval or generation time estimate
    group_by(variable) %>%
    summarise(n_samples=n(),
              meanx=mean(value),
              sdx=sd(value),
              lci=meanx-1.96*sdx,
              uci=meanx+1.96*sdx)
  
  names_df_lnorm_mean_sd[i]
  summary_table_si$study <- names_df_lnorm_mean_sd[i]
  collection_si2 <- rbind(collection_si2, summary_table_si)
  ######## transmission relative to symptoms plus bootstrap
  TRS <- sims1$SI-sims$IP
  sims2 <- data.frame(TRS=TRS)
  sims2x <- as.data.frame(sims2)
  sims2x$study <- names_df_lnorm_mean_sd[i]
  name_csv <- paste("TRS", names_df_lnorm_mean_sd[i], i, ".csv", sep="_")
  name_csv
  #setwd("~/Documents/covid_19/r_files/code/sim_b_output/trs")
  write.csv(sims2x, name_csv)
  sum(TRS<0, na.rm=T)
  length(a$TRS)
  length(a$TRS[!is.na(a$TRS)])
  strap_trs <- data.frame()
  
  for(l in 1:n_bootstrap_runs){
    a <- sample_n(sims2, n_bootstrap_samples, replace=TRUE)
    summary(a)
    Ax = data.frame(meanx=mean(a$TRS, na.rm=T), 
                    sdx=sd(a$TRS, na.rm=T),
                    X2.5=quantile(a$TRS, probs=0.025, na.rm=T),
                    X25=quantile(a$TRS, probs=0.25, na.rm=T),
                    X50=quantile(a$TRS, probs=0.5, na.rm=T),
                    X75=quantile(a$TRS, probs=0.75, na.rm=T),
                    X97.5=quantile(a$TRS, probs=0.975, na.rm=T),
                    n_samples=length(a$TRS[!is.na(a$TRS)]),
                    n_trs_neg=sum(a$TRS<0, na.rm=T))
    Ax$prop_trs_neg <- Ax$n_trs_neg/Ax$n_samples
    strap_trs=bind_rows(strap_trs, Ax)
  }
  name_csv <- paste("strap_TRS", names_df_lnorm_mean_sd[i], i, ".csv", sep="_")
  name_csv
#  setwd("~/Documents/covid_19/r_files/code/sim_b_output/strap_trs")
  write.csv(strap_trs, name_csv)
  melt_strap_trs <- melt(strap_trs)
  summary_table_trs <- melt_strap_trs  %>%  
    group_by(variable) %>%
    summarise(n_samples=n(),
              meanx=mean(value),
              sdx=sd(value),
              lci=meanx-1.96*sdx,
              uci=meanx+1.96*sdx)
  
  summary_table_trs$study <- names_df_lnorm_mean_sd[i]
  
  collection_trs2 <- rbind(collection_trs2, summary_table_trs)
  
}

##### Weibull mean and sd
#mean <- 4.8
#sd <- 2.3
#wp <- weibullpar(mean, sd)
#wp$shape
#wp$scale
#wp$loc

###### Weibull distribution with mean and standard deviation - Nishiura 18

summary(as.factor(d$ci_category))

df_weibull_mean_sd <- d[d$ci_category=="mean_sd_weibull",]
names_df_weibull_mean_sd <- unique(df_weibull_mean_sd$study_r)
collection_si3 <- data.frame()
collection_trs3 <- data.frame()

for (i in 1:length(names_df_weibull_mean_sd))
{
  df1 <- d[d$study_r==names_df_weibull_mean_sd[i],]
  mean_data <- df1$central_measure[df1$dist_par_r=="mean"]
  mean_data_lci <- df1$lower_ci[df1$dist_par_r=="mean"]
  mean_data_uci <- df1$upper_ci[df1$dist_par_r=="mean"]
  sd_data <- df1$central_measure[df1$dist_par_r=="sd"]
  sd_data_lci <- df1$lower_ci[df1$dist_par_r=="sd"]
  sd_data_uci <- df1$upper_ci[df1$dist_par_r=="sd"]
  
  sims1 <- data.frame ()
  for (j in 1:n_primary_runs){
    meanx <- rnorm(1, mean_data, (mean_data_uci-mean_data_lci)/2*1.96)
    meanx
    mean_data
    ####### try gamma distribution for standard deviation
    #####gamma shape=mean^2/sd^2  gamma rate=mean/sd^2
    sd_sd <- (sd_data_uci-sd_data_lci)/(2*1.96)
    sd_mean <- sd_data
    sd_shape <- sd_mean^2/(sd_sd^2)
    sd_rate <- sd_mean/(sd_sd^2)
    sdx <- rgamma(1, sd_shape, sd_rate)
    sdx
    sd_data
    wp <- weibullpar(meanx, sdx)
    shapex <- wp$shape
    shapex
    scalex <- wp$scale
    scalex
    B=data.frame(SI=rweibull(n_primary_samples, shapex, scalex))
    summary(B)
    sims1 <- rbind(sims1, B)
  }
  summary(sims1)
  sims1x <- as.data.frame(sims1)
  sims1x$study <- names_df_weibull_mean_sd[i]
  name_csv <- paste("SI", names_df_weibull_mean_sd[i], i, ".csv", sep="_")
  name_csv
 # setwd("~/Documents/covid_19/r_files/code/sim_b_output/si")
  write.csv(sims1x, name_csv)
  strap_si <- data.frame()
  for(k in 1:n_bootstrap_runs){
    a <- sample_n(sims1, n_bootstrap_samples, replace=TRUE)
    summary(a)
    Ax = data.frame(mean=mean(a$SI, na.rm=T), 
                    sd=sd(a$SI, na.rm=T),
                    X2.5=quantile(a$SI, probs=0.025, na.rm=T),
                    X25=quantile(a$SI, probs=0.25, na.rm=T),
                    X50=quantile(a$SI, probs=0.5, na.rm=T),
                    X75=quantile(a$SI, probs=0.75, na.rm=T),
                    X97.5=quantile(a$SI, probs=0.975, na.rm=T))
    strap_si=bind_rows(strap_si, Ax)
  }
  
  melt_strap_si <- melt(strap_si)
  summary_table_si <- melt_strap_si  %>%  #### summary at estimate level - that is based on each serial interval or generation time estimate
    group_by(variable) %>%
    summarise(n_samples=n(),
              meanx=mean(value),
              sdx=sd(value),
              lci=meanx-1.96*sdx,
              uci=meanx+1.96*sdx)
  
  names_df_weibull_mean_sd[i]
  summary_table_si$study <- names_df_weibull_mean_sd[i]
  collection_si3 <- rbind(collection_si3, summary_table_si)
  ######## transmission relative to symptoms plus bootstrap
  TRS <- sims1$SI-sims$IP
  sims2 <- data.frame(TRS=TRS)
  sims2x <- as.data.frame(sims2)
  sims2x$study <- names_df_weibull_mean_sd[i]
  name_csv <- paste("TRS", names_df_weibull_mean_sd[i], i, ".csv", sep="_")
  name_csv
  #setwd("~/Documents/covid_19/r_files/code/sim_b_output/trs")
  write.csv(sims2x, name_csv)
  sum(TRS<0, na.rm=T)
  length(a$TRS)
  length(a$TRS[!is.na(a$TRS)])
  strap_trs <- data.frame()
  
  for(l in 1:n_bootstrap_runs){
    a <- sample_n(sims2, n_bootstrap_samples, replace=TRUE)
    summary(a)
    Ax = data.frame(meanx=mean(a$TRS, na.rm=T), 
                    sdx=sd(a$TRS, na.rm=T),
                    X2.5=quantile(a$TRS, probs=0.025, na.rm=T),
                    X25=quantile(a$TRS, probs=0.25, na.rm=T),
                    X50=quantile(a$TRS, probs=0.5, na.rm=T),
                    X75=quantile(a$TRS, probs=0.75, na.rm=T),
                    X97.5=quantile(a$TRS, probs=0.975, na.rm=T),
                    n_samples=length(a$TRS[!is.na(a$TRS)]),
                    n_trs_neg=sum(a$TRS<0, na.rm=T))
    Ax$prop_trs_neg <- Ax$n_trs_neg/Ax$n_samples
    strap_trs=bind_rows(strap_trs, Ax)
  }
  name_csv <- paste("strap_TRS", names_df_weibull_mean_sd[i], i, ".csv", sep="_")
  name_csv
  #setwd("~/Documents/covid_19/r_files/code/sim_b_output/strap_trs")
  write.csv(strap_trs, name_csv)
  
  melt_strap_trs <- melt(strap_trs)
  summary_table_trs <- melt_strap_trs  %>%  
    group_by(variable) %>%
    summarise(n_samples=n(),
              meanx=mean(value),
              sdx=sd(value),
              lci=meanx-1.96*sdx,
              uci=meanx+1.96*sdx)
  
  summary_table_trs$study <- names_df_weibull_mean_sd[i]
  
  collection_trs3 <- rbind(collection_trs3, summary_table_trs)
}



si_df <- rbind(collection_si, collection_si1, collection_si2, collection_si3)
trs_df <- rbind(collection_trs, collection_trs1, collection_trs2, collection_trs3)
#setwd("~/Documents/covid_19/r_files/code")
write.csv(si_df, "si_df.csv")
write.csv(trs_df, "trs_df.csv")








