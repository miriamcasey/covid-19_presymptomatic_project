### load data and estimate transmission time relative to symptoms

#### load incubation period from McAloon et al metanalysis:
## https://www.medrxiv.org/content/10.1101/2020.04.24.20073957v1

d[d$study_r=="mca",]
### mean log
mca_ip_ml <- d$central_measure[d$dist_par_r=="ln_ml" & d$study_r=="mca"]
mca_ip_ml
## sd log
mca_ip_sdl <- d$central_measure[d$dist_par_r=="ln_sdl" & d$study_r=="mca"]
mca_ip_sdl
### sample distribution McAloon IP
mca_ip <- rlnorm(n_sample, mca_ip_ml, mca_ip_sdl)
summary(mca_ip)

############## loop for all gamma distributed serial intervals or generation times
colnames(d)
summary(as.factor(d$distribution))
gamma_df <- d[d$distribution=="Gamma",]
gamma_names <- unique(gamma_df$study_r)
gamma_names

gamma_df1 <- data.frame()### make a dataframe to collect data from loop

for(i in 1:length(gamma_names))
{
  df1 <- d[d$study_r==gamma_names[i],]  ### get data for a particular study
  g_sh <- df1$central_measure[df1$dist_par_r=="g_sh"] ### get gamma shape parameter for serial interval or generation time
  g_r <- df1$central_measure[df1$dist_par_r=="g_r"] ### get gamma rate parameter for serial interval or generation time
  si_or_gt1 <- rgamma(n_sample,g_sh, g_r) ##### generate distribution for serial interval or generation time 
  trs1 <- si_or_gt1-mca_ip #### subtract incubation period from serial interval or generation time to estimate transmission time relative to symptom onset
  df2 <- data.frame(study=rep(gamma_names[i], n_sample)) ### collect data in a dataframe
  df2$transmission_parameter <- rep(unique(df1$transmission_parameter), n_sample)
  df2$n_gt_or_si <- rep(unique(df1$n), n_sample)
  df2$si_or_gt <- si_or_gt1
  df2$trs <- trs1 
  df2$si_or_gt_dist <- "gamma"
  gamma_df1<- rbind(gamma_df1, df2)
}

############## loop for all normally distributed serial intervals or generation times
colnames(d)
summary(as.factor(d$distribution))
normal_df <- d[d$distribution=="Normal",]
normal_names <- unique(normal_df$study_r)
normal_names

normal_df1 <- data.frame()### make a dataframe to collect data from loop

for(i in 1:length(normal_names))
{
  df1 <- d[d$study_r==normal_names[i],]  ### get data for a particular study
  n_m <- df1$central_measure[df1$dist_par_r=="n_m"] ### get normal mean p for serial interval or generation time
  n_sd <- df1$central_measure[df1$dist_par_r=="n_sd"] ### get normal rstandard deviation for serial interval or generation time
  si_or_gt1 <- rnorm(n_sample,n_m, n_sd) ##### generate distribution for serial interval or generation time 
  trs1 <- si_or_gt1-mca_ip #### subtract incubation period from serial interval or generation time to estimate transmission time relative to symptom onset
  df2 <- data.frame(study=rep(normal_names[i], n_sample)) ### collect data in a dataframe
  df2$transmission_parameter <- rep(unique(df1$transmission_parameter), n_sample)
  df2$n_gt_or_si <- rep(unique(df1$n), n_sample)
  df2$si_or_gt <- si_or_gt1
  df2$trs <- trs1 
  df2$si_or_gt_dist <- "normal"
  normal_df1<- rbind(normal_df1, df2)
}

################
############## loop for all weibull distributed serial intervals or generation times Ferretti and Nishiura 18
colnames(d)
summary(as.factor(d$distribution))
weibull_df <- d[d$distribution=="Weibull",]
weibull_names <- unique(weibull_df$study_r)
weibull_names

weibull_df1 <- data.frame()### make a dataframe to collect data from loop

colnames(df2)
for(i in 1:length(weibull_names))
{
  df1 <- d[d$study_r==weibull_names[i],] 
  df1### get data for a particular study
  w_sh <- df1$central_measure[df1$dist_par_r=="w_sh"] ### get weibull shape parameter for serial interval or generation time
  w_sc <- df1$central_measure[df1$dist_par_r=="w_sc"] ### get weibull scale parameter for serial interval or generation time
  si_or_gt1 <- rweibull(n_sample,w_sh, w_sc) ##### generate distribution for serial interval or generation time 
  trs1 <- si_or_gt1-mca_ip #### subtract incubation period from serial interval or generation time to estimate transmission time relative to symptom onset
  df2 <- data.frame(study=rep(weibull_names[i], n_sample)) ### collect data in a dataframe
  df2$n_gt_or_si <- rep(unique(df1$n), n_sample)
  df2$transmission_parameter <- rep(unique(df1$transmission_parameter), n_sample)
  df2$si_or_gt <- si_or_gt1
  df2$trs <- trs1 
  df2$si_or_gt_dist <- "weibull"
  weibull_df1<- rbind(weibull_df1, df2)
}

############## lognormal distributed serial intervals or generation times - only Nishiura 28
colnames(d)
summary(as.factor(d$distribution))
lognorm_df <- d[d$distribution=="Lognormal",]
lognorm_names <- unique(lognorm_df$study_r[lognorm_df$study_r!="mca"])
lognorm_names  #### only 1 - Nishiura 28 so don't need loop

df1 <- d[d$study_r==lognorm_names,]  ### get data for a particular study - this time only Nishiura 28
ln_ml <- df1$central_measure[df1$dist_par_r=="ln_ml"] ### get lognormal mean log parameter for serial interval or generation time
ln_sdl <- df1$central_measure[df1$dist_par_r=="ln_sdl"] ### get lognormal sd log  parameter for serial interval or generation time
si_or_gt1 <- rlnorm(n_sample,ln_ml, ln_sdl) ##### generate distribution for serial interval or generation time 
trs1 <- si_or_gt1-mca_ip #### subtract incubation period from serial interval or generation time to estimate transmission time relative to symptom onset
df2 <- data.frame(study=rep("nishiura_28", n_sample)) ### collect data in a dataframe
df2$transmission_parameter <- rep("SI", n_sample)
df2$n_gt_or_si <- rep(unique(df1$n), n_sample)
df2$si_or_gt <- si_or_gt1
df2$trs <- trs1 
df2$si_or_gt_dist <- "lognormal"
df_lognorm <- df2

######### now for shifted gamma (7.6% negative SI's) of He et al
##### he SI - fitted gamma with 7.6% negative serial intervals
#### got shape and rate and how they shifted distributions from github shared code : https://github.com/ehylau/COVID-19
df1 <- d[d$study_r=="he",]
#### gamma shape
he_g_sh <- d$central_measure[d$dist_par_r=="g_sh" & d$study_r=="he"]
he_g_sh
#### gamma rate 
he_g_r <- d$central_measure[d$dist_par_r=="g_r" & d$study_r=="he"]
he_g_r
### sample distribution he SI
##### find value to subtract
##### what point is 7.6% into the distribution
qgammavalue <- qgamma(0.076, he_g_sh, he_g_r)
qgammavalue
he_si <- rgamma(n_sample, he_g_sh, he_g_r) - qgammavalue
hist(he_si)
plot(density(he_si), xlim=c(-5, 22))
summary(he_si)##### 
si_or_gt1 <- he_si ##### generate distribution for serial interval or generation time 
trs1 <- si_or_gt1-mca_ip #### subtract incubation period from serial interval or generation time to estimate transmission time relative to symptom onset
df2 <- data.frame(study=rep("he", n_sample)) ### collect data in a dataframe
df2$transmission_parameter <- rep("SI", n_sample)
df2$n_gt_or_si <- rep(unique(df1$n), n_sample)
df2$si_or_gt <- si_or_gt1
df2$trs <- trs1 
df2$si_or_gt_dist <- "shifted gamma (7.6% negative SI)"
df_he <- df2

###### load three generation times (Ganyani Singapore, Ganyani Tianjin and Ferretti) and estimate SI
### Ganyani described how to do this in paper and code:
### https://www.eurosurveillance.org/content/10.2807/1560-7917.ES.2020.25.17.2000257
## ganyani singapore_base GT
##### singapore mean gt = 5.2, sd gt = 1.72
##### mean of ip = 5.2, sd = 2.8  from Zhang https://www.thelancet.com/journals/laninf/article/PIIS1473-3099(20)30230-9/fulltext

ganyani_names <- c("ganyani_singapore_base", "ganyani_singapore_neg" , "ganyani_tianjin_base",
                   "ganyani_tianjin_neg" )

ganyani_df1 <- data.frame()### make a dataframe to collect data from loop

for(i in 1:length(ganyani_names))
{
  df1 <- d[d$study_r==ganyani_names[i],]  ### get data for a particular study
  g_sh <- df1$central_measure[df1$dist_par_r=="g_sh"] ### get ganyani shape parameter for generation time
  g_r <- df1$central_measure[df1$dist_par_r=="g_r"] ### get ganyani rate parameter for generation time
  gt1 <- rgamma(n_sample,g_sh, g_r)
  ###### Ganyani incubation period from Zhang https://www.thelancet.com/journals/laninf/article/PIIS1473-3099(20)30230-9/fulltext
  #####gamma shape=mean^2/sd^2  gamma rate=mean/sd^2
  inc1=rgamma(n_sample,shape=5.2^2/(2.8^2),rate = 5.2/(2.8^2))
  inc2=rgamma(n_sample,shape=5.2^2/(2.8^2),rate = 5.2/(2.8^2))
  si1 <- gt1 + inc1 - inc2
  ##### now estimate transmission time relative to symptoms from estimated Ganyani SI using same aproach as for other papers
  trs1 <- si1-mca_ip #### subtract incubation period from serial interval to estimate transmission time relative to symptom onset
  df2 <- data.frame(study=rep(ganyani_names[i], n_sample)) ### collect data in a dataframe
  df2$transmission_parameter <- rep("SI", n_sample)
  n_sample
  nrow(df2)
  unique(df1$n)
  df2$n_gt_or_si <- rep(unique(df1$n), n_sample)
  df2$si_or_gt <- si1
  df2$trs <- trs1 
  df2$si_or_gt_dist <- "not directly reported - estimated from GT and IP"
  ganyani_df1<- rbind(ganyani_df1, df2)
  rm(df2)
}

####### estimate serial interval from Ferretti's generation time using same approach as Ganyani
#### incubation period that Ferretti used is from the preprint of Lauer:
##### https://www.medrxiv.org/content/10.1101/2020.02.02.20020016v1
### final lauer: https://annals.org/aim/fullarticle/2762808/incubation-period-coronavirus-disease-2019-covid-19-from-publicly-reported
# (Lauer IP is a little different in preprint and final paper - taking parameters quoted from Ferretti paper)

df1 <- d[d$study_r=="ferretti",]  
df1### get data for a particular study
w_sh <- df1$central_measure[df1$dist_par_r=="w_sh"] ### get ferretti shape parameter for generation time
w_sc <- df1$central_measure[df1$dist_par_r=="w_sc"] ### get ferretti scale parameter for  generation time
gt1 <- rweibull(n_sample,w_sh, w_sc)
###### ferretti: IP is lognormal from Lauer, GT is weibull
#### The distribution of incubation period is lognormal with mean 5.5 days, median 5.2 days and standard deviation 2.1 days
inc1=rlnorm(n_sample,1.644, 0.363)
inc2=rlnorm(n_sample,1.644, 0.363)
si1 <- gt1 + inc1 - inc2
summary(si1)
summary(gt1)
sd(si1)
sd(gt1)
##### now estimate transmission time relative to symptoms from estimated Ganyani SI using same aproach as for other papers
trs1 <- si1-mca_ip #### subtract incubation period from serial interval to estimate transmission time relative to symptom onset
df2 <- data.frame(study=rep("ferretti", n_sample)) ### collect data in a dataframe
df2$transmission_parameter <- rep("SI", n_sample)
df2$n_gt_or_si <- rep(unique(df1$n), n_sample)
df2$si_or_gt <- si1
df2$trs <- trs1 
df2$si_or_gt_dist <- "not directly reported - estimated from GT and IP"
ferretti_df<- df2

########################################
##### merge dataframes
colnames(gamma_df1)
colnames(normal_df1)
colnames(weibull_df1)
colnames(df_lognorm)
colnames(ganyani_df1)
colnames(ferretti_df)


dfx <- rbind(gamma_df1, normal_df1, weibull_df1, df_lognorm, ganyani_df1, df_he, ferretti_df)