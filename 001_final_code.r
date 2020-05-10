rm(list=ls())
library(tidyverse)
library(janitor)
library(epitrix)
library(fitdistrplus)
library(rriskDistributions)


#### load data
d <- read_csv("final_data_5_may.csv")
#### remove liao pending clarification on distribution 
d <- d[d$study_r!="liao",]
##### remove cheng as cheng does not allow for pre-symptomatic transmission
d <- d[d$study_r!="cheng",]
d$central_measure

### sample numbers
n_sample <- 100000  #### reduce to 10,000 for code to run more quickly

#### load data and generate transmission relative to symptom onset distributions.
source("002_load_data_final.r")

#################
summary(dfx)
colnames(dfx)
dfx$si_or_gt <- as.numeric(dfx$si_or_gt)
hist(dfx$si_or_gt)
summary(dfx$si_or_gt)
is.numeric(dfx$si_or_gt)
mean(dfx$si_or_gt)
dfx$study1 <- paste(dfx$study, dfx$transmission_parameter)
colnames(dfx)

summary_table <- dfx  %>%  #### summary at estimate level - that is based on each serial interval or generation time estimate
  group_by(study1) %>%
  summarise(n_samples=n(),
            n_si_less_than_0=sum(si_or_gt<0), #### negative serial intervals for relevant studies
            prop_neg_si=n_si_less_than_0/n_samples, ##### proportion negative serial intervals for relevant studies
            N_trs_neg=sum(trs<0), ### negative transmission time relative to symptom onset
            prop_trs_neg=N_trs_neg/n_samples, ### proportion negative transmission time relative to symptom onset
            trans_par_si_or_gt=unique(transmission_parameter),
            n_from_si_or_gt_study=unique(n_gt_or_si),
            dist_si_or_gt=unique(si_or_gt_dist), ### distribution of serial interval or generation time
            mean_si_or_gt=mean(si_or_gt),
            sd_si_or_gt=sd(si_or_gt),
            si_or_gt_2.5_percentile=quantile(si_or_gt, 0.025),
            si_or_gt_25_percentile=quantile(si_or_gt, 0.25),
            si_or_gt_median=quantile(si_or_gt, 0.5),
            si_or_gt_75_percentile=quantile(si_or_gt, 0.75),
            si_or_gt_97.5_percentile=quantile(si_or_gt, 0.975),
            mean_trs=mean(trs),
            sd_trs=sd(trs),
            trs_2.5_percentile=quantile(trs, 0.025),
            trs_25_percentile=quantile(trs, 0.25),
            median_trs=quantile(trs, 0.5),
            trs_75_percentile=quantile(trs, 0.75),
            trs_97.5_percentile=quantile(trs, 0.975))

colnames(dfx)
dfx$transmission_parameter
dfx_si <- dfx[dfx$transmission_parameter=="SI",]  ### simple pooling of serial interval based estimates
warnings()
summary_table1 <- dfx_si  %>%
  summarise(n_samples=n(),
            N_trs_neg=sum(trs<0),
            prop_trs_neg=N_trs_neg/n_samples,
            mean_trs=mean(trs),
            median_trs=median(trs),
            sd_trs=sd(trs),
            quant_2.5_trs=quantile(trs, 0.025),
            quant_25_trs=quantile(trs, 0.25),
            quant_75_trs=quantile(trs, 0.75),
            quant_97.5_trs=quantile(trs, 0.975))
summary_table1

write.csv(summary_table1,"summary_table1_si.csv")

dfx_gt <- dfx[dfx$transmission_parameter=="GT",] ##### simple pooling of generation time based estimates
warnings()
summary_table1 <- dfx_gt  %>%
  summarise(n_samples=n(),
            N_trs_neg=sum(trs<0),
            prop_trs_neg=N_trs_neg/n_samples,
            mean_trs=mean(trs),
            median_trs=median(trs),
            sd_trs=sd(trs),
            quant_2.5_trs=quantile(trs, 0.025),
            quant_25_trs=quantile(trs, 0.25),
            quant_75_trs=quantile(trs, 0.75),
            quant_97.5_trs=quantile(trs, 0.975))
summary_table1

write.csv(summary_table1,"summary_table1_gt.csv")

################## boxplot
### order by median transmission relative to symptoms

summary_table$studyx[summary_table$study1=="bi SI"] <- "Bi et al."
summary_table$studyx[summary_table$study1=="cereda SI"] <- "Cereda et al."
summary_table$studyx[summary_table$study1=="cheng SI"] <- "Cheng et al."
summary_table$studyx[summary_table$study1=="du SI"] <- "Du et al. (a)"
summary_table$studyx[summary_table$study1=="du_a SI"] <- "Du et al. (b)"
summary_table$studyx[summary_table$study1=="ferretti GT"] <- "Ferretti et al. gt"
summary_table$studyx[summary_table$study1=="ferretti SI"] <- "Ferretti et al."
summary_table$studyx[summary_table$study1=="ganyani_singapore_base GT"] <- "Ganyani et al. (S) gt"
summary_table$studyx[summary_table$study1=="ganyani_singapore_base SI"] <- "Ganyani et al. (S)"
summary_table$studyx[summary_table$study1=="ganyani_singapore_neg GT"] <- "Ganyani et al. (Sn) gt"
summary_table$studyx[summary_table$study1=="ganyani_singapore_neg SI"] <- "Ganyani et al. (Sn)"
summary_table$studyx[summary_table$study1=="ganyani_tianjin_base GT"] <- "Ganyani et al. (T) gt"
summary_table$studyx[summary_table$study1=="ganyani_tianjin_base SI"] <- "Ganyani et al. (T)"
summary_table$studyx[summary_table$study1=="ganyani_tianjin_neg GT"] <- "Ganyani et al. (Tn) gt"
summary_table$studyx[summary_table$study1=="ganyani_tianjin_neg SI"] <- "Ganyani et al. (Tn)"
summary_table$studyx[summary_table$study1=="he SI"] <- "He et al."
summary_table$studyx[summary_table$study1=="lavezzo SI"] <- "Lavezzo et al."
summary_table$studyx[summary_table$study1=="li SI"] <- "Li et al."
summary_table$studyx[summary_table$study1=="ma SI"] <- "Ma et al."
summary_table$studyx[summary_table$study1=="nishiura_18 SI"] <- "Nishiura et al. (n=18)"
summary_table$studyx[summary_table$study1=="nishiura_28 SI"] <- "Nishiura et al. (n=28)"
summary_table$studyx[summary_table$study1=="tindale_singapore SI"] <- "Tindale et al. (S)"
summary_table$studyx[summary_table$study1=="tindale_tianjin SI"] <- "Tindale et al. (T)"
summary_table$studyx[summary_table$study1=="wang_teunis SI"] <- "Wang and Teunis"
summary_table$studyx[summary_table$study1=="wu SI"] <- "Wu et al."
summary_table$studyx[summary_table$study1=="zhang SI"] <- "Zhang et al.(b)"
summary_table$studyx[summary_table$study1=="zhao_a_12 SI"] <- "Zhao et al.(a) (n=12)"
summary_table$studyx[summary_table$study1=="zhao_a_21 SI"] <- "Zhao et al.(a) (n=21)"
summary_table$studyx[summary_table$study1=="zhao_b SI"] <- "Zhao et al.(b)"
write.csv(summary_table, "trial_summary.csv")

###### order by median transmission relative to symptom onset for plot
colnames(dfx)
unique(dfx$study)
summary_table$order <- as.numeric(as.factor(summary_table$median_trs))
summary_table$order
summary_table$order1 <- as.factor(summary_table$order)
summary_table$median_trs_mca
summary_table$Study <- factor(summary_table$studyx, 
                              levels = 
                                summary_table$studyx[order(summary_table$order)]) #### order by median transmission relative to symptoms
summary_table$Study
colnames(dfx)
dfx1 <- merge(dfx, summary_table, by.x="study1", by.y="study1", all=T)
colnames(dfx1)
colnames(summary_table)
summary_table
summary_table$transmission_parameter <- summary_table$trans_par_si_or_gt


############################Plot transmission time relative to symptom onset
colnames(dfx1)

p <- ggplot(dfx1, aes(x=Study, y=trs)) + 
  geom_boxplot(outlier.shape=NA) 
p
p <- p + coord_flip() + ylim(-20, 20)
p
p <- p+ ylab("Transmission time in days relative to symptom onset")
p <- p + geom_hline(yintercept=0, col="red")
p <- p + geom_point( data=summary_table,
                     aes(x=Study, y=mean_trs), col="purple", shape=17)
p
p <- p+ facet_grid(vars(rows=transmission_parameter),scales="free", space="free", drop=TRUE)
p

