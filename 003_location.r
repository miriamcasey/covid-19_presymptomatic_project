##### pooling, summarising and plotting by location where contact tracing data that serial interval or generation time data are based on were collected.
##### source("001_final_code.r") #### - this is to be run after "001_final_code.r"

##### pooling, summarising and plotting by location where contact tracing data that serial interval or generation time data are based on were collected.
##### source("001_final_code.r") #### - this is to be run after "001_final_code.r"

summary(dfx)

loc <- read.csv("location_and_reference.csv")

colnames(loc)
summary(loc$Location)
loc$loc1 <- as.character(loc$Location)
loc$loc1[loc$Location=="Hong Kong and Shenzen"] <- "HKSZ" ## shorten name
loc$Location <- loc$loc1
colnames(loc)
colnames(dfx)
dfx$study1x <- paste(dfx$study, dfx$transmission_parameter) #### unique id for merging
dfx_loc <- merge(dfx, loc, by.x="study1x", by.y="study1x")
summary(dfx_loc)

dfx_loc1 <- dfx_loc[dfx_loc$transmission_parameter=="SI",]  #### only analyse for SI


summary_table_loc <- dfx_loc1  %>%
  group_by(Location) %>%
  summarise(n_samples=n(),
            n_si_less_than_0=sum(si_or_gt<0), #### negative serial intervals for relevant studies
            prop_neg_si=n_si_less_than_0/n_samples, ##### proportion negative serial intervals for relevant studies
            N_trs_neg=sum(trs<0), ### negative transmission time relative to symptom onset
            prop_trs_neg=N_trs_neg/n_samples, ### proportion negative transmission time relative to symptom onset### distribution of serial interval or generation time
            mean_si_or_gt=mean(si_or_gt),
            sd_si_or_gt=sd(si_or_gt),
            si_or_gt_2.5_percentile=quantile(si_or_gt, 0.025),
            si_or_gt_25_percentile=quantile(si_or_gt, 0.25),
            si_or_gt_median=quantile(si_or_gt, 0.5),
            si_or_gt_75_percentile=quantile(si_or_gt, 0.75),
            si_or_gt_97.5_percentile=quantile(si_or_gt, 0.975),  ### best fit distribution for transmission time relative to symptom onset
            mean_trs=mean(trs),
            sd_trs=sd(trs),
            trs_2.5_percentile=quantile(trs, 0.025),
            trs_25_percentile=quantile(trs, 0.25),
            median_trs=quantile(trs, 0.5),
            trs_75_percentile=quantile(trs, 0.75),
            trs_97.5_percentile=quantile(trs, 0.975))

#write.csv(summary_table_loc, "summary_table_loc.csv")
colnames(dfx_loc1)
summary(dfx_loc1$Location)
summary_table_loc$order <- as.numeric(as.factor(summary_table_loc$median_trs))
summary_table_loc$Location <- factor(summary_table_loc$Location, 
                                     levels = 
                                       summary_table_loc$Location[order(summary_table_loc$order)])
summary_table_loc$Location1 <- summary_table_loc$Location

#### order by median transmission relative to symptoms
summary_table_loc$Location
colnames(dfx_loc1)
dfx_loc2 <- merge(dfx_loc1, summary_table_loc, by.x="Location", by.y="Location", all=T)

summary(dfx_loc2$Location1)

p <- ggplot(dfx_loc2, aes(x=Location1, y=trs)) + 
  geom_boxplot(outlier.shape=NA) 
p
p <- p + coord_flip() + ylim(-20, 20)
p <- p + xlab("Location")
#p <- p + scale_x_discrete(labels=nx)
p <- p+ ylab("Transmission time in days relative to symptom onset")
p <- p + geom_hline(yintercept=0, col="red")
p <- p + geom_point( data=summary_table_loc,
                     aes(x=Location1, y=mean_trs), col="purple", shape=17, size=8)
p

#pdf("Figure_3.pdf", width=20, height=22)
#p   + theme_for_font 
#+ ggtitle("Figure 4") +
# theme(plot.title = element_text(family = "Helvetica", face = "bold", size = (90)))
#dev.off()










