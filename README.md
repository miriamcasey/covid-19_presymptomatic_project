# covid-19_presymptomatic_project
data_and_code _to_estimate_presymptomatic_transmission_of_covid-19
###### MS ID#: MEDRXIV/2020/094870
#####  MS TITLE: Estimating pre-symptomatic transmission of COVID-19: a secondary analysis using published data
#### Background: Understanding the extent of SARS-CoV-2 transmission that can occur before symptom onset is required for targeting control measures against the COVID-19 pandemic. 
#Objective: Estimation of the proportion of pre-symptomatic transmission of COVID-19 that can occur and timing of transmission 
### relative to symptom onset.
#### Design: Secondary analysis of published data ####
#Methods: Simulations were generated of incubation period and of serial interval or generation time. From these, transmission times #relative to symptom onset were calculated and the proportion of pre-symptomatic transmission was estimated.
#### We hope the data and code uploaded here will save others effort identifying and parameterising incubation period, generation time and # serial interval distributions

#### "001_final_code.r" sources "load_data_final.r" and then summarises and plots data.

#### "002_load_data_final.r" simulates samples from a distribution for incubation period based on our separately published metanalysis: 
##### https://www.medrxiv.org/content/10.1101/2020.04.24.20073957v1
##### it then simulates samples from distributions for serial interval (n=23) or generation time (n=5) from published data. These were
#### sourced through our rapid systematic review - just uploaded on Medrxiv (Griffin et al.)
##### Incubation period samples are subtracted from each samples representing each serial interval or generation time to generate transmission times relative to symptoms "trs" and collected into a dataframe - this is summarised and plotted in "001_final_code.r
##### "location.r" merges location data for the different serial interval or generation time studies with teh data frame generated by "002_load_data_final.r". The simulation output is then summarised and plotted at location level.
