
############################################################################# 
## Programmer: Megan Wang
## Purposes: Develop R code to analyze the data collected from 374 clinical trials published on ClinicalTrials.gov
## Data:
## In order to understand how the phase I clinical trials are conducted and 
## how the design parameters are being used, I performed a literature review of 
## the published phase I trials on ClinicalTrials.gov and compiled a list of 374 Phase I clinical trials 
## This program is used to analyze the collected data for Figure 1 in the manuscript.
############################################################################# 

## Import data

mydata <- read.csv("/Users/Megan Wang/Documents/Research/dataset3.csv")

## Explore the data

head(mydata)
summary(mydata)
dim(mydata)
dimnames(mydata)
attach(mydata)

## [1] "NCT.Number"             "Study.Title"            "Study.URL"              "Study.Status"          
## [5] "Conditions"             "Interventions"          "Sponsor"                "Collaborators"         
## [9] "Study.Type"             "X"                      "single_combination"     "Dose_finding"          
## [13] "Design"                 "other_design"           "more_than_one_drug"     "n_finding"             
## [17] "dose.level"             "comment_finding"        "Dose.expansion"         "n_expansion"           
## [21] "toxicity.stopping.rule" "comment_expansion"      "phase2"                 "n_phase2"              
## [25] "comment_phase2" 

# Summarize data - entire data cohort

table(single_combination)
prop.table(table(single_combination))

table(Dose_finding)
round (prop.table(table(Dose_finding)), 3)

table(Design)
prop.table(table(Design))

table(other_design)

table(more_than_one_drug)
prop.table(table(more_than_one_drug))

table(Dose_expansion)
prop.table(table(Dose_expansion))

table(toxicity_stopping_rule)
prop.table(table(toxicity_stopping_rule))

table(phase2)
prop.table(table(phase2))

summary(n_finding)
mean(n_finding,na.rm=T)
sd(n_finding,na.rm=T)

summary(dose_level)
mean(dose.level,na.rm=T)
sd(dose.level,na.rm=T)

summary(n_expansion)
mean(n_expansion,na.rm=T)
sd(n_expansion,na.rm=T)

summary(n_phase2)
mean(n_phase2,na.rm=T)
sd(n_phase2,na.rm=T)


######################

# define the subgroup for the dose finding trials with at least 3 dose levels

mydata_sub<-mydata[which(mydata$dose_level>2 & Dose_finding=='Yes'),]

# Data summary for the subgroup

table(mydata_sub$single_combination)
prop.table(table(mydata_sub$single_combination))

table(mydata_sub$Dose_finding)
prop.table(table(mydata_sub$Dose_finding))

table(mydata_sub$Design)
prop.table(table(mydata_sub$Design))

table(mydata_sub$other_design)

table(mydata_sub$more_than_one_drug)
prop.table(table(mydata_sub$more_than_one_drug))

table(mydata_sub$Dose_expansion)
prop.table(table(mydata_sub$Dose_expansion))

table(mydata_sub$phase2)
prop.table(table(mydata_sub$phase2))

table(mydata_sub$phase2,mydata_sub$Dose_expansion)




