
#required packages
library("openxlsx")
library("stringr")
library(ComplexHeatmap)
library(circlize)
library(colorspace)
library(GetoptLong)
library(ggplot2)
library(gplots)
library(survival)
require("survival")
# install.packages("svglite")
library(survminer)
# library(reticulate)


#set the project path
project_path = "./project/"

#load patient data (include survival info)
entire_cohort = read.xlsx(paste(project_path, "data/tcga/entire_cohort.xlsx",sep=""), sheet = 1, colNames = TRUE)

length(rownames(entire_cohort))

individual = entire_cohort[entire_cohort$type=='BRCA',] # change the cancer type to results for others
length(rownames(individual))

fit = coxph(Surv(OS, Censor) ~ 

ENSG00000187185.4+
ENSG00000259641.4+
ENSG00000218510.5 + 
ENSG00000257989.1+
ENSG00000206567.8
, data=individual)

# summary(fit)

HR <- exp(coef(fit))
CI <- exp(confint(fit))
P <- round(coef(summary(fit))[,5],8)
COEF <- coef(fit)

colnames(CI) <- c("Lower", "Higher")

table2 <- as.data.frame(cbind(HR, CI,COEF, P))
table2

#calculate the risk score
individual$risk_lnc = (
    individual$ENSG00000187185.4   * table2[1,4] +
    individual$ENSG00000259641.4 * table2[2,4]+
    individual$ENSG00000218510.5 * table2[3,4] +
    individual$ENSG00000257989.1 * table2[4,4] +
    individual$ENSG00000206567.8 * table2[5,4] )

#calculate the time-dependent roc 
library(timeROC)
ROC<-timeROC(T=individual$OS,# survival time

delta=individual$Censor,# results, alive or death

marker=individual$risk_lnc,# variables

cause=1,

weighting="marginal",

times=c(365*3, 365*5),# time, survival 3,5 years

ROC = TRUE,

iid = TRUE)
ROC
