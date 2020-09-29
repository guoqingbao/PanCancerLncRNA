
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
library(reticulate)


#set the project path
project_path = "./project/"

#load patient data (include survival info)
study_cohort = read.xlsx(paste(project_path, "data/tcga/study_cohort.xlsx",sep=""), sheet = 1, colNames = TRUE)
entire_cohort = read.xlsx(paste(project_path, "data/tcga/entire_cohort.xlsx",sep=""), sheet = 1, colNames = TRUE)

length(rownames(study_cohort))

length(rownames(entire_cohort))

median(entire_cohort$ENSG00000206567.8)

#To build a risk score system based on lncRNA signature, we first estimate the coef of the signature
fit = coxph(Surv(OS, Censor) ~ 

ENSG00000187185.4+
ENSG00000259641.4+
ENSG00000218510.5 + 
ENSG00000257989.1+
ENSG00000206567.8
, data=study_cohort)

# summary(fit)

HR <- exp(coef(fit))
CI <- exp(confint(fit))
P <- round(coef(summary(fit))[,5],8)
COEF <- coef(fit)

colnames(CI) <- c("Lower", "Higher")

table2 <- as.data.frame(cbind(HR, CI,COEF, P))
table2

#calculate the risk score based on expression and coef of the signature
#use coef (estimated in the study cohort) in above table
entire_cohort$risk_lnc = (
    entire_cohort$ENSG00000187185.4   * table2[1,4] +
    entire_cohort$ENSG00000259641.4 * table2[2,4]+
    entire_cohort$ENSG00000218510.5 * table2[3,4] +
    entire_cohort$ENSG00000257989.1 * table2[4,4] +
    entire_cohort$ENSG00000206567.8 * table2[5,4] )

#calculate hazard ratio
fit1 = coxph(Surv(OS, Censor) ~ 

risk_lnc
, data=entire_cohort)

# summary(fit)

HR <- exp(coef(fit1))
CI <- exp(confint(fit1))
P <- coef(summary(fit1))[,5]
COEF <- coef(fit1)

colnames(CI) <- c("Lower", "Higher")

table22 <- as.data.frame(cbind(HR, CI,COEF, P))
table22

#calculate the time-dependent roc (lncrna signature on entire tcga cohort)
library(timeROC)
ROC<-timeROC(T=entire_cohort$OS,# survival time

delta=entire_cohort$Censor,# results, alive or death

marker=entire_cohort$risk_lnc,# variables

cause=1,

weighting="marginal",

times=c(365*3, 365*5, 365*10),# time, survival 3,5,10 years

ROC = TRUE,

iid = TRUE)
ROC

#the time-dependent roc/auc shown above result and confidence interval shown below
#the lncrna signature performance on entire TCGA dataset
confint(ROC)$CB_AUC

#To build a risk score system based on lncRNA signature + 3 clinical factors, 
#we first estimate the coef of the signature and 3 clinical factors

fit = coxph(Surv(OS, Censor) ~ stage + gender + age_at_diagnosis +

ENSG00000187185.4+
ENSG00000259641.4+
ENSG00000218510.5 + 
ENSG00000257989.1+
ENSG00000206567.8
, data=study_cohort)


# summary(fit)

HR <- exp(coef(fit))
CI <- exp(confint(fit))
P <- round(coef(summary(fit))[,5],8)
COEF <- coef(fit)

colnames(CI) <- c("Lower", "Higher")

table2 <- as.data.frame(cbind(HR, CI,COEF, P))
table2

#calculate the risk score based on expression and coef of the signature and 3 clinical factors
#use coef (estimated in the study cohort) in above table
entire_cohort$risk_lncClinical = (
    entire_cohort$stage   *         table2[1,4]+
    entire_cohort$gender    *       table2[2,4]+
    entire_cohort$age_at_diagnosis   * table2[3,4] +
    entire_cohort$ENSG00000187185.4   * table2[4,4] +
    entire_cohort$ENSG00000259641.4 * table2[5,4]+
    entire_cohort$ENSG00000218510.5 * table2[6,4] +
    entire_cohort$ENSG00000257989.1 * table2[7,4] +
    entire_cohort$ENSG00000206567.8 * table2[8,4] )

#calculate the time-dependent roc (lncrna signature + 3 clinical factors on entire cohort)
library(timeROC)
ROCAll<-timeROC(T=entire_cohort$OS,# survival time

delta=entire_cohort$Censor,# results, alive or death

marker=entire_cohort$risk_lncClinical,

cause=1,# positive results

weighting="marginal",

times=c(365*3, 365*5, 365*10),# time, survival 3,5,10 years

ROC = TRUE,

iid = TRUE)
ROCAll

#the time-dependent roc/auc is shown above and confidence interval shown below
#the lncrna signature + 3 clinical factors performance on entire dataset
confint(ROCAll)$CB_AUC

title = paste(project_path,"results/tcga_time-dependent-roc", ".pdf", sep="")

# pdf(title)

#(95% CI: 0.6328-0.6835), lncrna, 5 years
#(95% CI: 0.6771-0.7596), lncrna + 3 clinical factors, 10 years

#plot the time-dependent roc
plot(ROC,time=365*3,col = "red", lty=2,add =FALSE, title="")
plot(ROC,time=365*5,col = "red",lwd=2, add =TRUE)
plot(ROC,time=365*10,col = "red",lty=3, add =TRUE)

plot(ROCAll,time=365*3,col = "blue",lty=2,add =TRUE)
plot(ROCAll,time=365*5,col = "blue",lty=3,add =TRUE)
plot(ROCAll,time=365*10,col = "blue",lwd=2,add =TRUE)
legend("bottomright", c("5-lncRNA at 3 years (AUC = 0.64)","5-lncRNA at 5 years (AUC = 0.66)","5-lncRNA at 10 years (AUC = 0.63)", "5-lncRNA + 3 clinical factors at 3 years (AUC= 0.68)","5-lncRNA + 3 clinical factors at 5 years (AUC= 0.67)", "5-lncRNA + 3 clinical factors 10 years (AUC= 0.72)"), lty=c(2,1,3,2,3,1), lwd=c(1,2,1,1,1,2), col = c("red","red","red","blue","blue", "blue"), bty="n")

# dev.off()

#calculate the median risk score (lncrna signature)
mediam_risk_score = median(as.numeric(as.vector(entire_cohort$risk_lnc)))
entire_cohort$risk_cluster <- 1
entire_cohort$risk_cluster[which(entire_cohort$risk_lnc <= mediam_risk_score)] <- 0  
mediam_risk_score

#calculate the surviva different of the high- and low-risk score groups (lncrna signature)
survdiff(Surv(OS, Censor) ~risk_cluster, data = entire_cohort)

#plot the survival difference
library(survival)
require("survival")
# install.packages("svglite")
library(survminer)

dif = survdiff(Surv(OS, Censor) ~risk_cluster, data = entire_cohort)
pv = pchisq(dif$chisq, length(dif$n)-1, lower.tail = FALSE)
pv = formatC(pv, format = "e", digits = 2) 
ret = c("Log-rank\r\np = ",pv)
pvalue=paste(ret[1],ret[2], sep="")

fit<- survfit(Surv(OS, Censor) ~ risk_cluster, data = entire_cohort)

p <- ggsurvplot(fit, data=entire_cohort,risk.table = TRUE,tables.theme = theme_cleantable(),palette = c("#E7B890", "#2E9FCF"), ggtheme = theme_bw(),surv.median.line = "hv",legend.title = "Risk Score", legend.labs = c("Low", "High"), xlab='Time (days)', conf.int = TRUE,pval = pvalue,pval.method=TRUE,conf.int.style='ribbon', conf.int.alpha=0.1)
 
p = p + theme_survminer( font.main = c(16, "bold", "darkblue"), font.submain = c(15, "bold", "purple"), font.caption = c(14, "plain", "orange"), font.x = c(14, "bold"), font.y = c(14, "bold"), font.tickslab = c(12, "plain") )

p

#save the results
title = paste(project_path, "results/RiskScore-Survival-TCGA.tiff", sep="")
tiff(title, width=480*7, height=480*7, units="px", res=96*6, compression = "lzw")
print(p, newpage = FALSE)
dev.off()

# let's see the median survival time between high- and low-risk groups
# risk_cluster=0, low-risk score group (low risk), longer survival time
#risk_cluster=1, high-risk score group (high risk), shorter survival time
print(fit)

#we need to drop patients without age_at_diagnosis or OS info
validPatients = entire_cohort[!is.na(entire_cohort$OS),]
validPatients = validPatients[!is.na(validPatients$age_at_diagnosis),]
mediam_risk_score_lncClinical = median(as.numeric(as.vector(validPatients$risk_lncClinical)))
#stratified into two groups based on median risk score
validPatients$risk_cluster1 <- 1
validPatients$risk_cluster1[which(validPatients$risk_lncClinical <= mediam_risk_score_lncClinical)] <- 0 
mediam_risk_score_lncClinical

#calculate the survival difference
#risk_cluster1=0, low-risk score group (low risk), longer survival time
#risk_cluster1=1, high-risk score group (high risk), shorter survival time
fit<- survfit(Surv(OS, Censor) ~ risk_cluster1, data = validPatients)
fit

dif = survdiff(Surv(OS, Censor) ~risk_cluster1, data = validPatients)
pv = pchisq(dif$chisq, length(dif$n)-1, lower.tail = FALSE)
pv = formatC(pv, format = "e", digits = 2) 
pv # p-value
