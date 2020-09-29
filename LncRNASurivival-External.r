
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
# study_cohort = read.xlsx(paste(project_path, "data/study_cohort.xlsx",sep=""), sheet = 1, colNames = TRUE)
cptac_cohort = read.xlsx(paste(project_path, "data/external/cptac_cohort.xlsx",sep=""), sheet = 1, colNames = TRUE)
target_cohort = read.xlsx(paste(project_path, "data/external/target_cohort.xlsx",sep=""), sheet = 1, colNames = TRUE)

length(rownames(target_cohort))

length(rownames(cptac_cohort))

#calculate the risk score based on expression and coef of the signature
#use coef estimated in TCGA study cohort 


target_cohort$risk_lnc = (
    target_cohort$ENSG00000187185.4   * 0.11098351	+
    target_cohort$ENSG00000259641.4 * -0.76093818+
    target_cohort$ENSG00000218510.5 * -0.03033684 +
    target_cohort$ENSG00000257989.1 * -0.12374097 +
    target_cohort$ENSG00000206567.8 * -0.14095109 )


cptac_cohort$risk_lnc = (
    cptac_cohort$ENSG00000187185.4   * 0.11098351+
    cptac_cohort$ENSG00000259641.4 * -0.76093818+
    cptac_cohort$ENSG00000218510.5 * -0.03033684 +
    cptac_cohort$ENSG00000257989.1 * -0.12374097 +
    cptac_cohort$ENSG00000206567.8 * -0.14095109 )

#calculate the time-dependent roc (lncrna signature on external cohort)
library(timeROC)
targetROC<-timeROC(T=target_cohort$OS,# survival time

delta=target_cohort$Censor,# results, alive or death

marker=target_cohort$risk_lnc,# variables

cause=1,

weighting="marginal",

times=c(365*5,365*3),# time, survival at 3, 5 years

ROC = TRUE,

iid = TRUE)
targetROC

library(timeROC)
cptacROC<-timeROC(T=cptac_cohort$OS,# survival time

delta=cptac_cohort$Censor,# results, alive or death

marker=cptac_cohort$risk_lnc,# variables

cause=1,

weighting="marginal",

times=c(365*5,365*3),# time, survival at 3, 5 years

ROC = TRUE,

iid = TRUE)
cptacROC

#calculate hazard ratio
fit1 = coxph(Surv(OS, Censor) ~ 

risk_lnc
, data=target_cohort)

# summary(fit)

HR <- exp(coef(fit1))
CI <- exp(confint(fit1))
P <- coef(summary(fit1))[,5]
COEF <- coef(fit1)

colnames(CI) <- c("Lower", "Higher")

table22 <- as.data.frame(cbind(HR, CI,COEF, P))
table22

#calculate hazard ratio
fit1 = coxph(Surv(OS, Censor) ~ 

risk_lnc
, data=cptac_cohort)

# summary(fit)

HR <- exp(coef(fit1))
CI <- exp(confint(fit1))
P <- coef(summary(fit1))[,5]
COEF <- coef(fit1)

colnames(CI) <- c("Lower", "Higher")

table22 <- as.data.frame(cbind(HR, CI,COEF, P))
table22

title = paste(project_path,"results/target_cptac_time-dependent-roc", ".pdf", sep="")

#plot the time-dependent roc
plot(targetROC,time=365*3,col = "red",lwd=2,add =FALSE, title="")

plot(targetROC,time=365*5,col = "blue",lwd=2,add =TRUE)

plot(cptacROC,time=365*3,col = "green",lwd=2,add =TRUE)

legend("bottomright", c( "TARGET: ROC at 3 years (AUC = 0.60)", "TARGET: ROC at 5 years (AUC = 0.59)", "CPTAC: ROC at 3 years (AUC = 0.75)"), lty=c(1,1,1), lwd=c(2,2,2),col = c("red", "blue","green"), bty="n")


# dev.off()

#calculate the median risk score (lncrna signature)
mediam_risk_score = median(as.numeric(as.vector(target_cohort$risk_lnc)))
# mediam_risk_score = -0.221
target_cohort$risk_cluster <- 1
target_cohort$risk_cluster[which(target_cohort$risk_lnc <= mediam_risk_score)] <- 0  
# mediam_risk_score

#plot the survival difference
library(survival)
require("survival")
# install.packages("svglite")
library(survminer)

dif = survdiff(Surv(OS, Censor) ~risk_cluster, data = target_cohort)
pv = pchisq(dif$chisq, length(dif$n)-1, lower.tail = FALSE)
pv = formatC(pv, format = "e", digits = 2) 
ret = c("Log-rank\r\np = ",pv)
pvalue=paste(ret[1],ret[2], sep="")

fit<- survfit(Surv(OS, Censor) ~ risk_cluster, data = target_cohort)

p <- ggsurvplot(fit, data=target_cohort,risk.table = TRUE,tables.theme = theme_cleantable(),palette = c("#E7B890", "#2E9FCF"), ggtheme = theme_bw(),surv.median.line = "hv",legend.title = "Risk Score", legend.labs = c("Low", "High"), xlab='Time (days)', conf.int = TRUE,pval = pvalue,pval.method=TRUE,conf.int.style='ribbon', conf.int.alpha=0.1)
 
p = p + theme_survminer( font.main = c(16, "bold", "darkblue"), font.submain = c(15, "bold", "purple"), font.caption = c(14, "plain", "orange"), font.x = c(14, "bold"), font.y = c(14, "bold"), font.tickslab = c(12, "plain") )

p

# save the results
title = paste(project_path, "results/RiskScore-Survival-Target.tiff",sep="")
tiff(title, width=480*7, height=480*7, units="px", res=96*6, compression = "lzw")
print(p, newpage = FALSE)
dev.off()

mediam_risk_score = median(as.numeric(as.vector(cptac_cohort$risk_lnc)))
# mediam_risk_score = -0.221
cptac_cohort$risk_cluster <- 1
cptac_cohort$risk_cluster[which(cptac_cohort$risk_lnc <= mediam_risk_score)] <- 0  
# mediam_risk_score


dif = survdiff(Surv(OS, Censor) ~risk_cluster, data = cptac_cohort)
pv = pchisq(dif$chisq, length(dif$n)-1, lower.tail = FALSE)
pv = formatC(pv, format = "e", digits = 2) 
ret = c("Log-rank\r\np = ",pv)
pvalue=paste(ret[1],ret[2], sep="")

fit<- survfit(Surv(OS, Censor) ~ risk_cluster, data = cptac_cohort)

p <- ggsurvplot(fit, data=cptac_cohort,risk.table = TRUE,tables.theme = theme_cleantable(),palette = c("#E7B890", "#2E9FCF"), ggtheme = theme_bw(),surv.median.line = "hv",legend.title = "Risk Score", legend.labs = c("Low", "High"), xlab='Time (days)', conf.int = TRUE,pval = pvalue,pval.method=TRUE,conf.int.style='ribbon', conf.int.alpha=0.1)
 
p = p + theme_survminer( font.main = c(16, "bold", "darkblue"), font.submain = c(15, "bold", "purple"), font.caption = c(14, "plain", "orange"), font.x = c(14, "bold"), font.y = c(14, "bold"), font.tickslab = c(12, "plain") )

p

title = paste(project_path, "results/RiskScore-Survival-CPTAC.tiff",sep="")
tiff(title, width=480*7, height=480*7, units="px", res=96*6, compression = "lzw")
print(p, newpage = FALSE)
dev.off()
