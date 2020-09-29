
library("openxlsx")
library("stringr")
library(ComplexHeatmap)
library(circlize)
library(colorspace)
library(GetoptLong)
library(ggplot2)
library(gplots)
# library(reticulate)
# install.packages("svglite")
library(survival)
require("survival")
library(survminer)

#set the project path
project_path = "./project/"

#load patient data (include survival info)
entire_cohort = read.xlsx(paste(project_path, "data/tcga/entire_cohort.xlsx",sep=""), sheet = 1, colNames = TRUE)
entire_cohort = data.frame(entire_cohort)
entire_cohort$OS = round(entire_cohort$OS/365,2)

#stratify the entire cohort into different groups based on lncRNA expression
xl = entire_cohort[,"ENSG00000206567.8"]
entire_cohort$cluster1 <- 1
entire_cohort$cluster1[which(xl <= median(as.numeric(as.vector(xl))))] <- 0                        
                        
x2 = entire_cohort[,"ENSG00000187185.4"]
entire_cohort$cluster2 <- 1
entire_cohort$cluster2[which(x2 <= median(as.numeric(as.vector(x2))))] <- 0   

x3 = entire_cohort[,"ENSG00000259641.4"]
entire_cohort$cluster3 <- 1
entire_cohort$cluster3[which(x3 <= median(as.numeric(as.vector(x3))))] <- 0   

x4 = entire_cohort[,'ENSG00000218510.5']
entire_cohort$cluster4 <- 1
entire_cohort$cluster4[which(x4 <= median(as.numeric(as.vector(x4))))] <- 0   

x5 = entire_cohort[,"ENSG00000257989.1"]
entire_cohort$cluster5 <- 1
entire_cohort$cluster5[which(x5 <= median(as.numeric(as.vector(x5))))] <- 0   





dif = survdiff(Surv(OS, Censor) ~cluster1, data = entire_cohort)
pv = pchisq(dif$chisq, length(dif$n)-1, lower.tail = FALSE)
pv = formatC(pv, format = "e", digits = 2) 
ret = c("Log-rank\r\np = ",pv)
pvalue=paste(ret[1],ret[2], sep="")

fit<- survfit(Surv(OS, Censor) ~ cluster1, data = entire_cohort)
p <- ggsurvplot(fit, data=entire_cohort,risk.table = TRUE,tables.theme = theme_cleantable(),palette = c("#E7B800", "#2E9FDF"), ggtheme = theme_bw(),surv.median.line = "hv",legend.title = "Expression", legend.labs = c("Low", "High"), xlab='Time (years)', conf.int = TRUE,pval = pvalue,pval.method=TRUE,conf.int.style='ribbon', conf.int.alpha=0.1)
 
p = p + theme_survminer( font.main = c(16, "bold", "darkblue"), font.submain = c(15, "bold", "purple"), font.caption = c(14, "plain", "orange"), font.x = c(14, "bold"), font.y = c(14, "bold"), font.tickslab = c(12, "plain") )

p

#save figure
title = paste(project_path, "results/ENSG00000206567.8-Survival-AllData.tiff", sep="")
tiff(title, width=480*7, height=480*7, units="px", res=96*6, compression = "lzw")
print(p, newpage = FALSE)
dev.off()

dif = survdiff(Surv(OS, Censor) ~cluster2, data = entire_cohort)
pv = pchisq(dif$chisq, length(dif$n)-1, lower.tail = FALSE)
pv = formatC(pv, format = "e", digits = 2) 
ret = c("Log-rank\r\np = ",pv)
pvalue=paste(ret[1],ret[2], sep="")

fit<- survfit(Surv(OS, Censor) ~ cluster2, data = entire_cohort)
p <- ggsurvplot(fit, data=entire_cohort,risk.table = TRUE,tables.theme = theme_cleantable(),palette = c("#E7B800", "#2E9FDF"), ggtheme = theme_bw(),surv.median.line = "hv",legend.title = "Expression", legend.labs = c("Low", "High"), xlab='Time (years)', conf.int = TRUE,pval = pvalue,pval.method=TRUE,conf.int.style='ribbon', conf.int.alpha=0.1)
 
p = p + theme_survminer( font.main = c(16, "bold", "darkblue"), font.submain = c(15, "bold", "purple"), font.caption = c(14, "plain", "orange"), font.x = c(14, "bold"), font.y = c(14, "bold"), font.tickslab = c(12, "plain") )

p

#save figure
title = paste(project_path, "results/ENSG00000187185.4-Survival-AllData.tiff", sep="")
tiff(title, width=480*7, height=480*7, units="px", res=96*6, compression = "lzw")
print(p, newpage = FALSE)
dev.off()

dif = survdiff(Surv(OS, Censor) ~cluster3, data = entire_cohort)
pv = pchisq(dif$chisq, length(dif$n)-1, lower.tail = FALSE)
pv = formatC(pv, format = "e", digits = 2) 
ret = c("Log-rank\r\np = ",pv)
pvalue=paste(ret[1],ret[2], sep="")

fit<- survfit(Surv(OS, Censor) ~ cluster3, data = entire_cohort)
# title = paste(selectedLncRNA , " Expression")
p <- ggsurvplot(fit, data=entire_cohort,risk.table = TRUE,tables.theme = theme_cleantable(),palette = c("#E7B800", "#2E9FDF"), ggtheme = theme_bw(),surv.median.line = "hv",legend.title = "Expression", legend.labs = c("Low", "High"), xlab='Time (years)', conf.int = TRUE,pval = pvalue,pval.method=TRUE,conf.int.style='ribbon', conf.int.alpha=0.1)
 
p = p + theme_survminer( font.main = c(16, "bold", "darkblue"), font.submain = c(15, "bold", "purple"), font.caption = c(14, "plain", "orange"), font.x = c(14, "bold"), font.y = c(14, "bold"), font.tickslab = c(12, "plain") )

p

#save figure
title = paste(project_path, "results/ENSG00000259641.4-Survival-AllData.tiff", sep="")
tiff(title, width=480*7, height=480*7, units="px", res=96*6, compression = "lzw")
print(p, newpage = FALSE)
dev.off()

dif = survdiff(Surv(OS, Censor) ~cluster4, data = entire_cohort)
pv = pchisq(dif$chisq, length(dif$n)-1, lower.tail = FALSE)
pv = formatC(pv, format = "e", digits = 2) 
ret = c("Log-rank\r\np = ",pv)
pvalue=paste(ret[1],ret[2], sep="")

fit<- survfit(Surv(OS, Censor) ~ cluster4, data = entire_cohort)
# title = paste(selectedLncRNA , " Expression")
p <- ggsurvplot(fit, data=entire_cohort,risk.table = TRUE,tables.theme = theme_cleantable(),palette = c("#E7B800", "#2E9FDF"), ggtheme = theme_bw(),surv.median.line = "hv",legend.title = "Expression", legend.labs = c("Low", "High"), xlab='Time (years)', conf.int = TRUE,pval = pvalue,pval.method=TRUE,conf.int.style='ribbon', conf.int.alpha=0.1)
 
p = p + theme_survminer( font.main = c(16, "bold", "darkblue"), font.submain = c(15, "bold", "purple"), font.caption = c(14, "plain", "orange"), font.x = c(14, "bold"), font.y = c(14, "bold"), font.tickslab = c(12, "plain") )

p

#save figure
title = paste(project_path, "results/ENSG00000218510.5-Survival-AllData.tiff", sep="")
tiff(title, width=480*7, height=480*7, units="px", res=96*6, compression = "lzw")
print(p, newpage = FALSE)
dev.off()

dif = survdiff(Surv(OS, Censor) ~cluster5, data = entire_cohort)
pv = pchisq(dif$chisq, length(dif$n)-1, lower.tail = FALSE)
pv = formatC(pv, format = "e", digits = 2) 
ret = c("Log-rank\r\np = ",pv)
pvalue=paste(ret[1],ret[2], sep="")

fit<- survfit(Surv(OS, Censor) ~ cluster5, data = entire_cohort)
p <- ggsurvplot(fit, data=entire_cohort,risk.table = TRUE,tables.theme = theme_cleantable(),palette = c("#E7B800", "#2E9FDF"), ggtheme = theme_bw(),surv.median.line = "hv",legend.title = "Expression", legend.labs = c("Low", "High"), xlab='Time (years)', conf.int = TRUE,pval = pvalue,pval.method=TRUE,conf.int.style='ribbon', conf.int.alpha=0.1)
 
p = p + theme_survminer( font.main = c(16, "bold", "darkblue"), font.submain = c(15, "bold", "purple"), font.caption = c(14, "plain", "orange"), font.x = c(14, "bold"), font.y = c(14, "bold"), font.tickslab = c(12, "plain") )

p

#save figure
title = paste(project_path, "results/ENSG00000257989.1-Survival-AllData.tiff", sep="")
tiff(title, width=480*7, height=480*7, units="px", res=96*6, compression = "lzw")
print(p, newpage = FALSE)
dev.off()
