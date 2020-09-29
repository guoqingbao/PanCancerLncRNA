
library("openxlsx")
library("stringr")
library(ComplexHeatmap)
library(circlize)
library(colorspace)
library(GetoptLong)
library(ggplot2)
library(gplots)
library(reticulate)
# install.packages("svglite")

#set the project path
project_path = "./project/"

#load patient data (include survival info)
training_data = read.xlsx(paste(project_path, "data/tcga/training_data.xlsx",sep=""), sheet = 1, colNames = TRUE)

# study_cohort$OS = round(study_cohort$OS/30,2)

#convert expression value into expression levels
for(name in colnames(training_data[1,2:(length(colnames(training_data))-5)]))
{
    xl = training_data[,name[1]]
    training_data[,name[1]] <- 1
    training_data[,name[1]][which(xl <= median(as.numeric(as.vector(xl))))] <- 0   
}     

colnames(training_data)

#study_cohort$ENSG00000228623.3

library(survival)

tables <- NULL
for(name in colnames(training_data[1,2:(length(colnames(training_data))-5)]))
{
   fit = coxph(Surv(OS, Censor) ~  get(name), data=training_data)

summary(fit)

HR <- exp(coef(fit))
CI <- exp(confint(fit))
P <- coef(summary(fit))[,5]
    
colnames(CI) <- c("Lower", "Higher")

tmp <- as.data.frame(cbind(HR, CI, P),row.names=name)
    
    
    if(is.null(tables))
    {
    tables <- tmp
    }  
    else
    {
    tables <- rbind(tables, tmp)
    }
    
}
tables

length(rownames(tables))

tables[tables$P >0.05,] # 3 lncRNAs filtered in univariate cox analysis

remained = tables[tables$P <0.05,] # 23 lncRNAs remained after filtering

rownames(remained)

tables$HR = formatC(tables$HR , format = "f", digits = 3) 
tables$Lower = formatC(tables$Lower , format = "f", digits = 3) 
tables$Higher = formatC(tables$Higher , format = "f", digits = 3) 

tables$P = formatC(tables$P , format = "e", digits = 3) 
tables$CI = paste(tables$Lower, tables$Higher, sep="-")
tables = subset(tables, select=-c(Lower,Higher))

write.xlsx(tables, paste(project_path,"results/UnivariateCox.xlsx", sep=""),col.names=TRUE, row.names=TRUE)

library(survival)

fit = coxph(Surv(OS, Censor) ~ stage + gender + age_at_diagnosis +
            
ENSG00000214548.13 + ENSG00000251615.3 + ENSG00000226380.6 + ENSG00000225746.7 + ENSG00000257151.1 + ENSG00000206567.8 +
            ENSG00000253522.3 + ENSG00000218510.5 + ENSG00000172965.13 + ENSG00000229214.1 + ENSG00000227467.3 + 
            ENSG00000259641.4 + ENSG00000187185.4 + ENSG00000111206.11 + ENSG00000264247.1 + 
            ENSG00000254842.5 + ENSG00000182057.4 + ENSG00000272430.1 + ENSG00000257989.1 + ENSG00000212766.8 +
            ENSG00000183250.10 + ENSG00000253669.3 + ENSG00000146648.14
  , data=training_data)

HR <- exp(coef(fit))
CI <- exp(confint(fit))
P <- coef(summary(fit))[,5]

colnames(CI) <- c("Lower", "Higher")

table1 <- as.data.frame(cbind(HR, CI, P))
table1

lncTable=table1[4:length(rownames(table1)),]

lncTable[abs(lncTable$HR-1) <0.1,] # 11 lncRNAs filtered in MultivariateCox cox analysis

remained = lncTable[abs(lncTable$HR-1) >=0.1,] # 12 lncRNAs remained
rownames(remained)

remained

table1$HR = formatC(table1$HR , format = "f", digits = 3) 
table1$Lower = formatC(table1$Lower , format = "f", digits = 3) 
table1$Higher = formatC(table1$Higher , format = "f", digits = 3) 

table1$P = formatC(table1$P , format = "e", digits = 3) 
table1$CI = paste(table1$Lower, table1$Higher, sep="-")
table1 = subset(table1, select=-c(Lower,Higher))
write.xlsx(table1, paste(project_path,"results/MultivariateCox.xlsx", sep=""),col.names=TRUE, row.names=TRUE)

table2 <- NULL
for(name in rownames(remained))
{

dif = survdiff(Surv(OS, Censor) ~ get(name), data=training_data)
pv = pchisq(dif$chisq, length(dif$n)-1, lower.tail = FALSE)
    
P1 <- pv


tmp <- as.data.frame(cbind(P1),row.names=name)
    
    if(is.null(table2))
    {
    table2 <- tmp
    }  
    else
    {
    table2 <- rbind(table2, tmp)
    }
    
}
table2

table2[table2$P1 >=0.001,0] # 3 lncRNAs filtered

remained = table2[table2$P1 < 0.001,0]
rownames(remained) # 9 lncRNAs remained

table2$P1 = formatC(table2$P1 , format = "e", digits = 3) 
write.xlsx(table2, paste(project_path,"results/KaplanMeierAnalysis.xlsx", sep=""),col.names=TRUE, row.names=TRUE)
