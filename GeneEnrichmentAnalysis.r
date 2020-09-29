
library("openxlsx")
library("stringr")
library(ComplexHeatmap)
library(circlize)
library(colorspace)
library(GetoptLong)
library(ggplot2)
library(gplots)
# install.packages('OpasnetUtils')
library(reticulate)

memory.limit(22000) # be sure you have enough memory to conduct the enrichment analysis

#set the project path
project_path = "./project/"

# be sure geneframe_id.pkl file is accessible
genome_width_dataset_file = paste(project_path,"data/geneframe_id.pkl", sep="")

xlgene = py_load_object(genome_width_dataset_file, pickle = "pickle")


#specify which lncRNA you want to analysis
# you can specify one of the five lncRNAs: ENSG00000206567.8, ENSG00000187185.4, ENSG00000259641.4, ENSG00000218510.5， ENSG00000257989.1
selectedLncRNA = 'ENSG00000206567.8' # by default
# selectedLncRNA = 'ENSG00000218510.5' 

index = grep(selectedLncRNA, rownames(xlgene)) # the index of specified lncRNA in genome-width dataset
index

clinical_entire_cohort = read.xlsx(paste(project_path, "data/tcga/clinical_entire_cohort.xlsx",sep=""), sheet = 1, colNames = TRUE)
# clinical_entire_cohort

#find correlated genes (co-expressed genes), it may takes 10 mins!
rows = length(rownames(xlgene))
x = xlgene[index,]
realY = t(xlgene[-index,])

cormat = round(cor(t(x),realY, method="spearman"),2) #make correlation matrix


length(colnames(cormat)) #total
length(which(cormat >= 0.4)) #number of positively related genes
length(which(cormat <= -0.3))#number of negatively related genes

library(scales)

#get positive related genes and negative related genes separately using cutoff 0.4 and -0.4
df_positive = realY[,which(cormat >= 0.4)]
df_negative = realY[,which(cormat <= -0.3)]

#filter genes by p-values, remove genes which have correlations but not statistic meaning
cols_positive = c()
for ( i in 1: length(colnames(df_positive)))
    {
     if (cor.test(df_positive[,i],t(x)[,1])$p.value < 0.05)
     {
        cols_positive <- c(cols_positive,i)
     }   
}
 
length(cols_positive)

# we have a lncRNA that don't have negatively correlated genes, remember to ommit this part when dealing with that lncRNA
cols_negative = c()
for ( i in 1: length(colnames(df_negative)))
    {
     if (cor.test(df_negative[,i],t(x)[,1])$p.value < 0.05)
     {
        cols_negative <- c(cols_negative,i)
     }   
}


length(cols_negative)

df_sel_positive = df_positive[,cols_positive] #get final selected positively related genes

#get final selected negatively related genes, one of the lncRNAs have no negatively correlated genes,
#remember ommit this sentence when dealing with it
df_sel_negative = df_negative[,cols_negative] 


#rotate to set cluster (positive 1, negative 2), after that we can combine them
frame_positive = data.frame(t(df_sel_positive))
frame_negative = data.frame(t(df_sel_negative))
frame_positive$cluster <- 1
frame_negative$cluster <- 2

frame_sel_all = rbind(frame_positive, frame_negative) #combine selected genes, we have cluster to identify them


                            
df = frame_sel_all[,1:(length(colnames(frame_sel_all))-1)] #dataframe without cluster column

df = apply(df, 1, rescale, to=c(-10,10)) #apply scale


df = data.frame(t(df))
df$cluster = frame_sel_all$cluster #add cluster back to scaled dataset

dftemp = data.frame(t(df)) #a temp

dftemp$GeneA <- c(t(x)[,1],-100)#add GeneA expression, used for sorting

clinical_entire_cohort$tumor_stage[is.na(clinical_entire_cohort$tumor_stage)] = "unknown"

write.xlsx(frame_negative[,1:10], paste(project_path,"results/negative5.xlsx", sep=""),col.names=TRUE, row.names=TRUE) 

empty <- list()
for (age in clinical_entire_cohort$age_at_diagnosis)
{
    if(!is.na(age))
    {
         if(age > 65*365)
        {
            empty <- c(empty, 2)
        }   
       else if(age >45*365)
       {
         empty <- c(empty, 1)
       }
      else 
      {
         empty <- c(empty, 0)
       }
    }
    else
        {
        empty <- c(empty, 3)
    }

}

library(OpasnetUtils)
#add patient info
dftemp$gender = c(sapply(clinical_entire_cohort$gender,switch,"male"=1,"female"=0, "unknown"=2),0)
dftemp$System = c(clinical_entire_cohort$System,0)

dftemp$age_at_diagnosis = c(empty,0)

dftemp$stage = c(sapply(clinical_entire_cohort$tumor_stage,switch ,"ivb"=2,"i/ii nos"=1,"ia"=1,"x"=2,"iiic"=2,"ii"=1,"iib"=1,"iia"=1,"iic"=1,"iii"=2,"iva"=2,"ib"=1,"iiib"=2,"i"=1,"iiia"=2,"ivc"=2,"iv"=2,"0"=1,"is"=1, "unknown"=0), 0)

#female cancer heatmap
dftemp$cancer_type = c(sapply(clinical_entire_cohort$project_id,switch,"OV"=0,"MESO"=5,"LAML"=5,"LUAD"=5,"LGG"=5,"UCEC"=1,"GBM"=5,"CHOL"=5,"READ"=5,"BLCA"=5,"TGCT"=5,"PRAD"=5,"THYM"=5,"UCS"=5,"CESC"=2,"KIRC"=5,"SKCM"=5,"DLBC"=5,"THCA"=5,"PAAD"=5,"COAD"=5,"HNSC"=5,"KIRP"=5,"ACC"=5,"LUSC"=5,"STAD"=5,"BRCA"=3,"SARC"=5,"UVM"=5,"KICH"=5,"PCPG"=5,"ESCA"=5,"LIHC"=5),0)


length(rownames(dftemp))

dftemp = dftemp[dftemp$gender!=2,]
length(rownames(dftemp)) # in total of 4231 studies, another one is cluster for enrichment analysis

length(colnames(dftemp))

dfsort = dftemp[order(dftemp$GeneA),] #sort gene expression, the patients are reordered
frame = data.frame(t(dfsort)) # rotated for heatmap plotting
patient_frame = frame[(length(rownames(frame))-5):(length(rownames(frame))),]

clear_frame = frame[1:(length(rownames(frame))-6),]
pt = data.frame(t(patient_frame[,2:length(colnames(frame))]))
pt <- data.frame(sapply(pt, function(x) as.numeric(as.character(x))))

pt

length(rownames(pt))

length(rownames(pt[pt$age_at_diagnosis==3,]))

clear_frame <- data.frame(sapply(clear_frame, function(x) as.numeric(as.character(x))))
median(pt$GeneA[pt$GeneA > median(pt$GeneA)])

ret = paste0(c("Gene Group"), clear_frame$cluster) # two parts of the heatmap

median(pt$GeneA)

min(pt$GeneA)

max(pt$GeneA)

length(colnames(clear_frame[,2:length(colnames(clear_frame))]))

#show the heatmapp, the first column is cluster (positive or negative related genes), which is useless when plotting heatmap
#but, we can use cluster to split the heatmap into two parts, the top part is positively related gene expression, the lower part is negatively related gene expression
#we set :cluster_rows = FALSE, cluster_columns = FALSE, we do not use any of clustering methods provided by Heatmap function, we just plot the sorted gene expressions
ha1 = HeatmapAnnotation(df = pt,
    col = list(
               gender = c("1"="black","0"="pink"),
               GeneA=colorRamp2(c(min(pt$GeneA),median(pt$GeneA), median(pt$GeneA[pt$GeneA > median(pt$GeneA)]), max(pt$GeneA)), c("blue","#ebebeb", "#FF5C5C", "red")),
               stage = c("1"="#1e90ff","2"="orange", "0"="grey"),
               age_at_diagnosis = c("0" = "#C0C0C0", "1" = "#778899","2"= "#696969" ,  "3"="#2F4F4F"),   
        cancer_type = c("0"="#f24433","1"="#960c14","2"="#d82623","3"="#bc1419","5"="blue"),
        System = c("1"="#1764ab","2"="#bb1419","3"="#157e3a","4"="#404040","5"="#61409b","6"="#f3701b","7"="#4294c3","8"="#99017b","9"="#cccc00")

    ),
                       
     annotation_legend_param = list(
      gender = list(title = "Gender", at = c("1", "0"), labels = c("Male", "Female")),
         cancer_type = list(title = "FRC VS. Others", at = c("0", "1", "2", "3", "5"), labels = c("OV", "UCEC", "CESC", "BRCA", "Other 29 Cancers")),
         System = list(title = "System", at = c("1", "2", "3", "4", "5","6","7","8","9"), labels = c("Respiratory system", "Reproductive system", "Digestive system", "Endocrine system"
                                                                                                     , "Urinary system","Central nervous system","Immune system","Skin","Others")),
         stage = list(title = "Tumor stage", at = c("2","1", "0"), labels = c("I&II","III&IV", "Not reported")),
        age_at_diagnosis = list(title = "Age at diagnosis", at = c("3","2","1", "0"), labels = c("Unknown","Age > 65", "45 < Age ≤65", "Age ≤ 45")),
         # replace median value with above analysis results
         GeneA = list(title = paste("\n\n\n\n\n\nExpression", sep=""))

      )
   )

ht = Heatmap(clear_frame[,2:length(colnames(clear_frame))], name = "Expression", split = sapply(ret,switch, "Gene Group2"="C2 - Negative", "Gene Group1"="C1 - Positive"),
        show_row_dend = FALSE,show_column_dend = FALSE, cluster_rows = FALSE,
        cluster_columns = FALSE, show_row_names = FALSE, show_column_names = FALSE,
        show_heatmap_legend = FALSE,
       ,top_annotation=ha1)



p = draw(ht, heatmap_legend_side = "top")
decorate_annotation("GeneA", {grid.text("Expression", unit(1, "npc") + unit(2, "mm"), 0.5, default.units = "npc", just = "left" , gp = gpar(fontsize = 10))})

decorate_annotation("gender", {grid.text("Gender", unit(1, "npc") + unit(2, "mm"), 0.5, default.units = "npc", just = "left", gp = gpar(fontsize = 10))})
decorate_annotation("System", {grid.text("System", unit(1, "npc") + unit(2, "mm"), 0.5, default.units = "npc", just = "left", gp = gpar(fontsize = 10))})
decorate_annotation("stage", {grid.text("Tumor Stage", unit(1, "npc") + unit(2, "mm"), 0.5, default.units = "npc", just = "left", gp = gpar(fontsize = 10))})

decorate_annotation("cancer_type", {grid.text("FRC VS. Others", unit(1, "npc") + unit(2, "mm"), 0.5, default.units = "npc", just = "left", gp = gpar(fontsize = 10))})
decorate_annotation("age_at_diagnosis", {grid.text("Age at diagnosis", unit(1, "npc") + unit(2, "mm"), 0.5, default.units = "npc", just = "left", , gp = gpar(fontsize = 10) )})


# save results
title = paste(project_path, "results/", selectedLncRNA, "-Heatmap.tiff", sep="")
tiff(title, width=861*3, height=840*3, units="px", res=96*3, compression = "lzw")
p = draw(ht, heatmap_legend_side = "top")
decorate_annotation("GeneA", {grid.text("Expression", unit(1, "npc") + unit(2, "mm"), 0.5, default.units = "npc", just = "left" , gp = gpar(fontsize = 10))})

decorate_annotation("gender", {grid.text("Gender", unit(1, "npc") + unit(2, "mm"), 0.5, default.units = "npc", just = "left", gp = gpar(fontsize = 10))})
decorate_annotation("System", {grid.text("System", unit(1, "npc") + unit(2, "mm"), 0.5, default.units = "npc", just = "left", gp = gpar(fontsize = 10))})
decorate_annotation("stage", {grid.text("Tumor Stage", unit(1, "npc") + unit(2, "mm"), 0.5, default.units = "npc", just = "left", gp = gpar(fontsize = 10))})

decorate_annotation("cancer_type", {grid.text("FRC VS. Others", unit(1, "npc") + unit(2, "mm"), 0.5, default.units = "npc", just = "left", gp = gpar(fontsize = 10))})
decorate_annotation("age_at_diagnosis", {grid.text("Age at diagnosis", unit(1, "npc") + unit(2, "mm"), 0.5, default.units = "npc", just = "left", , gp = gpar(fontsize = 10) )})

dev.off()

