# PanCancerLncRNA
An intuitive and streamlined pan-cancer analysis framework for identification of lncRNAs associated with pan-cancer prognosis


### Please cite our work if you found it is useful for your research.

Guoqing Bao and Ran Xu et al., “Identification of lncRNA Signature Associated With Pan-cancer Prognosis”, IEEE Journal of Biomedical and Health Informatics, 2020, In press.

## Prerequisites

The following libraries are required:

Python: IPython, matplotlib, pandas, scikit-learn, plotly
R: survival, ggplot, ComplexHeatmap, svglite, reticulate

To run jupyter notebook (ipynb files), you also need to install: jupyter notebook/lab (https://jupyter.org/), R kernel for jupyter notebook (https://github.com/IRkernel/IRkernel), and R-3.5 (https://cran.r-project.org/bin/)

Jupyter notebook files require Python 3 Kernel: LncRNAPreprocessing.ipynb, RFStudyCohort.ipynb, SankeyDiagram.ipynb

Jupyter notebook files require R 3.5 Kernel: CoxKaplanAnalysis.ipynb, SurvivalPlotAllData.ipynb, LncRNASurivival.ipynb, GeneEnrichmentAnalysis.ipynb


#### All data used in this study can be found in folder "data".


## The whole project includes the following major components:

### 1. Data preprocessing
### (LncRNAPreprocessing.py)

To run the preprocessing code, you need to first download GDC transfer tool (https://gdc.cancer.gov/access-data/gdc-data-transfer-tool) and then use the provided metadata (under data/gdc_manifest.2018-09-20.txt) and GDC transfer tool to download RNA-Seq data (FPKM gene files).

You also need to download corresponding clinical data from the TCGA data portal (https://portal.gdc.cancer.gov/) or use the data we cached before (raw_clinical.xlsx, was converted from raw_clinical.tsv).

The preprocessing code will combine FPKM files into a single genome-wide dataset (data/geneframe_id.pkl). The code will also drop invalid data records and save full valid clinical data into data/clinical_entire_cohort.xlsx. The study cohort is then identified and corresponding clinical data is saved into data/clinical_study_cohort.xlsx.

### 2. Feature selection by Cox and Kaplan-Meier survival analysis
### (CoxKaplanAnalysis.R)

In this stage, we select prognostic related lncRNAs using R (version 3.5). You can either use R code provided. Using the training data, in this stage, we first performed Univariate and multivariate cox regression analysis, then we used Kaplan-Meier survival analysis to identify the most valuable lncRNAs based on outputs of Machine learning feature selection.

### 3. The machine learning model for prognosis prediction
### (LncRnaProject.py)

We carried out cross-validation and identified 5 lncRNAs as a prognostic signature. We trained an RF-based model and validated its performance for prognosis prediction on test data.

### 4. Prognostic validation of lncRNA signature on entire tcga cohort and two external cohorts
### (SurvivalPlotAllData.R/LncRNASurivival.R/LncRNASurivival-External.R)

In this stage, we validated the prognostic value of the signature on the entire TCGA cohort and two external cohorts. We first validated each of the LncRNA in the signature by carrying out Kaplan-Meier survival analysis on the entire TCGA cohort (SurvivalPlotAllData.R). Then, we performed time-dependent ROC analysis, created and validated a risk score system using the lncRNA signature on TCGA, TARGET, and CPTAC.  We also identified the median survival time between two risk groups that stratified by our signature (median risk score).

### 5. Caner prognosis indexing system 
### (SankeyDiagram.py)

In this stage, we built an indexing system based on lncRNA expression levels and corresponding patients' prognosis outcomes (reported by TCGA). To run the code, you need to install plotly (https://plot.ly/) and orca (https://github.com/plotly/orca). 

### 6. Gene enrichment and pathway analysis 
### (Metascape and GeneEnrichmentAnalysis.R)

In this stage, we identified positively- and negatively-correlated genes and analyzed the relationship between the correlated genes and corresponding clinical factors. We plotted the correlated genes as a heatmap (gene enrichment heatmap) and clinical factors are annotated accordingly. You can also use Metascape (http://metascape.org/) to analysis the pathways and biological processes related to the signature by using the correlated genes identified here.



