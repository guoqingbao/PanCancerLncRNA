
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import warnings;
warnings.filterwarnings('ignore');
import xlrd
import matplotlib
from matplotlib import pyplot as plt
import csv
import json
import os
import pickle
import seaborn as sns

from collections import defaultdict
from sklearn.metrics import cohen_kappa_score
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import RFECV, RFE
from sklearn.metrics import f1_score,confusion_matrix, classification_report, accuracy_score
# from sklearn import cross_validation, metrics
from sklearn.metrics import roc_curve, auc, precision_recall_fscore_support, roc_auc_score
from sklearn.model_selection import StratifiedKFold,StratifiedShuffleSplit
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import cross_val_score, GridSearchCV


# ## 1. Load preprocessed data

# In[2]:


path = './project/'


# In[3]:


#preprocessed clinical dataset (4231 patient studies)
frame = pd.read_excel(path + 'data/clinical_entire_cohort.xlsx')
frame = frame[frame.gender==frame.gender]


#lncrnas reported by TCGA for pan-cancer analysis
lncRnaSet = pd.read_csv(path + 'data/gene_set_lincRNA.2018-09-20.tsv',delimiter = '\t')

#preprocessed genome-width dataset
geneframe = pd.read_pickle(path + 'data/geneframe_id.pkl')
validFrame = frame.loc[frame.id.isin(list(geneframe.columns)),:]
validFrame = validFrame.reset_index(drop=True)
lncFrame = geneframe.ix[(index for index, row in geneframe.iterrows() if row.name[:row.name.find('.')] in list(lncRnaSet.id)), :]
len(frame)


# In[4]:


def stage(s):
    if s.find('iv') !=-1:
        return 4
    elif s.find('iii') !=-1:
        return 3
    elif s.find('ii') !=-1:
        return 2
    elif s.find('i') !=-1 or s.find('x') !=-1:
        return 1
    else:
        return 0
    
#the study cohort identified from entire cohort    
study_cohort_clinical = pd.read_excel(path + 'data/clinical_study_cohort.xlsx', index_col=0)
len(study_cohort_clinical)


# ## 3. Combine gene dataset with corresponding clinical dataset

# In[5]:


study_cohort = lncFrame.T.loc[(row['id'] for index, row in study_cohort_clinical.iterrows() ), :]
gender = {'male': 1,'female': 0} 
race = {'unknown': 0,'white': 1, 'asian':2, 'black or african american':3, 'american indian or alaska native':4, 'native hawaiian or other pacific islander':5} 
study_cohort.loc[:,'gender'] = [gender[item] for item in study_cohort_clinical.gender] 
study_cohort.loc[:,'race'] = [race[item] for item in study_cohort_clinical.race] 
study_cohort.loc[:,'age_at_diagnosis']  = list(round(study_cohort_clinical.age_at_diagnosis/365,0))
study_cohort.loc[:,'stage']  = [stage(item) for item in study_cohort_clinical.tumor_stage]
types = list(set(study_cohort_clinical.project_id))
global_random_state = 5

def getindex(tp):
    for i in range(len(types)):
        if types[i] == tp:
            return i + 1
    return 0
study_cohort.loc[:,'type'] = [getindex(item) for item in study_cohort_clinical.project_id]

study_cohort.loc[:,'threehalf'] = list(study_cohort_clinical.threehalf)

data = study_cohort.sample(frac=1, random_state = global_random_state)
data.fillna(0, inplace=True)


# ## 4. Normalize and separate train/cross-validation and test datasets.

# In[6]:


data.ix[:,:-1]= (data.ix[:,:-1] - data.ix[:,:-1].mean())/data.ix[:,:-1].std()
x_train, x_test, y_train, y_test = train_test_split(data.ix[:,:-1], data.ix[:,-1], stratify=data.ix[:,-1], test_size=0.15, shuffle=True, random_state=global_random_state)


# In[7]:


train_backup = x_train.copy()
test_backup = x_test.copy()

x_train.reset_index(drop=True,inplace=True)
x_test.reset_index(drop=True,inplace=True)
y_train.reset_index(drop=True,inplace=True)
y_test.reset_index(drop=True,inplace=True)
print("Training/validation cases {} ; Testing cases {}".format(len(x_train),len(x_test)))


# In[8]:


training_data_in_study_cohort = pd.read_excel(path + 'data/study_cohort.xlsx', index_col=0).loc[train_backup.index]


# In[9]:


# training_data_in_study_cohort.to_excel(path + "results/rf1/train_data_in_study.xlsx")


# ## 5. Use identifed lncRNAs and/or clinical features

# In[9]:


lncfeatures = [
    'ENSG00000206567.8',
    'ENSG00000259641.4',
'ENSG00000187185.4',
    'ENSG00000257989.1',
    'ENSG00000218510.5'
]

clinicalfeatures = ['age_at_diagnosis','stage', 'gender']
lncClinicFeatures = []
lncClinicFeatures.extend(lncfeatures)
lncClinicFeatures.extend(clinicalfeatures)


# ## 6. Tune hyperparameters under nested cross-validation

# In[18]:


p_grid = {"n_estimators": [500, 1000],
          "max_depth": [2, 4, 8],
         "min_samples_leaf": [4, 8],
         "min_samples_split": [15, 25]}

clf = GridSearchCV(estimator=RandomForestClassifier(), param_grid=p_grid, scoring = "accuracy", cv=5, n_jobs=-1)
clf.fit(x_train.ix[:, lncClinicFeatures], y_train)

# Nested CV with parameter optimization
nested_score = cross_val_score(clf, X=x_train.ix[:, lncClinicFeatures], y=y_train, cv=5, n_jobs=-1)


# In[19]:


print(np.mean(nested_score))


# In[20]:


clf.best_params_

#{'max_depth': 8,
# 'min_samples_leaf': 4,
# 'min_samples_split': 25,
# 'n_estimators': 1000}


# # 7. Testing with tuned hyperparameters

# In[21]:


#lncrna+clinical performance on test set
classifier_rf = RandomForestClassifier(n_estimators=1000, max_depth=8, oob_score=True, min_samples_split=25, random_state=global_random_state,
                           min_samples_leaf=4)  
clf = classifier_rf.fit(x_train.ix[:, lncClinicFeatures], y_train)
lncClinicTestProb = clf.predict_proba(x_test.ix[:, lncClinicFeatures])


# In[22]:


#lncrna performance on test set
classifier_rf = RandomForestClassifier(n_estimators=1000, max_depth=8, oob_score=True, min_samples_split=25, random_state=global_random_state,
                           min_samples_leaf=4)  
clf = classifier_rf.fit(x_train.ix[:, lncfeatures], y_train)
lncTestProb = clf.predict_proba(x_test.ix[:, lncfeatures])


# In[23]:


#clinical performance on test set
classifier_rf = RandomForestClassifier(n_estimators=1000, max_depth=8, oob_score=True, min_samples_split=25, random_state=global_random_state,
                           min_samples_leaf=4)  
clf = classifier_rf.fit(x_train.ix[:, clinicalfeatures], y_train)
clinicTestProb = clf.predict_proba(x_test.ix[:, clinicalfeatures])


# In[24]:


testProbs = [lncClinicTestProb, lncTestProb, clinicTestProb]


# In[25]:


import numpy as np
from scipy import interp
import matplotlib.pyplot as plt
from itertools import cycle

from sklearn import svm, datasets


plt.rcParams['font.sans-serif']=['Arial']
plt.rcParams['axes.unicode_minus']=False 
ax = plt.gca()

fig = plt.gcf()
fig.set_size_inches( 7, 6)
from sklearn.metrics import precision_recall_fscore_support

colors = ['green', 'blue', 'grey']
names = ['LncRNA+Clinical','LncRNA','Clinical']
for i in range(len(testProbs)):
    fpr, tpr, thresholds = roc_curve(y_test, testProbs[i][:, 1])
    roc_auc = auc(fpr, tpr)
    plt.plot(fpr, tpr, lw=1.5, alpha=0.9,color=colors[i],
             label=names[i]+' ROC (AUC = %0.3f)' % (roc_auc))
      
plt.plot([0, 1], [0, 1], linestyle='--', lw=2, color='black')
plt.xlim([0, 1.0])
plt.ylim([0, 1.0])

plt.xlabel('1-Specificity',fontsize=12,fontweight='bold')
plt.ylabel('Sensitivity',fontsize=12,fontweight='bold')
plt.legend(loc="lower right")


ax.xaxis.set_tick_params(labelsize=12)
ax.yaxis.set_tick_params(labelsize=12)

plt.tight_layout()

plt.savefig(path + 'results/rf/ROC_curve1.svg',format='svg')

plt.savefig(path + 'results/rf/ROC_curve1.tiff',format='tiff', dpi=300, compression = "tiff_lzw")


plt.show()

