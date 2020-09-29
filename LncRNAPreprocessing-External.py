
# coding: utf-8

# In[62]:


import numpy as np
import pandas as pd
import xlrd
import matplotlib
from matplotlib import pyplot as plt
import csv
import json
import os


# In[63]:


path = './project/'


# In[64]:


# the path where FPKM files downloaded (need to download)
gene_raw_data_path = "./project/data/externalsets/"


# In[65]:


#the algorithm process FPKM file and generate the pkl file (geneframe_external.pkl)
pkl_path = './project/data/pkl/'


# ## Load lncRNA metadata and clinical data

# In[66]:


# corresponding full clinical data
clinical = pd.read_excel(path + 'data/raw/external_raw_clinical.xlsx') 
# the lncrnas that reported by TCGA for pan-cancer analysis
lncRnaSet = pd.read_csv(path + 'data/raw/gene_set_lincRNA.2018-09-20.tsv',delimiter = '\t')
#lncrna metada, used for download the FPKM files (with the metadata, you can use GDC transfer tool to redownload them)
lncRnaFiles = pd.read_csv(path + 'data/raw/external_gdc_manifest.2020-03-03.txt',delimiter = '\t')
#for mapping clinical data with gene expression data
with open(path + 'data/raw/external_files.2020-03-04.json') as json_file:
    jsonFile = json.load(json_file)


# In[67]:


lncRnaSet[:10] #let's see some lncrna ids


# In[11]:


lncRnaFiles[:10] #let's see some lncrna metadata


# In[68]:


len(clinical)#let's see some clinical data


# In[72]:


set(clinical.primary_diagnosis)


# ## Binding clinical data with lncRNA metadata file

# In[13]:


from collections import defaultdict
dictID = defaultdict()
for item in jsonFile:
    dictID[item['file_name']] = item['cases'][0]['case_id']


# In[14]:


caseids = [dictID[row['filename']] for index, row in lncRnaFiles.iterrows()]


# In[15]:


lncRnaFullFiles = lncRnaFiles.copy()
lncRnaFullFiles['case_id'] = caseids


# In[16]:


len(set(lncRnaFullFiles.case_id))


# In[17]:


lncRnaFullFiles = lncRnaFullFiles.drop_duplicates(subset='case_id')
lncRnaFullFiles.reset_index(inplace=True)
clinical_full = pd.merge(clinical, lncRnaFullFiles, how='left', on='case_id')
len(clinical_full)


# In[21]:


clinical_full = clinical_full.dropna(axis=1, how='all')
clinical_full = clinical_full.dropna(axis=0, how='all')
frame = clinical_full.drop(['md5','size','state_y'], axis=1)


# ## Unzip FPKM files and merge them as a genome-width dataframe 

# In[54]:


problems = []
emptys = []
startIndex = 0
import gzip
from io import StringIO
curPercent = 0


# In[ ]:


remain_frames = []
for index, row in frame.iterrows():
    if index < startIndex:
        continue
    percent = round(index*1.0/frame.shape[0]*100,0)
    if percent != curPercent:
        curPercent = percent
        print("progress: {}， curIndex {}".format(percent,index))
        
    if row['id'] != row['id']:
        emptys.append(row['case_id'])     
    else:
        try:           
            f=gzip.open(gene_raw_data_path + row['id'] + str("/") + row['filename'],'r')
            tmp = pd.read_csv(StringIO(f.read().decode("utf-8")), delimiter = '\t',header=None).set_index(0)
            tmp.rename(columns={1:row['id']}, inplace=True)
            remain_frames.append(tmp)
        except KeyboardInterrupt:
            remain_frame = pd.concat(remain_frames,axis=1)
            remain_frame.to_pickle(pkl_path + 'geneframe_external.pkl')
#             print("progress: {}， curIndex {}".format(percent,index))
            print("interrupted!")
            break
        except:
            print(row['case_id'], " problem:", row['id'] + str("/") + row['filename'])
            problems.append(row['id'] + str("/") + row['filename'])
            pass
remain_frame = pd.concat(remain_frames,axis=1)
# The lncRNA dataframe, contains genome-width FPKM expression
remain_frame.to_pickle(pkl_path + 'geneframe_external.pkl')


# In[26]:


geneframe = pd.read_pickle(pkl_path +'geneframe_external.pkl')


# ## Drop invalid patient data

# In[27]:


validFrame = frame.loc[frame.id.isin(list(geneframe.columns)),:]
validFrame = validFrame.reset_index(drop=True)


# In[29]:


#save clinical info (entire cohort)
validFrame.tumor_stage.fillna('unknown',inplace=True)
validFrame.race.fillna('unknown', inplace=True)
validFrame.gender.fillna('unknown', inplace=True)
validFrame.tumor_stage = validFrame.tumor_stage.str.replace("stage ","")
# validFrame.to_excel(path + 'data/external/clinical_external.xlsx')


# In[33]:


f1 = validFrame[(validFrame.vital_status=='Alive')|(validFrame.vital_status=='Dead')]


# In[34]:


f1.reset_index(inplace=True)


# In[35]:


f2 = f1[(f1.days_to_death!='--')|(f1.days_to_last_follow_up!='--')]


# In[36]:


validFrame = f2.drop(['level_0'],axis=1).reset_index().drop(['level_0'],axis=1)


# ## We extract LncRNA frame from genome-width dataset

# In[37]:


lncFrame = geneframe.ix[(index for index, row in geneframe.iterrows() if row.name[:row.name.find('.')] in list(lncRnaSet.id)), :]


# In[43]:


# validFrame.head()


# In[45]:


items = pd.read_excel(path + "data/tcga/clinical_study_cohort.xlsx").columns[:-1]


# In[48]:


valid_clinical = validFrame.ix[:,items]
valid_clinical.reset_index(drop=True,inplace=True)


# In[47]:


# lncFrame.T.loc[valid_clinical.id].to_excel(path + "data/external/external_cohort.xlsx")


# In[50]:


set(valid_clinical.project_id)


# In[55]:


gender = {'male': 1,'female': 0} 
censor = {'Dead': 1,'Alive': 0} 


# In[61]:


for dataset in ['TARGET','CPTAC']:
    ex_clinical = valid_clinical[valid_clinical.project_id.str.find(dataset)==0]
    ex_clinical.reset_index(drop=True, inplace=True)
    lncSubFrame = lncFrame.T.loc[ex_clinical.id]
    lncSubFrame.loc[:,'gender'] = [gender[item] for item in ex_clinical.gender] 
    lncSubFrame.loc[:,'age_at_diagnosis']  = list(round(ex_clinical.age_at_diagnosis/365,0))
    # lncSubFrame.loc[:,'stage']  = [stage(item) for item in valid_clinical.tumor_stage]

    censors = [censor[status] for status in ex_clinical.vital_status]
    oss = []
    for index, row in ex_clinical.iterrows():
        if row.vital_status == 'Dead' and row.days_to_death==row.days_to_death:
    #         print(row.days_to_death)
            oss.append(row.days_to_death)
        else:
            oss.append(row.days_to_last_follow_up)
    lncSubFrame.loc[:, 'OS'] = oss
    lncSubFrame.loc[:, 'Censor'] = censors
    lncSubFrame.to_excel(path + 'data/external/'+ dataset.lower() + '_cohort.xlsx')  
    ex_clinical.to_excel(path + 'data/external/'+dataset.lower()+'_clinical.xlsx')

