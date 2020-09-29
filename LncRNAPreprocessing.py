
# coding: utf-8

# In[3]:


import numpy as np
import pandas as pd
import xlrd
import matplotlib
from matplotlib import pyplot as plt
import csv
import json
import os


# In[ ]:


path = './project/'


# In[2]:


# the path where FPKM files downloaded (need to download)
gene_raw_data_path = "./project/data/sets/"


# In[5]:


#the algorithm process FPKM file and generate the pkl file (geneframe_external.pkl)
pkl_path = './project/data/pkl/'


# ## Load lncRNA metadata and clinical data

# In[3]:


# corresponding full clinical data
clinical = pd.read_excel(path + 'data/raw/tcga_raw_clinical.xlsx')  # you can convert tcga_raw_clinical.tsv to tcga_raw_clinical.xlsx using excel
# the lncrnas that reported by TCGA for pan-cancer analysis
lncRnaSet = pd.read_csv(path + 'data/raw/gene_set_lincRNA.2018-09-20.tsv',delimiter = '\t')
#lncrna metada, used for download the FPKM files (with the metadata, you can use GDC transfer tool to redownload them)
lncRnaFiles = pd.read_csv(path + 'data/raw/tcga_gdc_manifest.2018-09-20.txt',delimiter = '\t')
#for mapping clinical data with gene expression data
with open(path + 'data/raw/tcga_files.2018-09-20.json') as json_file:
    jsonFile = json.load(json_file)


# In[4]:


lncRnaSet[:10] #let's see some lncrna ids


# In[5]:


lncRnaFiles[:10] #let's see some lncrna metadata


# In[6]:


clinical[:5] #let's see some clinical data


# ## Binding clinical data with lncRNA metadata file

# In[7]:


from collections import defaultdict
dictID = defaultdict()
for item in jsonFile:
    dictID[item['file_name']] = item['cases'][0]['case_id']


# In[8]:


caseids = [dictID[row['filename']] for index, row in lncRnaFiles.iterrows()]


# In[9]:


lncRnaFullFiles = lncRnaFiles.copy()
lncRnaFullFiles['case_id'] = caseids


# In[10]:


lncRnaFullFiles[:10]


# In[11]:


clinical_full = pd.merge(clinical, lncRnaFullFiles, how='left', on='case_id')


# In[12]:


clinical_full[:10]


# In[13]:


clinical_full = clinical_full.dropna(axis=1, how='all')


# In[14]:


clinical_full = clinical_full.dropna(axis=0, how='all')


# In[15]:


frame = clinical_full.drop(['md5','size','state'], axis=1)


# In[16]:


lst = frame[frame.vital_status == 'dead'].year_of_death.dropna()
np.max(list(lst))


# In[17]:


frame.project_id = frame.project_id.str.replace("TCGA-","")


# In[18]:


frame #let's see some clinical data, you can see they have corresponding FPKM file info


# ## Unzip FPKM files and merge them as a genome-width dataframe 

# In[19]:


problems = []
emptys = []
startIndex = 0
import gzip
from io import StringIO
curPercent = 0


# In[21]:


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
            f=gzip.open(gene_path + row['id'] + str("/") + row['filename'],'r')
            tmp = pd.read_csv(StringIO(f.read().decode("utf-8")), delimiter = '\t',header=None).set_index(0)
            tmp.rename(columns={1:row['id']}, inplace=True)
            remain_frames.append(tmp)
        except KeyboardInterrupt:
            remain_frame = pd.concat(remain_frames,axis=1)
            remain_frame.to_pickle('./gene/geneframe_id.pkl')
            print("progress: {}， curIndex {}".format(percent,index))
            print("interrupted!")
            break
        except:
            print(row['case_id'], " problem:", row['id'] + str("/") + row['filename'])
            problems.append(row['id'] + str("/") + row['filename'])
            pass
remain_frame = pd.concat(remain_frames,axis=1)
# The lncRNA dataframe, contains genome-width FPKM expression
remain_frame.to_pickle(pkl_path + 'geneframe_id.pkl')


# In[20]:


geneframe = pd.read_pickle(pkl_path + 'geneframe_id.pkl')


# ## Drop invalid patient data

# In[51]:


validFrame = frame.loc[frame.id.isin(list(geneframe.columns)),:]
    


# In[52]:


validFrame = validFrame.reset_index(drop=True)


# In[53]:


validFrame # we have 4235 patient data that have corresponding gene expressions, however, there are 4 of them without clinical info (the last four rows)


# In[ ]:


#save clinical info (entire cohort)
validFrame.tumor_stage.fillna('unknown',inplace=True)
validFrame.race.fillna('unknown', inplace=True)
validFrame.gender.fillna('unknown', inplace=True)
validFrame.tumor_stage = validFrame.tumor_stage.str.replace("stage ","")
validFrame.to_excel(path + 'data/tcga/clinical_entire_cohort.xlsx')


# ## We extract LncRNA frame from genome-width dataset

# In[48]:


lncFrame = geneframe.ix[(index for index, row in geneframe.iterrows() if row.name[:row.name.find('.')] in list(lncRnaSet.id)), :]


# ## Extract study cohort from the entire cohort

# In[54]:


#we identify the balanced prognosis cases, with the help of target prognosis, we can identify most significant prognostic lncrnas
above_cutoff = np.where(((validFrame.days_to_last_follow_up > 3.5*365) & (validFrame.vital_status=='alive')) | ((validFrame.days_to_death > 3.5*365) & (validFrame.vital_status!='alive')))


# In[55]:


above_cutoff = list(above_cutoff[0])


# In[57]:


below_cutoff = np.where((validFrame.days_to_death==validFrame.days_to_death) & (validFrame.days_to_death < 3.5*365) & (validFrame.vital_status=='dead') )


# In[58]:


below_cutoff = list(below_cutoff[0])


# In[45]:


candidate = above_cutoff
candidate.extend(below_cutoff)


# In[118]:



import numpy as np
import matplotlib.pyplot as plt
plt.subplots(1,1, figsize=(8,5))

N = 3
days_to_death_above_cut = len(np.where((validFrame.days_to_death > 3.5*365) & (validFrame.vital_status!='alive'))[0]) 
days_to_death_below_cut = len(np.where((validFrame.days_to_death==validFrame.days_to_death) & (validFrame.days_to_death < 3.5*365) & (validFrame.vital_status=='dead'))[0])
            
days_to_last_follow_up_above_cut = len(np.where((validFrame.days_to_last_follow_up > 3.5*365) & (validFrame.vital_status=='alive'))[0])
days_to_last_follow_up_below_cut = len(np.where((validFrame.days_to_last_follow_up <= 3.5*365) & (validFrame.vital_status=='alive'))[0])
                                       

days_to_death = ( days_to_death_below_cut, days_to_last_follow_up_above_cut )
days_to_last_follow_up = (days_to_death_above_cut, days_to_last_follow_up_below_cut)

ind = np.arange(N)    # the x locations for the groups
width = 0.35       # the width of the bars: can also be len(x) sequence

p1 = plt.bar(ind, (days_to_death_below_cut, days_to_last_follow_up_below_cut,days_to_death_below_cut), width)
p2 = plt.bar(ind, (days_to_death_above_cut, days_to_last_follow_up_above_cut,days_to_death_above_cut+days_to_last_follow_up_above_cut), width,
             bottom=(days_to_death_below_cut, days_to_last_follow_up_below_cut,days_to_death_below_cut))


plt.ylabel('Count')
plt.title('Clinical Survival Data and Cutoff Selection')
plt.xticks(ind, ('days_to_death', 'days_to_last_follow_up', 'Merged'))
# plt.yticks(np.arange(0, 81, 10))
plt.legend((p2[0], p1[0]), ('Above Cutoff (3.5 years)', 'Below Cutoff (3.5 years)'))
for r1, r2 in zip(p1, p2):
    h1 = r1.get_height()
    h2 = r2.get_height()
    plt.text(r1.get_x() + r1.get_width() / 2., h1 / 2., "%d" % h1, ha="center", va="bottom", color="white", fontsize=16, fontweight="bold")
    plt.text(r2.get_x() + r2.get_width() / 2., h1 + h2 / 2. -70, "%d" % h2, ha="center", va="bottom", color="white", fontsize=16, fontweight="bold")
    
plt.show()


# In[46]:


#the dataframe contains balanced prognosis studies (use the above cutoff)
candidate_frame = validFrame.ix[candidate]


# In[47]:


candidate_frame.loc[fiveyears,'threehalf'] = 1
candidate_frame.loc[lessfiveyears,'threehalf'] = 0


# In[48]:


candidate_frame


# In[49]:


#save study cohort
candidate_frame.reset_index(drop=True, inplace=True)
candidate_frame.loc[candidate_frame.project_id == 'GBM', 'tumor_stage'] = 'stage iv'
candidate_frame.loc[candidate_frame.project_id == 'LGG', 'tumor_stage'] = 'stage iic'
candidate_frame.tumor_stage.fillna('unknown', inplace=True)
candidate_frame.race.fillna('unknown', inplace=True)
candidate_frame.gender.fillna('unknown', inplace=True)

candidate_frame.tumor_stage = study_cohort.tumor_stage.str.replace("stage ","")

candidate_frame.to_excel(path + 'data/tcga/clinical_study_cohort.xlsx')

