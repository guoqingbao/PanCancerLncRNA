
# coding: utf-8

# In[1]:


import numpy as np
from PIL import Image
import sqlite3
get_ipython().run_line_magic('matplotlib', 'inline')
import matplotlib.pyplot as plt
import pandas as pd
import os
from plotly import __version__
from plotly.offline import download_plotlyjs, init_notebook_mode, iplot
from plotly.graph_objs import Scatter
get_ipython().run_line_magic('matplotlib', 'inline')
init_notebook_mode(connected=True)
print("plotly version:", __version__)


# ## Load clinical and LncRNA data

# In[2]:


#You may change the path
path = './project/'


# In[3]:


frame = pd.read_excel(path + 'data/clinical_entire_cohort.xlsx')
frame.project_id = frame.project_id.str.replace("TCGA-","")
labels = list(set(frame.project_id))


# In[4]:


study_cohort_clinical = pd.read_excel(path + "data/clinical_study_cohort.xlsx", index_col=0)
study_cohort = pd.read_excel(path + "data/study_cohort.xlsx", index_col=0)


# In[5]:


lncStudyCohort = study_cohort.copy()
lncStudyCohort.ix[:,'threehalf'] = list(study_cohort_clinical.threehalf)
lncStudyCohort.ix[:,'project_id'] = list(study_cohort_clinical.project_id)


# In[6]:


lncStudyCohort


# ## Mapping lncRNA expression levels with cancers and prognostic outcomes

# In[7]:


Amean = np.median(lncStudyCohort["ENSG00000206567.8"])
Bmean = np.median(lncStudyCohort["ENSG00000259641.4"])
Cmean = np.median(lncStudyCohort["ENSG00000257989.1"])
Dmean = np.median(lncStudyCohort["ENSG00000187185.4"])
Emean = np.median(lncStudyCohort["ENSG00000218510.5"])


# In[8]:


lncStudyCohort.loc[:,"A"] = "L"
lncStudyCohort.loc[lncStudyCohort["ENSG00000206567.8"]>Amean,"A"] = "H"

lncStudyCohort.loc[:,"B"] = "L"
lncStudyCohort.loc[lncStudyCohort["ENSG00000259641.4"]>Bmean,"B"] = "H"

lncStudyCohort.loc[:,"C"] = "L"
lncStudyCohort.loc[lncStudyCohort["ENSG00000257989.1"]>Cmean,"C"] = "H"


lncStudyCohort.loc[:,"D"] = "L"
lncStudyCohort.loc[lncStudyCohort["ENSG00000187185.4"]>Dmean,"D"] = "H"

lncStudyCohort.loc[:,"E"] = "L"
lncStudyCohort.loc[lncStudyCohort["ENSG00000218510.5"]>Emean,"E"] = "H"


# In[9]:


combination = list(lncStudyCohort.A + lncStudyCohort.B + lncStudyCohort.C + lncStudyCohort.D + lncStudyCohort.E)
lncStudyCohort.loc[:,"combination"] = combination
uniqCom = set(combination)


# In[10]:


labels.extend(list(uniqCom))


# In[11]:


labels.append("> 3.5 Years")
labels.append("≤ 3.5 Years")


# In[12]:


labels


# In[13]:


len(labels)


# In[14]:


import seaborn as sns
colors = sns.color_palette('hls', len(labels)).as_hex()
from random import shuffle
shuffle(colors)
sns.palplot(colors)
labelFrame = pd.DataFrame(np.array(colors), index=labels, columns=['cl'])


# In[15]:


labelFrame.describe()


# In[16]:


color_frame = pd.read_excel(path + "data/color_frame.xlsx", index_col=0)
distribution = pd.read_excel(path + "data/Distribution.xlsx", index_col=0)


# In[17]:


colors = sns.color_palette('hls', 34).as_hex()
from random import shuffle
shuffle(colors)
sns.palplot(colors)


# In[18]:


color_frame


# In[19]:


labelFrame


# In[20]:


index = 0
index1 = 0
for item in list(labelFrame.index):
    if len(color_frame[color_frame.Gene==item].Color)>0:
        strr = list(color_frame[color_frame.Gene==item].Color)[0]
        strr = strr.replace("[","")
        strr = strr.replace("]","")
        labelFrame.cl[index]= "rgba"+str(tuple([float(it) for it in strr.strip().split(" ") if it!=""]))
        index = index +1
    else:
        labelFrame.cl[index + index1] = colors[index1]
        index1 = index1 +1


# In[21]:


len(labelFrame)


# In[22]:


from collections import defaultdict
index = 0
indexList = []
for item in labelFrame.iterrows():
    sub = lncStudyCohort[lncStudyCohort.project_id == item[0]]
    
    targets = defaultdict(int)

    for row in sub.iterrows():
        pos = np.where(np.array(labels)==row[1][-1])[0][0]
        targets[pos] = targets[pos] + 1
    for (key, value) in targets.items():   
        indexList.append([index, key , value, labelFrame.ix[index,"cl"]])
    index = index + 1
    if index >= 33:
        break


# In[23]:


indexList


# In[24]:


sub = lncStudyCohort[lncStudyCohort.age_at_diagnosis >45*356]
for row in sub.iterrows():
    print(row[1][745])
    break


# In[25]:


sourceAbove = defaultdict(int)
sub = lncStudyCohort[lncStudyCohort.threehalf>0]
for row in sub.iterrows():
    pos = np.where(np.array(labels)==row[1][-1])[0][0]
    sourceAbove[pos] = sourceAbove[pos] + 1  
    
for (key, value) in sourceAbove.items():   
    indexList.append([key, len(labelFrame)-2 , value, labelFrame.ix[key,"cl"]])
    
sourceBelow = defaultdict(int)
sub = lncStudyCohort[lncStudyCohort.threehalf==0]
for row in sub.iterrows():
    pos = np.where(np.array(labels)==row[1][-1])[0][0]
    sourceBelow[pos] = sourceBelow[pos] + 1  
    
for (key, value) in sourceBelow.items():   
    indexList.append([key, len(labelFrame)-1 , value, labelFrame.ix[key,"cl"]])


# ## Plot the sankey diagram

# In[26]:


# indexList
def hextorgba(hex):
    if hex.find("rgba")!=-1:
        return hex.replace("1.0)", "0.6)")
    h = hex.lstrip('#')
    tp = tuple(int(h[i:i+2], 16) for i in (0, 2 ,4))
#     print(tp[0]

    s = 'rgba({},{},{},{})'.format(tp[0] ,tp[1],tp[2],0.7)
    return s


# In[27]:


hextorgba('rgba(0.44, 0.0, 0.46, 1.0)')


# In[28]:


indexList[0]


# In[29]:


# the sankey diagram script (used by plotly)
data_trace = dict(
    type='sankey',
    domain = dict(
      x =  [0,1],
      y =  [0,1]
    ),
    orientation = "h",
    valueformat = ".0f",
    node = dict(
      pad = 10,
      thickness = 30,
      line = dict(
        color = "black",
        width = 0
      ),
      label =  labelFrame.index,
      color = labelFrame['cl']
        
    ),
    link = dict(
      source = np.array(indexList)[:, 0],
      target = np.array(indexList)[:, 1],
      value = np.array(indexList)[:, 2],
      color = [hextorgba(item) for item in np.array(indexList)[:, 3]],
  )
)

layout =  dict(

    height = 720,
    font = dict(
      size = 11,
      family = "Times New Roman",
    ),    

)


# In[30]:


#make sure you have installed plotly and orca
import plotly
import plotly.io as pio
plotly.io.orca.config.executable = '/home/gbao5100/anaconda3/orca.sh'
plotly.io.orca.config.save() 


# In[31]:


get_ipython().run_cell_magic('time', '', "fig = dict(data=[data_trace], layout=layout)\nimg_bytes = pio.to_image(fig, format='pdf', width=1000, height=750)\npio.write_image(fig, path + 'results/sankey_diagram1.pdf', width=1000, height=750)\nimg_bytes = pio.to_image(fig, format='png', width=1000, height=750)\n\nfrom IPython.display import SVG,Image, display\ndisplay(Image(img_bytes))")

