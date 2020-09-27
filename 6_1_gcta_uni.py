#!/usr/bin/env python
# coding: utf-8

# jupyter nbconvert 6_1_gcta_uni.ipynb --to script
# 
# for i in {00..101};do python 6_1_gcta_uni.py $i;done
# 
# for i in {00..10};do python 4_association.py $i;done
# 
# 
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6292650/

# In[1]:


import os
import sys

import logging
import re
import pathlib

import datetime

import numpy as np
import pandas as pd
from scipy.stats import pearsonr

import matplotlib.pyplot as plt

from pyplink import PyPlink

from basic_tools import *


# # load plink, aa and check integrity

# In[2]:


plink_KCHIP_HLA_AA_SNP_1000G=PyPlink(plink_KCHIP_HLA_AA_SNP_1000G_path)
plink_KCHIP_HLA_AA_SNP_1000G_fam=plink_KCHIP_HLA_AA_SNP_1000G.get_fam().astype({'fid':str,'iid':str}).rename(columns={'fid':'FID','iid':'IID'})
plink_KCHIP_HLA_AA_SNP_1000G_bim=plink_KCHIP_HLA_AA_SNP_1000G.get_bim()


# In[3]:


grm_path='data/genotype/4_merge/KCHIP_HLA_AA_SNP_1000G.grm'


# In[4]:


#final_plink_aa_grm_path


# # load phenotype and check integrity

# In[5]:


phenotypes=pd.read_csv(pheno_all_file_path,sep='\t')
'  '.join(phenotypes.columns)


# In[6]:


phenotypes=phenotypes.set_index('ID').loc[plink_KCHIP_HLA_AA_SNP_1000G_fam['IID']]
phenotypes.shape


# In[7]:


#np.all(phenotypes['ALP'].isnull())


# In[8]:


assert (phenotypes.index!=plink_KCHIP_HLA_AA_SNP_1000G_fam['IID']).sum()==0


# In[9]:


binary_continuous_traits=phenotypes.columns.difference(['age','sex','cohort','diabetes'])
binary_continuous_traits,len(binary_continuous_traits)


# # parse parameter

# In[10]:


if 'ipykernel' in sys.argv[0]:
    ipykernel=True
    phenotype_name='diabetes'
    #phenotype_name='height'
else:
    ipykernel=False
    phenotype_name=sys.argv[1]
    
if phenotype_name.isdigit():
    phenotype_name=int(phenotype_name)
    phenotype_name=binary_continuous_traits[phenotype_name]      


# In[10]:


pheno=pd.read_csv(data_out_pheno_path+phenotype_name+'.phe',sep='\t',header=None,names=['FID','IID','pheno'])
phenotype_type='binary' if len(pheno['pheno'][pheno['pheno']!=-9].value_counts())<3 else 'continuous'
phenotype_type


# In[ ]:





# In[16]:


for a,i in enumerate(binary_continuous_traits):
    if not os.path.exists(data_out_gcta_path+i+'.HEreg'):
        print(a)
        #print(a)
        print(i,os.path.exists(data_out_gcta_path+i+'.HEreg'))


# In[14]:


log = logging.getLogger('logger')
log.setLevel(logging.DEBUG)

log_file_name=datetime.datetime.now().strftime('%Y%m%d_%H%M%S')+'.log'
log_file_path=data_out_gcta_path+log_file_name
fileHandler = logging.FileHandler(log_file_path)
streamHandler = logging.StreamHandler()

formatter = logging.Formatter(' %(asctime)s [%(levelname)s] %(lineno)d > %(message)s')
fileHandler.setFormatter(formatter)
streamHandler.setFormatter(formatter)

log.addHandler(fileHandler)
log.addHandler(streamHandler)


# In[15]:


log.info("phenotype_name: {}, phenotype_type:{}".format(phenotype_name,phenotype_type))


# In[16]:


pheno[pheno['pheno']!=-9].to_csv(data_out_gcta_path+phenotype_name+'.phe',sep='\t',index=None,header=None)


# In[29]:


if phenotype_type=='binary':
    with open(data_out_pheno_path+phenotype_name+'.phe'+'.prev','r') as f:
        prev=float(f.read())
    print('prev',prev)


# In[30]:


log.info("#########################################  Run GCTA  #########################################")
#Run omnibus association test
command='gcta64 --HEreg --grm {} --pheno {} --out {} --thread-num 40'.format(grm_path,
                                                                             data_out_gcta_path+phenotype_name+'.phe',
                                                                             data_out_gcta_path+phenotype_name,                                         
                                                                            )
log.info(command)
stdout,stderr=run_subprocess(command,dry=False)
log.info(stdout)
log.error(stderr)    


# In[20]:





# In[48]:


#pheno0=pd.read_csv(data_out_assoc_path+phenotype_list[0]+'/'+'phenotype.phe',header=None,sep='\t',names=['FID','IID','pheno0'])
#pheno1=pd.read_csv(data_out_assoc_path+phenotype_list[1]+'/'+'phenotype.phe',header=None,sep='\t',names=['FID','IID','pheno1']);pheno0['pheno1']=pheno1['pheno1']
#pheno2=pd.read_csv(data_out_assoc_path+phenotype_list[1]+'/'+'phenotype.phe',header=None,sep='\t',names=['FID','IID','pheno2']);pheno0['pheno2']=pheno2['pheno2']


# command='gcta64 --HEreg-bivar 1 2 --grm {} --pheno {} --out {} --thread-num 40'.format(final_plink_aa_grm_path,
#                                                                             'temp.phe',
#                                                                              'testout2'  
#                                                                             )
# command

# gcta64 --HEreg-bivar 1 2 --grm data/genotype/4_merge/grm --pheno temp2.phe --out data/out_assoc/height/HE2 --thread-num 40

# (pheno0['pheno0']!=-9).sum(),(pheno0['pheno1']!=-9).sum(),(pheno0['pheno2']!=-9).sum()

# In[29]:


#!cat data/out_gcta/ALP.phe


# In[ ]:




