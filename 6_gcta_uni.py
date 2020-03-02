#!/usr/bin/env python
# coding: utf-8

# jupyter nbconvert 6_gcta_uni.ipynb --to script
# 
# for i in {00..101};do python 6_gcta_uni.py $i;done
# 
# for i in {00..10};do python 4_association.py $i;done
# 
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


# In[2]:


final_aa_path,final_plink_path,final_plink_aa_path


# # load plink, aa and check integrity

# In[27]:


plink_path=final_plink_path
plink_aa_path=final_plink_aa_path
grm_path='data/genotype/4_merge/KCHIP_SNP_1000G_merged.grm'
aa_path=final_aa_path


# In[31]:


#final_plink_aa_grm_path


# In[ ]:


plink=PyPlink(plink_path)
fam=plink.get_fam().astype({'fid':str,'iid':str}).rename(columns={'fid':'FID','iid':'IID'})
bim=plink.get_bim()


# In[5]:


plink_aa=PyPlink(plink_aa_path)
fam_aa=plink_aa.get_fam().astype({'fid':str,'iid':str}).rename(columns={'fid':'FID','iid':'IID'})
bim_aa=plink_aa.get_bim()


# In[6]:


assert (fam['IID']!=fam_aa['IID']).sum()==0


# In[7]:


f=open(aa_path,'r');aa_ind=f.readline().strip().split(' ')[2:];f.close()


# In[8]:


aa_ind_1=[aa_ind[i] for i in range(0,len(aa_ind),2)]
aa_ind_2=[aa_ind[i+1] for i in range(0,len(aa_ind),2)]


# In[9]:


assert (fam['IID']!=aa_ind_1).sum()==0
assert (fam['IID']!=aa_ind_2).sum()==0


# # load phenotype and check integrity

# In[10]:


phenotypes=pd.read_csv(pheno_all_file_path,sep='\t')
'  '.join(phenotypes.columns)


# In[11]:


phenotypes=phenotypes.set_index('ID').loc[fam.IID]
phenotypes.shape


# In[12]:


assert (phenotypes.index!=fam['IID']).sum()==0


# In[13]:


binary_continuous_traits=sorted(phenotypes.columns[~phenotypes.columns.str.contains('x_ray')])
binary_continuous_traits


# # parse parameter

# In[14]:


if 'ipykernel' in sys.argv[0]:
    ipykernel=True
    phenotype_name='0'
    #phenotype_name='height'
else:
    ipykernel=False
    phenotype_name=sys.argv[1]
if phenotype_name.isdigit():
    phenotype_name=int(phenotype_name)
    phenotype_name=binary_continuous_traits[phenotype_name]      


# In[18]:


pheno=pd.read_csv(data_out_pheno_path+phenotype_name+'.phe',sep='\t',header=None,names=['FID','IID','pheno'])
phenotype_type='binary' if len(pheno['pheno'].value_counts())<3 else 'continuous'


# In[ ]:





# In[20]:


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


# In[21]:


log.info("phenotype_name: {}, phenotype_type:{}".format(phenotype_name,phenotype_type))


# In[24]:


pheno[pheno['pheno']!=-9].to_csv(data_out_gcta_path+phenotype_name+'.phe',sep='\t',index=None,header=None)


# In[ ]:


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


# In[48]:


#pheno0=pd.read_csv(data_out_assoc_path+phenotype_list[0]+'/'+'phenotype.phe',header=None,sep='\t',names=['FID','IID','pheno0'])
#pheno1=pd.read_csv(data_out_assoc_path+phenotype_list[1]+'/'+'phenotype.phe',header=None,sep='\t',names=['FID','IID','pheno1']);pheno0['pheno1']=pheno1['pheno1']
#pheno2=pd.read_csv(data_out_assoc_path+phenotype_list[1]+'/'+'phenotype.phe',header=None,sep='\t',names=['FID','IID','pheno2']);pheno0['pheno2']=pheno2['pheno2']


# In[50]:


#pheno0.to_csv('temp.phe',header=None,index=None,sep='\t')


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




