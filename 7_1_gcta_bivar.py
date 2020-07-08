#!/usr/bin/env python
# coding: utf-8

# jupyter nbconvert 7_1_gcta_bivar.ipynb --to script
# 
# 
# for i in {00..20};do python 7_1_gcta_bivar.py $i;done
# for i in {21..40};do python 7_1_gcta_bivar.py $i;done
# for i in {41..60};do python 7_1_gcta_bivar.py $i;done
# for i in {61..80};do python 7_1_gcta_bivar.py $i;done
# for i in {81..101};do python 7_1_gcta_bivar.py $i;done
# 
# 
# 
# python 7_1_gcta_bivar.py $WINDOW
# 
# for i in {00..101};do python 6_gcta_uni.py $i;done
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


# In[2]:


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


# In[ ]:





# In[ ]:





# In[ ]:





# # load plink, aa and check integrity

# In[3]:


plink_KCHIP_HLA_AA_SNP_1000G=PyPlink(plink_KCHIP_HLA_AA_SNP_1000G_path)
plink_KCHIP_HLA_AA_SNP_1000G_fam=plink_KCHIP_HLA_AA_SNP_1000G.get_fam().astype({'fid':str,'iid':str}).rename(columns={'fid':'FID','iid':'IID'})
plink_KCHIP_HLA_AA_SNP_1000G_bim=plink_KCHIP_HLA_AA_SNP_1000G.get_bim()


# In[4]:


grm_path='data/genotype/4_merge/KCHIP_HLA_AA_SNP_1000G.grm'


# In[5]:


#final_plink_aa_grm_path


# # load phenotype and check integrity

# In[6]:


phenotypes=pd.read_csv(pheno_all_file_path,sep='\t')
'  '.join(phenotypes.columns)


# In[7]:


phenotypes=phenotypes.set_index('ID').loc[plink_KCHIP_HLA_AA_SNP_1000G_fam['IID']]
phenotypes.shape


# In[8]:


#np.all(phenotypes['ALP'].isnull())


# In[9]:


assert (phenotypes.index!=plink_KCHIP_HLA_AA_SNP_1000G_fam['IID']).sum()==0


# In[10]:


binary_continuous_traits=phenotypes.columns.difference(['age','sex','cohort'])
binary_continuous_traits,len(binary_continuous_traits)


# # parse parameter

# In[11]:


if 'ipykernel' in sys.argv[0]:
    ipykernel=True
    i=0
    #phenotype_name='height'
else:
    ipykernel=False
    i=int(sys.argv[1])  
    
log.info(i)


# In[ ]:





# In[14]:


#i=0


# In[15]:


#(pheno['pheno_x']!=-9).sum(),(pheno['pheno_y']!=-9).sum(),((pheno['pheno_x']!=-9) & (pheno['pheno_y']!=-9)).sum()


# In[16]:


#pd.read_csv(data_out_gcta_path+phenotype_name1+'-'+phenotype_name2+'.phe')


# In[15]:


for j in range(i+1,len(binary_continuous_traits)):
    phenotype_name1=binary_continuous_traits[i]
    phenotype_name2=binary_continuous_traits[j]
    
    if os.path.exists(data_out_gcta_path+phenotype_name1+'-'+phenotype_name2+'.HEreg') or os.path.exists(data_out_gcta_path+phenotype_name2+'-'+phenotype_name1+'.HEreg'):
        continue



    pheno1=pd.read_csv(data_out_pheno_path+phenotype_name1+'.phe',sep='\t',header=None,names=['FID','IID','pheno'])
    pheno2=pd.read_csv(data_out_pheno_path+phenotype_name2+'.phe',sep='\t',header=None,names=['FID','IID','pheno'])
    pheno=pheno1.merge(right=pheno2,left_on=['FID','IID'],right_on=['FID','IID'])

    pheno_filter=pheno[(pheno['pheno_x']!=-9) & (pheno['pheno_y']!=-9)]
    pheno_filter.to_csv(data_out_gcta_path+phenotype_name1+'-'+phenotype_name2+'.phe',sep='\t',index=None,header=None)

    log.info("phenotype_name1: {}, phenotype_name2:{}".format(phenotype_name1,phenotype_name2))
    log.info('pheno1 mising {}'.format((pheno['pheno_x']!=-9).sum()))
    log.info('pheno2 mising {}'.format((pheno['pheno_y']!=-9).sum()))

    if os.path.exists(data_out_gcta_path+phenotype_name1+'-'+phenotype_name2+'.HEreg'):
        log.info(data_out_gcta_path+phenotype_name1+'-'+phenotype_name2+'.HEreg exists')
        break
    
    log.info("#########################################  Run GCTA  #########################################")
    #Run omnibus association test
    command='gcta64 --HEreg-bivar 1 2 --grm {} --pheno {} --out {} --thread-num 10'.format(grm_path,
                                                                                 data_out_gcta_path+phenotype_name1+'-'+phenotype_name2+'.phe',
                                                                                 data_out_gcta_path+phenotype_name1+'-'+phenotype_name2,                                         
                                                                                )
    log.info(command)
    stdout,stderr=run_subprocess(command,dry=False)
    log.info(stdout)
    log.error(stderr)          

