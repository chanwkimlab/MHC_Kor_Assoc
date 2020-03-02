#!/usr/bin/env python
# coding: utf-8

# jupyter nbconvert 7_1_gcta_bivar.ipynb --to script
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


# In[3]:


result_uni=pd.read_csv(data_out_gcta_path+'result_uni.tsv',sep='\t',index_col=0)


# In[4]:


grm_path='data/genotype/4_merge/KCHIP_HLA_AA_SNP_1000G.grm'


# In[ ]:


i=int(sys.argv[1])


# In[5]:


if 'ipykernel' in sys.argv[0]:
    ipykernel=True
    i=0
    #phenotype_name='height'
else:
    ipykernel=False
    i=int(sys.argv[1])  
    
log.info(i)


# In[6]:


for j in range(i+1,len(result_uni)):

    phenotype_name1=result_uni.iloc[i].name
    phenotype_name2=result_uni.iloc[j].name

    pheno1=pd.read_csv(data_out_pheno_path+phenotype_name1+'.phe',sep='\t',header=None,names=['FID','IID','pheno'])
    pheno2=pd.read_csv(data_out_pheno_path+phenotype_name2+'.phe',sep='\t',header=None,names=['FID','IID','pheno'])
    pheno=pheno1.merge(right=pheno2,left_on=['FID','IID'],right_on=['FID','IID'])

    pheno_filter=pheno[(pheno['pheno_x']!=-9) & (pheno['pheno_y']!=-9)]
    pheno_filter.to_csv(data_out_gcta_path+phenotype_name1+'-'+phenotype_name2+'.phe',sep='\t',index=None,header=None)

    log.info("phenotype_name1: {}, phenotype_name2:{}".format(phenotype_name1,phenotype_name2))
    log.info('pheno1 mising {}'.format((pheno['pheno_x']!=-9).sum()))
    log.info('pheno2 mising {}'.format((pheno['pheno_y']!=-9).sum()))

    log.info("#########################################  Run GCTA  #########################################")
    #Run omnibus association test
    command='gcta64 --HEreg-bivar 1 2 --grm {} --pheno {} --out {} --thread-num 1'.format(grm_path,
                                                                                 data_out_gcta_path+phenotype_name1+'-'+phenotype_name2+'.phe',
                                                                                 data_out_gcta_path+phenotype_name1+'-'+phenotype_name2,                                         
                                                                                )
    log.info(command)
    stdout,stderr=run_subprocess(command,dry=False)
    log.info(stdout)
    log.error(stderr)          


# In[ ]:




