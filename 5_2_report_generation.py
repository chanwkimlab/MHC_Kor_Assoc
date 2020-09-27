#!/usr/bin/env python
# coding: utf-8

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
import seaborn as sns

import statsmodels.api as sm

from basic_tools import *

"""
jupyter nbconvert 5_2_report_generation.ipynb --to script

for i in {00..20};do python 5_2_report_generation.py $i;done
for i in {21..40};do python 5_2_report_generation.py $i;done
for i in {41..60};do python 5_2_report_generation.py $i;done
for i in {61..80};do python 5_2_report_generation.py $i;done
for i in {81..98};do python 5_2_report_generation.py $i;done
"""


# In[2]:


plink_KCHIP_HLA_AA_SNP_1000G=PyPlink(plink_KCHIP_HLA_AA_SNP_1000G_path)
plink_KCHIP_HLA_AA_SNP_1000G_fam=plink_KCHIP_HLA_AA_SNP_1000G.get_fam().astype({'fid':str,'iid':str}).rename(columns={'fid':'FID','iid':'IID'})
plink_KCHIP_HLA_AA_SNP_1000G_bim=plink_KCHIP_HLA_AA_SNP_1000G.get_bim()


# In[3]:


phenotypes=pd.read_csv(pheno_all_file_path,sep='\t')
phenotypes=phenotypes.set_index('ID').loc[plink_KCHIP_HLA_AA_SNP_1000G_fam['IID']]

binary_continuous_traits=phenotypes.columns.difference(['age','sex','cohort','diabetes'])
print(len(binary_continuous_traits))

if 'ipykernel' in sys.argv[0]:
    ipykernel=True
    phenotype_name='thyroid_disease'
    #phenotype_name='height'
else:
    ipykernel=False
    phenotype_name=sys.argv[1]
    
if phenotype_name.isdigit():
    phenotype_name=int(phenotype_name)
    phenotype_name=binary_continuous_traits[phenotype_name]      


# In[4]:


data_out_assoc_phenotype_path=data_out_assoc_path+phenotype_name+'/'
pathlib.Path(data_out_assoc_phenotype_path).mkdir(parents=True, exist_ok=True)


# In[5]:


pheno=pd.read_csv(data_out_pheno_path+phenotype_name+'.phe',sep='\t',names=['FID','IID','pheno'])
phenotype_type='binary' if len(pheno['pheno'][pheno['pheno']!=-9].value_counts())<3 else 'continuous'
phenotype_type


# In[6]:


log = logging.getLogger('logger')
log.setLevel(logging.DEBUG)

log_file_name=datetime.datetime.now().strftime('%Y%m%d_%H%M%S')+'.log'
log_file_path=data_out_assoc_phenotype_path+log_file_name
fileHandler = logging.FileHandler(log_file_path)
streamHandler = logging.StreamHandler()

formatter = logging.Formatter('%(message)s')
fileHandler.setFormatter(formatter)
streamHandler.setFormatter(formatter)

log.addHandler(fileHandler)
log.addHandler(streamHandler)


# In[7]:


log.info_head=lambda x: log.info('-'*int((100-len(x))/2)+x+'-'*int((100-len(x))/2))


# In[8]:


def get_a1_freq_case_control(row):
    if 'bialle' in row['note']:
        dosage=plink_KCHIP_HLA_AA_SNP_1000G.get_geno_marker(row['marker_name']).astype(float)
        dosage[dosage==-1]=np.nan
        if plink_KCHIP_HLA_AA_SNP_1000G_bim.loc[row['marker_name']]['a1']==row['A1'] and plink_KCHIP_HLA_AA_SNP_1000G_bim.loc[row['marker_name']]['a2']==row['A2']:
            pass
        elif plink_KCHIP_HLA_AA_SNP_1000G_bim.loc[row['marker_name']]['a1']==row['A2'] and plink_KCHIP_HLA_AA_SNP_1000G_bim.loc[row['marker_name']]['a2']==row['A1']:
            dosage=2-dosage

        dosage_control=dosage[pheno['pheno']==1]
        dosage_case=dosage[pheno['pheno']==2]

        a1_freq_case=(2*(dosage_case==2).sum()+1*(dosage_case==1).sum())/(2*(dosage_case!=-1).sum())
        a1_freq_control=(2*(dosage_control==2).sum()+1*(dosage_control==1).sum())/(2*(dosage_control!=-1).sum())

        return a1_freq_case,a1_freq_control
    else:
        return np.nan,np.nan


# In[9]:


def name_to_pos(name):
    if name[:3]=='AA_':
        return int(name.split('_')[3])
    elif name[:4]=='HLA_':
        return plink_KCHIP_HLA_AA_SNP_1000G_bim.loc[plink_KCHIP_HLA_AA_SNP_1000G_bim.index.str.contains('HLA_'+name.split('*')[0].split('_')[1]+'*',regex=False)].iloc[0]['pos']
    else:
        raise


# In[ ]:


for step_idx_sub in range(1,100):
    log.info_head("phenotype_name: {}, phenotype_type:{} , Step : {} ".format(phenotype_name,phenotype_type,step_idx_sub))
    
    if os.path.exists(data_out_assoc_phenotype_path+'step_{:02d}.cond.stop'.format(step_idx_sub)):
        log.info("step_{}.cond.stop->stops ".format(step_idx_sub))
        sys.exit()
    
    result_merge=pd.read_csv(data_out_assoc_phenotype_path+'step_{:02d}.result.tsv'.format(step_idx_sub),sep='\t')                                        
    result_merge['note']=result_merge['note'].str.replace('phased bialleic','phased biallelic')
    
    result_merge['POS'][result_merge['POS'].isnull()]=result_merge[result_merge['POS'].isnull()]['marker_name'].apply(name_to_pos)
    assert len(result_merge[result_merge['POS'].isnull()])==0

    result_merge['phenotype_name']=phenotype_name
    result_merge['step']=step_idx_sub    

    with open(data_out_assoc_phenotype_path+'step_{:02d}.cond'.format(step_idx_sub)) as f:
        result_merge['condition']=','.join(f.read().split())      
        
        
    if phenotype_type=='binary':
        assert set(pheno['pheno'].unique())=={-9,1,2}
        if step_idx_sub==1:
            a1_freq_case_control=result_merge.apply(get_a1_freq_case_control,axis=1)
        else:
            assert np.all(result_merge_save['phenotype_name']==phenotype_name)
            assert np.all(result_merge_save['marker_name']==result_merge['marker_name'])
            assert np.all(result_merge_save['A1'].fillna(0)==result_merge['A1'].fillna(0))
            
        result_merge['A1_freq_case']=[i[0] for i in a1_freq_case_control]
        result_merge['A1_freq_control']=[i[1] for i in a1_freq_case_control]

        result_merge['samples(case/control)']='{}/{}'.format((pheno['pheno']==2).sum(),(pheno['pheno']==1).sum())
    else:
        result_merge['A1_freq_case'],result_merge['A1_freq_control']=np.nan,np.nan

        result_merge['samples(case/control)']='{}'.format((pheno['pheno']!=-9).sum())
        
        
    result_merge_save=result_merge[['phenotype_name','samples(case/control)','step','condition','marker_name','note','term','POS','A1','A2','A1_freq_case','A1_freq_control','multi_allele','Z','coef','std','chisq','df','nobs','P']]
    result_merge_save.sort_values('POS').to_csv(data_out_assoc_phenotype_path+'step_{:02d}.merge.result.tsv'.format(step_idx_sub),sep='\t',index=None)        

