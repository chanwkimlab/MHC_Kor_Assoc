#!/usr/bin/env python
# coding: utf-8

# jupyter nbconvert 4_association-Copy1.ipynb --to script
# 
# for i in {00..102};do python 4_association-Copy1.py $i;done
# 
# cd /data/ch6845/MHC*;screen -S assoc;
# 
# for i in {0..10};do python 4_association.py $i;done
# 
# 
# Rscript OmnibusTest_BHv5_modified.R logistic data/out_assoc/allergic_disease/step_02.omnibus data/genotype/4_merge/KCHIP_HLA_SNP_1000G_merged.fam data/out_assoc/allergic_disease/step_02.aa data/out_assoc/allergic_disease/phenotype.pheomnibus pheno data/out_assoc/allergic_disease/step_02.omnibus.covar header0 NA
# 
# 5*3=15
# 
# for i in {00..06};do python 4_association.py $i;done
# for i in {07..13};do python 4_association.py $i;done
# for i in {14..20};do python 4_association.py $i;done
# 
# for i in {21..27};do python 4_association.py $i;done
# for i in {28..34};do python 4_association.py $i;done
# for i in {35..41};do python 4_association.py $i;done
# 
# for i in {42..48};do python 4_association.py $i;done
# for i in {49..55};do python 4_association.py $i;done
# for i in {56..62};do python 4_association.py $i;done
# 
# for i in {63..69};do python 4_association.py $i;done
# for i in {70..76};do python 4_association.py $i;done
# for i in {77..83};do python 4_association.py $i;done
# 
# for i in {84..90};do python 4_association.py $i;done
# for i in {91..97};do python 4_association.py $i;done
# for i in {98..101};do python 4_association.py $i;done
# 
# 
# 
# for i in {21..40};do python 4_association.py $i;done 
# for i in {41..60};do python 4_association.py $i;done 
# for i in {61..80};do python 4_association.py $i;done 
# for i in {81..101};do python 4_association.py $i;done
# 
# Rscript OmnibusTest_BHv5_modified-Copy1.R linear data/out_assoc/age/step_01.omnibus data/genotype/4_merge/KCHIP_HLA_SNP_1000G_merged.fam data/out_assoc/age/step_01.aa data/out_assoc/age/phenotype.phe data/out_assoc/age/step_01.omnibus.covar header0,header1,header2,header3 NA
# 
# Rscript OmnibusTest_BHv5_modified-Copy1.R linear data/out_assoc/grip_strength/step_01.omnibus data/genotype/4_merge/KCHIP_HLA_SNP_1000G_merged.fam data/out_assoc/grip_strength/step_01.aa data/out_assoc/grip_strength/phenotype.phe data/out_assoc/grip_strength/step_01.omnibus.covar header0,header1,header2,header3 NA
# 
# Rscript OmnibusTest_BHv5_modified-Copy1.R linear data/out_assoc/cervical_cancer/step_01.omnibus data/genotype/4_merge/KCHIP_HLA_SNP_1000G_merged.fam data/out_assoc/cervical_cancer/step_01.aa data/out_assoc/cervical_cancer/phenotype.phe data/out_assoc/cervical_cancer/step_01.omnibus.covar header0,header1,header2,header3 NA
# 
# Rscript OmnibusTest_BHv5_modified.R linear data/out_assoc/grip_strength/step_01.omnibus data/genotype/4_merge/KCHIP_HLA_SNP_1000G_merged.fam data/out_assoc/grip_strength/step_01.aa data/out_assoc/grip_strength/phenotype.phe data/out_assoc/grip_strength/step_01.omnibus.covar header0,header1,header2,header3 NA
# 
# Rscript OmnibusTest_BHv5_modified-Copy1.R logistic data/out_assoc/asthma/step_01.omnibus data/genotype/4_merge/KCHIP_HLA_SNP_1000G_merged.fam data/out_assoc/asthma/step_01.aa data/out_assoc/asthma/phenotype.phe data/out_assoc/asthma/step_01.omnibus.covar header0,header1,header2,header3 NA
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

# In[3]:


plink_path=final_plink_path
plink_aa_path=final_plink_aa_path
aa_path=final_aa_path


# In[4]:


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


binary_traits=phenotypes.columns[phenotypes.apply(lambda x: (not 'x_ray' in x.name) & (len(x.value_counts())<3),axis=0)]
binary_traits,len(binary_traits)


# In[14]:


continuous_traits=phenotypes.columns[phenotypes.apply(lambda x: (not 'x_ray' in x.name) & (len(x.value_counts())>=3),axis=0)]
continuous_traits


# In[29]:


binary_continuous_traits=sorted(binary_traits.union(continuous_traits))
binary_continuous_traits


# # parse parameter

# In[16]:


if 'ipykernel' in sys.argv[0]:
    ipykernel=True
    phenotype_name='protein_in_blood'
    #phenotype_name='height'
    
else:
    ipykernel=False
    phenotype_name=sys.argv[1]
    if phenotype_name.isdigit():
        phenotype_name=int(phenotype_name)
        phenotype_name=binary_continuous_traits[phenotype_name]

if phenotype_name in binary_traits:
    phenotype_type='binary'
elif phenotype_name in continuous_traits:
    phenotype_type='continuous'        


# In[17]:


data_out_assoc_phenotype_path=data_out_assoc_path+phenotype_name+'/'
pathlib.Path(data_out_assoc_phenotype_path).mkdir(parents=True, exist_ok=True)


# In[18]:


log = logging.getLogger('logger')
log.setLevel(logging.DEBUG)

log_file_name='test.log'
log_file_path=log_file_name
fileHandler = logging.FileHandler(log_file_path)
streamHandler = logging.StreamHandler()

formatter = logging.Formatter(' %(asctime)s [%(levelname)s] %(lineno)d > %(message)s')
fileHandler.setFormatter(formatter)
streamHandler.setFormatter(formatter)

log.addHandler(fileHandler)
log.addHandler(streamHandler)


# In[19]:


log.info("phenotype_name: {}, phenotype_type:{}".format(phenotype_name,phenotype_type))


# In[20]:


phenotype_define=np.full(len(phenotypes.index),np.nan)

if phenotype_type=='binary':
    for cohort in sorted(phenotypes['cohort'].unique()):
        log.info('------------------per cohort---------------------------')
        cohort_check=(phenotypes['cohort']==cohort)
        cohort_case_check=(phenotypes['cohort']==cohort)&(phenotypes[phenotype_name]==2)
        cohort_control_check=(phenotypes['cohort']==cohort)&(~(phenotypes[phenotype_name]==2))
        log.info('cohort: {}, {:5d}/{:5d} ({:.3f}%)'.format(cohort,cohort_case_check.sum(),cohort_check.sum(),100*cohort_case_check.sum()/cohort_check.sum()))

        if cohort_case_check.sum()>0:
            phenotype_define[cohort_case_check]=2
            phenotype_define[cohort_control_check]=1
        elif np.isnan(cohort):
            raise
            cohort_check_temp=(phenotypes['cohort'].isnull())
            phenotype_define[cohort_check_temp]=-9
            log.info('missing individuals founded: {}'.format(cohort_check_temp.sum()))
        else:
            log.info('cohort {} ignored. it may due to nonexistence of questionnaire'.format(cohort))

        log.info('Total case:'+str((phenotype_define==2).sum()))
    log.info("phenotype defined\n"+str(pd.Series(phenotype_define).value_counts()))
    
elif phenotype_type=='continuous':
    for cohort in sorted(phenotypes['cohort'].unique()):
        log.info('------------------per cohort---------------------------')
        cohort_check=(phenotypes['cohort']==cohort)
        cohort_notnull_check=(phenotypes['cohort']==cohort)&(~(phenotypes[phenotype_name].isnull()))
        #print(cohort_notnull_check,type(cohort_notnull_check))
        log.info('cohort: {}, {:5d}/{:5d} ({:.3f}%)'.format(cohort,cohort_notnull_check.sum(),cohort_check.sum(),100*cohort_notnull_check.sum()/cohort_check.sum()))

        if cohort_notnull_check.sum()>0:
            phenotype_define[cohort_notnull_check]=phenotypes[phenotype_name][cohort_notnull_check]
        elif np.isnan(cohort):
            raise
            cohort_check_temp=(phenotypes['cohort'].isnull())
            phenotype_define[cohort_check_temp]=-9
            log.info('missing individuals founded: {}'.format(cohort_check_temp.sum()))
        else:
            log.info('cohort {} ignored. it may due to nonexistence of questionnaire'.format(cohort))

        log.info('Total values: {}'.format((~np.isnan(phenotype_define)).sum()))
        
    log.info("median:{:.3f}, mean: {:.3f}, std: {:.3f}, max: {:.3f}, min: {:.3f}".format(pd.Series(phenotype_define).median(),
                                                         pd.Series(phenotype_define).mean(),
                                                         pd.Series(phenotype_define).std(),
                                                         pd.Series(phenotype_define).max(),
                                                         pd.Series(phenotype_define).min()
                                                        )
         )
    log.info(">mean+3std:{}, <mean-3std:{}".format((phenotype_define>pd.Series(phenotype_define).mean()+3*pd.Series(phenotype_define).std()).sum(),
                                                (phenotype_define<pd.Series(phenotype_define).mean()-3*pd.Series(phenotype_define).std()).sum()
                                               )
         )
    phenotype_define[phenotype_define>pd.Series(phenotype_define).mean()+3*pd.Series(phenotype_define).std()]=np.nan
    phenotype_define[phenotype_define<pd.Series(phenotype_define).mean()-3*pd.Series(phenotype_define).std()]=np.nan
    
    log.info('Total values: {}'.format((~np.isnan(phenotype_define)).sum()))                                                  
    pd.Series(phenotype_define).hist()


# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6708789/
#     
# age, gender, race, diabetes, hyperlipidemia, hypertension, and all significant alleles.
# 
# 
# 
# phenotypes[]
# union 
# * diabetes
# * hyperlipidemia
# * hypertension
# * allergic_disease
# * colon polyps
# * rheumatoid_arthritis
# 
# -> unhealthy individuals -> if overlap with case-> set as missing

# In[21]:


if phenotype_type=='binary' and phenotype_name!='sex':
    unhealthy_individuals=(phenotypes['diabetes']==2)|                            (phenotypes['hyperlipidemia']==2)|                            (phenotypes['hypertension']==2)|                            (phenotypes['allergic_disease']==2)|                            (phenotypes['colon_polyps']==2)|                            (phenotypes['rheumatoid_arthritis']==2)
    log.info("unhealthy individuals: {}".format(unhealthy_individuals.sum()))
    log.info("unhealthy individuals among control removed: {}".format(((phenotype_define==1) & (unhealthy_individuals)).sum()))
    phenotype_define[(phenotype_define==1) & (unhealthy_individuals)]=np.nan
    ## change to np.nan and test!!
    log.info("phenotype defined\n"+str(pd.Series(phenotype_define).value_counts()))


# In[22]:


if phenotype_type=='binary':
    if phenotype_name=='breast_cancer' or phenotype_name=='cervical_cancer':
        log.info('exclude men: {}'.format(((~np.isnan(phenotype_define))&(phenotypes['sex']==1)).sum()))
        phenotype_define[(~np.isnan(phenotype_define))&(phenotypes['sex']==1)]=np.nan
        log.info("phenotype defined\n"+str(pd.Series(phenotype_define).value_counts()))
    elif phenotype_name=='prostate_cancer':
        log.info('exclude women: {}'.format(((~np.isnan(phenotype_define))&(phenotypes['sex']==2)).sum()))
        phenotype_define[(~np.isnan(phenotype_define))&(phenotypes['sex']==2)]=np.nan
        log.info("phenotype defined\n"+str(pd.Series(phenotype_define).value_counts()))


# In[23]:


phenotype_define_df=pd.DataFrame(phenotype_define,index=phenotypes.index)

phenotype_define_df=phenotype_define_df.loc[fam['IID']].fillna(-9)

if phenotype_name in binary_traits:
    phenotype_define_df=phenotype_define_df.astype(int)
    
phenotype_define_df_noindex=phenotype_define_df.reset_index().rename(columns={0:'pheno'})

phenotype_define_df_noindex[[phenotype_define_df_noindex.columns[0],phenotype_define_df_noindex.columns[0],phenotype_define_df_noindex.columns[1]]].to_csv(data_out_assoc_phenotype_path+'phenotype.phe',index=None,header=None,sep='\t')
#phenotype_define_df_noindex[[phenotype_define_df_noindex.columns[0],phenotype_define_df_noindex.columns[0],phenotype_define_df_noindex.columns[1]]].to_csv(data_out_assoc_phenotype_path+'phenotype.pheomnibus',index=None,sep='\t')


# In[24]:


phenotype_define_df_noindex.shape


# In[25]:


assert (phenotype_define_df_noindex['pheno']==np.nan).sum()==0
if phenotype_type=='binary':
    log.info(phenotype_define_df_noindex['pheno'].value_counts())
elif phenotype_type=='continuous':
    log.info(phenotype_define_df_noindex['pheno'].value_counts().iloc[:5])
    log.info('...')
    log.info(phenotype_define_df_noindex['pheno'].value_counts().iloc[-5:])
else:
    raise

