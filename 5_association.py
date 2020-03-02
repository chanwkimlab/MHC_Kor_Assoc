#!/usr/bin/env python
# coding: utf-8

# jupyter nbconvert 5_association.ipynb --to script
# 
# cd /data/ch6845/MHC*;screen -S assoc;
# 
# for i in {0..10};do python 5_association.py $i 1;done
# 
# 
# 5*3=15
# 
# for i in {00..101};do python 5_association.py $i;done
# 
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6292650/

# In[ ]:


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

import statsmodels.api as sm

from basic_tools import *


# plink_KCHIP_HLA_AA_SNP_1000G_path=data_genotype_merge_path+'KCHIP_HLA_AA_SNP_1000G'
# plink_KCHIP_HLA_SNP_1000G_path=data_genotype_merge_path+'KCHIP_HLA_SNP_1000G' # for association
# plink_1000G_path=data_genotype_merge_path+'1000G'

# In[ ]:


#plink_KCHIP_HLA_SNP_1000G=PyPlink(plink_KCHIP_HLA_SNP_1000G_path)
#plink_KCHIP_HLA_SNP_1000G_fam=plink_KCHIP_HLA_SNP_1000G.get_fam().astype({'fid':str,'iid':str}).rename(columns={'fid':'FID','iid':'IID'})
#plink_KCHIP_HLA_SNP_1000G_bim=plink_KCHIP_HLA_SNP_1000G.get_bim()


# In[ ]:


plink_KCHIP_HLA_AA_SNP_1000G=PyPlink(plink_KCHIP_HLA_AA_SNP_1000G_path)
plink_KCHIP_HLA_AA_SNP_1000G_fam=plink_KCHIP_HLA_AA_SNP_1000G.get_fam().astype({'fid':str,'iid':str}).rename(columns={'fid':'FID','iid':'IID'})
plink_KCHIP_HLA_AA_SNP_1000G_bim=plink_KCHIP_HLA_AA_SNP_1000G.get_bim()


# In[ ]:


#plink_1000G=PyPlink(plink_1000G_path)
#plink_1000G_fam=plink_1000G.get_fam().astype({'fid':str,'iid':str}).rename(columns={'fid':'FID','iid':'IID'})
#plink_1000G_bim=plink_1000G.get_bim()
# plink_KCHIP_HLA_SNP_1000G_bim -> for parameter


# In[ ]:


phased_KCHIP_HLA_AA_path
phased_KCHIP_HLA_AA_SNP_path


# # parse parameter

# In[ ]:


phenotypes=pd.read_csv(pheno_all_file_path,sep='\t')
phenotypes=phenotypes.set_index('ID').loc[plink_KCHIP_HLA_AA_SNP_1000G_fam['IID']]

binary_continuous_traits=sorted(phenotypes.columns[~phenotypes.columns.str.contains('x_ray')])
#binary_continuous_traits

if 'ipykernel' in sys.argv[0]:
    ipykernel=True
    phenotype_name='ALP'
    step_idx=1
    #phenotype_name='height'
else:
    ipykernel=False
    phenotype_name=sys.argv[1]
    step_idx=int(sys.argv[2])
    
if phenotype_name.isdigit():
    phenotype_name=int(phenotype_name)
    phenotype_name=binary_continuous_traits[phenotype_name]      


# In[ ]:


pheno=pd.read_csv(data_out_pheno_path+phenotype_name+'.phe',sep='\t',names=['FID','IID','pheno'])
phenotype_type='binary' if len(pheno['pheno'].value_counts())<3 else 'continuous'


# In[ ]:


data_out_assoc_phenotype_path=data_out_assoc_path+phenotype_name+'/'
pathlib.Path(data_out_assoc_phenotype_path).mkdir(parents=True, exist_ok=True)


# In[ ]:


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


# In[ ]:


log.info_head=lambda x: log.info('-'*int((100-len(x))/2)+x+'-'*int((100-len(x))/2))


# In[ ]:


log.info_head("phenotype_name: {}, phenotype_type:{}".format(phenotype_name,phenotype_type))


# # Association

# plink_KCHIP_HLA_AA_SNP_1000G_HLAassign=plink_KCHIP_HLA_AA_SNP_1000G_bim[['pos']]
# HLA_names=np.unique([i[0].split('_')[1] for i in plink_KCHIP_HLA_AA_SNP_1000G_bim[plink_KCHIP_HLA_AA_SNP_1000G_bim.index.str.contains('HLA_')].index.str.split('*')])
# 
# for HLA_name in HLA_names:
#     plink_KCHIP_HLA_AA_SNP_1000G_HLAassign[HLA_name]=0
#     gene_select=plink_KCHIP_HLA_AA_SNP_1000G_HLAassign[plink_KCHIP_HLA_AA_SNP_1000G_HLAassign.index.str.contains('HLA_'+HLA_name)|plink_KCHIP_HLA_AA_SNP_1000G_HLAassign.index.str.contains('SNPS_'+HLA_name)|plink_KCHIP_HLA_AA_SNP_1000G_HLAassign.index.str.contains('AA_'+HLA_name)]#print(gene_select.sort_values('pos').iloc[0],gene_select.sort_values('pos').iloc[-1])
#     plink_KCHIP_HLA_AA_SNP_1000G_HLAassign[HLA_name][(plink_KCHIP_HLA_AA_SNP_1000G_HLAassign['pos']>=plink_KCHIP_HLA_AA_SNP_1000G_HLAassign['pos'].min())&(plink_KCHIP_HLA_AA_SNP_1000G_HLAassign['pos']<=plink_KCHIP_HLA_AA_SNP_1000G_HLAassign['pos'].max())]=1
# plink_KCHIP_HLA_AA_SNP_1000G_HLAassign

# In[ ]:


#for i in gene_assign[gene_assign.A==1].index:
#    print(i,gene_assign.loc[i])


# In[ ]:


plink_KCHIP_HLA_AA_SNP_1000G_HLAassign=plink_KCHIP_HLA_AA_SNP_1000G_bim[['pos']]
HLA_names=np.unique([i[0].split('_')[1] for i in plink_KCHIP_HLA_AA_SNP_1000G_bim[plink_KCHIP_HLA_AA_SNP_1000G_bim.index.str.contains('HLA_')].index.str.split('*')])

for HLA_name in HLA_names:
    plink_KCHIP_HLA_AA_SNP_1000G_HLAassign[HLA_name]=0
    gene_select=plink_KCHIP_HLA_AA_SNP_1000G_HLAassign[plink_KCHIP_HLA_AA_SNP_1000G_HLAassign.index.str.contains('HLA_'+HLA_name)|plink_KCHIP_HLA_AA_SNP_1000G_HLAassign.index.str.contains('SNPS_'+HLA_name)|plink_KCHIP_HLA_AA_SNP_1000G_HLAassign.index.str.contains('AA_'+HLA_name)]#print(gene_select.sort_values('pos').iloc[0],gene_select.sort_values('pos').iloc[-1])
    plink_KCHIP_HLA_AA_SNP_1000G_HLAassign[HLA_name][(plink_KCHIP_HLA_AA_SNP_1000G_HLAassign['pos']>=gene_select['pos'].min())&(plink_KCHIP_HLA_AA_SNP_1000G_HLAassign['pos']<=gene_select['pos'].max())]=1
plink_KCHIP_HLA_AA_SNP_1000G_HLAassign


# In[ ]:


plink_KCHIP_HLA_AA_SNP_1000G_r2_list_dict={}


# In[ ]:


#conditional_variant_list=['AA_A_-1','SNPS_A_28_30018337_exon1','6:30165273_C/T']


# In[ ]:


phased_KCHIP_HLA_AA_path


# debug=True
# #debug=True
# 
# if debug:
#     arg_split='\
# --assoc linear \
# --out sample_output \
# --bgl-phased /data/ch6845/MHC_phewas_testbench/data/genotype/4_merge/KCHIP_HLA_AA_SNP.bgl.phased \
# --bfile /data/ch6845/MHC_phewas_testbench/data/genotype/4_merge/KCHIP_HLA_SNP_1000G \
# --multialleic (?P<name>HLA_[0-9A-Z]*)\*(?P<allele>[0-9:]*) \
# --multialleic-always (?P<name>AA_[A-Z0-9]*_[\-0-9]*_[0-9]*_exon[0-9]*)_*(?P<allele>[A-Z]*) \
# --pheno /data/ch6845/MHC_phewas_testbench/data/out_pheno/FEV_predicted.phe \
# --covar /data/ch6845/MHC_phewas_testbench/data/out_assoc/FEV_predicted/step_01.plink.covar \
# --condition-list /data/ch6845/MHC_phewas_testbench/data/out_assoc/FEV_predicted/step_01.plink.cond\
# '.split(' ')
#     args=parser.parse_args(arg_split)
# else:
#     args=parser.parse_args()
#     
# if args.bfile is None and args.bgl_phased is None:
#     raise argparse.ArgumentTypeError("either --bfile or --bgl-phased parameter is needed")    

# In[ ]:


plink_KCHIP_HLA_SNP_1000G_path


# In[ ]:


conditional_list=[]
pd.Series(conditional_list).to_csv(data_out_assoc_phenotype_path+'step_{:02d}.cond'.format(step_idx),index=None,sep='\t',header=False)


# In[ ]:


covariate_df=pd.read_csv(PC_path,sep='\t').set_index('ID').loc[plink_KCHIP_HLA_AA_SNP_1000G_fam['IID']]
covariate_df['age']=phenotypes['age']
covariate_df['sex']=phenotypes['sex']-1
covariate_df['AS']=phenotypes['cohort'].replace(1,1).replace(2,0).replace(3,0)
covariate_df['NC']=phenotypes['cohort'].replace(1,0).replace(2,0).replace(3,1)


# In[ ]:


plink_KCHIP_HLA_AA_SNP_1000G_fam.iloc[:,:2].merge(right=covariate_df,left_on='IID',right_index=True).fillna(-9).to_csv(data_out_assoc_phenotype_path+'step_{:02d}.covar'.format(step_idx),index=None,sep='\t')


# In[ ]:


log.info("######################################### step {:02d} Association  #########################################".format(step_idx))
  


# In[ ]:


command="python Generic_Association_Tool/GAT.py --assoc {assoc_mode} --out {out} --bgl-phased {bgl_phased} --bfile {bfile} --pheno {pheno} --covar {covar} --condition-list {cond} --multialleic (?P<name>HLA_[0-9A-Z]*)\*(?P<allele>[0-9:]*) --multialleic-always (?P<name>AA_[A-Z0-9]*_[\-0-9]*_[0-9]*_exon[0-9]*)_*(?P<allele>[A-Z]*)".format(
assoc_mode='logistic' if phenotype_type=='binary' else 'linear',
out=data_out_assoc_phenotype_path+'step_{:02d}'.format(step_idx),
bgl_phased='data/genotype/4_merge/KCHIP_HLA_AA_SNP.bgl.phased',
bfile=plink_KCHIP_HLA_SNP_1000G_path,  
pheno=data_out_pheno_path+phenotype_name+'.phe',
covar=data_out_assoc_phenotype_path+'step_{:02d}.covar'.format(step_idx),   
cond=data_out_assoc_phenotype_path+'step_{:02d}.cond'.format(step_idx)  
)    
    
log.info(command)
stdout,stderr=run_subprocess(command,dry=False)
log.info(stdout)
log.error(stderr)     


# In[ ]:




