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

import statsmodels.api as sm

from basic_tools import *

"""

jupyter nbconvert 5_1_association.ipynb --to script

python 5_1_association.py height 1 0,1,2




for i in {00..98};do python 5_1_association.py $i 0 0,1,2;done


for i in {00..10};do python 5_1_association.py $i 1 0;done
for i in {11..20};do python 5_1_association.py $i 1 0;done
for i in {21..30};do python 5_1_association.py $i 1 0;done
for i in {31..40};do python 5_1_association.py $i 1 0;done
for i in {41..50};do python 5_1_association.py $i 1 0;done
for i in {51..60};do python 5_1_association.py $i 1 0;done
for i in {61..70};do python 5_1_association.py $i 1 0;done
for i in {71..80};do python 5_1_association.py $i 1 0;done
for i in {81..90};do python 5_1_association.py $i 1 0;done
for i in {91..97};do python 5_1_association.py $i 1 0;done
"""


# In[2]:


plink_KCHIP_HLA_AA_SNP_1000G=PyPlink(plink_KCHIP_HLA_AA_SNP_1000G_path)
plink_KCHIP_HLA_AA_SNP_1000G_fam=plink_KCHIP_HLA_AA_SNP_1000G.get_fam().astype({'fid':str,'iid':str}).rename(columns={'fid':'FID','iid':'IID'})
plink_KCHIP_HLA_AA_SNP_1000G_bim=plink_KCHIP_HLA_AA_SNP_1000G.get_bim()


# In[3]:


#len(binary_continuous_traits)


# In[4]:


phenotypes=pd.read_csv(pheno_all_file_path,sep='\t')
phenotypes=phenotypes.set_index('ID').loc[plink_KCHIP_HLA_AA_SNP_1000G_fam['IID']]

binary_continuous_traits=phenotypes.columns.difference(['age','sex','cohort'])
#binary_continuous_traits

if 'ipykernel' in sys.argv[0]:
    ipykernel=True
    phenotype_name='blood_in_urine'
    step_idx=1
    mode_list=[0]
    #phenotype_name='height'
else:
    ipykernel=False
    phenotype_name=sys.argv[1]
    step_idx=int(sys.argv[2])
    mode_list=[int(i) for i in sys.argv[3].strip().split(',')]
    
if phenotype_name.isdigit():
    phenotype_name=int(phenotype_name)
    phenotype_name=binary_continuous_traits[phenotype_name]      


# In[5]:


binary_continuous_traits.tolist().index('blood_in_urine'),binary_continuous_traits.tolist().index('glucose_in_blood'),binary_continuous_traits.tolist().index('t2_diabetes'),binary_continuous_traits.tolist().index('diabetes'),


# In[6]:


data_out_assoc_phenotype_path=data_out_assoc_path+phenotype_name+'/'
pathlib.Path(data_out_assoc_phenotype_path).mkdir(parents=True, exist_ok=True)


# In[7]:


#for i in binary_continuous_traits:
#    if not os.path.exists('data/out_assoc/{}/step_01.plink.PHENO2.glm.linear'.format(i)) and not os.path.exists('data/out_assoc/{}/step_01.plink.PHENO2.glm.logistic'.format(i)):
#        print(i)


# In[8]:


pheno=pd.read_csv(data_out_pheno_path+phenotype_name+'.phe',sep='\t',names=['FID','IID','pheno'])
phenotype_type='binary' if len(pheno['pheno'][pheno['pheno']!=-9].value_counts())<3 else 'continuous'
phenotype_type


# In[9]:


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


# In[10]:


for step_idx_sub in range(1,step_idx):
    if os.path.exists(data_out_assoc_phenotype_path+'step_{:02d}.cond.stop'.format(step_idx_sub)):
        log.info("step_{}.cond.stop->stops ".format(step_idx_sub))
        sys.exit()


# In[11]:


log.info_head=lambda x: log.info('-'*int((100-len(x))/2)+x+'-'*int((100-len(x))/2))


# In[12]:


log.info_head("phenotype_name: {}, phenotype_type:{} , Step : {} ".format(phenotype_name,phenotype_type,step_idx))


# In[13]:


#os.path.exists(data_out_assoc_phenotype_path+'step_{:02d}.cond.stop'.format(step_idx))


# In[14]:


for step_idx_sub in range(1,step_idx+1):

    if os.path.exists(data_out_assoc_phenotype_path+'step_{:02d}.cond.stop'.format(step_idx_sub)):
        log.info("stops cond.stop")
        sys.exit()


# In[15]:


#gene_bed['name2'][gene_bed['name2'].str.contains('HLA')]


# In[16]:


#gene_assign.shape


# In[17]:


gene_bed_path='data/mart_export.txt'
gene_bed=pd.read_csv(gene_bed_path,sep='\t')
gene_bed=gene_bed.drop(columns='Exon stable ID')
gene_bed=gene_bed[(gene_bed['Gene start (bp)']>=plink_KCHIP_HLA_AA_SNP_1000G_bim.pos.min())&(gene_bed['Gene end (bp)']<=plink_KCHIP_HLA_AA_SNP_1000G_bim.pos.max())]
gene_bed=gene_bed[(gene_bed['Transcript type']=='protein_coding')]
gene_bed=gene_bed[~gene_bed.duplicated(['Gene name','Gene start (bp)','Gene end (bp)'])]
print(gene_bed.shape)
gene_bed=gene_bed[~gene_bed.duplicated(['Gene name'])]
print(gene_bed.shape)

gene_assign=plink_KCHIP_HLA_AA_SNP_1000G_bim[['pos']]

for idx,row in gene_bed.iterrows():
    gene_assign[row['Gene name']]=0
    
for idx,row in gene_bed.iterrows():    
    gene_assign[row['Gene name']][(gene_assign['pos']>=row['Gene start (bp)'])&(gene_assign['pos']<=row['Gene end (bp)'])]=1

gene_assign.columns=gene_assign.columns.str.replace('HLA-','HLA_')        
    
HLA_names=np.unique([i[0].split('_')[1] for i in plink_KCHIP_HLA_AA_SNP_1000G_bim[plink_KCHIP_HLA_AA_SNP_1000G_bim.index.str.contains('HLA_')].index.str.split('*')])

for HLA_name in HLA_names:
    gene_select=gene_assign[gene_assign.index.str.contains('HLA_'+HLA_name)|gene_assign.index.str.contains('SNPS_'+HLA_name)|gene_assign.index.str.contains('AA_'+HLA_name)]#print(gene_select.sort_values('pos').iloc[0],gene_select.sort_values('pos').iloc[-1])
    HLA_name='HLA_{}'.format(HLA_name)
    gene_assign[HLA_name][(gene_assign['pos']>=gene_select['pos'].min())&(gene_assign['pos']<=gene_select['pos'].max())]=1 


# In[18]:


"""
gene_bed_path='data/known_genes_chr6.hg19.txt'
gene_bed=pd.read_csv(gene_bed_path,sep='\t')
gene_bed=gene_bed[(gene_bed['txStart']>=plink_KCHIP_HLA_AA_SNP_1000G_bim.pos.min())&(gene_bed['txEnd']<=plink_KCHIP_HLA_AA_SNP_1000G_bim.pos.max())]
gene_bed=gene_bed[~gene_bed.duplicated(['name2'])]

gene_assign=plink_KCHIP_HLA_AA_SNP_1000G_bim[['pos']]

for idx,row in gene_bed.iterrows():
    gene_assign[row['name2']]=0
    
for idx,row in gene_bed.iterrows():    
    gene_assign[row['name2']][(gene_assign['pos']>=row['txStart'])&(gene_assign['pos']<=row['txEnd'])]=1

gene_assign.columns=gene_assign.columns.str.replace('HLA-','HLA_')        
    
HLA_names=np.unique([i[0].split('_')[1] for i in plink_KCHIP_HLA_AA_SNP_1000G_bim[plink_KCHIP_HLA_AA_SNP_1000G_bim.index.str.contains('HLA_')].index.str.split('*')])

for HLA_name in HLA_names:
    gene_select=gene_assign[gene_assign.index.str.contains('HLA_'+HLA_name)|gene_assign.index.str.contains('SNPS_'+HLA_name)|gene_assign.index.str.contains('AA_'+HLA_name)]#print(gene_select.sort_values('pos').iloc[0],gene_select.sort_values('pos').iloc[-1])
    HLA_name='HLA_{}'.format(HLA_name)
    gene_assign[HLA_name][(gene_assign['pos']>=gene_select['pos'].min())&(gene_assign['pos']<=gene_select['pos'].max())]=1 
"""    


# In[19]:


covariate_df=pd.read_csv(PC_path,sep='\t').set_index('ID').loc[plink_KCHIP_HLA_AA_SNP_1000G_fam['IID']]
covariate_df['age']=phenotypes['age']
covariate_df['sex']=phenotypes['sex']-1
if np.all((pheno['pheno']==-9).values | (phenotypes['cohort']==1).values):
    pass
elif np.all((pheno['pheno']==-9).values | (phenotypes['cohort']==2).values):
    pass
elif np.all((pheno['pheno']==-9).values | (phenotypes['cohort']==3).values):
    pass
elif np.all((pheno['pheno']==-9).values | ((phenotypes['cohort']==1).values|(phenotypes['cohort']==2).values)):
    covariate_df['AS']=phenotypes['cohort'].replace(1,1).replace(2,0).replace(3,0)
elif np.all((pheno['pheno']==-9).values | ((phenotypes['cohort']==2)|(phenotypes['cohort']==3).values)):
    covariate_df['CT']=phenotypes['cohort'].replace(1,0).replace(2,1).replace(3,0)    
elif np.all((pheno['pheno']==-9).values | ((phenotypes['cohort']==1)|(phenotypes['cohort']==3).values)):
    covariate_df['AS']=phenotypes['cohort'].replace(1,1).replace(2,0).replace(3,0)        
else:
    covariate_df['AS']=phenotypes['cohort'].replace(1,1).replace(2,0).replace(3,0)
    covariate_df['CT']=phenotypes['cohort'].replace(1,0).replace(2,1).replace(3,0)
    
plink_KCHIP_HLA_AA_SNP_1000G_fam.iloc[:,:2].merge(right=covariate_df,left_on='IID',right_index=True).fillna(-9).to_csv(data_out_assoc_phenotype_path+'covar',index=None,sep='\t')


# In[ ]:





# In[22]:


#covariate_df[(covariate_df['CT']==0) &(covariate_df['AS']==0)]['age'].mean()


# In[39]:


#phenotypes.index[phenotypes['cohort']==1]


# In[40]:


if 0 in mode_list:
    if os.path.exists(data_out_assoc_phenotype_path+'step_{:02d}.cond'.format(step_idx)):
        log.warning("Tried to construct .cond but. cond already exists")
    else:
        if step_idx==1:
            conditional_list=[]
            pd.Series(conditional_list).to_csv(data_out_assoc_phenotype_path+'step_{:02d}.cond'.format(step_idx),index=None,sep='\t',header=False)                        
        else:
            conditional_list=[]
            for step_idx_sub in range(1,step_idx):
                log.info("step idx sub {} cond finding".format(step_idx_sub))
                GAT_result=pd.read_csv(data_out_assoc_phenotype_path+'step_{:02d}.GAT.result.tsv'.format(step_idx_sub),sep='\t')
                plink_result=pd.read_csv(data_out_assoc_phenotype_path+'step_{:02d}.plink.PHENO2.glm.{}'.format(step_idx_sub,'logistic' if phenotype_type=='binary' else 'linear'),sep='\t')

                plink_result_munge=plink_result[plink_result['TEST']=='ADD'].drop(columns='#CHROM').rename(columns={'TEST':'term','ID':'marker_name','BETA':'coef','SE':'std','T_STAT':'Z','Z_STAT':'Z','OBS_CT':'nobs'})
                plink_result_munge['A2']=plink_result_munge.apply(lambda x: x['ALT'] if x['A1']==x['REF'] else x['REF'],axis=1)
                plink_result_munge=plink_result_munge.drop(columns=['REF','ALT'])
                plink_result_munge['note']='unphased bialleic'  
                result_merge=pd.concat([GAT_result,plink_result_munge],sort=True)[['marker_name','note','term','POS','Z','coef','std','chisq','df','A1','A2','multi_allele','nobs','P']]
            
                result_merge_sorted=result_merge.astype({'P':float}).sort_values('P')   

                result_merge.to_csv(data_out_assoc_phenotype_path+'step_{:02d}.result.tsv'.format(step_idx_sub),index=None,sep='\t',header=True)                                        
                result_merge_sorted.to_csv(data_out_assoc_phenotype_path+'step_{:02d}.result_sorted.tsv'.format(step_idx_sub),index=None,sep='\t',header=True)                                        
                
                log.info(result_merge_sorted.iloc[:5])
                if np.isnan(result_merge_sorted.iloc[0]['P']) or result_merge_sorted.iloc[0]['P']>5e-8:
                    pd.Series(conditional_list).to_csv(data_out_assoc_phenotype_path+'step_{:02d}.cond.stop'.format(step_idx),index=None,sep='\t',header=False)                        
                    log.info("p value insignificant")
                    break

                marker_name=result_merge_sorted.iloc[0].marker_name


                if marker_name[:3]=='AA_':
                    conditional_list.append('HLA_'+marker_name.split('_')[1])
                elif marker_name[:5]=='SNPS_':
                    conditional_list.append('HLA_'+marker_name.split('_')[1])
                elif marker_name[:4]=='HLA_':
                    conditional_list.append('HLA_'+marker_name.split('_')[1].split('*')[0])
                elif marker_name[:9]=='INS_SNPS_':
                    conditional_list.append('HLA_'+marker_name.split('_')[2])
                else:
                    #conditional_list.append(marker_name)                    
                    r2_list=[]
                    for idx_bim,(SNP,row) in enumerate(plink_KCHIP_HLA_AA_SNP_1000G_bim.iterrows()):
                        r2=pearsonr(plink_KCHIP_HLA_AA_SNP_1000G.get_geno_marker(marker_name),plink_KCHIP_HLA_AA_SNP_1000G.get_geno_marker(SNP))[0]**2
                        r2_list.append(r2)
                    r2_df=pd.DataFrame(r2_list,index=plink_KCHIP_HLA_AA_SNP_1000G_bim.index)                

                    if gene_assign[r2_df[0]>0.95][['HLA_A', 'HLA_B', 'HLA_C', 'HLA_DPA1', 'HLA_DPB1', 'HLA_DQA1', 'HLA_DQB1', 'HLA_DRB1']].sum().sum()==0:
                        conditional_list.append(marker_name)                     
                        log.info("{} not in HLA polymorphism-> added".format(marker_name))
                        
                    HLA_count=gene_assign[r2_df[0]>0.7][['HLA_A', 'HLA_B', 'HLA_C', 'HLA_DPA1', 'HLA_DPB1', 'HLA_DQA1', 'HLA_DQB1', 'HLA_DRB1']].sum(axis=0).sort_values(ascending=False)
                    print(HLA_count[HLA_count>0].index)
                    for i in HLA_count[HLA_count>0].index:
                        print(i)
                        conditional_list.append(i)
                        log.info("{} is strong LD with {}".format(marker_name,i))
                        break
                
                log.info('step idx {} conditional_list: {}'.format(step_idx_sub,conditional_list))
                conditional_list=np.unique(conditional_list).tolist()
                if step_idx_sub==step_idx-1:
                    pd.Series(conditional_list).to_csv(data_out_assoc_phenotype_path+'step_{:02d}.cond'.format(step_idx),index=None,sep='\t',header=False)                        
                
        log.info('conditional_list: {}'.format(conditional_list))
                        


# In[54]:


if 1 in mode_list:
    if not os.path.exists(data_out_assoc_phenotype_path+'step_{:02d}.cond'.format(step_idx)):
        log.warning("cond not existing... stop GAT...") 
        
        
    elif os.path.exists(data_out_assoc_phenotype_path+'step_{:02d}.GAT.result.tsv'.format(step_idx)):
        log.warning("GAT result already exits") 
        
    else:
        log.info("######################################### step {:02d} Phased Association  #########################################".format(step_idx))


        command='python Generic_Association_Tool/GAT.py         --assoc {assoc_mode}         --out {out}         --bfile {bfile}         --bgl-phased {bgl_phased}         --pheno {pheno}         --covar {covar}         --condition-list {cond}         --skip "(?P<name>6:[0-9]*_[A-Z]*/[\<\>A-Z\:0-9]*),(?P<name>AX\-[0-9]*),(?P<name>AFFX\-SP\-[0-9]*),(?P<name>SNPS_.*),(?P<name>INS_SNPS_.*)"         --multiallelic "(?P<name>HLA_[0-9A-Z]*)\*(?P<allele>[0-9:]*)"         --multiallelic-always "(?P<name>AA_[A-Z0-9]*_[\-0-9]*_[0-9]*_exon[0-9]*)_*(?P<allele>[A-Z]*)"'.format(
        assoc_mode='logistic' if phenotype_type=='binary' else 'linear',
        out=data_out_assoc_phenotype_path+'step_{:02d}.GAT'.format(step_idx),
        bfile=plink_1000G_path,
        bgl_phased=phased_KCHIP_HLA_AA_SNP_path,
        pheno=data_out_pheno_path+phenotype_name+'.phe',
        covar=data_out_assoc_phenotype_path+'covar',   
        cond=data_out_assoc_phenotype_path+'step_{:02d}.cond'.format(step_idx)  
        )    

        log.info(command)
        stdout,stderr=run_subprocess(command,dry=False)
        log.info(stdout)
        log.error(stderr) 



if 2 in mode_list:
    if not os.path.exists(data_out_assoc_phenotype_path+'step_{:02d}.cond'.format(step_idx)):
        log.warning("cond not existing... stop plink...") 

    elif os.path.exists(data_out_assoc_phenotype_path+'step_{:02d}.plink.PHENO2.glm.{}'.format(step_idx,'logistic' if phenotype_type=='binary' else 'linear')):
        log.warning("plink result already exits")         
        
    else:    
        log.info("######################################### step {:02d} Unphased Association  #########################################".format(step_idx))

        command='plink2         --bfile {bfile}         {assoc_mode}         --pheno {pheno}         --covar {covar}         --out {out}         --covar-variance-standardize         --threads 40'.format(
        bfile=plink_KCHIP_SNP_1000G_path,
        assoc_mode='--logistic' if phenotype_type=='binary' else '--linear',
        pheno=data_out_pheno_path+phenotype_name+'.phe',
        covar=data_out_assoc_phenotype_path+'step_{:02d}.GAT.covar_unphased.tsv'.format(step_idx),
        out=data_out_assoc_phenotype_path+'step_{:02d}.plink'.format(step_idx)                                                                         
        )

        log.info(command)
        stdout,stderr=run_subprocess(command,dry=False)
        log.info(stdout)
        log.error(stderr)  


# # Deprecated

# for i in range(0,101):
#     phenotype_name=binary_continuous_traits[i]      
#     pheno=pd.read_csv(data_out_pheno_path+phenotype_name+'.phe',sep='\t',names=['FID','IID','pheno'])
#     
#     if len(pheno['pheno'].unique())<4:
#         continue
#     pheno['pheno'][pheno['pheno']!=-9].hist(bins=50)
#     print(phenotype_name)
#     plt.show()
# 

# r2_list=[]
# for idx_bim,(SNP,row) in enumerate(plink_KCHIP_HLA_AA_SNP_1000G_bim.iterrows()):
#     r2=pearsonr(plink_KCHIP_HLA_AA_SNP_1000G.get_geno_marker('6:30165273_C/T'),plink_KCHIP_HLA_AA_SNP_1000G.get_geno_marker(SNP))[0]**2
#     r2_list.append(r2)
# r2_df=pd.DataFrame(r2_list,index=plink_KCHIP_HLA_AA_SNP_1000G_bim.index)  

# r2_df.sort_values(0)

# #pheno['pheno']

# #pheno['pheno'].loc[phenotypes.index[phenotypes['cohort']==1]].value_counts()

# pheno.set_index('IID').loc[phenotypes.index[phenotypes['cohort']==1]]['pheno'].value_counts()

# #pheno.set_index('IID').loc[phenotypes.index[phenotypes['cohort']==2]]['pheno'].value_counts()

# pheno.set_index('IID').loc[phenotypes.index[phenotypes['cohort']==3]]['pheno'].value_counts()

# """
# #r2_df.loc[gene_assign.index[gene_assign['HLA_DPB1']==1]].sort_values(0)
# r2_list=[]
# for idx_bim,(SNP,row) in enumerate(plink_KCHIP_HLA_AA_SNP_1000G_bim.iterrows()):
#     r2=pearsonr(plink_KCHIP_HLA_AA_SNP_1000G.get_geno_marker('6:31367865_G/A'),plink_KCHIP_HLA_AA_SNP_1000G.get_geno_marker(SNP))[0]**2
#     r2_list.append(r2)
# r2_df=pd.DataFrame(r2_list,index=plink_KCHIP_HLA_AA_SNP_1000G_bim.index)  
# #r2_df.loc['6:32665115_T/C']
# 
# r2_df.loc[gene_assign.index[gene_assign['HLA_B']==1]].sort_values(0)
# """
# 

# #pheno['pheno'].loc[pheno['pheno']!=-9]
# #covariate_df['age']
# #pheno['pheno']

# plt.scatter(covariate_df['age'].loc[pheno.set_index('IID')['pheno']!=-9],pheno.set_index('IID')['pheno'].loc[pheno.set_index('IID')['pheno']!=-9])

# pheno

# #phenotypes['cohort']==2

# (plink_KCHIP_HLA_AA_SNP_1000G.get_geno_marker(marker_name1)==0)[phenotypes['cohort']==2].sum(),(plink_KCHIP_HLA_AA_SNP_1000G.get_geno_marker(marker_name1)==1)[phenotypes['cohort']==2].sum(),(plink_KCHIP_HLA_AA_SNP_1000G.get_geno_marker(marker_name1)==2)[phenotypes['cohort']==2].sum()

# marker_name1='6:30165273_C/T'
# (plink_KCHIP_HLA_AA_SNP_1000G.get_geno_marker(marker_name1)==0).sum(),(plink_KCHIP_HLA_AA_SNP_1000G.get_geno_marker(marker_name1)==1).sum(),(plink_KCHIP_HLA_AA_SNP_1000G.get_geno_marker(marker_name1)==2).sum()

# ((plink_KCHIP_HLA_AA_SNP_1000G.get_geno_marker(marker_name1)==1).sum()+2*(plink_KCHIP_HLA_AA_SNP_1000G.get_geno_marker(marker_name1)==2).sum())/\
# (2*len(((plink_KCHIP_HLA_AA_SNP_1000G.get_geno_marker(marker_name1)))))

# marker_name1='rs2523942'
# (plink_KCHIP_HLA_AA_SNP_1000G.get_geno_marker(marker_name1)==0).sum(),(plink_KCHIP_HLA_AA_SNP_1000G.get_geno_marker(marker_name1)==1).sum(),(plink_KCHIP_HLA_AA_SNP_1000G.get_geno_marker(marker_name1)==2).sum()

# ((plink_KCHIP_HLA_AA_SNP_1000G.get_geno_marker(marker_name1)==1).sum()+2*(plink_KCHIP_HLA_AA_SNP_1000G.get_geno_marker(marker_name1)==2).sum())/\
# (2*len(((plink_KCHIP_HLA_AA_SNP_1000G.get_geno_marker(marker_name1)))))

# 

# marker_name1='6:29925086_T/G'
# (plink_KCHIP_HLA_AA_SNP_1000G.get_geno_marker(marker_name1)==0).sum(),(plink_KCHIP_HLA_AA_SNP_1000G.get_geno_marker(marker_name1)==1).sum(),(plink_KCHIP_HLA_AA_SNP_1000G.get_geno_marker(marker_name1)==2).sum()

# ((plink_KCHIP_HLA_AA_SNP_1000G.get_geno_marker(marker_name1)==1).sum()+2*(plink_KCHIP_HLA_AA_SNP_1000G.get_geno_marker(marker_name1)==2).sum())/\
# (2*len(((plink_KCHIP_HLA_AA_SNP_1000G.get_geno_marker(marker_name1)))))

# (plink_KCHIP_HLA_AA_SNP_1000G.get_geno_marker('6:30165273_C/T')==-1).sum()

# 2*np.sqrt(4036)*np.sqrt(86072)

# phenotypes['cohort']==1

# pheno['pheno'][pheno['pheno']>0].mean(),pheno['pheno'][pheno['pheno']>0].std()

# 
# marker_name1='6:30165273_C/T'
# marker_name2='6:30167835_C/T'

# pearsonr(plink_KCHIP_HLA_AA_SNP_1000G.get_geno_marker(marker_name1),plink_KCHIP_HLA_AA_SNP_1000G.get_geno_marker(marker_name2))

# ((plink_KCHIP_HLA_AA_SNP_1000G.get_geno_marker(marker_name1)==2)&(plink_KCHIP_HLA_AA_SNP_1000G.get_geno_marker(marker_name2)==2)).sum()

# pheno['pheno'][plink_KCHIP_HLA_AA_SNP_1000G.get_geno_marker(marker_name2)==0].replace(-9,np.nan).mean(),pheno['pheno'][plink_KCHIP_HLA_AA_SNP_1000G.get_geno_marker(marker_name2)==0].replace(-9,np.nan).hist(bins=50)
# 

# pheno['pheno'][plink_KCHIP_HLA_AA_SNP_1000G.get_geno_marker(marker_name2)==1].replace(-9,np.nan).mean(),pheno['pheno'][plink_KCHIP_HLA_AA_SNP_1000G.get_geno_marker(marker_name2)==1].replace(-9,np.nan).hist(bins=50)

# pheno['pheno'][plink_KCHIP_HLA_AA_SNP_1000G.get_geno_marker(marker_name2)==2].replace(-9,np.nan).mean(),pheno['pheno'][plink_KCHIP_HLA_AA_SNP_1000G.get_geno_marker(marker_name2)==2].replace(-9,np.nan).hist(bins=50)
# 

# pheno['pheno'].replace(-9,np.nan).hist(bins=50)

# pheno['pheno'][plink_KCHIP_HLA_AA_SNP_1000G.get_geno_marker(marker_name1)==0].replace(-9,np.nan).mean(),pheno['pheno'][plink_KCHIP_HLA_AA_SNP_1000G.get_geno_marker(marker_name1)==0].replace(-9,np.nan).hist(bins=50)

# pheno['pheno'][plink_KCHIP_HLA_AA_SNP_1000G.get_geno_marker(marker_name1)==1].replace(-9,np.nan).mean(),pheno['pheno'][plink_KCHIP_HLA_AA_SNP_1000G.get_geno_marker(marker_name1)==1].replace(-9,np.nan).hist(bins=50)

# pheno['pheno'][plink_KCHIP_HLA_AA_SNP_1000G.get_geno_marker(marker_name1)==2].replace(-9,np.nan).mean(),pheno['pheno'][plink_KCHIP_HLA_AA_SNP_1000G.get_geno_marker(marker_name1)==2].replace(-9,np.nan).hist(bins=50)
# 
# 

# phenotype_raw_CT=pd.read_csv(phenotype_raw_CT_path,sep='\t',index_col='ID')

# #phenotype_raw_CT['CT1_ALP'][(phenotype_raw_CT['CT1_ALP']!=66666)&(phenotype_raw_CT['CT1_ALP']!=99999)&(phenotype_raw_CT['CT1_ALP']<5000)].hist(bins=100)
