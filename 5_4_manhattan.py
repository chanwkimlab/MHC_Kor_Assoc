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

import matplotlib.patches as patches
import matplotlib

from basic_tools import *

"""

jupyter nbconvert 5_4_manhattan.ipynb --to script

for i in {00..20};do python 5_4_manhattan.py $i;done
for i in {21..40};do python 5_4_manhattan.py $i;done
for i in {41..60};do python 5_4_manhattan.py $i;done
for i in {61..80};do python 5_4_manhattan.py $i;done
for i in {81..98};do python 5_4_manhattan.py $i;done

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
    phenotype_name='allergic_disease'
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

#gene_assign.columns=gene_assign.columns.str.replace('HLA-','HLA_')        
    
HLA_names=np.unique([i[0].split('_')[1] for i in plink_KCHIP_HLA_AA_SNP_1000G_bim[plink_KCHIP_HLA_AA_SNP_1000G_bim.index.str.contains('HLA_')].index.str.split('*')])

for HLA_name in HLA_names:
    gene_select=gene_assign[gene_assign.index.str.contains('HLA_'+HLA_name)|gene_assign.index.str.contains('SNPS_'+HLA_name)|gene_assign.index.str.contains('AA_'+HLA_name)]#print(gene_select.sort_values('pos').iloc[0],gene_select.sort_values('pos').iloc[-1])
    HLA_name='HLA-{}'.format(HLA_name)
    gene_assign[HLA_name][(gene_assign['pos']>=gene_select['pos'].min())&(gene_assign['pos']<=gene_select['pos'].max())]=1 


# In[7]:


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


# In[8]:


log.info_head=lambda x: log.info('-'*int((100-len(x))/2)+x+'-'*int((100-len(x))/2))


# In[11]:


plt.rcParams["figure.figsize"] = (20,5)
plt.rcParams["font.size"] = 15
plt.rcParams['font.family']='Arial'


for step_idx_sub in range(1,100):
    log.info_head("phenotype_name: {}, phenotype_type:{} , Step : {} ".format(phenotype_name,phenotype_type,step_idx_sub))
    plt.clf()
    #print(step_idx_sub)
    if os.path.exists(data_out_assoc_phenotype_path+'step_{:02d}.cond.stop'.format(step_idx_sub)):
        print('meets end',step_idx_sub)
        break
    result_merge=pd.read_csv(data_out_assoc_phenotype_path+'step_{:02d}.merge.result.tsv'.format(step_idx_sub),sep='\t')



    ##########################################
    #############Plotting#####################
    ##########################################

    data=result_merge[['marker_name','POS','P']].copy()



    data['-log10_P']=-np.log10(data['P'])
    data['size']=-np.log10(data['P'])*5
    data['size']=350*(data['size']-data['size'].min())/(data['size'].max()-data['size'].min())+0.001

    POS_MB_min=27.9
    POS_MB_max=35.1


    tick_pad=20
    gene_annot_pad=data['-log10_P'].max()/30
    gene_annot_height=data['-log10_P'].max()/30

    data['POS_MB']=data['POS']/1000000
    data['check']=False

    data=data.sort_values('P',ascending=False).reset_index()

    ax = plt.subplot(111)


    data_select_list=[]

    for gene_name in ['HLA-A','HLA-B','HLA-C']:
        pos_min=gene_assign[gene_assign[gene_name]==1]['pos'].min()/1000000
        pos_max=gene_assign[gene_assign[gene_name]==1]['pos'].max()/1000000
        color='#fa3c14'

        rect = patches.Rectangle((pos_min,-gene_annot_height-gene_annot_pad),pos_max-pos_min,gene_annot_height,linewidth=1,edgecolor=None,facecolor=color)
        ax.add_patch(rect)
        #plt.annotate(gene_name.replace('_','-'),xy=((pos_min+pos_max)/2,-5),xytext=((pos_min+pos_max)/2,-10),fontsize=3,ha='center',arrowprops=dict(arrowstyle="- >",connectionstyle="arc3,rad=0",ls='dashed'))#,bbox= dict(boxstyle="round,pad=0.3", fc=(1,1,1,0.5), ec="black", lw=0.3),size=15)

        select=(data['POS_MB']>pos_min)&(data['POS_MB']<pos_max)&(data['check']==False)
        data.loc[select,'check']=True
        data_select=data.loc[select].copy()   
        data_select['color']=color

        data_select_list.append(data_select.sort_values('-log10_P'))



    for gene_name in ['HLA-DPA1', 'HLA-DPB1','HLA-DQA1', 'HLA-DQB1','HLA-DRB1','HLA-DRA']:
        pos_min=gene_assign[gene_assign[gene_name]==1]['pos'].min()/1000000
        pos_max=gene_assign[gene_assign[gene_name]==1]['pos'].max()/1000000
        color='#2850c8'

        rect = patches.Rectangle((pos_min,-gene_annot_height-gene_annot_pad),pos_max-pos_min,gene_annot_height,linewidth=1,edgecolor=None,facecolor=color)
        ax.add_patch(rect)
        #plt.annotate(gene_name.replace('_','-'),xy=((pos_min+pos_max)/2,-5),xytext=((pos_min+pos_max)/2,-10),fontsize=3,ha='center',arrowprops=dict(arrowstyle="- >",connectionstyle="arc3,rad=0",ls='dashed'))#,bbox= dict(boxstyle="round,pad=0.3", fc=(1,1,1,0.5), ec="black", lw=0.3),size=15)

        select=(data['POS_MB']>pos_min)&(data['POS_MB']<pos_max)&(data['check']==False)
        data.loc[select,'check']=True
        data_select=data.loc[select].copy()       
        data_select['color']=color
        data_select_list.append(data_select.sort_values('-log10_P'))


    #for gene_name in gene_assign.columns[gene_assign.columns.str.contains('HLA')].difference(['HLA_A','HLA_B','HLA_C','HLA_DPA1','HLA_DPB1','HLA_DQA1', 'HLA_DQB1', 'HLA_DRB1','HLA-DRA']):#['HLA_A29.1']):
    for gene_name in ['HLA-F', 'HLA-G', 'HLA-E', 'HLA-DRB5', 'HLA-DQA2', 'HLA-DQB2', 'HLA-DOB', 'HLA-DMB', 'HLA-DMA', 'HLA-DOA', 'MICA','MICB','TAP2','TAP1']:
    #for gene_name in ['HLA-F', 'HLA-G',  'HLA-E', 'HLA-DRB5','HLA-DMB', 'HLA-DMA', 'HLA-DOA', 'HLA-DOB','HLA-DQA2', 'HLA-DQB2']:
        pos_min=gene_assign[gene_assign[gene_name]==1]['pos'].min()/1000000
        pos_max=gene_assign[gene_assign[gene_name]==1]['pos'].max()/1000000
        color='#28c828'

        rect = patches.Rectangle((pos_min,-gene_annot_height-gene_annot_pad),pos_max-pos_min,gene_annot_height,linewidth=1,edgecolor=None,facecolor=color)
        ax.add_patch(rect)
        #plt.annotate(gene_name.replace('_','-'),xy=((pos_min+pos_max)/2,-5),xytext=((pos_min+pos_max)/2,-10),fontsize=3,ha='center',arrowprops=dict(arrowstyle="- >",connectionstyle="arc3,rad=0",ls='dashed'))#,bbox= dict(boxstyle="round,pad=0.3", fc=(1,1,1,0.5), ec="black", lw=0.3),size=15)


        select=(data['POS_MB']>pos_min)&(data['POS_MB']<pos_max)&(data['check']==False)
        data.loc[select,'check']=True
        data_select=data.loc[select].copy()       
        data_select['color']=color
        data_select_list.append(data_select.sort_values('-log10_P'))  


    for gene_name in gene_assign.columns[~gene_assign.columns.str.contains('HLA')]:
        #print(gene_name)
        pos_min=gene_assign[gene_assign[gene_name]==1]['pos'].min()/1000000
        pos_max=gene_assign[gene_assign[gene_name]==1]['pos'].max()/1000000
        color='#fafac8'
        #color=matplotlib.colors.rgb2hex(0.5*np.array(plt.cm.rainbow(np.random.rand())[:3]))
        rect = patches.Rectangle((pos_min,-gene_annot_height-gene_annot_pad),pos_max-pos_min,gene_annot_height,linewidth=1,edgecolor='grey',facecolor=color,alpha=0.5)
        #rect = patches.Rectangle((pos_min,-5.5),pos_max-pos_min,3,linewidth=1,edgecolor='grey',facecolor=color,alpha=0.5)
        ax.add_patch(rect)
        #plt.annotate(gene_name.replace('_','-'),xy=((pos_min+pos_max)/2,-5),xytext=((pos_min+pos_max)/2,-10),fontsize=3,ha='center',arrowprops=dict(arrowstyle="- >",connectionstyle="arc3,rad=0",ls='dashed'))#,bbox= dict(boxstyle="round,pad=0.3", fc=(1,1,1,0.5), ec="black", lw=0.3),size=15)


        select=(data['POS_MB']>pos_min)&(data['POS_MB']<pos_max)&(data['check']==False)
        data.loc[select,'check']=True
        data_select=data.loc[select].copy()       
        data_select['color']=color
        data_select_list.append(data_select.sort_values('-log10_P'))     


    color='#fafac8'
    select=(data['check']==False)
    data.loc[select,'check']=True
    data_select=data.loc[select].copy()    
    data_select['color']=color
    data_select_list.append(data_select.sort_values('-log10_P'))  
    #data_select=data.loc[select].sort_values('-log10_P')
    #plt.scatter(x=data_select['POS_MB'],y=data_select['-log10_P'],sizes=data_select['size'],marker='d',linewidth=0.3,edgecolor='black',color=color)   

    for data_select in data_select_list[::-1]:
        plt.scatter(x=data_select['POS_MB'],y=data_select['-log10_P'],sizes=data_select['size'],marker='d',linewidth=0.3,edgecolor='black',color=data_select['color'])   



    plt.plot(np.linspace(POS_MB_min,POS_MB_max,100),-np.log10(5e-8)*np.ones_like(np.linspace(POS_MB_min,POS_MB_max,100)),'--',color='black',alpha=0.7,linewidth=0.7)    
    #ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)  
    ax.spines['bottom'].set_position('zero')#ax.spines['bottom'].set_visible(False) 

    ax.tick_params(axis='x', which='major', pad=tick_pad)

    plt.xlim(POS_MB_min,POS_MB_max)
    #ax.set_ylim(bottom=0.)
    #ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    #plt.savefig('savefig_polygenicity.png',bbox_inches="tight")        


    ##########################################
    #############Plotting#####################
    ##########################################
    #plt.savefig('savefig_polygenicity.png',bbox_inches="tight") 
    #plt.show()
    plt.savefig(data_out_assoc_phenotype_path+'step_{:02d}.merge.manhattan.png'.format(step_idx_sub),bbox_inches="tight")
    #plt.show()
    #break
#break


# tar -czvf s_gwas_all.tar data/out_assoc/*/*merge*.png
# tar -czvf s_figure_manhattan.tar data/out_assoc/*/*merge*
#for HLA_name in gene_assign.columns[gene_assign.columns.str.contains('HLA')].difference(['HLA_A29.1']):
#    gene_assign[gene_assign[HLA_name]==1]['pos'].min(),gene_assign[gene_assign[HLA_name]==1]['pos'].max()
# In[95]:


#ax.tick_params?

