#!/usr/bin/env python
# coding: utf-8

# jupyter nbconvert omnibus_test.ipynb --to script
# 
# 
# 
# example
# omnibus_test.py    --assoc linear \
#     --out data/out_assoc/hba1c/step_02.omnibus \
#     --pheno data/out_assoc/hba1c/phenotype.phe \
#     --fam data/genotype/4_merge/KCHIP_HLA_SNP_1000G_merged.fam \
#     --covar data/out_assoc/hba1c/step_02.omnibus.covar.temp \
#     --aa data/out_assoc/hba1c/step_02.aa \
#     --condition-list data/out_assoc/hba1c/step_02.omnibus.cond\
#     
# numpy              1.17.4 
# pandas             0.25.3   
# scipy              1.3.2 
# 
# statsmodels        0.10.1

# In[1]:


import os
import argparse
import logging
import time

import pandas as pd
import numpy as np
from scipy.stats import chi2

import statsmodels.api as sm


# In[2]:


def dir_path(path):
    if os.path.exists(os.path.dirname(path)):
        return path
    else:
        raise argparse.ArgumentTypeError(f"readable_dir:{path} is not a valid path")

def file_path(path):
    if os.path.isfile(path):
        return path
    else:
        raise argparse.ArgumentTypeError(f"readable_dir:{path} is not a valid path")


parser = argparse.ArgumentParser(description='Omnibus test')

#required mode
parser.add_argument('--assoc',choices=['linear','logistic'],required=True)

#required output file
parser.add_argument('--out', type=dir_path,required=True,help='output file prefix. (prefix.log, prefix.assoc will be generated)')

#required input files
parser.add_argument('--fam', type=file_path,required=True,help='plink fam file')
parser.add_argument('--aa', type=file_path,required=True,help='AA call text file. For more information, refer to documentation. (like plink ped)')
parser.add_argument('--pheno', type=file_path,required=True,help='format is the same as plink. Tab-delimited file without header of which the first and second columns is family and within-family IDs respectively and the third columns are pheotype')

#optional
parser.add_argument('--covar', type=file_path,help='format is the same as plink')
parser.add_argument('--condition-list',type=file_path,help='format is the same as plink')


# In[3]:


#debug=False
debug=True

if debug:
    args=parser.parse_args('--assoc linear --out data/out_assoc/hba1c/step_02.omnibus --pheno data/out_assoc/hba1c/phenotype.phe --fam data/genotype/4_merge/KCHIP_HLA_SNP_1000G_merged.fam --covar data/out_assoc/hba1c/step_02.omnibus.covar.temp --aa data/out_assoc/hba1c/step_02.aa --condition-list data/out_assoc/hba1c/step_02.omnibus.cond'.split(' '))
else:
    args=parser.parse_args()


# In[4]:


log = logging.getLogger('logger')
log.setLevel(logging.DEBUG)

log_file_path=args.out+'.log'
fileHandler = logging.FileHandler(log_file_path)
streamHandler = logging.StreamHandler()

formatter = logging.Formatter('%(message)s')
fileHandler.setFormatter(formatter)
streamHandler.setFormatter(formatter)

log.addHandler(fileHandler)
log.addHandler(streamHandler)


# In[10]:


log.info_head=lambda x: log.info('-'*int((100-len(x))/2)+x+'-'*int((100-len(x))/2))


# In[11]:


log.info_head("Amino Acid Association Test (Omnibus Test)")


# In[12]:


log.info('Parameters\n'+'\n'.join(['--{} {}'.format(key,value) for key,value in vars(args).items()]))


# In[13]:


log.info('Working directory: '+os.getcwd())


# In[14]:


log.info("Start time: "+time.strftime('%c', time.localtime(time.time())))


# In[15]:


assoc=args.assoc
out=args.out


# In[18]:


log.info_head("Data Loading")


# # parse input files

# In[19]:


fam=pd.read_csv(args.fam,header=None,sep=' ',names=['FID','IID','fID','mID','sex','pheno']).astype({'FID':str,'IID':str,'fID':str,'mID':str,'sex':int,'pheno':float})


# In[20]:


log.info("{} samples ({} males, {} females) loaded from {}".format(fam.shape[0],(fam['sex']==1).sum(),(fam['sex']==2).sum(),args.fam))


# In[28]:


aa_IID_list=None
aa_sex_list=None
aa_marker_name_list=[]
aa_marker_data_list=[]
with open(args.aa,'r') as f:
    line_cnt=0
    while True:
        line=f.readline()            
        if not line:
            break        
        line_cnt+=1
        line_split=line.strip().split(' ')
        line_type,line_id,line_data=line_split[0],line_split[1],line_split[2:]
        if line_type=='P':
            pass
        elif line_type=='fID':
            pass
        elif line_type=='mID':
            pass        
        elif line_type=='I':        
            aa_IID_list1=np.array([line_data[i] for i in range(0,len(line_data),2)])
            aa_IID_list2=np.array([line_data[i+1] for i in range(0,len(line_data),2)])
            if np.all(aa_IID_list1==aa_IID_list2):
                aa_IID_list=aa_IID_list1
            else:
                raise
        elif line_type=='C':
            aa_sex_list1=np.array([line_data[i] for i in range(0,len(line_data),2)])
            aa_sex_list2=np.array([line_data[i+1] for i in range(0,len(line_data),2)])
            if np.all(aa_sex_list1==aa_sex_list2):
                aa_sex_list=aa_sex_list1.astype(int)
            else:
                raise
        elif line_type=='M':
            aa_marker_name_list.append(line_id)
            line_data=np.array([np.nan if x=='NA' else x for x in line_data])
            aa_marker_data_list.append(line_data)
        else:
            print(line_type)
            raise 
            
aa_marker_name_list_aaonly=pd.Series(aa_marker_name_list)[pd.Series(aa_marker_name_list).str.slice(stop=3)=='AA_'].values            


# In[29]:


log.info("{} variants ({} AAs) loaded from {}".format(len(aa_marker_name_list),len(aa_marker_name_list_aaonly),args.aa))


# In[30]:


pheno=pd.read_csv(args.pheno,header=None,sep='\t',names=['FID','IID','pheno'])
pheno['pheno']=pheno['pheno'].replace(-9,np.nan)


# In[31]:


log.info("pheotype loaded from {}".format(args.pheno))


# In[32]:


if args.assoc=='linear':
    assert len(pheno['pheno'].unique())>2
else:
    assert np.all(np.isnan(pheno['pheno'])|(pheno['pheno']==1)|(pheno['pheno']==2))
    pheno['pheno']=pheno['pheno']-1


# In[33]:


log.info("{} pheotype loaded from {}".format(pheno.shape[0],args.pheno))
log.info("Among them, valid: {}, missing: {}".format((~pheno['pheno'].isnull()).sum(),pheno['pheno'].isnull().sum()))
if assoc=='linear':
    log.info("mean={:4f} std={:4f} median={:4f} min={:4f} max={:4f}".format(pheno['pheno'].mean(),pheno['pheno'].std(),pheno['pheno'].median(),pheno['pheno'].min(),pheno['pheno'].max()))
else:
    log.info("case: {} / control: {}".format((pheno['pheno']==1).sum(),(pheno['pheno']==0).sum()))


# # parse optional input files

# In[34]:


if args.covar is None:
    covar=fam.iloc[:,:2]
else:
    covar=pd.read_csv(args.covar,sep='\t')
    covar.columns=['FID','IID']+covar.columns[2:].tolist()
    covar=covar.astype({'FID':str,'IID':str})
    
    covar.iloc[:,2:]=covar.iloc[:,2:].astype(float)
    covar.iloc[:,2:]=covar.iloc[:,2:].replace(-9,np.nan)
    
    log.info("{} covariates loaded from {}".format(len(covar.columns[2:]),args.covar))


# In[35]:


if args.condition_list is None:
    condition_list=[]
else:
    with open(args.condition_list,'r') as f:
        condition_list=f.read().strip().split('\n')
    log.info("{} covariates loaded from {}".format(len(condition_list),args.condition_list))
    log.info(', '.join(condition_list))


# # check idx integrity

# In[20]:


assert np.all(fam['IID']==aa_IID_list)
assert np.all(fam['sex']==aa_sex_list)
assert len(aa_marker_name_list)==len(aa_marker_data_list)

assert np.all(fam['IID']==covar['IID'])
assert len(set(condition_list).difference(set(aa_marker_name_list)))==0


# # Run regression

# In[36]:


log.info_head("Regression")


# In[74]:


def marker_data_to_onehot(aa_marker_data):
    
    aa_marker_data_unique=np.unique(aa_marker_data)
    aa_marker_data_unique_nonan=aa_marker_data_unique[aa_marker_data_unique!='nan'].tolist()
    aa_marker_data_unique_nan=aa_marker_data_unique_nonan+['nan']
    
    aa_marker_data_int_nan=list(map(lambda x: aa_marker_data_unique_nan.index(x),aa_marker_data))
    
    aa_marker_data_onehot_nan=np.zeros((len(aa_marker_data),len(aa_marker_data_unique_nan)))
    aa_marker_data_onehot_nan[np.arange(len(aa_marker_data)),aa_marker_data_int_nan]=1
    aa_marker_data_onehot_nonan=aa_marker_data_onehot_nan[:,:-1]
    
    return aa_marker_data_unique_nonan,aa_marker_data_onehot_nonan

def prepare_onehot(aa_marker_data_unique_nonan,aa_marker_data_onehot_nonan,set_nan=True,cut_mostfrequent=True):
    assert len(aa_marker_data_unique_nonan)==aa_marker_data_onehot_nonan.shape[1]
    assert np.isnan(aa_marker_data_onehot_nonan).sum()==0
    assert (np.nan not in aa_marker_data_unique_nonan) and ('nan' not in aa_marker_data_unique_nonan)
    
    aa_marker_data_onehot_nonan_sumrow=np.sum(aa_marker_data_onehot_nonan,axis=1)
    aa_marker_data_onehot_nonan_sumcol=np.sum(aa_marker_data_onehot_nonan,axis=0)
    #print(aa_marker_data_onehot_nonan_sumcol,np.argmax(aa_marker_data_onehot_nonan_sumcol))
    #print(aa_marker_data_onehot_nonan_sumrow.shape,aa_marker_data_onehot_nonan_sumcol.shape)
    #print(np.argmax(aa_marker_data_onehot_nonan_sumcol))
    if set_nan:
        aa_marker_data_onehot_nonan[aa_marker_data_onehot_nonan_sumrow!=1,:]=np.nan
    if cut_mostfrequent:
        aa_marker_data_unique_nonan=np.delete(aa_marker_data_unique_nonan, np.argmax(aa_marker_data_onehot_nonan_sumcol))
        aa_marker_data_onehot_nonan=np.delete(aa_marker_data_onehot_nonan, np.argmax(aa_marker_data_onehot_nonan_sumcol), axis=1)
    return aa_marker_data_unique_nonan,aa_marker_data_onehot_nonan


# ## y data

# In[75]:


y_data_pheno=np.repeat(pheno['pheno'].values,2)


# In[76]:


y_data_pheno_nan=np.isnan(y_data_pheno)


# In[77]:


if np.var(y_data_pheno[~y_data_pheno_nan])==0:
    log.error("No variance in y_data")
    raise


# In[85]:


log.info("missing values in y_data_pheno: {}".format(y_data_pheno_nan.sum()))


# ## x data

# In[78]:


x_data_intercept=np.array([np.ones(2*fam.shape[0])]).transpose()


# In[86]:


x_data_covar=np.repeat(covar.iloc[:,2:].values,2,axis=0)


# In[87]:


x_data_covar_nan=np.any(np.isnan(x_data_covar),axis=1)
log.info("missing values in x_data_covar: {}".format(x_data_covar_nan.sum()))


# In[81]:


x_data_covar_varcheck=np.array([np.var(covar) for covar in x_data_covar[~x_data_covar_nan].transpose()])==0

if np.any(x_data_covar_varcheck):
    log.info("Removed covariates of no variance"+", ".join(covar.columns[2:][x_data_covar_varcheck].tolist()))
    x_data_covar=np.delete(x_data_covar,np.arange(len(x_data_covar_varcheck))[x_data_covar_varcheck],axis=1)


# In[88]:


if (x_data_covar.size!=0) and (np.linalg.matrix_rank(x_data_covar)<x_data_covar.shape[1]):
    log.info("duplicated covariates were found. (rank(covariates)< # of covarirates))")
    #raise
    #x_data_covar=np.delete(x_data_covar,np.arange(len(x_data_covar_varcheck))[x_data_covar_varcheck],axis=1)


# In[ ]:





# In[89]:


temp_list=[x_data_intercept[:,0:0]]
for aa_marker_name in condition_list:
    aa_marker_data=aa_marker_data_list[aa_marker_name_list.index(aa_marker_name)]

    aa_marker_data_unique,aa_marker_data_onehot=marker_data_to_onehot(aa_marker_data)
    aa_marker_data_unique_cut,aa_marker_data_onehot_cut=prepare_onehot(aa_marker_data_unique,aa_marker_data_onehot,set_nan=True,cut_mostfrequent=True)        
    temp_list.append(aa_marker_data_onehot_cut)
    
x_data_condition=np.concatenate(temp_list,axis=1)


# In[136]:


log.info("missing values in x_data_condition: {}".format(x_data_condition_nan.sum()))


# In[135]:


"""
x_data_condition=np.delete(x_data_condition,list(todel),axis=1)
todel=set()
for i in range(x_data_condition.shape[1]):
    for j in range(i+1,x_data_condition.shape[1]):
        if np.corrcoef(x_data_condition[~x_data_condition_nan][:,i],x_data_condition[~x_data_condition_nan][:,j])[1,0]>0.9:
            todel.add(i)
            todel.add(j)
x_data_condition_nan=np.any(np.isnan(x_data_condition),axis=1)  
np.linalg.matrix_rank(x_data_condition[~x_data_condition_nan]),x_data_condition.shape[1]
"""            


# In[ ]:


log.info('[{:3d}/{:3d}] {:10s} {:15s} {:5s} {:.5s}({}) {:.5s}'.format(
                            0,
                            len(aa_marker_name_list_aaonly),
                            'ID',
                            'residues',
                            'n_obs',
                            'chisq',
                            'df',
                            'P'
                        ))
assoc_result_list=[]

for idx,aa_marker_name in enumerate(aa_marker_name_list_aaonly):
    aa_marker_data=aa_marker_data_list[idx]
    aa_marker_data_unique,aa_marker_data_onehot=marker_data_to_onehot(aa_marker_data)
    aa_marker_data_unique_cut,aa_marker_data_onehot_cut=prepare_onehot(aa_marker_data_unique,aa_marker_data_onehot,set_nan=True,cut_mostfrequent=True)
    
    x_data_aa_marker=aa_marker_data_onehot_cut
    x_data_aa_marker_nan=np.any(np.isnan(x_data_aa_marker),axis=1)
    
    x_data_nan=np.logical_or.reduce([x_data_covar_nan,x_data_condition_nan,x_data_aa_marker_nan])
    
    y_data=y_data_pheno
    y_data_nan=y_data_pheno_nan    
    x_y_data_nan=(x_data_nan)|(y_data_nan)


    x_data_null=np.concatenate([x_data_intercept,x_data_covar,x_data_condition,x_data_aa_marker],axis=1)[~x_y_data_nan]
    x_data_alt=np.concatenate([x_data_intercept,x_data_covar,x_data_condition],axis=1)[~x_y_data_nan]
    y_data=y_data[~x_y_data_nan]
    
    
    family=(sm.families.Gaussian() if assoc=='linear' else sm.families.Binomial())
    model_null = sm.GLM(y_data,x_data_null, family=family,missing='raise')
    model_alt = sm.GLM(y_data,x_data_alt, family=family,missing='raise')
    
    try:
        result_alt = model_alt.fit()
        result_null = model_null.fit()
    except sm.tools.sm_exceptions.PerfectSeparationError as e:
        nobs=np.nan
        chisq_diff=np.nan
        df_diff=np.nan
        p_value=np.nan
    else:
        assert result_alt.nobs==result_null.nobs
        nobs=result_alt.nobs
        chisq_diff=2*(result_null.llf-result_alt.llf)
        df_diff=result_null.df_model-result_alt.df_model
        p_value=1 - chi2.cdf(chisq_diff,df_diff)

    #print(result_alt.summary())
    assoc_result={'idx':idx+1,
                  'ID':aa_marker_name,
                  'residues':','.join(aa_marker_data_unique),
                  'n_obs':nobs,
                  'chisq':chisq_diff,
                  'df':df_diff,
                'P':p_value}
    
    
    assoc_result_list.append(assoc_result)
    
    log.info('[{:3d}/{:3d}] {:10s} {:15s} {:5d} {:.5f}({}) {:.5f}'.format(
                                assoc_result['idx'],
                                len(aa_marker_name_list_aaonly),
                                assoc_result['ID'],
                                assoc_result['residues'],
                                assoc_result['n_obs'],
                                assoc_result['chisq'],
                                assoc_result['df'],
                                assoc_result['P']
                            ))
    
    if idx==120:
        break


# In[174]:


pd.DataFrame(assoc_result_list).to_csv(out+'.result',sep='\t',index=None)

