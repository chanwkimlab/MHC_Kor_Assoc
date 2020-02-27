import subprocess
import re
import pandas as pd
import numpy as np

from path_configure import *

def run_subprocess(command,quiet=False,dry=False):
    print("------{}-----".format("RUN"))
    print(command)
    if dry:
        return
    process = subprocess.Popen(command.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    if quiet==False:
        print("------{}-----".format("ERROR"))
        print(stderr.decode())
        print("------{}-----".format("OUTPUT"))
        print(stdout.decode())
    return stdout.decode(),stderr.decode()


def replace_characters(str_input,character_dict,catch_exception=True):
    for to_replace,value in character_dict.items():
        if catch_exception:
            try:
                str_input=str_input.replace(to_replace,value)
            except:
                pass
        else:
            str_input=str_input.replace(to_replace,value)
    return str_input  

def parse_plink_assoc(file_name):
    with open(file_name,'r') as f:
        lines=f.readlines()
        #lines=[line.strip().split(' ') for line in lines]
        header=re.split('\s+',lines[0].strip())
        lines=[re.split('\s+',line.strip()) for line in lines[1:]]
    ret=pd.DataFrame(lines,columns=header).replace('NA',np.nan)
    if file_name.split('.')[-1]=='assoc':
        return ret.astype({'BP':int,'A1':str,'A2':str,'F_A':float,'F_U':float,'P':float,'OR':float})
    elif file_name.split('.')[-1]=='qassoc':
        return ret.astype({'CHR':int,'SNP':str,'BP':int,'NMISS':int,'BETA':float,'SE':float,'R2':float,'T':float,'P':float})     
    elif file_name.split('.')[-1]=='logistic':
        return ret.astype({'CHR':int,'SNP':str,'BP':int,'A1':str,'NMISS':int,'OR':float,'STAT':float,'P':float})        
    elif file_name.split('.')[-1]=='linear':
        return ret.astype({'CHR':int,'SNP':str,'BP':int,'A1':str,'NMISS':int,'BETA':float,'STAT':float,'P':float})        
    else:
        raise