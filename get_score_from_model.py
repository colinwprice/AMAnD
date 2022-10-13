# -*- coding: utf-8 -*-

import argparse
import os
import json
from tensorflow import keras
import pandas as pd
import glob
import itertools
import numpy as np
from datetime import datetime


def is_dir(string):
    if os.path.isdir(string):
        return string
    else:
        raise NotADirectoryError(string)
        
def is_file(string):
	if os.path.isfile(string):
		return string
	else:
		raise ValueError("Not a file")        

def run_pangia(paired, dir_to_score):    
    if paired: #Paired end reads
        # PanGIA control
        for i in os.listdir(dir_to_score):
            if i.endswith(('R1.fa', 'R1.fasta', 'R1.fastq', 'R1.fq')):
                split_str = i.split('R1.')
                out_str = split_str[0]+'_CONTROL_PanGIA_out'
                
                R1 = i
                R2 = split_str[0]+'_R2.'+split_str[1]
                
                command_string = "python tools/pangia/pangia.py -db tools/pangia_db/*.fa -i "+str(dir_to_score)+str(R1)+" "+str(dir_to_score)+str(R2)+" -t "+str(t)+" -rm minimap2 -st illumina -da -o "+str(dir_to_score)+"/"+str(out_str) 
                os.system(command_string)
    
    
    else:
        for i in os.listdir(dir_to_score):
            if i.endswith(('.fa', '.fasta', '.fastq', '.fq')):
                split_str = i.split('.')
                out_str = split_str[0]+'_CONTROL_PanGIA_out'
                
                command_string = "python tools/pangia/pangia.py -db tools/pangia_db/*.fa -i "+str(dir_to_score)+"/"+str(i)+" -t "+str(t)+" -rm minimap2 -st illumina -da -se -o "+str(dir_to_score)+"/"+str(out_str) 
                os.system(command_string)

def run_jellyfish(kmer_size, dir_to_score):
    for i in os.listdir(dir_to_score):
            if i.endswith(('.fa', '.fasta', '.fastq', '.fq')):
                split_str = i.split('.')
                out_str = split_str[0]+'_TEST_kmer_out'
                
                command_string = "jellyfish count -m "+str(kmer_size)+" -s 200M -t "+str(t)+" -C "+str(dir_to_score)+"/"+str(i)+" -o "+str(dir_to_score)+"/"+str(out_str)+".jf"
                os.system(command_string)
                command_string = "jellyfish dump "+str(dir_to_score)+"/"+str(out_str)+".jf -c -t -o "+str(dir_to_score)+"/"+str(out_str)+".tsv"
                os.system(command_string)
                os.system("rm "+str(dir_to_score)+"/"+str(out_str)+".jf")

def merge_kmer_type(d):
    tsv_dfs = []
    for file in os.listdir(d):
        if file.endswith('kmer_out.tsv'):
            tsv_df = pd.read_table(d+'/'+file, index_col=0, header=0, names=[file])
            tsv_df_T = tsv_df.T
            tsv_dfs.append(tsv_df_T)
    
    master_df = pd.concat(tsv_dfs, copy=False).fillna(0)
    return master_df    
    


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Load a model and score some samples")
    
    parser.add_argument('-t', '--threads', type=int, default=4, required=False, help='number of threads to employ where applicable (default=4)')
    parser.add_argument('-d', '--deepsvdd_model', default='./models/PanGIA_svdd_model', help='which model to test on (default=./models/PanGIA_svdd_model)')
    parser.add_argument('-m', '--metadata', type=is_file, default='./models/PanGIA_svdd_model_metadata.json', help='metadata file that acompanies name of model (default=./models/PanGIA_svdd_model_metadata.json)')
    parser.add_argument('-s', '--dir_to_score', type=is_dir, default='./test_samples', help='in which directory are the files to score stored (default=test_samples)')
    parser.add_argument('-p', '--paired', type = bool, required=False, default=False, help='Are samples paired? Paired samples must have "R1" and "R2" before the .fastq/.fastq/.fa/.fq (default=False)')
    # infer length of truncation from metadata file
    # can infer kmer for jellyfish from metadata size
    parser.parse_args()
    
    args = parser.parse_args()
    t = args.threads
    deepsvdd_model = args.deepsvdd_model
    metadata_file = args.metadata
    paired = args.paired
    dir_to_score = args.dir_to_score
    
    
    with open(metadata_file) as json_file:
        metadata = json.load(json_file)
    
    model_type = metadata['type']
    truncation = metadata['truncation']
    
    df_master=None
    inds_to_write=None
    
    if model_type == 'kmer':
        kmer_size = metadata['k_size']
        run_jellyfish(kmer_size, dir_to_score)
        df_master = merge_kmer_type(dir_to_score)
        inds_to_write=list(df_master.index)
        inds_to_write = [x.split('_TEST_kmer_out')[0] for x in inds_to_write]
        df_master = df_master.to_numpy(dtype='float32')
        
        
    if model_type == 'pangia':    
        run_pangia(paired, dir_to_score)
        files = glob.glob(dir_to_score+'/*/*.tsv')
    
        metric_names=['GEN_ri', 'GEN_rnr', 'SPE_ri', 'SPE_rnr', 'STR_ri', 'STR_rnr']
        with open('tools/species_taxids_from_pangia_db.txt') as nf:
            species_names = nf.readlines()
            species_names = [name.rstrip() for name in species_names]
        
        master_pangia_df_header = [n[0]+"-"+n[1] for n in list(itertools.product(species_names, metric_names))]
        master_pangia_df = pd.DataFrame(columns=master_pangia_df_header)  
        
        samples_df = pd.DataFrame()
        
        for file in files:
            this_df = pd.read_table(file, usecols=['LEVEL','TAXID','GEN_ri', 'GEN_rnr', 'SPE_ri', 'SPE_rnr', 'STR_ri', 'STR_rnr'])
            taxids_this_df = list(this_df['TAXID'])
            this_df.set_index('TAXID', inplace=True)
            this_df = this_df[this_df['LEVEL']=='species']
            this_df = this_df.drop('LEVEL', axis=1)
            this_df_dict = this_df.to_dict()
            
            cat_df = pd.DataFrame(columns=['NAME']+list(itertools.product(taxids_this_df, metric_names)))
            cat_df['NAME'] = file.split('/')[-1].split('.')[0]
            
            cat_dict = {}
            
            for k, v in this_df_dict.items():
                for k_prime, v_prime in v.items():
                    cat_dict[str(int(k_prime))+'-'+str(k)] = [float(v_prime)]
                    
            cat_df = pd.DataFrame.from_dict(cat_dict)
            cat_df['NAME'] = file.split('/')[-1].split('.')[0]
            cat_df.set_index('NAME', inplace=True, drop=True)
           
            samples_df = pd.concat([samples_df, cat_df])
        
        df_master = pd.concat([master_pangia_df, samples_df])
        
        rank_dict={}
        
        for i, fname in enumerate(df_master.index):
            list_b = df_master.iloc[i,:].tolist()
            list_b = np.array([float(x) for x in list_b])
            top_k_indices = np.sort(np.argpartition(-list_b, truncation)[:truncation])
            rank_dict[fname] = list_b[top_k_indices]
    
    
        rank_ordered_pangia=pd.DataFrame.from_dict(rank_dict, orient="index")
        rank_ordered_pangia.reset_index(drop=True, inplace=True)
        rank_ordered_pangia = rank_ordered_pangia.to_numpy(dtype='float32')
        inds_to_write=list(df_master.index)
        df_master = rank_ordered_pangia
        
    
    
    
    model = keras.models.load_model(deepsvdd_model)
    predictions = model.predict(df_master)
    
    print(inds_to_write)
    print(predictions)
    
    d_to_write = {'Index': inds_to_write, 'Predictions': predictions}
    df_to_write = pd.DataFrame.from_dict(d_to_write)
    model_name = deepsvdd_model.split('/')[-1]
    now=datetime.now()
    now_str = now.strftime("%m-%d_%H-%M-%S")
    df_to_write.to_csv('./test_sample_scores/'+now_str+'_'+model_name+'scores.csv', sep=',', index=False)
       
    
    

    
    
    
    
    
    