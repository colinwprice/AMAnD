# -*- coding: utf-8 -*-

import argparse
import os
import glob
import pandas as pd
import itertools

def is_file(string):
	if os.path.isfile(string):
		return string
	else:
		raise ValueError("Not a file")

def merge_kmer_type(label, d):
    tsv_dfs = []
    for file in os.listdir(d):
        if file.endswith('kmer_out.tsv'):
            tsv_df = pd.read_table(d+'/'+file, index_col=0, header=0, names=[file])
            tsv_df_T = tsv_df.T
            tsv_dfs.append(tsv_df_T)
    
    master_df = pd.concat(tsv_dfs, copy=False).fillna(0)
    master_df['LABEL'] = label
    return master_df

# file structure to parse:
# the .tsvs of kmer counts will just be inside the control_dir and anomaly_dir
# just give them a 'label' column and be on with it 0=control or 1=anomaly

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Second step in the pipeline is to assign the binary label and merge/compress all outputs")
    
    parser.add_argument('-a', '--anomaly_dir', type = str, required=False, default='anomaly_dir', help='Directory of anomalous sample outputs (default=anomaly_dir)')
    parser.add_argument('-c', '--control_dir', type = str, required=False, default='control_dir', help='Directory of typical, normal, control sample outputs (default=control_dir)')
    parser.add_argument('-ok', '--output_kmer', type = str, required=False, default='all_kmer_outputs', help='The name of output kmer file (default=all_kmer_outputs)')
    parser.add_argument('-op', '--output_pangia', type = str, required=False, default='all_pangia_outputs', help='The name of the output pangia file (default=all_pangia_outputs)')
    parser.add_argument('-n', '--names_file', type=is_file, default='tools/species_taxids_from_pangia_db.txt', help='file containing names of the pangia database taxa (default=tools/species_taxids_from_pangia_db.txt)')
    parser.add_argument('-m', '--make_names_file', type=bool, default=False, help='remake names file? (should not be updated between runs and only if PanGIA DB has changed) (default=False)')
    
    args = parser.parse_args()
    
    anomaly_dir = args.anomaly_dir
    control_dir = args.control_dir
    out_kmer = args.output_kmer
    out_pangia = args.output_pangia
    names_file = args.names_file
    remake = args.make_names_file
    
    
    print('Attempting merge of kmers.')
    if not len(os.listdir(anomaly_dir))==0:
        kmer_df_anom = merge_kmer_type(1, anomaly_dir)
        
    kmer_df_con = merge_kmer_type(0, control_dir)
    
    if not len(os.listdir(anomaly_dir))==0:
        kmer_df_master = pd.concat([kmer_df_con, kmer_df_anom], copy=False)
    else:
        kmer_df_master=kmer_df_con
    kmer_df_master.to_pickle(out_kmer+".pkl")
    
    
    # from a names file - create n names by 6 different features of pangia
    # this creates the list of all columns which will hopefully make the merge/cat operation easier
    
    if remake:   
        for file in os.listdir('tools/pangia_db/'):
            if file.endswith('.headers'):
                with open('tools/pangia_db/'+file) as headers_file:
                    lines = headers_file.readlines()
                    lines = [line.split('|')[2].split('.')[0] for line in lines]
                    lines = list(set(lines))
                    
                    writefile = open('tools/species_taxids_from_pangia_db.txt', 'w')
                    for line in lines:
                        writefile.write(f"{line}\n")
                    writefile.close()
        
        
    files = glob.glob('*_dir/*/*.tsv')
    
    metric_names=['GEN_ri', 'GEN_rnr', 'SPE_ri', 'SPE_rnr', 'STR_ri', 'STR_rnr']
    with open(names_file) as nf:
        species_names = nf.readlines()
        species_names = [name.rstrip() for name in species_names]
    
    master_pangia_df_header = [n[0]+"-"+n[1] for n in list(itertools.product(species_names, metric_names))]
    master_pangia_df = pd.DataFrame(columns=master_pangia_df_header)  
    
    samples_df = pd.DataFrame()
    
    print('Attempting merge of pangia reports.')
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
        
        if file.startswith('control_dir'):
            cat_df['LABEL']=0
        elif file.startswith('anomaly_dir'):
            cat_df['LABEL']=1
        else:
            print('WARNING: Ambiguous file detected, unable to assign label')
            cat_df['LABEL']='NA'
       
        samples_df = pd.concat([samples_df, cat_df])
    
     
    master_pangia_df = pd.concat([master_pangia_df, samples_df])
    master_pangia_df.fillna(0, inplace=True)
    kmer_df_master.to_pickle(out_pangia+".pkl")








