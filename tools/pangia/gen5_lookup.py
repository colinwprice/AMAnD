import os
import numpy as np
from scipy.fftpack import fft
from Bio import SeqIO
import pandas as pd
from annoy import AnnoyIndex
import argparse

parser = argparse.ArgumentParser(description="Annoy Look Up Table for BSAT targets")
parser.add_argument('-a', required=True, help="File Location of Annoy Index")
parser.add_argument('-n',type=int, default=40,help="N number of neighbors to query.")
parser.add_argument('-t', required=True, help="Taxid dictionary for Annoy Index (pickled)")
parser.add_argument('-r', required=True, help="Fasta file of reads")
parser.add_argument('-o', default=None, help="Output file name (Default: Inherit from reads file name.)")

myargs=parser.parse_args()
index_file = myargs.a
reads_file = myargs.r
taxid_file = myargs.t
neighbors = myargs.n

reads = [s for s in SeqIO.parse(reads_file,'fasta') if len(s.seq) >= 124]

#integer encoder for vectorization
nucl_dict = {'A':0,'C':1,'G':1,'T':0}

#encoding function
def encode(x):
    '''this should readily lend itself to cythonization for encoding speed up'''
    v = np.array([nucl_dict[n] for n in x])
    #this could be replaced with equivalent call to pyFFTW
    yf = fft(v)
    return 2.0/len(x) * np.abs(yf[0:len(x)//2])

def tally(x_a,x_b):
    '''function to evaluate nearest neighbors'''
    table = pd.DataFrame({'taxid':x_a,'distance':x_b}).pivot_table(values='distance',index='taxid',aggfunc=np.min)
    if table[table['distance'] == 0.0].empty or table[table['distance'] == 0.0].shape[0] > 1:
        return 'Non-specific'
    else:
        return table[table['distance'] == 0.0].index.unique()[0]

#dependent on min read length 124
vec_size = 62
t = AnnoyIndex(vec_size, 'euclidean')

#prefault flag for MEM MAPPED index
t.load(index_file,prefault=True)
taxids = pd.read_pickle(taxid_file)

#TODO: Multiprocess this look up
results_taxid = []
results_distance = []
for s in reads:
    a,b = t.get_nns_by_vector(encode(str(s.seq)[:124]), neighbors, search_k=-1, include_distances=True)
    results_taxid.append(taxids.iloc[a].values)
    results_distance.append(b)

results_tally = [tally(a,b) for a,b in zip(results_taxid,results_distance)]

output = pd.DataFrame({'taxid':results_tally,'read_id':[r.id for r in reads]})

print(output['taxid'].value_counts())

if myargs.o:
    output_file = myargs.o
else:
    output_file = reads_file.rsplit('.',maxsplit=1)[0] + '_classifications.tsv'

output.to_csv(output_file,sep='\t',index=False)

