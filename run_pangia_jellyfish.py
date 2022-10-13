# -*- coding: utf-8 -*-

import argparse
import os

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="First step in the pipeline is to run PanGIA on all the provided reads")
    
    parser.add_argument('-a', '--anomaly_dir', type = str, required=False, default='anomaly_dir', help='Directory of anomalous samples. Can be left empty if none exist (default=anomaly_dir)')
    parser.add_argument('-c', '--control_dir', type = str, required=False, default='control_dir', help='Directory of typical, normal, control samples (default=control_dir)')
    parser.add_argument('-t', '--threads', type = int, required=False, default=4, help='Threads to use in processing (default=4)')
    parser.add_argument('-p', '--paired', type = bool, required=False, default=False, help='Are samples paired? Paired samples must have "R1" and "R2" before the .fastq/.fastq/.fa/.fq (default=False)')
    parser.add_argument('-k', '--kmer_size', type=int, required=False, default=4, help='Specify size of kmers to count (default=4)')
    
    args = parser.parse_args()
    
    anomaly_dir = args.anomaly_dir
    control_dir = args.control_dir
    t = args.threads
    paired = args.paired
    kmer_size = args.kmer_size
    
    
    if paired: #Paired end reads
        # PanGIA control
        for i in os.listdir(control_dir):
            if i.endswith(('R1.fa', 'R1.fasta', 'R1.fastq', 'R1.fq')):
                split_str = i.split('R1.')
                out_str = split_str[0]+'_CONTROL_PanGIA_out'
                
                R1 = i
                R2 = split_str[0]+'_R2.'+split_str[1]
                
                command_string = "python tools/pangia/pangia.py -db tools/pangia_db/*.fa -i "+str(control_dir)+str(R1)+" "+str(control_dir)+str(R2)+" -t "+str(t)+" -rm minimap2 -st illumina -da -o "+str(control_dir)+"/"+str(out_str) 
                os.system(command_string)
    
        # PanGIA anomaly
        for i in os.listdir(anomaly_dir):
            if i.endswith(('R1.fa', 'R1.fasta', 'R1.fastq', 'R1.fq')):
                split_str = i.split('R1.')
                out_str = split_str[0]+'_ANOMALY_PanGIA_out'
                
                R1 = i
                R2 = split_str[0]+'_R2.'+split_str[1]
                
                command_string = "python tools/pangia/pangia.py -db tools/pangia_db/*.fa -i "+str(anomaly_dir)+str(R1)+" "+str(anomaly_dir)+str(R2)+" -t "+str(t)+" -rm minimap2 -st illumina -da -o "+str(anomaly_dir)+"/"+str(out_str) 
                os.system(command_string)
    
    
    
    else:
        for i in os.listdir(control_dir):
            if i.endswith(('.fa', '.fasta', '.fastq', '.fq')):
                split_str = i.split('.')
                out_str = split_str[0]+'_CONTROL_PanGIA_out'
                
                command_string = "python tools/pangia/pangia.py -db tools/pangia_db/*.fa -i "+str(control_dir)+"/"+str(i)+" -t "+str(t)+" -rm minimap2 -st illumina -da -se -o "+str(control_dir)+"/"+str(out_str) 
                os.system(command_string)
    
    
        for i in os.listdir(anomaly_dir):
            if i.endswith(('.fa', '.fasta', '.fastq', '.fq')):
                split_str = i.split('.')
                out_str = split_str[0]+'_ANOMALY_PanGIA_out'
                
                
                command_string = "python tools/pangia/pangia.py -db tools/pangia_db/*.fa -i "+str(anomaly_dir)+"/"+str(i)+" -t "+str(t)+" -rm minimap2 -st illumina -da -se -o "+str(anomaly_dir)+"/"+str(out_str) 
                os.system(command_string)
                
     
        
    # kmer counting can be done independently of pairedness
    for i in os.listdir(control_dir):
        if i.endswith(('.fa', '.fasta', '.fastq', '.fq')):
            split_str = i.split('.')
            out_str = split_str[0]+'_CONTROL_kmer_out'
            
            command_string = "jellyfish count -m "+str(kmer_size)+" -s 200M -t "+str(t)+" -C "+str(control_dir)+"/"+str(i)+" -o "+str(control_dir)+"/"+str(out_str)+".jf"
            os.system(command_string)
            command_string = "jellyfish dump "+str(control_dir)+"/"+str(out_str)+".jf -c -t -o "+str(control_dir)+"/"+str(out_str)+".tsv"
            os.system(command_string)
            os.system("rm "+str(control_dir)+"/"+str(out_str)+".jf")
    
    
    for i in os.listdir(anomaly_dir):
        if i.endswith(('.fa', '.fasta', '.fastq', '.fq')):
            split_str = i.split('.')
            out_str = split_str[0]+'_ANOMALY_kmer_out'
            
            command_string = "jellyfish count -m "+str(kmer_size)+" -s 200M -t "+str(t)+" -C "+str(anomaly_dir)+"/"+str(i)+" -o "+str(anomaly_dir)+"/"+str(out_str)+".jf"
            os.system(command_string)
            command_string = "jellyfish dump "+str(anomaly_dir)+"/"+str(out_str)+".jf -c -t -o "+str(anomaly_dir)+"/"+str(out_str)+".tsv"
            os.system(command_string)
            os.system("rm "+str(anomaly_dir)+"/"+str(out_str)+".jf")
                
            
    
    
    
            
    
