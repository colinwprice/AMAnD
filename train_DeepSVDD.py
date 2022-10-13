# -*- coding: utf-8 -*-

import pandas as pd
from pyod.models.deep_svdd import DeepSVDD
import json
import os
import numpy as np
import argparse
import sys

def is_file(string):
	if os.path.isfile(string):
		return string
	else:
		raise ValueError("Not a file")

def train_SVDD(params):
    os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2' 
    i = params[0]
    hidden_layer = params[1]
    activation_func = params[2]
    output_func = params[3]
    #loss = tf.keras.losses.get(params[4])
    batch_size = params[5]
    dropout = params[6]
    l2 = params[7]
    use_AE = params[8]
    control_set = params[9]
    truncation = params[10]
    anomaly_set = params[11]
    epochs_to_train = params[12]
    clf2 = DeepSVDD(hidden_activation = activation_func,
                    output_activation = output_func,
                    epochs=epochs_to_train,
                    batch_size=batch_size,
                    dropout_rate = dropout,
                    l2_regularizer=l2,
                    preprocessing=True,
                    verbose=1,
                    contamination=.00000001,
                    hidden_neurons=hidden_layer,
                    random_state=42,
                    use_ae=use_AE)
    
    clf2.fit(control_set)
    if anomaly_set is None:
        anomaly_scores=None
    else:
        anomaly_scores = clf2.decision_function(anomaly_set)
    control_scores = clf2.decision_function(control_set)
    history = clf2.history_
    
    print('Finished Model #'+str(i))
    
    return (anomaly_scores, control_scores, history, clf2.model_, truncation)
    
    

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description="Train a PanGIA and a kmer model from the .pkl data objects")
    
    parser.add_argument('-k', '--kmer_pkl', type=is_file, required=False, default='all_kmer_outputs.pkl', help='name of compiled compressed kmer object (default=all_kmer_outputs.pkl)')
    parser.add_argument('-p', '--pangia_pkl', type=is_file, required=False, default='all_pangia_outputs.pkl', help='name of compiled compressed PanGIA object (default=all_pangia_outputs.pkl)')
    parser.add_argument('-e', '--epochs_to_train', type=int, required=False, default=200, help='how many epochs should each model parameterization be trained (default=200)')
    parser.add_argument('-a', '--anomaly_present', type=bool, required=False, default=None, help='set to True if anomaly data provided (default=None)')
    # THREADS ARG THREADS ARG THREADS ARG??????#
    args = parser.parse_args()
    kmer_file = args.kmer_pkl
    pangia_file = args.pangia_pkl
    epochs_to_train = args.epochs_to_train
    has_anomaly = args.anomaly_present
    
    kmer_pickle = pd.read_pickle(kmer_file)
    pangia_pickle = pd.read_pickle(pangia_file)
    
    if has_anomaly:
        kmer_anomaly=kmer_pickle.loc[kmer_pickle.LABEL==1]
        kmer_anomaly.drop('LABEL', axis=1, inplace=True)
    else:
        kmer_anomaly=None
    kmer_control=kmer_pickle.loc[kmer_pickle.LABEL==0]
    
    if has_anomaly:
        pangia_anomaly=pangia_pickle.loc[pangia_pickle.LABEL==1]
        pangia_anomaly.drop('LABEL', axis=1, inplace=True)
    else:
        pangia_anomaly=None
    pangia_control=pangia_pickle.loc[pangia_pickle.LABEL==0]
    
    
    kmer_control.drop('LABEL', axis=1, inplace=True)
    pangia_control.drop('LABEL', axis=1, inplace=True)
    
    
    tasks = []
    for truncation in [40, 60, 120]:
        rank_dict = {}
        
        for i, fname in enumerate(pangia_control.index):
            list_b = pangia_control.iloc[i,:].tolist()
            list_b = np.array([float(x) for x in list_b])
            top_k_indices = np.sort(np.argpartition(-list_b, truncation)[:truncation])
            rank_dict[fname] = list_b[top_k_indices]
    
    
        rank_ordered_pangia_control=pd.DataFrame.from_dict(rank_dict, orient="index")
        rank_ordered_pangia_control.reset_index(drop=True, inplace=True)
        rank_ordered_pangia_control = rank_ordered_pangia_control.to_numpy(dtype='float32')
    
        if has_anomaly:
            rank_dict = {}
            
            for i, fname in enumerate(pangia_anomaly.index):
                list_b = pangia_anomaly.iloc[i,:].tolist()
                list_b = np.array([float(x) for x in list_b])
                top_k_indices = np.sort(np.argpartition(-list_b, truncation)[:truncation])
                rank_dict[fname] = list_b[top_k_indices]
                
            rank_ordered_pangia_anomaly=pd.DataFrame.from_dict(rank_dict, orient="index")
            
            rank_ordered_pangia_anomaly.reset_index(drop=True, inplace=True)
            rank_ordered_pangia_anomaly = rank_ordered_pangia_anomaly.to_numpy(dtype='float32')
        else:
            rank_ordered_pangia_anomaly=None
    
        i=0
        for hl in [[10, 5], [20, 10], [30, 15], [50, 20]]:
            for af in ['relu']:
                for of in ['sigmoid']:
                    for bs in [1]:
                        for do in [0,.1,.2]:
                            for l2 in [0,.1,.2]:
                                for ua in [True, False]:
                                    for ls in ['mse']:
                                        tasks.append((i, hl, af, of, ls, bs, do, l2, ua, rank_ordered_pangia_control, truncation, rank_ordered_pangia_anomaly, epochs_to_train))    
                                        i+=1
                                        
                                        
    best_training_convergence_pangia=sys.maxsize
    best_model_pangia=None
    best_ascores_pangia=None
    best_cscores_pangia=None
    best_threshold_pangia=None
    best_truncation_pangia=None
    
    for task in tasks:
        result = train_SVDD(task)
        # write the best result at the end to a model file
        a_scores = result[0] #(anomaly_scores, control_scores, history, clf2.model_)
        c_scores = result[1]
        history = result[2]
        model = result[3] # model.save('path to location')
        truncation = result[4]
        # best training convergence is average delta between last 20 or so epochs
        val_loss = history['val_loss']
        train_loss = history['loss']
        
        cutoff = min(20, int(epochs_to_train))
        val_loss_cutoff = np.array(val_loss[-cutoff:])
        train_loss_cutoff = np.array(train_loss[-cutoff:])
        
        final_state = np.subtract(val_loss_cutoff, train_loss_cutoff)
        final_state_avg = np.mean(final_state)
        
        if final_state_avg < best_training_convergence_pangia:
            best_training_convergence_pangia = final_state_avg
            best_model_pangia=model
            best_ascores_pangia=a_scores
            best_cscores_pangia=c_scores
            best_truncation_pangia=truncation
            
    pangia_dict_to_dump = {}
    if has_anomaly:
        pangia_dict_to_dump['anomaly_scores'] = best_ascores_pangia.tolist()
    pangia_dict_to_dump['control_scores'] = best_cscores_pangia.tolist()
    pangia_dict_to_dump['truncation'] = best_truncation_pangia
    pangia_dict_to_dump['type'] = 'pangia'
    
    # write model file
    best_model_pangia.save('./models/PanGIA_svdd_model')
    # write metadata file
    with open("./models/PanGIA_svdd_model_metadata.json", "w") as outfile:
        json.dump(pangia_dict_to_dump, outfile)
    
    tasks = []
    i=0
    truncation = kmer_pickle.shape[1]
    for hl in [[10, 5], [20, 10], [30, 15]]:
        for af in ['relu']:
            for of in ['sigmoid']:
                for bs in [1]:
                    for do in [0,.1,.2]:
                        for l2 in [0,.1,.2]:
                            for ua in [True, False]:
                                for ls in ['mse']:
                                    tasks.append((i, hl, af, of, ls, bs, do, l2, ua, kmer_control, truncation, kmer_anomaly, epochs_to_train))
                                    i+=1
    
    
    best_training_convergence_kmer=sys.maxsize
    best_model_kmer=None
    best_ascores_kmer=None
    best_cscores_kmer=None   
    best_threshold_kmer=None
    best_truncation_kmer=None
    
    for task in tasks:
        result = train_SVDD(task)
        # write the best result at the end to a model file
        a_scores = result[0] #(anomaly_scores, control_scores, history, clf2.model_)
        c_scores = result[1]
        history = result[2]
        model = result[3]
        truncation = result[4]
    
        val_loss = history['val_loss']
        train_loss = history['loss']
        
        cutoff = min(20, int(epochs_to_train))
        val_loss_cutoff = np.array(val_loss[-cutoff:])
        train_loss_cutoff = np.array(train_loss[-cutoff:])
        
        final_state = np.subtract(val_loss_cutoff, train_loss_cutoff)
        final_state_avg = np.mean(final_state)
        
        if final_state_avg < best_training_convergence_kmer:
            best_training_convergence_kmer = final_state_avg
            best_model_kmer=model
            best_ascores_kmer=a_scores
            best_cscores_kmer=c_scores
            best_truncation_kmer=truncation
    
    kmer_dict_to_dump = {}
    if has_anomaly:
        kmer_dict_to_dump['anomaly_scores'] = best_ascores_kmer.tolist()
    kmer_dict_to_dump['control_scores'] = best_cscores_kmer.tolist()
    kmer_dict_to_dump['truncation'] = best_truncation_kmer
    kmer_dict_to_dump['type'] = 'kmer'
    kmer_dict_to_dump['k_size'] = str(len(kmer_pickle.columns[0]))
    
    # write model file
    best_model_kmer.save('./models/kmer_svdd_model')
    # write metadata file
    with open("./models/kmer_svdd_model_metadata.json", "w") as outfile:
        json.dump(kmer_dict_to_dump, outfile)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    