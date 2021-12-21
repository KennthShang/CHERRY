import numpy as np
import pandas as pd
import os
import Bio
from Bio import SeqIO
import pandas as pd
import subprocess
import argparse
import re


#####################################################################
##########################  Input Params  ###########################
#####################################################################

parser = argparse.ArgumentParser(description='manual to this script')
parser.add_argument('--contigs', type=str, default = 'test_contigs.fa')
parser.add_argument('--len', type=int, default=8000)
parser.add_argument('--gpus', type=int, default = 0)
parser.add_argument('--mode', type=str, default='virus')
parser.add_argument('--model', type=str, default='pretrain')
parser.add_argument('--taxa',  type=str, default='Species')
parser.add_argument('--topk',  type=int, default=3)
inputs = parser.parse_args()

taxa_in = inputs.taxa
taxa_list = taxa_in.split(',')


def check_folder(file_name):
    if not os.path.exists(file_name):
        _ = os.makedirs(file_name)
    else:
        print("folder {0} exist... cleaning dictionary".format(file_name))
        if os.listdir(file_name):
            try:
                _ = subprocess.check_call(f"rm -rf {file_name}", shell=True)
                _ = os.makedirs(file_name)
                print("Dictionary cleaned")
            except:
                print("Cannot clean your folder... permission denied")
                exit(1)

check_folder("input")
check_folder("pred")
check_folder("Split_files")
check_folder("tmp_pred")
check_folder("train_phage/")



#####################################################################
#########################   processing    ###########################
#####################################################################

for record in SeqIO.parse('dataset/nucl.fasta', 'fasta'):
    _ = SeqIO.write(record, 'train_phage/'+record.id, 'fasta')

#####################################################################
#########################  Start Program  ###########################
#####################################################################
# split into sub files
cnt = 0
file_id = 0
records = []
for record in SeqIO.parse(inputs.contigs, 'fasta'):
    if cnt !=0 and cnt%2000 == 0:
        SeqIO.write(records, f"Split_files/contig_{file_id}.fasta","fasta") 
        records = []
        file_id+=1
        cnt = 0
    seq = str(record.seq)
    seq = seq.upper()
    if len(record.seq) > inputs.len:
        records.append(record)
        cnt+=1

SeqIO.write(records, f"Split_files/contig_{file_id}.fasta","fasta")
file_id+=1

# run sub files
for i in range(file_id):
    cmd = f"mv Split_files/contig_{i}.fasta input/"
    try:
        out = subprocess.check_call(cmd, shell=True)
    except:
        print(f"Moving file Error for file contig_{i}")
        exit()


    cmd = "python edge_virus_virus.py"
    try:
        out = subprocess.check_call(cmd, shell=True)
    except:
        print(f"phage_phage Error for file contig_{i}")
        cmd = "rm input/*"
        out = subprocess.check_call(cmd, shell=True)
        exit()

    cmd = f"python edge_virus_prokaryote.py --mode {inputs.mode}"
    try:
        out = subprocess.check_call(cmd, shell=True)
    except:
        print(f"phage_host Error for file contig_{i}")
        cmd = "rm input/*"
        out = subprocess.check_call(cmd, shell=True)
        exit()

    cmd = f"python create_feature.py --mode {inputs.mode}"
    try:
        check_folder("node_feature")
        out = subprocess.check_call(cmd, shell=True)
    except:
        print(f"Pre-trained CNN Error for file contig_{i}")
        cmd = "rm input/*"
        out = subprocess.check_call(cmd, shell=True)
        exit()


    for taxa in taxa_list:
        cmd = f"python multimodal_graph.py --mode {inputs.mode} --taxa {taxa}"
        try:
            out = subprocess.check_call(cmd, shell=True)
        except:
            print(f"multimodal Graph Error for file contig_{i}")
            cmd = "rm input/*"
            out = subprocess.check_call(cmd, shell=True)
            exit()

        cmd = f"python run_Cherry.py --mode {inputs.mode} --model {inputs.model} --gpus {inputs.gpus} --taxa {taxa} --topk {inputs.topk}"
        try:
            out = subprocess.check_call(cmd, shell=True)
        except:
            print("GCN Error for file contig_{i}")
            cmd = "rm input/*"
            out = subprocess.check_call(cmd, shell=True)
            exit()

    # Clean files
    cmd = "rm input/*"
    out = subprocess.check_call(cmd, shell=True)



    # load prediction
    if inputs.mode == 'virus':
        tmp_pred = pd.read_csv(f'tmp_pred/predict_{taxa}.csv')
        name_list = pd.read_csv("name_list.csv")
        prediction = tmp_pred.rename(columns={'contig_names':'idx'})
        contig_to_pred = pd.merge(name_list, prediction, on='idx')
        contig_to_pred.to_csv(f"pred/file_{i}_{taxa}.csv", index = None)

        cmd = "rm name_list.csv"
        out = subprocess.check_call(cmd, shell=True)

        cmd = "rm tmp_pred/*"
        out = subprocess.check_call(cmd, shell=True)
    
    elif inputs.mode == 'prokaryote':
        cmd = f"mv tmp_pred/predict.csv pred/file_{i}.csv"
        out = subprocess.check_call(cmd, shell=True)

if inputs.mode == 'virus':
    prediction_df = []
    for taxa in taxa_list:
        for i in range(file_id):
            prediction_df.append(pd.read_csv(f'pred/file_{i}_{taxa}.csv'))
            prediction_df = pd.concat(prediction_df)
            prediction_df = prediction_df.drop(columns=['idx'])
            prediction_df.to_csv(f'final_prediction_{taxa}.csv', index = None)


elif inputs.mode == 'prokaryote':
    for i in range(file_id):
        prediction_df.append(pd.read_csv(f'pred/file_{i}.csv'))
        prediction_df = pd.concat(prediction_df)
        prediction_df.to_csv(f'final_prediction.csv', index = None)


