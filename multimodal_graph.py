import os
import sys
import Bio
import logging
import argparse
import subprocess
import scipy as sp
import numpy as np
import pandas as pd
import pickle as pkl
import networkx as nx
import scipy.stats as stats
import scipy.sparse as sparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord



################################################################################
###############################  create folder   ###############################
################################################################################

phage_phage_ntw = "out/phage_phage.ntw"
phage_host_ntw = "out/phage_host.ntw"
parser = argparse.ArgumentParser(description='manual to this script')
parser.add_argument('--mode', type=str, default = 'virus')
inputs = parser.parse_args()

def check_folder(file_name):
    if not os.path.exists(file_name):
        _ = os.makedirs(file_name)
    else:
        print("folder {0} exist... cleaning dictionary".format(file_name))
        if os.listdir(file_name):
            try:
                _ = subprocess.check_call("rm -rf {0}".format(file_name), shell=True)
                _ = os.makedirs(file_name)
                print("Dictionary cleaned")
            except:
                print("Cannot clean your folder... permission denied")
                exit(1)

check_folder("GCN_data")

################################################################################
############################  Edge construction  ###############################
################################################################################

# Add virus-virus edges
G = nx.Graph()
with open(phage_phage_ntw) as file_in:
    for line in file_in.readlines():
        tmp = line[:-1].split(",")
        node1 = tmp[0].split('.')[0]
        node2 = tmp[1].split('.')[0]
        G.add_edge(node1, node2, weight = 1)


# Add blastn edges
with open(phage_host_ntw) as file_in:
    for line in file_in.readlines():
        tmp = line[:-1].split(",")
        node1 = tmp[0].split('.')[0]
        node2 = tmp[1].split('.')[0]
        G.add_edge(node1, node2, weight = 1)



bacteria_df = pd.read_csv('dataset/prokaryote.csv')
virus_df = pd.read_csv('dataset/virus.csv')

bacteria_list = os.listdir('prokaryote/')
bacteria_list = [name.split('.')[0] for name in bacteria_list]

# add crispr edges
species2bacteria = {bacteria_df[bacteria_df['Accession'] == item]['Species'].values[0]: item for item in bacteria_list}
crispr_pred = pkl.load(open('out/crispr_pred.dict', 'rb'))
for virus, host in crispr_pred.items():
    if host in species2bacteria:
        G.add_edge(virus, species2bacteria[host])

# add dataset edges
for bacteria in bacteria_list:
    species = bacteria_df[bacteria_df['Accession'] == bacteria]['Species'].values[0]
    phage_list = virus_df[virus_df['Species'] == species]['Accession'].values
    for phage in phage_list:
        if phage in G.nodes():
            G.add_edge(bacteria, phage, weight = 1)

################################################################################
############################  Nodes construction ###############################
################################################################################

virus2id      = pkl.load(open("node_feature/virus.dict",'rb'))
virusF        = pkl.load(open("node_feature/virus.F",'rb'))
prokaryote2id = pkl.load(open("node_feature/prokaryote.dict",'rb'))
prokaryoteF   = pkl.load(open("node_feature/prokaryote.F",'rb'))
 
if inputs.mode == 'virus':
    test_virus2id      = pkl.load(open("node_feature/test_virus.dict",'rb'))
    test_virusF        = pkl.load(open("node_feature/test_virus.F",'rb'))
    test_prokaryote2id = {}
elif inputs.mode == 'prokaryote':
    test_prokaryote2id = pkl.load(open("node_feature/test_prokaryote.dict",'rb'))
    test_prokaryoteF   = pkl.load(open("node_feature/test_prokaryote.F",'rb'))
    test_virus2id      = pkl.load(open("node_feature/test_virus.dict",'rb'))
    test_virusF        = pkl.load(open("node_feature/test_virus.F",'rb'))
    test_virus2id = {}
else:
    print("mode error")


node_feature = []
for node in G.nodes():
    # if prokaryote node
    if node in prokaryote2id.keys():
        node_feature.append(prokaryoteF[prokaryote2id[node]])
    # if virus node
    elif node in virus2id.keys():
        node_feature.append(virusF[virus2id[node]])
    # if test virus node
    elif node in test_virus2id.keys():
        node_feature.append(test_virusF[test_virus2id[node]])
    # if test prokaryote node
    elif node in test_prokaryote2id.keys():
        node_feature.append(test_prokaryoteF[test_prokaryote2id[node]])
    else:
        print(f"node error {node}")
        exit()

node_feature = np.array(node_feature)

################################################################################
############################  Label construction ###############################
################################################################################

crispr_pred = pkl.load(open('out/crispr_pred.dict', 'rb'))
virus_pred = pkl.load(open('out/virus_pred.dict', 'rb'))
virus_df = pd.read_csv("dataset/virus.csv")
prokaryote_df = pd.read_csv("dataset/prokaryote.csv")


idx = 0
test_id = {}
node2label = {}
cnt = 0
for node in G.nodes():
    # if test virus node
    if "cherry" in node:
        neighbor_label = []
        for _, neighbor in G.edges(node):
            if neighbor in virus2id.keys():
                virus_label = virus_df[virus_df['Accession'] == neighbor]['Species'].values[0]
                neighbor_label.append(virus_label)
            elif neighbor in prokaryote2id.keys():
                prokaryote_label = prokaryote_df[prokaryote_df['Accession'] == neighbor]['Species'].values[0]
                neighbor_label.append(prokaryote_label)
        # subgraph
        if len(set(neighbor_label)) == 1:
            node2label[node] = neighbor_label[0]
            test_id[node] = 1
        # CRISPR
        elif node in crispr_pred:
            node2label[node] = prokaryote_df[prokaryote_df['Accession'] == crispr_pred[node]]['Species'].values[0]
            test_id[node] = 1
        elif node in virus_pred:
            node2label[node] = virus_df[virus_df['Accession'] == virus_pred[node]]['Species'].values[0]
            test_id[node] = 1
        # unlabelled
        else:
            node2label[node] = 'unknown'
            test_id[node] = 2
    # if phage or host node
    elif node in prokaryote2id.keys():
        prokaryote_label = prokaryote_df[prokaryote_df['Accession'] == node]['Species'].values[0]
        node2label[node] = prokaryote_label
        test_id[node] = 0
    elif node in test_prokaryote2id.keys():
        prokaryote_label = prokaryote_df[prokaryote_df['Accession'] == node]['Species'].values[0]
        node2label[node] = prokaryote_label
        test_id[node] = 0
    elif node in virus2id.keys():
        virus_label = virus_df[virus_df['Accession'] == node]['Species'].values[0]
        node2label[node] = virus_label
        test_id[node] = 0
    else: 
        print("Error: " + node)
    idx += 1


# check subgraph situation 1
for sub in nx.connected_components(G):
    flag = 0
    for node in sub:
        if "cherry" not in node:
            flag = 1
    # use CRISPR
    if not flag:
        CRISPR_label = ""
        CRISPR_cnt = 0
        for node in sub:
            if node in crispr_pred:
                CRISPR_cnt+=1
                CRISPR_label = crispr_pred[node]
        if CRISPR_cnt == 1:
            for node in sub:
                node2label[node] = CRISPR_label

# check subgraph situation 2
for sub in nx.connected_components(G):
    sub_label = []
    for node in sub:
        if node in virus2id.keys():
            virus_label = virus_df[virus_df['Accession'] == neighbor]['Species'].values[0]
            sub_label.append(virus_label)
        elif neighbor in prokaryote2id.keys():
            prokaryote_label = prokaryote_df[prokaryote_df['Accession'] == neighbor]['Species'].values[0]
            sub_label.append(prokaryote_label)
    if set(sub_label) == 1:
        for node in sub:
            node2label[node] = sub_label[0]



id2node = {idx: node for idx, node in enumerate(G.nodes())}
node2id = {node: idx for idx, node in enumerate(G.nodes())}

adj = nx.adjacency_matrix(G)
pkl.dump(adj,          open("GCN_data/graph.list", "wb" ))
pkl.dump(node_feature, open("GCN_data/feature.list", "wb" ))
pkl.dump(node2label,   open("GCN_data/node2label.dict", "wb" ))
pkl.dump(id2node,      open("GCN_data/id2node.dict", "wb" ))
pkl.dump(node2id,      open("GCN_data/node2id.dict", "wb" ))
pkl.dump(test_id,      open("GCN_data/test_id.dict", "wb" ))


