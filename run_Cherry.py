import  torch
from    torch import nn
from    torch import optim
from    torch.nn import functional as F
import torch.utils.data as Data
 
import os
import  numpy as np
from    data import load_data, preprocess_features, preprocess_adj, sample_mask
import  model
from    config import  args
from    utils import masked_loss, masked_acc, masked_ECE
import  pickle as pkl
import  scipy.sparse as sp
import argparse
from scipy.special import softmax
from sklearn.metrics import classification_report
from sklearn.metrics import accuracy_score
from collections import Counter
from importlib import reload
import networkx as nx
import pandas as pd
import random
 



seed = 123
np.random.seed(seed)
torch.random.manual_seed(seed)
inputs = args.parse_args()

# Check whether cuda is enable
if torch.cuda.is_available():
    torch.cuda.set_device(inputs.gpus)
    device = torch.device('cuda')
else:
    print("Running with cpu")
    device = torch.device('cpu')


# load data
adj         = pkl.load(open("GCN_data/graph.list",'rb'))
features    = pkl.load(open("GCN_data/feature.list",'rb'))
id2node     = pkl.load(open("GCN_data/id2node.dict",'rb'))
node2id     = pkl.load(open("GCN_data/node2id.dict", "rb" ))
idx_test    = pkl.load(open("GCN_data/test_id.dict", 'rb'))
node2label  = pkl.load(open("GCN_data/node2label.dict",'rb'))
crispr_pred = pkl.load(open('out/crispr_pred.dict', 'rb'))
bacteria_df = pd.read_csv('dataset/prokaryote.csv')



# prokaryotes in the training set
trainable_host = []
for file in os.listdir('prokaryote/'):
    trainable_host.append(file.rsplit('.', 1)[0])



host2id = {}
label2hostid =  {}
trainable_host_idx = []
trainable_label = []
for idx, node in id2node.items():
    # if prokaryote
    if node in trainable_host:
        host2id[node] = idx
        trainable_host_idx.append(idx)
        trainable_label.append(node2label[node])
        label2hostid[node2label[node]] = idx




# pre-processing
features = sp.csc_matrix(features)
print('adj:', adj.shape)
print('features:', features.shape)


# convert to torch tensor
features = preprocess_features(features)
supports = preprocess_adj(adj)
num_classes = len(set(list(node2label.values())))+1
# graph
i = torch.from_numpy(features[0]).long().to(device)
v = torch.from_numpy(features[1]).to(device)
feature = torch.sparse.FloatTensor(i.t(), v, features[2]).float().to(device)
feature = feature.to_dense()
i = torch.from_numpy(supports[0]).long().to(device)
v = torch.from_numpy(supports[1]).to(device)
support = torch.sparse.FloatTensor(i.t(), v, supports[2]).float().to(device)
support = support.to_dense()


print('x :', feature)
print('sp:', support)
feat_dim = adj.shape[0]
node_dim = feature.shape[1]


# Definition of the model
net = model.encoder(feat_dim, node_dim, node_dim, 0)
decoder = model.decoder(node_dim, 128, 32)


# Load pre-trained model
encoder_dict = torch.load(f"dataset/pkl/Encoder_Species.pkl", map_location='cpu')
decoder_dict = torch.load(f"dataset/pkl/Decoder_Species.pkl", map_location='cpu')
net.load_state_dict(encoder_dict)
decoder.load_state_dict(decoder_dict)

net.to(device)
decoder.to(device)

# end-to-end training
params = list(net.parameters()) + list(decoder.parameters())
optimizer = optim.Adam(params, lr=0.001)#args.learning_rate
loss_func = nn.BCEWithLogitsLoss()


#################################################################
#####################  evaluation metrics #######################
#################################################################

def train_accuracy():
    with torch.no_grad():
        total = 0
        correct = 0
        for i in range(len(encode)):
            if idx_test[id2node[i]] != 0 or i in trainable_host_idx:
                continue
            virus_feature = encode[i]
            max_pred = 0
            pred_label = ""
            for label in trainable_label:
                prokaryote_feature = encode[label2hostid[label]]
                pred = decoder(virus_feature - prokaryote_feature)
                if pred > max_pred:
                    max_pred = pred
                    pred_label = label
            if pred_label == node2label[id2node[i]]:
                correct+=1
            total += 1
    return correct/total


def test_accuracy():
    with torch.no_grad():
        total = 0
        correct = 0
        for i in range(len(encode)):
            if idx_test[id2node[i]] != 0 or i in trainable_host_idx:
                continue
            if 'NC_' not in id2node[i]:
                continue
            virus_feature = encode[i]
            max_pred = 0
            pred_label = ""
            for label in trainable_label:
                prokaryote_feature = encode[label2hostid[label]]
                pred = decoder(virus_feature - prokaryote_feature)
                if pred > max_pred:
                    max_pred = pred
                    pred_label = label
            if pred_label == node2label[id2node[i]]:
                correct+=1
            total += 1
    return correct/total


def train_topk_accuracy(k):
    with torch.no_grad():
        total = 0
        correct = 0
        for i in range(len(encode)):
            if idx_test[id2node[i]] != 0 or i in trainable_host_idx:
                continue
            virus_feature = encode[i]
            max_pred = 0
            pred_label = []
            for label in trainable_label:
                prokaryote_feature = encode[label2hostid[label]]
                pred = decoder(virus_feature - prokaryote_feature)
                pred_label.append(pred)
            pred_label = sorted(pred_label, reverse=True)
            real_pred = decoder(virus_feature - encode[label2hostid[node2label[id2node[i]]]])
            if real_pred in pred_label[:k]:
                correct += 1
            total += 1
    return correct/total


#################################################################
##########################  Training  ###########################
#################################################################

if inputs.model == 'retrain':
    _ = net.train()
    _ = decoder.train()
    for epoch in range(120):
        encode = net((feature, support))
        loss = 0
        for label in trainable_label:
            virus_idx_list = [idx for idx in range(len(encode)) if idx not in trainable_host_idx and node2label[id2node[idx]] == label]
            prokaryote_feature = encode[label2hostid[label]]
            for virus_idx in virus_idx_list:
                # neg sampling loss
                virus_feature = encode[virus_idx]
                pred = decoder(virus_feature - prokaryote_feature)
                loss += loss_func(pred, torch.ones([1]).to(device))
                cnt=0
                while cnt < 1:
                    fake_label = random.choice(trainable_label)
                    if fake_label != label:
                        fake_prokaryote_feature = encode[label2hostid[fake_label]]
                        pred = decoder(virus_feature - fake_prokaryote_feature)
                        loss += loss_func(pred, torch.zeros([1]).cuda())
                        cnt+=1
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        #print(loss.cpu().detach().numpy())
        if epoch % 50 == 0:
            train_acc = train_accuracy()
            print("Train acc: {:.2f} ".format(train_acc))
            if train_acc > 0.7:
                torch.save(net.state_dict(), f'Custom_Encoder_params.pkl')
                torch.save(decoder.state_dict(), f'Custom_Decoder_params.pkl')
                break
    torch.save(net.state_dict(), f'Custom_Encoder_params.pkl')
    torch.save(decoder.state_dict(), f'Custom_Decoder_params.pkl')



#################################################################
#########################  Prediction  ##########################
#################################################################

# predicting host
if inputs.mode == 'virus':
    node2pred = {}
    with torch.no_grad():
        encode = net((feature, support))
        for i in range(len(encode)):
            confident_label = 'unknown'
            if idx_test[id2node[i]] == 0:
                continue
            if idx_test[id2node[i]] == 1:
                confident_label = node2label[id2node[i]]
            virus_feature = encode[i]
            pred_label_score = []
            for label in set(trainable_label):
                if label == confident_label:
                    pred_label_score.append((label, 1))
                    continue
                prokaryote_feature = encode[label2hostid[label]]
                pred = decoder(virus_feature - prokaryote_feature)
                pred_label_score.append((label, torch.sigmoid(pred).detach().cpu().numpy()[0]))
            node2pred[id2node[i]] = sorted(pred_label_score, key=lambda tup: tup[1], reverse=True)
        for virus in crispr_pred:
            if virus not in node2pred:
                pred = prokaryote_df[prokaryote_df['Accession'] == crispr_pred[virus]]['Species'].values[0]
                node2pred[virus] = (pred, 1)
        # dump the prediction
        with open(f"tmp_pred/predict.csv", 'w') as file_out:
            file_out.write('contig,')
            for i in range(inputs.topk):
                file_out.write(f'Top_{i+1}_label,Score_{i+1},')
            file_out.write('\n')
            for contig in node2pred:
                file_out.write(f'{contig},')
                cnt = 1
                for label, score in node2pred[contig]:
                    if cnt > inputs.topk:
                        break
                    cnt+=1
                    file_out.write(f'{label},{score:.2f},')
                file_out.write('\n')




# predicting virus
if inputs.mode == 'prokaryote':
    candidate_host = []
    for file in os.listdir('new_prokaryote/'):
        candidate_host.append(file.rsplit('.', 1)[0])
    candidateidx = []
    for host in candidate_host:
        if host in node2id:
            candidateidx.append(node2id[host])
    host2pred = {}
    with torch.no_grad():
        encode = net((feature, support))
        for host in candidate_host:
            if host not in node2id:
                host2pred[host] = "unknown"
            else:
                prokaryote_feature = encode[node2id[host]]
                for i in range(len(encode)):
                    if i in trainable_host_idx or i in candidateidx:
                        continue
                    virus_feature = encode[i]
                    logit = torch.sigmoid(decoder(virus_feature - prokaryote_feature))
                    if logit > inputs.t:
                        try:
                            host2pred[host].append(id2node[i])
                        except:
                            host2pred[host] = [id2node[i]]
    with open('tmp_pred/predict.csv', 'w') as file:
        file.write('prokaryote,virus\n')
        for prokaryote in host2pred:
            file.write(prokaryote+','+ "|".join(host2pred[prokaryote])+'\n')













