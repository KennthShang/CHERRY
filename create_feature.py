import os
import numpy as np
import Bio
from Bio import SeqIO
import pickle as pkl
import argparse



parser = argparse.ArgumentParser(description='manual to this script')
parser.add_argument('--mode', type=str, default = 'virus')
inputs = parser.parse_args()


def return_4mer(file_in_fn):
    # alphbet
    k_list = ["A", "C", "G", "T"]
    nucl_list = ["A", "C", "G", "T"]
    for i in range(3):
        tmp = []
        for item in nucl_list:
            for nucl in k_list:
                tmp.append(nucl+item)
        k_list = tmp
    # dictionary
    mer2dict = {mer: idx for idx, mer in enumerate(k_list)}
    # search files
    file_list = os.listdir(file_in_fn)
    num_file = len(file_list)
    file2idx = {}
    # convert to words
    feature = np.zeros((num_file, 256))
    for idx, file in enumerate(file_list):
        file2idx[file.rsplit('.', 1)[0]] = idx
        for record in SeqIO.parse(file_in_fn + file, 'fasta'):
            seq = str(record.seq)
            seq = seq.upper()
            for pos in range(len(seq)-3):
                try:
                    feature[idx][mer2dict[seq[pos:pos+4]]] += 1
                except:
                    #print(seq[pos:pos+4])
                    pass
    # nomarlization
    norm_feature = np.zeros((num_file, 256))
    for i in range(len(feature)):
        norm_feature[i] = (feature[i] - np.min(feature[i]))/(np.max(feature[i]) - np.min(feature[i]))
    return norm_feature, file2idx


virus, virus2id = return_4mer('train_phage/')
pkl.dump(virus2id, open('node_feature/virus.dict', 'wb'))
pkl.dump(virus, open('node_feature/virus.F', 'wb'))

prokaryote, prokaryote2id = return_4mer('prokaryote/')
pkl.dump(prokaryote2id, open('node_feature/prokaryote.dict', 'wb'))
pkl.dump(prokaryote, open('node_feature/prokaryote.F', 'wb'))

if inputs.mode == 'virus':
    test_virus, test_virus2id = return_4mer('single_contig/')
    pkl.dump(test_virus2id, open('node_feature/test_virus.dict', 'wb'))
    pkl.dump(test_virus, open('node_feature/test_virus.F', 'wb'))
elif inputs.mode == 'prokaryote':
    test_prokaryote, test_prokaryote2id = return_4mer('new_prokaryote/')
    pkl.dump(test_prokaryote2id, open('node_feature/test_prokaryote.dict', 'wb'))
    pkl.dump(test_prokaryote, open('node_feature/test_prokaryote.F', 'wb'))
    test_virus, test_virus2id = return_4mer('single_contig/')
    pkl.dump(test_virus2id, open('node_feature/test_virus.dict', 'wb'))
    pkl.dump(test_virus, open('node_feature/test_virus.F', 'wb'))
else:
    print('wrong parameters')
    exit()



