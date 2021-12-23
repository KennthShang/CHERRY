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
from Bio.Blast.Applications import NcbiblastnCommandline


parser = argparse.ArgumentParser(description='manual to this script')
parser.add_argument('--mode', type=str, default = 'virus')
inputs = parser.parse_args()


# Defined folder
prokaryote_fn= "prokaryote/"
new_prokaryote = "new_prokaryote/"
blast_database_out = "blast_db/"
new_blast_database_out = "new_blast_db/"
blast_tab_out = "blast_tab/"
new_blast_tab_out = "new_blast_tab/"
Knowledge_graph = "Cyber_data/"

################################################################################
############################  Check the folder #################################
################################################################################
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

check_folder(blast_database_out)
check_folder(blast_tab_out)
if inputs.mode != 'virus':
    check_folder(new_blast_tab_out)
    check_folder(new_blast_database_out)


#  combine phage file 
_ = subprocess.check_call("cat dataset/nucl.fasta single_contig/* > out/query.fa", shell=True)




################################################################################
###############################  Run CRISRP   ##################################
################################################################################

query_file = "out/test.fa"
db_host_crispr_prefix = "dataset/crispr_db/allCRISPRs"
output_file = "out/crispr_out.tab"
crispr_call = NcbiblastnCommandline(query=query_file,db=db_host_crispr_prefix,out=output_file,outfmt="6 qseqid sseqid evalue pident length slen", evalue=1,gapopen=10,penalty=-1,
                                  gapextend=2,word_size=7,dust='no',
                                 task='blastn-short',perc_identity=90,num_threads=16)
crispr_call()


crispr_pred = {}
with open(output_file) as file_out:
    for line in file_out.readlines():
        parse = line.replace("\n", "").split("\t")
        virus = parse[0]
        prokaryote = parse[1].split('|')[1]
        prokaryote = prokaryote.split('.')[0]
        ident = float(parse[-3])
        length = float(parse[-2])
        slen = float(parse[-1])
        if virus not in crispr_pred and length/slen > 0.9 and ident > 0.9:
            crispr_pred[virus] = prokaryote

pkl.dump(crispr_pred, open('out/crispr_pred.dict', 'wb'))




################################################################################
###############################  Run BLASTN(train)   ###########################
################################################################################


genome_list = os.listdir(prokaryote_fn)
for genome in genome_list:
    make_blast_cmd = 'makeblastdb -in '+ prokaryote_fn + genome +' -dbtype nucl -parse_seqids -out '+ blast_database_out + genome.split(".")[0]
    print("Creating blast database...")
    _ = subprocess.check_call(make_blast_cmd, shell=True)
    blast_cmd = 'blastn -query out/query.fa -db blast_db/'+genome.split(".")[0]+' -outfmt 6 -out '+ blast_tab_out + genome.split(".")[0]+'.tab -num_threads 16'
    print("Running blastn...")
    _ = subprocess.check_call(blast_cmd, shell=True)


################################################################################
###############################  Run BLASTN(test)   ############################
################################################################################

if inputs.mode != 'virus':
    genome_list = os.listdir(new_prokaryote)
    for genome in genome_list:
        make_blast_cmd = 'makeblastdb -in '+ new_prokaryote + genome +' -dbtype nucl -parse_seqids -out '+ new_blast_database_out + genome.split(".")[0]
        print("Creating blast database...")
        _ = subprocess.check_call(make_blast_cmd, shell=True)
        blast_cmd = 'blastn -query out/query.fa -db new_blast_db/'+genome.split(".")[0]+' -outfmt 6 -out '+ new_blast_tab_out + genome.split(".")[0]+'.tab -num_threads 16'
        print("Running blastn...")
        _ = subprocess.check_call(blast_cmd, shell=True)




################################################################################
###############################  virus-host   ##################################
################################################################################

# add connections between prokaryotes and viruses
tab_file_list = os.listdir(blast_tab_out)
prokaryote2virus = {}
for file in tab_file_list:
    prokaryote_id = file.split('.')[0]
    virus_id_list = []
    with open(blast_tab_out+file) as file_in:
        for line in file_in.readlines():
            tmp = line.split('\t')
            virus_id = tmp[0]
            try:
                prokaryote2virus[prokaryote_id].append(virus_id)
            except:
                prokaryote2virus[prokaryote_id] = [virus_id]

if inputs.mode != 'virus':
    tab_file_list = os.listdir(new_blast_tab_out)
    for file in tab_file_list:
        prokaryote_id = file.split('.')[0]
        virus_id_list = []
        with open(blast_tab_out+file) as file_in:
            for line in file_in.readlines():
                tmp = line.split('\t')
                virus_id = tmp[0]
                try:
                    prokaryote2virus[prokaryote_id].append(virus_id)
                except:
                    prokaryote2virus[prokaryote_id] = [virus_id]


# De-duplication
for key in prokaryote2virus:
    prokaryote2virus[key] = list(set(prokaryote2virus[key]))


# Save the virus-host graph
with open("out/phage_host.ntw", 'w') as file_out:
    for prokaryote in prokaryote2virus:
        for virus in prokaryote2virus[prokaryote]:
            _ = file_out.write(prokaryote + "," + virus + "\n")
