import os
import sys
import Bio
import argparse
import subprocess
import scipy as sp
import numpy as np
import pandas as pd
import pickle as pkl
import networkx as nx
import scipy.stats as stats
import scipy.sparse as sparse
import utils
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast.Applications import NcbiblastnCommandline
from utils import make_protein_clusters_mcl, load_mcl_clusters, build_clusters
from utils import build_pc_matrices, to_clusterer, create_network

# Defined folder
out_fn = "out/"
contig_in = "input/"
contig_out = "single_contig/"
file_in_fn = "single_contig/"
file_out_fnn = "all_proteins/"
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

check_folder(out_fn)
check_folder(file_in_fn)
check_folder(file_out_fnn)

################################################################################
############################  Rename the files  ################################
################################################################################

# Split contigs files into single contig per file.
file_list = sorted(os.listdir(contig_in))
seq = []
old_file_id = 0
contig_id = 0
with open("name_list.csv",'w') as list_out:
    list_out.write("contig_name,idx\n")
    for file_n in file_list:
        for record in SeqIO.parse(contig_in+file_n, "fasta"):
            name = "cherry_"+str(old_file_id) + "_" + str(contig_id)
            _ = list_out.write(record.id + "," + name + "\n")
            record.id = "cherry_"+str(old_file_id) + "_" + str(contig_id)
            _ = SeqIO.write(record, contig_out+name+".fasta", "fasta")
            contig_id += 1
        old_file_id += 1


################################################################################
###################### Translate contigs into 6 ORFs ###########################
################################################################################

file_list = os.listdir(file_in_fn)
for file in file_list:
    prodigal_cmd = 'prodigal -i ' + file_in_fn + file + ' -a '+ file_out_fnn + file +' -f gff -p meta'
    print("Running prodigal...")
    _ = subprocess.check_call(prodigal_cmd, shell=True)

all_protein_f = out_fn + "all_proteins.fa"
_ = subprocess.check_call("cat {0} {1} > {2}".format(file_out_fnn+"*", "dataset/protein.fasta", all_protein_f), shell=True)

################################################################################
############################## Run diamond BLASTp  #############################
################################################################################

print("\n\n" + "{:-^80}".format("Diamond BLASTp"))
print("Creating Diamond database and running Diamond...")

try:
    make_diamond_cmd = 'diamond makedb --threads 8 --in out/all_proteins.fa -d out/database.dmnd'
    print("Creating Diamond database...")
    _ = subprocess.check_call(make_diamond_cmd, shell=True)
    
    diamond_cmd = 'diamond blastp --threads 8 --sensitive -d out/database.dmnd -q out/all_proteins.fa -o out/database.self-diamond.tab'
    print("Running Diamond...")
    _ = subprocess.check_call(diamond_cmd, shell=True)
    diamond_out_fp = "out/database.self-diamond.tab"
    database_abc_fp = "out/database.self-diamond.tab.abc"
    _ = subprocess.check_call("awk '$1!=$2 {{print $1,$2,$11}}' {0} > {1}".format(diamond_out_fp, database_abc_fp), shell=True)
except:
    print("create database failed")
    exit(1)




# Generating gene-to-genome.csv: protein_id, contig_id, keywords
blastp = pd.read_csv(database_abc_fp, sep=' ', names=["contig", "ref", "e-value"])
protein_id = sorted(list(set(blastp["contig"].values)|set(blastp["ref"].values)))
contig_protein = [item for item in protein_id if "cherry" == item.split("_")[0]]
contig_id = [item.rsplit("_", 1)[0] for item in contig_protein]
description = ["hypothetical protein" for item in contig_protein]
gene2genome = pd.DataFrame({"protein_id": contig_protein, "contig_id": contig_id ,"keywords": description})
gene2genome.to_csv(out_fn+"contig_gene_to_genome.csv", index=None)


_ = subprocess.check_call("cat dataset/database_gene_to_genome.csv {0}contig_gene_to_genome.csv > {1}gene_to_genome.csv".format(out_fn, out_fn), shell=True)




################################################################################
################################## Run MCL #####################################
################################################################################
print("\n\n" + "{:-^80}".format("Protein clustering"))
print("Loading proteins...")
gene2genome_fp = out_fn+"gene_to_genome.csv"
gene2genome_df = pd.read_csv(gene2genome_fp, sep=',', header=0)


# Parameters for MCL
pc_overlap, pc_penalty, pc_haircut, pc_inflation = 0.8, 2.0, 0.1, 2.0
pcs_fp = make_protein_clusters_mcl(database_abc_fp, out_fn, pc_inflation)
print("Building the cluster and profiles (this may take some time...)")


# Dump MCL results
protein_df, clusters_df, profiles_df, contigs_df = build_clusters(pcs_fp, gene2genome_df)
print("Saving files")
dfs = [gene2genome_df, contigs_df, clusters_df]
names = ['proteins', 'contigs', 'pcs']
output_dir = out_fn

for name, df in zip(names, dfs):
    fn = "Cyber_{}.csv".format(name)
    fp = os.path.join(output_dir, fn)
    index_id = name.strip('s') + '_id'
    if not os.path.exists(fp):
        df.set_index(index_id).to_csv(fp)
    else:
        print("File {} exists and will be used. Use -f to overwrite.".format(fn))

profiles_fn = "Cyber_profiles.csv"
profiles_fp = os.path.join(out_fn, profiles_fn)
if not os.path.exists(profiles_fp):
    profiles_df.to_csv(profiles_fp, index=False)
else:
    print("File {} exists and will be used. Use -f to overwrite.".format(profiles_fn))




################################################################################
################################## Run MCL #####################################
################################################################################

# Loding dataset
contigs_df = pd.read_csv("out/Cyber_contigs.csv")
clusters_df = pd.read_csv("out/Cyber_pcs.csv")
profiles_df = pd.read_csv("out/Cyber_profiles.csv")

# Replace names
contigs_csv_df = contigs_df.copy()
print("Read {} entries from {}".format(len(contigs_csv_df), os.path.join(output_dir, '{}_contigs.csv'.format(name))))
contigs_csv_df.index.name = "pos"
contigs_csv_df.reset_index(inplace=True)

pcs_csv_df = clusters_df.copy()
profiles = profiles_df.copy()

# Filtering the PC profiles that appears only once
before_filter = len(profiles)
cont_by_pc = profiles.groupby("pc_id").count().contig_id.reset_index()

# get the number of contigs for each pcs and add it to the dataframe
cont_by_pc.columns = ["pc_id", "nb_proteins"]
pcs_csv_df = pd.merge(pcs_csv_df, cont_by_pc, left_on="pc_id", right_on="pc_id", how="left")
pcs_csv_df.fillna({"nb_proteins": 0}, inplace=True)

# Drop the pcs that <= 1 contig from the profiles.
pcs_csv_df = pcs_csv_df[pcs_csv_df['nb_proteins'] > 1]  # .query("nb_contigs>1")
at_least_a_cont = cont_by_pc[cont_by_pc['nb_proteins'] > 1]  # cont_by_pc.query("nb_contigs>1")
profiles = profiles[profiles['pc_id'].isin(at_least_a_cont.pc_id)]
print("Read {} entries (dropped {} singletons) from {}".format(len(profiles), (before_filter - len(profiles)), profiles_fp))
pcs_csv_df = pcs_csv_df.reset_index(drop=True)
pcs_csv_df.index.name = "pos"
pcs_csv_df = pcs_csv_df.reset_index()

matrix, singletons = build_pc_matrices(profiles, contigs_csv_df, pcs_csv_df)
profiles_csv = {"matrix": matrix, "singletons": singletons}
merged_df = contigs_csv_df
merged_fp = os.path.join(output_dir, 'merged_df.csv')
merged_df.to_csv(merged_fp)

ntw = create_network(matrix, singletons, thres=1, max_sig=300)
fi = to_clusterer(ntw, out_fn+"intermediate.ntw", merged_df.copy())

################################################################################
################################# Run BLASTN ###################################
################################################################################
_ = subprocess.check_call("cat single_contig/* > out/test.fa", shell=True)
query_file = "out/test.fa"
db_virus_prefix = "dataset/virus_db/allVIRUS"
output_file = "out/virus_out.tab"
virus_call = NcbiblastnCommandline(query=query_file,db=db_virus_prefix,out=output_file,outfmt="6 qseqid sseqid evalue pident length qlen", evalue=1,gapopen=10,penalty=-1,
                                  gapextend=2,word_size=7,dust='no',
                                 task='blastn-short',perc_identity=90,num_threads=16)
virus_call()


virus_pred = {}
with open(output_file) as file_out:
    for line in file_out.readlines():
        parse = line.replace("\n", "").split("\t")
        virus = parse[0]
        ref_virus = parse[1].split('|')[1]
        ref_virus = ref_virus.split('.')[0]
        ident = float(parse[-3])
        length = float(parse[-2])
        qlen = float(parse[-1])
        if virus not in virus_pred and length/qlen > 0.9 and ident > 0.9:
            virus_pred[virus] = ref_virus

pkl.dump(virus_pred, open('out/virus_pred.dict', 'wb'))


################################################################################
############################### Dump the graph #################################
################################################################################

G = nx.Graph()
# Create graph
with open(out_fn+"intermediate.ntw") as file_in:
    for line in file_in.readlines():
        tmp = line[:-1].split(" ")
        node1 = tmp[0]
        node2 = tmp[1]
        G.add_edge(node1, node2, weight = 1)


graph = "out/phage_phage.ntw"
with open(graph, 'w') as file_out:
    for node1 in G.nodes():
        for _,node2 in G.edges(node1):
            _ = file_out.write(node1+","+node2+"\n")
