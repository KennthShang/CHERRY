import os
import sys
import Bio
import subprocess
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast.Applications import NcbiblastnCommandline



prodigal_cmd = 'prodigal -i dataset/nucl.fasta -a dataset/protein.fasta -f gff -p meta'
print("Running prodigal...")
_ = subprocess.check_call(prodigal_cmd, shell=True)



proteins = []
contigs  = []
keywords = []
for record in SeqIO.parse('dataset/protein.fasta', 'fasta'):
    name = record.id
    contigs.append(name.rsplit("_", 1)[0])
    proteins.append(name)
    keywords.append('hypothetical protein')

gene2genome_df = pd.DataFrame({'protein_id': proteins, 'contig_id': contigs, 'keywords': keywords})
gene2genome_df.to_csv('dataset/database_gene_to_genome.csv', index=False)