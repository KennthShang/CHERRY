import os
import subprocess

accession_list = os.listdir("ncbi_dataset/data/")

for accession in accession_list:
    for name in os.listdir("ncbi_dataset/data/"+accession):
        if name.rsplit(".", 1)[-1] == "fna":
            try:
                _ = os.system("mv ncbi_dataset/data/"+accession+"/"+name+" bacteria_genome/"+accession.split(".")[0]+".fasta")
            except:
                print(accession+"error")


import os
import subprocess

accession_list = os.listdir("ncbi_dataset/data/")

for accession in accession_list:
    for name in os.listdir("ncbi_dataset/data/"+accession):
        if name.rsplit(".", 1)[-1] == "faa":
            try:
                _ = os.system("mv ncbi_dataset/data/"+accession+"/"+name+" bacteria_protein/"+accession.split(".")[0]+".fasta")
            except:
                print(accession+"error")


############################### Diamond #################################
import os
import subprocess
genome_list = os.listdir("bacteria_protein")
for genome in genome_list:
    make_diamond_cmd = 'diamond makedb --threads 8 --in bacteria_protein/'+ genome +' -d diamond/'+genome.split(".")[0]+'.dmnd'
    print("Creating Diamond database...")
    _ = subprocess.check_call(make_diamond_cmd, shell=True)
    diamond_cmd = 'diamond blastx --threads 8 --sensitive -d diamond/'+genome.split(".")[0]+'.dmnd -q query.fa -o tab/'+ genome.split(".")[0] + '.tab --evalue 1e-200'
    print("Running Diamond...")
    _ = subprocess.check_call(diamond_cmd, shell=True)
    diamond_out_fp = "tab/" + genome.split(".")[0] + '.tab'
    database_abc_fp = "abc/" + genome.split(".")[0] + '.abc'
    _ = subprocess.check_call("awk '$6 < 20 {{print $1,$2,$11}}' {0} > {1}".format(diamond_out_fp, database_abc_fp), shell=True)
   


############################### BLAST #################################
import os
import subprocess
genome_list = os.listdir("real")
for genome in genome_list:
    make_blast_cmd = 'makeblastdb -in real/'+ genome +' -dbtype nucl -parse_seqids -out bacteria_db/'+genome.split(".")[0]
    print("Creating blast database...")
    _ = subprocess.check_call(make_blast_cmd, shell=True)
    blast_cmd = 'blastn -query query.fa -db bacteria_db/'+genome.split(".")[0]+' -outfmt 6 -out blast_tab/'+genome.split(".")[0]+'.tab -num_threads 8'
    print("Running blastn...")
    _ = subprocess.check_call(blast_cmd, shell=True)
    break
    

############################### CRISPR ################################# 
import os
import subprocess
import Bio
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

genome_list = os.listdir("bacteria_protein")
for genome in genome_list:
    CRISPR_cmd = 'java -cp CRT1.2-CLI.jar crt bacteria_genome/'+ genome +' CRISPR/'+genome.split(".")[0]+'.cri'
    print("Capturing CRISPR")
    _ = subprocess.check_call(CRISPR_cmd, shell=True)


def special_match(strg, search=re.compile(r'[^ACGT]').search):
    return not bool(search(strg))


accession_file_list = os.listdir("CRISPR/")

for accession_file in accession_file_list:
    CRISPR_list = []
    with open('CRISPR/'+accession_file) as file_in:
        txt = list(file_in.readlines())
        for i in range(len(txt)):
            if 'No CRISPR elements were found.' in txt[i]:
                break
            try:
                CRISPR_list.append(txt[i].split('\t')[3])
            except:
                continue
    # remove nonCRISPR
    clean_CRISPR_list = []
    for seq in CRISPR_list:
        if special_match(seq) and seq != '':
            clean_CRISPR_list.append(seq)
    CRISPR_list = clean_CRISPR_list
    # save file
    if len(CRISPR_list) > 0:
        cnt = 1
        accession = accession_file.split('.')[0]
        record_list = []
        for CRISPR in CRISPR_list:
            rec = SeqRecord(Seq(CRISPR), id=accession+ "." + str(cnt), description="")
            record_list.append(rec)
            cnt += 1
        _ = SeqIO.write(record_list, "CRISPR_fasta/"+accession+".fasta", "fasta")


accession_file_list = os.listdir("CRISPR_fasta/")
for accession_file in accession_file_list:
    accession = accession_file.split('.')[0]
    blastn_cmd = 'blastn -query CRISPR_fasta/'+accession_file+' -db blastdb/query -outfmt 6 -out CRISPR_blast_tab/'+accession+'.tab -num_threads 8'
    _ = subprocess.check_call(blastn_cmd, shell=True)
    break



import os
import subprocess
genome_list = os.listdir("bacteria_genome")
for genome in genome_list:
    blast_cmd = 'blastn -query query.fa -db bacteria_db/'+genome.split(".")[0]+' -outfmt 6 -out o_3k/'+genome.split(".")[0]+'.tab -num_threads 8'
    print("Running blastn...")
    _ = subprocess.check_call(blast_cmd, shell=True)