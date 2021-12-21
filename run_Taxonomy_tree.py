import os
import sys
import Bio
import subprocess
import numpy as np
import pandas as pd



parser = argparse.ArgumentParser(description='manual to this script')
parser.add_argument('--k', type=int, default = 1)
inputs = parser.parse_args()

pred_df  = pd.read_csv('final_prediction.csv')
label_df = pd.read_csv('dataset/prokaryote.csv')

accession = pred_df['contig_name'].values
prediction = pred_df[f'Top_{inputs.k}_label'].values

tree_df = []
for label in prediction:
    tree_df.append(label_df[label_df['Species'] == label])

tree_df = pd.concat(tree_df)
tree_df['Accession'] = accession
tree_df = tree_df.rename(columns={'Accession':'Contig_name'})
tree_df.to_csv(f'Top_{inputs.k}_prediction_taxonomy.csv', index=False)