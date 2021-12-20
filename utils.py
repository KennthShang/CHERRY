import os
import sys
import Bio
import scipy as sp
import numpy as np
import pandas as pd
import pickle as pkl
import networkx as nx
import scipy.stats as stats
import scipy.sparse as sparse
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import  torch
from    torch import nn
from    torch.nn import functional as F


def make_diamond_db(cpu: int):
    diamond_db_bp = "database/database.dmnd"
    aa_fp = "database/protein.fasta"

    make_diamond_cmd = ['diamond', 'makedb', '--threads', str(cpu), '--in', aa_fp, '-d', diamond_db_bp]
    print("Creating Diamond database...")
    res = subprocess.run(make_diamond_cmd, check=True, stdout=subprocess.PIPE)
    if res.returncode != 0:
        print('Error creating Diamond database')
        exit(1)
    diamond_db_fp = diamond_db_bp + '.dmnd'
    return diamond_db_fp




def run_diamond(aa_fp, db_fp, cpu: int, diamond_out_fnn):
    # More sensitive as an option?
    diamond_cmd = ['diamond', 'blastp', '--threads', str(cpu), '--sensitive', '-d', db_fp, '-q', aa_fp,
                   '-o', diamond_out_fnn]
    print("Running Diamond...")
    res = subprocess.run(diamond_cmd, check=True, stdout=subprocess.PIPE)
    if res.returncode != 0:
        print('Error running Diamond')
        exit(1)
    return diamond_out_fnn


def make_protein_clusters_mcl(abc_fp, out_p, inflation=2):
    print("Running MCL...")
    abc_fn = "merged"
    mci_fn = '{}.mci'.format(abc_fn)
    mci_fp = os.path.join(out_p, mci_fn)
    mcxload_fn = '{}_mcxload.tab'.format(abc_fn)
    mcxload_fp = os.path.join(out_p, mcxload_fn)
    subprocess.check_call("mcxload -abc {0} --stream-mirror --stream-neg-log10 -stream-tf 'ceil(200)' -o {1}"
                          " -write-tab {2}".format(abc_fp, mci_fp, mcxload_fp), shell=True)
    mcl_clstr_fn = "{0}_mcl{1}.clusters".format(abc_fn, int(inflation*10))
    mcl_clstr_fp = os.path.join(out_p, mcl_clstr_fn)
    subprocess.check_call("mcl {0} -I {1} -use-tab {2} -o {3}".format(
        mci_fp, inflation, mcxload_fp, mcl_clstr_fp), shell=True)
    return mcl_clstr_fp



def load_mcl_clusters(fi):
    with open(fi) as f:
        c = [line.rstrip("\n").split("\t") for line in f]
    c = [x for x in c if len(c) > 1]
    nb_clusters = len(c)
    formatter = "PC_{{:>0{}}}".format(int(round(np.log10(nb_clusters))+1))
    name = [formatter.format(str(i)) for i in range(nb_clusters)]
    size = [len(i) for i in c]
    clusters_df = pd.DataFrame({"size": size, "pc_id": name}).set_index("pc_id")
    return clusters_df, name, c


def build_clusters(fp, gene2genome):
    # Read MCL
    clusters_df, name, c = load_mcl_clusters(fp)
    print("Using MCL to generate PCs.")
    # Assign each prot to its cluster
    gene2genome.set_index("protein_id", inplace=True)  # id, contig, keywords, cluster
    for prots, clust in zip(c, name):
        try:
            gene2genome.loc[prots, "cluster"] = clust
        except KeyError:
            prots_in = [p for p in prots if p in gene2genome.index]
            not_in = frozenset(prots) - frozenset(prots_in)
            print("{} protein(s) without contig: {}".format(len(not_in), not_in))
            gene2genome.loc[prots_in, "cluster"] = clust
    # Keys
    for clust, prots in gene2genome.groupby("cluster"):
        clusters_df.loc[clust, "annotated"] = prots.keywords.count()
        if prots.keywords.count():
            keys = ";".join(prots.keywords.dropna().values).split(";")
            key_count = {}
            for k in keys:
                k = k.strip()
                try:
                    key_count[k] += 1
                except KeyError:
                    key_count[k] = 1
            clusters_df.loc[clust, "keys"] = "; ".join(["{} ({})".format(x, y) for x, y in key_count.items()])
    gene2genome.reset_index(inplace=True)
    clusters_df.reset_index(inplace=True)
    profiles_df = gene2genome.loc[:, ["contig_id", "cluster"]].drop_duplicates()
    profiles_df.columns = ["contig_id", "pc_id"]
    contigs_df = pd.DataFrame(gene2genome.fillna(0).groupby("contig_id").count().protein_id)
    contigs_df.index.name = "contig_id"
    contigs_df.columns = ["proteins"]
    contigs_df.reset_index(inplace=True)
    return gene2genome, clusters_df, profiles_df, contigs_df


# Create P-edges
def build_pc_matrices(profiles, contigs, pcs):
    pc_by_cont = profiles.groupby("contig_id").count().pc_id
    pc_by_cont = pd.merge(contigs.sort_values("pos").loc[:, ["pos", "contig_id", "proteins"]], pc_by_cont.to_frame(), how="left",
                          left_on="contig_id", right_on="contig_id").fillna(0)
    singletons = (pc_by_cont.proteins - pc_by_cont.pc_id).values
    singletons = sparse.lil_matrix(singletons).transpose()
    # Matrix
    profiles.index.name = "pos"
    profiles.reset_index(inplace=True)
    profiles = pd.merge(profiles, pcs.loc[:, ["pc_id", "pos"]], left_on="pc_id", right_on="pc_id", how="inner",
                            suffixes=["", "_pc"])  # pos, contig_id, pc_id, id (pc), pos_pc
    profiles = pd.merge(profiles, contigs.loc[:, ["contig_id", "pos"]], left_on="contig_id", right_on="contig_id", how="inner",
                            suffixes=["", "_contig"])
    profiles = profiles.loc[:, ["pos_contig", "pos_pc"]]
    matrix = sparse.coo_matrix(([1]*len(profiles), (zip(*profiles.values))), shape=(len(contigs), len(pcs)),
                               dtype="bool")
    return matrix.tocsr(), singletons.tocsr()


def to_clusterer(matrix, fi, contigs=None,names=None):
    names = contigs if names is None else names
    names = names.set_index("pos").contig_id
    with open(fi, "wt") as f:
        matrix = sparse.dok_matrix(matrix)
        for r, c in zip(*matrix.nonzero()):
            f.write(" ".join([str(x) for x in (names[r], names[c], matrix[r, c])]))
            f.write("\n")
    print("Saving network in file {0} ({1} lines).".format(fi, matrix.getnnz()))
    return fi


def create_network(matrix, singletons, thres=1, max_sig=1000):
    contigs, pcs = matrix.shape
    pcs += singletons.sum()
    # Number of comparisons
    T = 0.5 * contigs * (contigs - 1)
    logT = np.log10(T)
    # Number of protein clusters in each contig
    # = # shared pcs + #singletons
    number_of_pc = matrix.sum(1) + singletons
    number_of_pc = number_of_pc.A1  # Transform into a flat array
    # Number of common protein clusters between two contigs, tuple + commons
    commons_pc = matrix.dot(sparse.csr_matrix(matrix.transpose(), dtype=int))
    S = sparse.lil_matrix((contigs, contigs))
    total_c = float(commons_pc.getnnz())
    i = 0  # Display
    for A, B in zip(*commons_pc.nonzero()):  # For A & B sharing contigs
        if A != B:
            a, b = sorted([number_of_pc[A], number_of_pc[B]])
            pval = stats.hypergeom.sf(commons_pc[A, B] - 1, pcs, a, b)
            sig = min(max_sig, np.nan_to_num(-np.log10(pval) - logT))
            if sig > thres:
                S[min(A, B), max(A, B)] = sig
            # Display
            i += 1
            if i % 1000 == 0:
                sys.stdout.write(".")
            if i % 10000 == 0:
                sys.stdout.write("{:6.2%} {}/{}\n".format(i / total_c, i, total_c))
    S += S.T
    S = S.tocsr()
    if len(S.data) != 0:
        print("Hypergeometric contig-similarity network:\n {0:10} contigs,\n {1:10} edges (min:{2:.2}"
                    "max: {3:.2}, threshold was {4})".format(contigs, S.getnnz(), S.data.min(), S.data.max(), thres))
    else:
        raise ValueError("No edge in the similarity network !") 
    return S

def masked_ECE(out, train_label, train_mask):
    device = torch.device('cuda')
    m = nn.Softmax(dim=1)
    prob_box = dict.fromkeys([0, 1, 2, 3, 4])
    pred_box = dict.fromkeys([0, 1, 2, 3, 4])
    label_box = dict.fromkeys([0, 1, 2, 3, 4])
    # calculate pred
    prob = torch.max(m(out), 1)[0]
    pred = torch.argmax(out, 1)
    # calculate box
    for i in range(len(prob)):
        if test_mask[i] == True:
            if prob[i] > 0.9:
                try:
                    prob_box[9].append(prob[i])
                    pred_box[9].append(pred[i])
                    label_box[9].append(test_label[i])
                except:
                    prob_box[9] = [prob[i]]
                    pred_box[9] = [pred[i]]
                    label_box[9]= [test_label[i]]
            elif prob[i] > 0.8:
                try:
                    prob_box[8].append(prob[i])
                    pred_box[8].append(pred[i])
                    label_box[8].append(test_label[i])
                except:
                    prob_box[8] = [prob[i]]
                    pred_box[8] = [pred[i]]
                    label_box[8]= [test_label[i]]
            elif prob[i] > 0.7:
                try:
                    prob_box[7].append(prob[i])
                    pred_box[7].append(pred[i])
                    label_box[7].append(test_label[i])
                except:
                    prob_box[7] = [prob[i]]
                    pred_box[7] = [pred[i]]
                    label_box[7]= [test_label[i]]
            elif prob[i] > 0.6:
                try:
                    prob_box[6].append(prob[i])
                    pred_box[6].append(pred[i])
                    label_box[6].append(test_label[i])
                except:
                    prob_box[6] = [prob[i]]
                    pred_box[6] = [pred[i]]
                    label_box[6]= [test_label[i]]
            elif prob[i] > 0.5:
                try:
                    prob_box[5].append(prob[i])
                    pred_box[5].append(pred[i])
                    label_box[5].append(test_label[i])
                except:
                    prob_box[5] = [prob[i]]
                    pred_box[5] = [pred[i]]
                    label_box[5]= [test_label[i]]
            elif prob[i] > 0.4:
                try:
                    prob_box[4].append(prob[i])
                    pred_box[4].append(pred[i])
                    label_box[4].append(test_label[i])
                except:
                    prob_box[4] = [prob[i]]
                    pred_box[4] = [pred[i]]
                    label_box[4]= [test_label[i]]
            elif prob[i] > 0.3:
                try:
                    prob_box[3].append(prob[i])
                    pred_box[3].append(pred[i])
                    label_box[3].append(test_label[i])
                except:
                    prob_box[3] = [prob[i]]
                    pred_box[3] = [pred[i]]
                    label_box[3]= [test_label[i]]
            elif prob[i] > 0.2:
                try:
                    prob_box[2].append(prob[i])
                    pred_box[2].append(pred[i])
                    label_box[2].append(test_label[i])
                except:
                    prob_box[2] = [prob[i]]
                    pred_box[2] = [pred[i]]
                    label_box[2]= [test_label[i]]
            elif prob[i] > 0.1:
                try:
                    prob_box[1].append(prob[i])
                    pred_box[1].append(pred[i])
                    label_box[1].append(test_label[i])
                except:
                    prob_box[1] = [prob[i]]
                    pred_box[1] = [pred[i]]
                    label_box[1]= [test_label[i]]
            else:
                try:
                    prob_box[0].append(prob[i])
                    pred_box[0].append(pred[i])
                    label_box[0].append(test_label[i])
                except:
                    prob_box[0] = [prob[i]]
                    pred_box[0] = [pred[i]]
                    label_box[0]= [test_label[i]]
    # calculate ECE
    ECE = 0
    for key in label_box.keys():
        if label_box[key] == None:
            continue
        accuarcy = torch.sum(torch.FloatTensor(label_box[key]) == torch.FloatTensor(pred_box[key]))/len(pred_box[key])
        confidence = torch.mean(torch.FloatTensor(prob_box[key]))
        gap = len(pred_box[key])*torch.abs(accuarcy-confidence)
        ECE+=gap
    return ECE.to(device)/len(out)


def masked_loss(out, label, mask):
    #if torch.cuda.is_available():
    #    w = torch.Tensor([3.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 1.0]).cuda()
    #else:
    #    w = torch.Tensor([3.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 1.0])
    #loss = F.cross_entropy(out, label, w, reduction='none')
    loss = F.cross_entropy(out, label, reduction='none')
    #all phage
    #w = torch.Tensor([3.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 1.0, 3.0, 2.0, 3.0]).cuda()
    #loss = F.cross_entropy(out, label, w, reduction='none')
    mask = mask.float()
    mask = mask / mask.mean()
    loss *= mask
    loss = loss.mean()
    return loss


def masked_acc(out, label, mask):
    # [node, f]
    pred = out.argmax(dim=1)
    correct = torch.eq(pred, label).float()
    mask = mask.float()
    mask = mask / mask.mean()
    correct *= mask
    acc = correct.mean()
    return acc



def sparse_dropout(x, rate, noise_shape):
    """

    :param x:
    :param rate:
    :param noise_shape: int scalar
    :return:
    """
    random_tensor = 1 - rate
    random_tensor += torch.rand(noise_shape).to(x.device)
    dropout_mask = torch.floor(random_tensor).byte()
    i = x._indices() # [2, 49216]
    v = x._values() # [49216]

    # [2, 4926] => [49216, 2] => [remained node, 2] => [2, remained node]
    i = i[:, dropout_mask]
    v = v[dropout_mask]

    out = torch.sparse.FloatTensor(i, v, x.shape).to(x.device)

    out = out * (1./ (1-rate))

    return out


def dot(x, y, sparse=False):
    if sparse:
        res = torch.sparse.mm(x, y)
    else:
        res = torch.mm(x, y)

    return res

