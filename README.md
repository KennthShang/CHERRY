![CHERRY](logo.png)
CHERRY is a python library for predicting the interactions between viral and prokaryotic genomes. CHERRY is based on a deep learning model, which consists of a graph convolutional encoder and a link prediction decoder.


## News !!!

1. If you want to use cherry on your own bacterial assemblies. Please visit [https://github.com/KennthShang/CHERRY_MAGs](https://github.com/KennthShang/CHERRY_MAGs/tree/main) to check. This MAGs version allows you to use your own bacterial assemblies and predict the interactions between your phages and your bacteria.

2. Only CRISPR version of CHERRY is available now. Please visit [https://github.com/KennthShang/CHERRY_crispr](https://github.com/KennthShang/CHERRY_crispr) to check. This lite version only uses the CRISPRs information for prediction (low number of predicted phages but high precision)

3. This folder will no longer be maintained. The program has been updated and moved to PhaBOX [https://github.com/KennthShang/PhaBOX], which is more user-friendly. Hope you will enjoy it.

4. Our web server for phage-related tasks (including phage identification, taxonomy classification, lifestyle prediction, and host prediction) is available! You can visit [https://phage.ee.cityu.edu.hk/] to use the GUI. We also provided more detailed intermediate files and visualization for further analysis. You can also find the local version in [https://github.com/KennthShang/PhaBOX]

# Overview
There are two kind of tasks that CHERRY can work:
1. Host prediction for virus
2. Identifying viruses that infect pathogenic bacteria


## Required Dependencies

### Easy way to install
*Note*: we suggest you to install all the package using conda (both miniconda and [Anaconda](https://anaconda.org/) are ok)

After cloning this respository, you can use anaconda to install the **CHERRY.yaml**. The command is: `conda env create -f CHERRY.yaml -n cherry`


### Prepare the database
Due to the limited size of the GitHub, we zip the database. Before using CHEERY, you need to unpack them using the following commands.

```
cd CHEERY
conda env create -f CHERRY.yaml -n cherry
conda activate cherry
cd dataset
bzip2 -d protein.fasta.bz2
bzip2 -d nucl.fasta.bz2
cd ../prokaryote
gunzip *
cd ..
```

You only need to activate your 'cherry' environment before using CHERRY in the next time.

```
conda activate cherry
```

## Usage
### 1 Predicting host for your viruses
The input should be a fasta file containing the viral sequences. We provide an example file named "test_contigs.fa". Then, the only command that you need to run is 

    python run_Speed_up.py [--contigs INPUT_FA] [--len MINIMUM_LEN] [--model MODEL] [--topk TOPK_PRED]
    
**Options**


      --contigs INPUT_FA
                            input fasta file
      --len MINIMUM_LEN
                            predict only for sequence >= len bp (default 8000)
      --model MODEL (pretrain or retrain)
                            predicting host with pretrained parameters or retrained paramters (default pretrain)
      --topk TOPK_PRED
                            The host prediction with topk score (default 1)
               

**Example**

Prediction on species level with pretrained paramters:

    python run_Speed_up.py --contigs test_contigs.fa --len 8000 --model pretrain --topk 1
    
*Note:* Commonly, you do not need to retrain the model, especially when you do not have gpu unit. 
    
**OUTPUT**

The format of the output file is a csv file ("final_prediction.csv") which contain the prediction of each virus. Column *contig_name* is the accession from the input. 

We will supply a script for you to convert the prediction into a complte taxonmoy tree. Use the following command to generate taxonomy tree:

    python run_Taxonomy_tree.py [--k TOPK_PRED]
    
Because there are k prediction in the "final_prediction.csv" file, you need to specify the k to generate the tree. The output of program is 'Top_k_prediction_taxonomy.csv'.


### Extension of the virus-prokaryote interactions database

If you know more virus-prokaryote interactions than our pre-trained model (given in Interactiondata), you can add them to train a custom model. Several steps you need to do to train your model:

1. Add your viral genomes into the *nucl.fasta* file and run the *python refresh.py* to generate new *protein.fasta* and *database_gene_to_genome.csv* files. They will replace the old one in the *dataset/* folder automatically. 
2. Add the entrys of host taxonomy information into *dataset/virus.csv*. The corresponding header of the entry is: Accession (of the virus), Superkingdom, Phylum, Class, Order, Family, Genus, Species. The required field is **Species**. You can left it blank if you do not know other fields. Also, the accession of the virus shall be the same as your fasta entry. 
3. Place your prokaryotic genomes into the the *prokaryote/* folder and add an entry in *dataset/prokaryote.csv*. The guideline is the same as the previous section.
4. Use **retrain** as the parameter for *--mode* option to run the program.




### 2 Predicting virus infecting prokaryote
If you want to predict candidate viruses that infect a set of given bacteria, you need to supply three kinds of inputs:
1. Place your prokaryotic genomes in *new_prokaryote/* folder.
2. A fasta file containing the virus squences.
3. Add the taxa information in 'database/prokaryote.csv'. (The example can be found in the *Extension of the parokaryotic genomes database*)
Then, the program will output which virus in your fasta file will infect the prkaryotes in the *new_prokaryote/* folder.

The command is simlar to the previous one but two more paramter is need:


    python run_Speed_up.py [--mode MODE] [--t THRESHOLD]
    
**Example**


    python run_Speed_up.py --contigs test_contigs.fa --mode prokaryote --t 0.98


**Options**


      --mode MODE (prokaryote or virus)
                            Switch mode for predicting virus or predicting host
      --t THRESHOLD
                            The confident threshold for predicting virus, the higier the threshold the higher the precision. (default 0.98)
                            
                   
**OUTPUT**

The format of the output file is a csv file which contain the prediction of each virus. Column *prokaryote* is the accession of your given prokaryotic genomes. Column *virus* is the list of viruses that might infect these genomes.





### References
The paper published in the *Briefings in Bioinformatics*: 
```
Jiayu Shang, Yanni Sun, CHERRY: a Computational metHod for accuratE pRediction of virus–pRokarYotic interactions using a graph encoder–decoder model, Briefings in Bioinformatics, 2022;, bbac182, https://doi.org/10.1093/bib/bbac182
```

The arXiv version can be found via: [CHERRY: a Computational metHod for accuratE pRediction of virus-pRokarYotic interactions using a graph encoder-decoder model](https://arxiv.org/abs/2201.01018)

### Contact
If you have any questions, please email us: jyshang2-c@my.cityu.edu.hk


### Notes
1. if the program output an error (which is caused by your machine):
`Error: mkl-service + Intel(R) MKL: MKL_THREADING_LAYER=INTEL is incompatible with libgomp.so.1 library.`
You can type in the command `export MKL_SERVICE_FORCE_INTEL=1` before runing *run_Speed_up.py*


