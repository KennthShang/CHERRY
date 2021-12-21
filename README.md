![CHERRY](logo.png)
CHERRY is a python library for predicting the interactions between viral and prokaryotic genomes. CHERRY is based on a deep learning model, which consists of a graph convolutional encoder and a link prediction decoder.

# Overview
There are two kind of tasks that CHERRY can work:
1. Host prediction for virus
2. Identifying viruses that infect pathogenic bacteria

Users can choose one of the task when running CHERRY. If you have any trouble installing or using CHERRY, please let us know by opening an issue on GitHub or emailing us (jyshang2-c@my.cityu.edu.hk).

## Required Dependencies
* Python 3.x
* Numpy
* Pytorch>1.8.0
* Networkx
* Pandas
* [Diamond](https://github.com/bbuchfink/diamond)
* BLAST
* MCL
* [Prodigal](https://github.com/hyattpd/Prodigal)

All these packages can be installed using Anaconda.

If you want to use the gpu to accelerate the program:
* cuda
* Pytorch-gpu

### An easiler way to install
We recommend you to install all the package with [Anaconda](https://anaconda.org/)

After cloning this respository, you can use anaconda to install the **CHERRY.yaml**. This will install all packages you need with gpu mode (make sure you have installed cuda on your system to use the gpu version. Othervise, it will run with cpu version). The command is: `conda env create -f CHERRY.yaml`

* For cpu version pytorch: `conda install pytorch torchvision torchaudio cpuonly -c pytorch`
* For gpu version pytorch: Search [pytorch](https://pytorch.org/) to find the correct cuda version according to your computer
*Note*: we suggest you to install all the package using conda (both miniconda and anaconda are ok). We supply a 

### Prepare the database
Due to the limited size of the GitHub, we zip the database. Before using CHEERY, you need to unpack them using the following commands.

```
cd CHEERY/dataset
bzip2 -d protein.fasta.bz2
bzip2 -d nucl.fasta.bz2
cd ../prokaryote
gunzip *
cd ..
```

## Usage
### 1 Predicting host for viruses
If you want to predict hosts for viruses, the input should be a fasta file containing the virual sequences. We support an example file named "test_contigs.fa" in the Github folder. Then, the only command that you need to run is 



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

    python run_Speed_up.py --contigs test_contigs.fa --len 8000 --model pretrain --topk 3
    
    
**OUTPUT**

The format of the output file is a csv file which contain the prediction of each virus. *contig_name* is the accession from the input. 


### 2 Predicting virus infecting prokaryote
If you want to predict hosts for viruses, you need to supply two kinds of inputs:
1. Place your prokaryotic genomes in *new_prokaryote* folder.
2. A fasta file containing the virus squences.
Then, the program will output which virus in your fasta file will infect the prkaryotes in the *new_prokaryote* folder.

The command is simlar to the previous one but one more paramter is need:


    python run_Speed_up.py [--mode PROKARYOTE]
    
**Example**


    python run_Speed_up.py --contigs test_contigs.fa --len 8000 --model pretrain --mode prokaryote

**OUTPUT**

The format of the output file is a csv file which contain the prediction of each virus. *prokaryote* is the accession of your given prokaryotic genomes. *virus* is the list of viruses that might infect these genomes.


# References
The paper is submitted to the *Nucleic Acids Research*.

## Contact
If you have any questions, please email us: jyshang2-c@my.cityu.edu.hk


## Notes
1. if the program output an error (which is caused by your machine):
`Error: mkl-service + Intel(R) MKL: MKL_THREADING_LAYER=INTEL is incompatible with libgomp.so.1 library.`
You can type in the command `export MKL_SERVICE_FORCE_INTEL=1` before runing *run_Speed_up.py*
2. Since the topk method is given, we cannot give the how taxaonmic tree for each prediction. However, we will supply a script for you to convert the prediction into a complte taxonmoy tree. (will be avaliable soon)

