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
* Pytorch>1.6.0
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

After cloning this respository, you can use anaconda to install the **CHERRY.yaml**. This will install all packages you need with gpu mode (make sure you have installed cuda on your system to use the gpu version. Othervise, it will run with cpu version). The command is: `conda env create -f HostG.yaml`

* For cpu version pytorch: `conda install pytorch torchvision torchaudio cpuonly -c pytorch`
* For gpu version pytorch: Search [pytorch](https://pytorch.org/) to find the correct cuda version according to your computer
*Note*: we suggest you to install all the package using conda (both miniconda and anaconda are ok). We supply a 

