## A flexible document organization process 

A tool for perform text clustering and cluster labelling.

#### Features 

* Clustering methods (FCM, PCM, PFCM) 
* Distance norms (Euclidian, Cosine, Jaccard)
* Cluster labelling/Descriptor extractor methods (SoftO-FDCL, PDCL, Mixed-PFDCL)
* Generate ARFF file after process, thus its possible just submit for classification benchmarks in WEKA tool
* Parameters selection
* Input data should be in programming contest problems format


#### Setup

    make all  

#### Run

Show help

    ./clustering --help 

Executing process

    ./clustering -x -m 1.2 -n 1.2 -k < X12.in 


#### Input format

**N** is the number of terms<br> 
**M** is the number of documents


    N M
    term1
    term2
    ...
    termN
    value11 value12 ... value1N
    value21 value22 ... value2N
    ...
    valueM1 valueM2 ... valueMN
    
#### Ploting results with Sammon's Mapping

    Rscript sammons.r samples/X12.frequencys X12 .
    
![X12 clustered](https://github.com/niltonvasques/fuzzy-text-clustering/blob/master/samples/X12.png)
