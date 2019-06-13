# Single Cell Analysis Applications in Python

This repository contains applications of [Scanpy](https://github.com/theislab/scanpy) used by the Jason Spence Lab. **basic_analysis_script.py** contains a basic outline to run a single cell RNA sequencing analysis. **scanpy_spence.py** contains our collection of functions used to run the analysis.

## Getting Started

To run this project, install Scanpy from [here](https://github.com/theislab/scanpy). Make sure to have a recent version of Python. As of June 2019, we are using Python version 3.7.3.

To understand the single cell analysis tools, check out the Scanpy and Seurat tutorials side-by-side [here](https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html) and [here](https://satijalab.org/seurat/v3.0/pbmc3k_tutorial.html)

## Data Storage

We keep track of our data using a meta-data table in a .tsv file. Each sample is represented in the table with the following format

```
<Sample ID Number> <path/from/storage_mount_point/to/raw_data_matrix.h5> <metadata fields>
```
**Example**
```
2182-1 01_RNAseq_RAW_Data/Run_2182_Czerwinski_HIO_fetal_intestine_scRNAseq/Data/Intensities/BaseCalls/1-HIO-Fresh/outs/filtered_gene_bc_matrices_h5.h5 age:30 tissue:HIO gel:Matrigel media:ENR sex:Male sampleName:1-HIO-Fresh
```


<!-- ### Prerequisites

What things you need to install the software and how to install them

```
Give examples
```

### Installing

A step by step series of examples that tell you how to get a development env running

Say what the step will be

```
Give the example
```

And repeat

```
until finished
```

End with an example of getting some data out of the system or using it for a little demo

## Running the tests

Explain how to run the automated tests for this system

### Break down into end to end tests

Explain what these tests test and why

```
Give an example
```

### And coding style tests

Explain what these tests test and why

```
Give an example
``` -->

## Authors

* **Joshua Wu** wujos@med.umich.edu
* **Mike Czerwinski** czerwmj@med.umich.edu

Any questions, comments and/or feedback are welcome!

## Acknowledgments

* [**Jason Spence Lab**](http://www.jasonspencelab.com/) for support from our colleagues
* [**Theis Lab**](https://github.com/theislab) and [**Satija Lab**](https://satijalab.org/) for their work on modern single cell analysis techniques
