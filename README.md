# Gene Regulatory Network construction with kernel CCA
--
Here is a condition-specific gene regulatory network (GRN) construction method from transcriptome data, utilizing kernel canonical correlation (kernel CCA). 

![workflow](https://github.com/DabinJeong/GRN_construction_with_kernelCCA/blob/master/workflow/figN_workflow.jpg?raw=true)

The novelty and advantages of our method are listed as follows.

* Complex relationship of multiple transcription factors (TFs) regulating multiple target genes (TGs) are modeled
* Different combinations of co-working TFs specific to certain conditions (e.g. time-series responses against heat stress) are modeled
*  Our method outputs GRN in a feasible size of sub-networks, where regulatory relations from co-working TFs and co-regulated TGs are suggested as a subnetwork.

## Usage
`snakemake --configfile config.yaml --cores 1`
Since the method is composed of multiple scripts, the method is implemented with workflow management system, Snakemake.
config.yaml contains hyperparameters and informations about datasets. You can run the method with your dataset, modifying config.yaml files. Input data should be stored in data folder. 

## Dependencies

**python 3.8.2**

* [snakemake 5.31.1](https://github.com/snakemake/snakemake)
* [pandas 1.0.4](https://pandas.pydata.org)
* [pyrcca](https://github.com/gallantlab/pyrcca)

<br>

**R 4.0.2**

* [dplyr 2.1.0](https://dplyr.tidyverse.org)
* [igraph](https://igraph.org/r/)

<br>

## References



