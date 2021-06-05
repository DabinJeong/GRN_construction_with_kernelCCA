# Gene Regulatory Network construction with kernel CCA
Here is a condition-specific gene regulatory network (GRN) construction method from transcriptome data, utilizing kernel canonical correlation (kernel CCA). 

![workflow](https://github.com/DabinJeong/GRN_construction_with_kernelCCA/blob/master/workflow/figN_workflow.jpg?raw=true)

The novelty and advantages of our method are listed as follows.

* Complex relationship of multiple transcription factors (TFs) regulating multiple target genes (TGs) are modeled
* Different combinations of co-working TFs specific to certain conditions (e.g. time-series responses against heat stress) are modeled
*  Our method outputs GRN in a feasible size of sub-networks, where regulatory relations from co-working TFs and co-regulated TGs are suggested as a subnetwork.

## Usage
For time-series transcriptome data of Homo Sapiens or Arabidopsis Thaliana, we provide you with one-line command for analysis. <br>
`snakemake --configfile config.yaml --cores 1` <br>
Since the method is composed of multiple scripts, the method is implemented with workflow management system, Snakemake.
config.yaml contains hyperparameters and informations about datasets. You can run the method with your dataset, modifying config.yaml files. Input data should be stored in data folder. 
<br><br>
[Input]

* `./data/{taxon}_{dataset identifier}_gene_exp_profile`
<br>: gene expression profile (gene X sample matrix) with all time points concatenated.
* `./data/{taxon}_{dataset identifier}_gene_exp_profile_tp{i}`
<br>: gene expression profile (gene X sample matrix)of each time point i
* `./data/{taxon}_{dataset identifier}_DEG_tp{i}`
<br>: list of genes that are identified as DEG each time point i

<br>
[config.yaml]

* parameters
	1. **nThreads**: number of threads for multi processing
	2. **reg_kernel**: reguluraization coefficient for regularized kernel CCA.
	3. **corr_thr**: edge threshold for condition-specific protein-protein interaction (PPI) network construction in step 1.
	4. **edge_thr**: edge threshold for filtering edges inferred using kernel CCA in step3.
	
* datasets
<br> : consists of multiple dataset analyzed at a time with a snakemake pipeline
	1. 	**GEO_id**: dataset identifier
	2. **taxon**: taxonomy of the dataset, "Human" or "Ara"
	3. **num_tp**: numer of time points in the dataset
	4. **num_comp**: number of canonical components, should be less than number of replicates (or sample size)
	5. **norm**: if True, minMaxScaler is applied to the embedding of gene retreived from kernel CCA.
<br><br>

[Output]

* `./TFTG_cluster/{taxon}_{dataset identifier}_community_TF.tsv`,`./TFTG_cluster/{taxon}_{dataset identifier}_community_TG.tsv` <br>
: Co-working TFs and co-regulated TGs are clustered (represneted as cricles in the figure above) and suggested as tsv files. First column represents genes and second column represents community id that the gene is assigned for.
* `./result/{taxon}_{dataset identifier}_{num_comp}_norm{norm}_GRN_{edge_thr}_tp{i}.pkl` <br>
: results are saved as pickle file, consists of dictionary where keys are (TF community id, TG community id) and values are edges list of GRN sub-networks. 

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

2021, Dabin Jeong et al. Frontiers in Genetics, <a href="https://www.frontiersin.org/articles/10.3389/fgene.2021.652623/abstract" style="text-decoration:none" hover="text_decoration:underline"><font size="+1"><b>Construction of Condition-Specific Gene Regulatory Network using Kernel Canonical Correlation Analysis</b></font></a>


