'''
preprocess -- 
Initial preprocessing of AnnData object for downstream analyses

Written by Joshua H Wu
09 March, 2021
'''
import scanpy as sc
import numpy as np

class preprocess:

	def __init__(self, 
				 gene_dict,
				 pp_params):
		'''
		gene_dict: Dictionary of genes_lists to plot and score
		Preprocessing Params --
		    combat: Run combat batch correction across the specified metadata field
		    min_mean: Low cutoff for feature (gene) means
		    max_mean: High cutoff for feature (gene) means
		    min_disp: Low cutoff for feature (gene) dispersions
		    regress_out: Variables to regress out, for example, percent_mito
		'''
		self.gene_dict = gene_dict
		self.combat = pp_params.combat
		self.min_mean = pp_params.min_mean
		self.max_mean = pp_params.max_mean
		self.min_disp = pp_params.min_disp
		self.regress_out = pp_params.regress_out
		self.adata_unscaled = None

	## Standardize and normalize the data set
	# Includes additional processing such as removing biased and uninformative data
	def run_preprocess(self,
					   adata):
		print("Preprocessing data")
		## Normalize the expression matrix to median reads per cell, so that counts become comparable among cells.
		# This corrects for differences in sequencing depth between cells and samples
		sc.pp.normalize_total(adata)#,target_sum=10000)

		## Log transform the data.
		sc.pp.log1p(adata)

		## Set the .raw attribute of AnnData object to the logarithmized raw gene expression for later use in differential testing and visualizations of gene expression.
		# We need to do this because the expression matrix will be rescaled and centered which flattens expression too much for some purposes
		adata.raw = adata.copy()

		## Find cell type score for each cell based on a predefined set of gene lists
		if any([self.gene_dict[key].cell_score_list for key in self.gene_dict.keys()]):
			adata_scaled = sc.pp.scale(adata, max_value=10, copy=True)
			for key in self.gene_dict.keys():
				if self.gene_dict[key].cell_score_list:
					adata.obs[key] = adata_scaled.X[:,adata_scaled.var_names.isin(self.gene_dict[key].markers)].mean(1)

		## Identify highly-variable genes based on dispersion relative to expression level
		sc.pp.highly_variable_genes(adata, min_mean=self.min_mean, max_mean=self.max_mean, min_disp=self.min_disp)

		## Filter the genes to remove non-variable genes since they are uninformative
		adata = adata[:, adata.var['highly_variable']].copy()

		## Regress out effects of total reads per cell and the percentage of mitochondrial genes expressed.
		sc.pp.regress_out(adata, self.regress_out)

		if not self.combat is 'None':
			print("Conducting combat batch correction")
			sc.pp.combat(adata, key=self.combat)

		self.adata_unscaled = adata.copy()

		print('\nDoing, final filtering...\nKeeping', len(adata.obs_names),'cells and', len(adata.var_names),'genes.\n')

		## Scale each gene to unit variance. Clip values exceeding standard deviation 10 to remove extreme outliers
		sc.pp.scale(adata, max_value=10)
		return adata.copy()