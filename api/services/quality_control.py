'''
quality_control -- 
Service that quality controls raw single-cell data and filters out unwanted data points

Written by Joshua H Wu
09 March, 2021
'''
import scanpy as sc
import numpy as np

class quality_control:
	'''
	Service that quality controls raw single-cell data and filters out unwanted data points
	'''
	def __init__(self, qc_params, species):
		'''
		Quality Control Params --
			min_cells: Filter out genes with a few number of cells
			min_genes: Filter out cells with fewer genes to remove dead cells
			max_genes: Filter out cells with more genes to remove most doublets
			max_counts: Filter out cells with more UMIs to catch a few remaining doublets
			max_mito: Filter out cells with high mitochondrial gene content
			doublet_detection: Run DoubletDetection by Jonathan Shor

		species: The species of samples used, examples: human, mouse
		'''
		self.min_cells = qc_params.min_cells
		self.min_genes = qc_params.min_genes
		self.max_genes = qc_params.max_genes
		self.max_counts = qc_params.max_counts
		self.max_mito = qc_params.max_mito
		self.doublet_detection = qc_params.doublet_detection
		self.species = species

		self.doublet_clf = None
		self.adata_preQC = None
		self.adata_doublet = None

	## Filters data based on certain parameters
	# Attempts to remove "bad" data such as dead cells, doublets, etc.
	def run_qc(self, adata):
		'''
		Removes cells expressing low to no genes, and genes expressed in few to no cells
		Filters out cells based on mitochondrial genes, UMI and max gene expression
		'''
		print("Conducting quality control")

		# Conducting DoubletDetector analysis by Jonathan Shor
		if self.doublet_detection:
			print("Starting doublet detection")
			import DoubletDetection as doubletdetection
			self.doublet_clf = doubletdetection.BoostClassifier(n_iters=50, use_phenograph=False, standard_scaling=True)
			adata.obs['doublet_label'] = self.doublet_clf.fit(adata.X).predict(p_thresh=1e-16, voter_thresh=0.5)
			adata.obs['doublet_score'] = self.doublet_clf.doublet_score()
			self.adata_doublet = adata.copy()
			if self.species=='human':
				self.adata_doublet.var["mito"] = self.adata_doublet.var_names.str.startswith('MT-')
			elif self.species=='mouse':
				self.adata_doublet.var["mito"] = self.adata_doublet.var_names.str.startswith('mt-')
			else:
				print("No valid species name - assuming human")
				self.adata_doublet.var["mito"] = self.adata_doublet.var_names.str.startswith('MT-')
			sc.pp.calculate_qc_metrics(self.adata_doublet, qc_vars=["mito"], inplace=True)	

		## Basic filtering to get rid of useless cells and unexpressed genes
		sc.pp.filter_genes(adata, min_cells=self.min_cells)
		sc.pp.filter_cells(adata, min_genes=self.min_genes)

		# Calculate the percent of genes derived from mito vs genome
		# the `.A1` is only necessary as X is sparse (to transform to a dense array after summing)
		if self.species=='human':
			adata.var["mito"] = adata.var_names.str.startswith('MT-')
			mito_genes = adata.var_names.str.startswith('MT-')
		elif self.species=='mouse':
			adata.var["mito"] = adata.var_names.str.startswith('mt-')
			mito_genes = adata.var_names.str.startswith('mt-')
		else:
			print("No valid species name - assuming human")
			adata.var["mito"] = adata.var_names.str.startswith('MT-')
			mito_genes = adata.var_names.str.startswith('MT-')

		# sc.pp.calculate_qc_metrics(adata, qc_vars=["mito"], inplace=True)

		try:
			adata.obs['percent_mito'] = np.sum(adata[:,mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
			# add the total counts per cell as observations-annotation to adata
			adata.obs['n_counts'] = adata.X.sum(axis=1).A1
		except:
			adata.obs['percent_mito'] = np.sum(adata[:,mito_genes].X, axis=1) / np.sum(adata.X, axis=1)

			# add the total counts per cell as observations-annotation to adata
			adata.obs['n_counts'] = adata.X.sum(axis=1)


		self.adata_preQC = adata.copy() # Saving pre-Filtered AnnData

		## Actually do the filtering.
		if self.doublet_detection:
			adata = adata[((adata.obs['doublet_label'] != 1)
						& (adata.obs['pct_counts_mito'] < self.max_mito))].copy()
		else:
			adata = adata[((adata.obs['n_genes'] < self.max_genes)   # Keep cells with less than __ genes to remove most doublets
						& (adata.obs['n_counts'] < self.max_counts)   # Keep cells with less than __ UMIs to catch a few remaining doublets
						& (adata.obs['percent_mito'] < self.max_mito))].copy()   # Keep cells with less than __ mito/genomic gene ratio

		return adata