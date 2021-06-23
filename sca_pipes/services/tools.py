'''
tools -- 
Main analysis tools for single-cell RNA seq data

Written by Joshua H Wu
09 March, 2021
'''
import scanpy as sc
import numpy as np

class tools:

	def __init__(self,
				 analysis_params):
		'''
		Analysis Params --
			n_neighbors: Size of the local neighborhood used for manifold approximation
			n_pcs: Number of principal components to use in construction of neighborhood graph
			spread: In combination with min_dist determines how clumped embedded points are
			min_dist: Minimum distance between points on the umap graph
			resolution: High resolution attempts to increases # of clusters identified
			do_bbknn: Run batch balanced k-nearest neighbors batch correction algorithm
			do_tSNE: Run tSNE dimensional reduction analysis
			clustering_choice: Run leiden clustering or louvain clustering based on user's choice. Default is leiden clustering 
			dpt: Run diffusion pseudotime analysis with input: ['metadata_cat','group']
		'''
		self.n_neighbors = analysis_params.n_neighbors
		self.n_pcs = analysis_params.n_pcs
		self.spread = analysis_params.spread
		self.min_dist = analysis_params.min_dist
		self.resolution = analysis_params.resolution
		self.do_bbknn = analysis_params.do_bbknn
		self.do_tSNE = analysis_params.do_tSNE
		self.clustering_choice = analysis_params.clustering_choice
		self.dpt = analysis_params.dpt

	## Run dimensional reduction analysis and clustering using KNN graph
	def run_tools(self,
				  adata):
		## Run PCA to compute the default number of components
		sc.tl.pca(adata, svd_solver='arpack')

		## Save the existing data to disk for later
		## self.adata_postPCA = adata.copy()

		## Remove batch effects
		# Note that doing this may override previous sc.pp.neighbors()
		if self.do_bbknn:
			import bbknn
			bbknn.bbknn(adata, batch_key='sampleName', copy=False)#, 
						# n_pcs=self.n_pcs, neighbors_within_batch=self.n_neighbors)
			#sc.pp.external.mnn_correct(adata,batch_key='sampleName') # Testing another algorithm
		else:
			## Compute nearest-neighbors
			sc.pp.neighbors(adata, n_neighbors=self.n_neighbors, n_pcs=self.n_pcs)

		## Run UMAP Dim reduction
		sc.tl.umap(adata, spread=self.spread, min_dist=self.min_dist)#, n_components=50) # Min_dist needs to be between 0.01 to 0.5

		## Run tSNE analysis
		if self.do_tSNE:
			sc.tl.tsne(adata, n_pcs=self.n_pcs)

		## Calculate cell clusters via the chosen clustering algorithm
		getattr(sc.tl, self.clustering_choice)(adata, resolution=self.resolution)

		## Run diffusion pseudotime analysis
		if self.dpt:
			adata.uns['iroot'] = np.flatnonzero(adata.obs[self.dpt[0]].isin(self.dpt[1]))[0]
			print(adata.uns['iroot'])
			sc.tl.diffmap(adata)
			sc.tl.dpt(adata, n_branchings=0, n_dcs=10)

		## Run PAGA to predict non-directional cluster-cluster relationships to infer possible developmental progressions
		sc.tl.paga(adata, groups=self.clustering_choice, model='v1.2')

		## Do Dendrogram analysis based on PCs
		sc.tl.dendrogram(adata, groupby=self.clustering_choice, n_pcs=self.n_pcs, linkage_method="median", use_raw=True)
		return adata