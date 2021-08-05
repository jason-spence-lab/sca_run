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
			umap_check: Run UMAP embedding if true
			spread: In combination with min_dist determines how clumped embedded points are
			min_dist: Minimum distance between points on the umap graph
			resolution: High resolution attempts to increases # of clusters identified
			do_bbknn: Run batch balanced k-nearest neighbors batch correction algorithm
			do_tSNE: Run tSNE dimensional reduction analysis
			clustering_choice: Run leiden clustering or louvain clustering based on user's choice. Default is leiden clustering 
			dpt: Run diffusion pseudotime analysis with input: ['metadata_cat',['group']]
			draw_force_atlas: Run force-directed graphing of data
			umap_init_pos: Basis from which to initiate umap embedding ('paga','spectra','random')
			phate: Run Potential of Heat-diffusion for Affinity-based Trajectory Embedding (PHATE)
		'''
		self.n_neighbors = analysis_params.n_neighbors
		self.n_pcs = analysis_params.n_pcs
		self.umap_check = analysis_params.umap_check
		self.spread = analysis_params.spread
		self.min_dist = analysis_params.min_dist
		self.resolution = analysis_params.resolution
		self.do_bbknn = analysis_params.do_bbknn
		self.do_tSNE = analysis_params.do_tSNE
		self.clustering_choice = analysis_params.clustering_choice
		self.dpt = analysis_params.dpt
		self.draw_force_atlas = analysis_params.draw_force_atlas
		self.umap_init_pos = analysis_params.umap_init_pos
		self.phate = analysis_params.phate

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
			bbknn.bbknn(adata, batch_key=self.do_bbknn, copy=False)#, 
						# n_pcs=self.n_pcs, neighbors_within_batch=self.n_neighbors)
			#sc.pp.external.mnn_correct(adata,batch_key='sampleName') # Testing another algorithm
		else:
			## Compute nearest-neighbors
			sc.pp.neighbors(adata, n_neighbors=self.n_neighbors, n_pcs=self.n_pcs)

		## Calculate cell clusters via the chosen clustering algorithm
		getattr(sc.tl, self.clustering_choice)(adata, resolution=self.resolution)


		# ## Run PAGA to predict non-directional cluster-cluster relationships to infer possible developmental progressions
		sc.tl.paga(adata, groups=self.clustering_choice, model='v1.2')

		# Set the thresholds and scaling factors for drawing the paga map/plot
		node_size_scale=1.25
		node_size_power=0.9
		edge_width_scale=1
		min_edge_width=0.035
		max_edge_width=2
		threshold=0.05
		sc.pl.paga(adata, layout='fr', threshold=threshold, node_size_scale=node_size_scale, 
			node_size_power=node_size_power, edge_width_scale=edge_width_scale,
			min_edge_width=min_edge_width, max_edge_width=max_edge_width, show=False, save=False,
			title='PAGA: Fruchterman Reingold',frameon=False)

		## Run UMAP Dim reduction
		if self.umap_check:
			sc.tl.umap(adata, spread=self.spread, min_dist=self.min_dist, init_pos=self.umap_init_pos)#, n_components=50) # Min_dist needs to be between 0.01 to 0.5

		## Run tSNE analysis
		if self.do_tSNE:
			sc.tl.tsne(adata, n_pcs=self.n_pcs)

		# Graph Force Atlas
		if self.draw_force_atlas:
			adata_fa = adata.raw.to_adata().copy()
			sc.pp.highly_variable_genes(adata_fa, min_mean=0.0125, max_mean=3, min_disp=0.5)
			adata_fa = adata_fa[:, adata_fa.var['highly_variable']].copy()

			if self.do_bbknn:
				import bbknn
				bbknn.bbknn(adata_fa, batch_key=self.do_bbknn, copy=False)#, 
							# n_pcs=self.n_pcs, neighbors_within_batch=self.n_neighbors)
				#sc.pp.external.mnn_correct(adata,batch_key='sampleName') # Testing another algorithm
			else:
				## Compute nearest-neighbors
				sc.pp.neighbors(adata_fa, n_neighbors=self.n_neighbors, n_pcs=self.n_pcs)

			adata_fa.uns['paga'] = adata.uns['paga']
			sc.tl.draw_graph(adata_fa, init_pos='paga')
			adata.obsm['X_draw_graph_fa'] = adata_fa.obsm['X_draw_graph_fa']
			adata.uns['draw_graph'] = adata_fa.uns['draw_graph']

		## Run diffusion pseudotime analysis
		if self.dpt:
			adata_dpt = adata.raw.to_adata().copy()
			# sc.pp.neighbors(adata_dpt, n_neighbors=self.n_neighbors, n_pcs=self.n_pcs, method='gauss')

			if self.do_bbknn:
				import bbknn
				bbknn.bbknn(adata_dpt, batch_key=self.do_bbknn, copy=False)#, 
							# n_pcs=self.n_pcs, neighbors_within_batch=self.n_neighbors)
				#sc.pp.external.mnn_correct(adata,batch_key='sampleName') # Testing another algorithm
			else:
				## Compute nearest-neighbors
				sc.pp.neighbors(adata_dpt, n_neighbors=self.n_neighbors, n_pcs=self.n_pcs, method='gauss')

			adata_dpt.uns['iroot'] = np.flatnonzero(adata.obs[self.dpt[0]].isin(self.dpt[1]))[0]
			print(adata_dpt.uns['iroot'])
			sc.tl.diffmap(adata_dpt, n_comps=20)
			sc.tl.dpt(adata_dpt)
			adata.uns['iroot'] = adata_dpt.uns['iroot']
			adata.uns['diffmap_evals'] = adata_dpt.uns['diffmap_evals']
			adata.obsm['X_diffmap'] = adata_dpt.obsm['X_diffmap']
			adata.obs['dpt_pseudotime'] = adata_dpt.obs['dpt_pseudotime']

		if self.phate:
			print("Running PHATE")
			import phate
			X_phate = phate.PHATE(n_jobs=-2, knn=3).fit_transform(adata)
			adata.obsm['X_phate'] = X_phate

		## Do Dendrogram analysis based on PCs
		sc.tl.dendrogram(adata, groupby=self.clustering_choice, n_pcs=self.n_pcs, linkage_method="median", use_raw=True)
		return adata