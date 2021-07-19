'''
Applications of single cell data analysis techniques
Written by Josh Wu, Mike Czerwinski, and Zhiwei Xiao
7 July, 2021

Relies heavily on the Scanpy Python module developed by the Theis Lab
Read more about Scanpy at https://scanpy.readthedocs.io/en/latest/index.html

'''

import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from pathlib import Path
import csv
import copy
import os
import json
from matplotlib import rcParams
from matplotlib import pyplot as plt
import matplotlib as mpl
import matplotlib.pyplot as plt
import random
import math
import _pickle as pickle
import classes.sca_params as sca_params
import services.load_data as load_data
import services.quality_control as quality_control
import services.preprocess as preprocess
import services.tools as tools
import services.plotting as plotting

class SCARunner:
	sc.settings.verbosity = 3 # verbosity: errors (0), warnings (1), info (2), hints (3)
	'''
	Organizes single-cell functions into distinct pipelines
	'''

	## Most basic pipeline - Input data and output all figures
	def pipe_basic(self,
				   sca_params,
				   figdir='./figures/',
				   adata_filtered=None,
				   adata_loaded=None,
				   load_save=None, 
				   new_save='adata_save.p',
				   remove_genes=None,
				   only_plot=False,
				   load_file_type='.h5'):
		'''
		sca_params: Class that handles all relevant parameters for setting up a SCARunner session
		figdir: The path for saving figures generated during the analysis
		adata_filtered: AnnData that's already been quality controlled 
		adata_loaded: AnnData that already exists, eg., other lab's AnnData downloaded online
		load_save: AnnData saved for the analysis to duplicate
		new_save: Pickle module that contains analysis information and relevant AnnData objects
		remove_genes: List of unnecessary genes that should not be used in analysis
		only_plot: Choose to only plot graphs with existing AnnData objects

		Uses the pickle module to save adata instances for easy access in future analyses
		*Note that pickled adata files are very large - min 4GB (make sure to clean out unused pickle files)
		'''
		if load_save: # See if there is already a save file for the analysis to duplicate
			run_save = pickle.load(open(''.join([figdir,load_save]),"rb"))

			adata = run_save.adata.copy()
			# sca_params.vmin_list = run_save.vmin_list
			# sca_params.vmax_list = run_save.vmax_list
			sca_params.adata_preQC = run_save.adata_preQC
			sca_params.adata_unscaled = run_save.adata_unscaled
			sca_params.adata_postQC = run_save.adata_postQC
			sca_params.annotation_dict = run_save.annotation_dict
			sca_params.initial_cell_count = run_save.initial_cell_count
			sca_params.initial_gene_count= run_save.initial_gene_count
		else:
			if not adata_filtered:
				if adata_loaded:
					adata = adata_loaded.copy()
				else:
					ld = load_data.load_data(storage_mount_point = sca_params.storage_mount_point,
											 sample_list = sca_params.sample_list,
											 remove_genes = remove_genes)
					adata = ld.load(file_type = load_file_type).copy()
					sca_params.initial_cell_count = ld.initial_cell_count
					sca_params.initial_gene_count = ld.initial_gene_count
					sca_params.annotation_dict = ld.annotation_dict

				## Filter and process data
				qc = quality_control.quality_control(sca_params.qc_params, sca_params.species)
				adata = qc.run_qc(adata).copy()
				sca_params.doublet_clf = qc.doublet_clf
				sca_params.adata_preQC = qc.adata_preQC.copy()
				
				if sca_params.qc_params.doublet_detection:
					adata_doublet = qc.adata_doublet
					pp = preprocess.preprocess(sca_params.gene_dict, sca_params.pp_params)
					adata_doublet = pp.run_preprocess(adata_doublet)
					tl = tools.tools(sca_params.analysis_params)
					adata_doublet = tl.run_tools(adata_doublet)
					sca_params.adata_doublet = adata_doublet.copy()
			else: 
				adata = adata_filtered.copy()

			sca_params.adata_postQC = adata.copy() # Save for later in case of necessary extraction
			pp = preprocess.preprocess(sca_params.gene_dict, sca_params.pp_params)
			adata = pp.run_preprocess(adata).copy()
			sca_params.adata_unscaled = pp.adata_unscaled.copy()



		if not only_plot:
			## Dimensional reduction and clustering - construction of the neighborhood graph
			tl = tools.tools(sca_params.analysis_params)
			adata = tl.run_tools(adata)



		## Plot figures
		scplt = plotting.plotting(sca_params.plot_params)
		adata = scplt.plot_sca(adata,sca_params,figdir=figdir)

		sca_params.adata = adata.copy()
		## Write a summary of the analysis to a text file including sample information and parameters
		sca_params.write_summary(figdir=figdir)

		## Save analysis information and relevant AnnData objects to the disk using the Pickle module
		if new_save:
			pickle.dump(sca_params,open(''.join([figdir,new_save]),"wb"),protocol=4)

		print("\nAll done!\n")

		return self

	## Pipeline for analysis in which you extract interesting clusters/observations after an initial run
	# Extracts clusters to an filtered but unprocessed AnnData object, then reprocesses and reclusters
	def pipe_ext(self, sca_params, extracted, figdir='./figures/', load_save=None, new_save='extracted_adata_save.p', label='',
				 no_preprocess=False):
		'''
		Allows loading of a saved pickle adata file, file must contained adata that has gone through a complete pipeline
		Otherwise will complete a new full analysis, and then extracts clusters

		sca_params: Class that handles all relevant parameters for setting up a SCARunner session
		extracted: List of the numbers of interested clusters
		figdir: The path for saving figures generated during the analysis
		load_save: AnnData saved for the analysis to duplicate
		new_save: Pickle module that contains analysis information and relevant AnnData objects from the extraction
		label: Name of the file folder that contains the output files generated in the extraction
		no_preprocess: Indicator of whether the AnnData object is processed or not
		'''
		## Ask user for input on which clusters to extract
		if extracted is None:
			print("Which clusters would you like to extract? (Separate your answer with commas) ")
			extracted = list(map(str,input().split(',')))

		if not label:
			label='_'.join(extracted)

		if not np.any([not cluster.isdigit() for cluster in extracted]):
			## Load pickle file if one is given
			run_save = pickle.load(open(''.join([figdir,load_save]),"rb"))

		
			adata = sca_params.adata = run_save.adata
			# self.vmin_list = run_save.vmin_list
			# self.vmax_list = run_save.vmax_list
			sca_params.adata_preQC = run_save.adata_preQC
			sca_params.adata_unscaled = run_save.adata_unscaled
			sca_params.adata_postQC = run_save.adata_postQC
			sca_params.annotation_dict = run_save.annotation_dict
			sca_params.initial_cell_count = run_save.initial_cell_count
			sca_params.initial_gene_count= run_save.initial_gene_count

			if not no_preprocess:
				## Create an unprocessed AnnData object with the desired clusters
				adata_ext = sca_params.adata_postQC[adata.obs[sca_params.analysis_params.clustering_choice].isin(extracted)].copy()
				sca_params.adata_postQC = adata_ext.copy()

				## Reprocess and recluster extracted cells
				pp = preprocess.preprocess(sca_params.gene_dict, sca_params.pp_params)
				adata_ext = pp.run_preprocess(adata_ext)
				sca_params.adata_unscaled = pp.adata_unscaled.copy()
			else:
				sca_params.adata_postQC = sca_params.adata_postQC[adata.obs[sca_params.analysis_params.clustering_choice].isin(extracted)].copy()
				# Need to add save sca_params.adata_unscaled
				adata_ext = sca_params.adata[adata.obs[sca_params.analysis_params.clustering_choice].isin(extracted)].copy()
		else:
			print(extracted, "is not a valid input")
			return

		tl = tools.tools(sca_params.analysis_params)
		adata_ext = tl.run_tools(adata_ext)

		## Plot figures
		scplt = plotting.plotting(sca_params.plot_params)
		adata_ext = scplt.plot_sca(adata_ext,sca_params,figdir = ''.join([figdir,'extracted/',label,'/']))

		sca_params.adata = adata_ext.copy()
		## Write a summary of the analysis to a text file including sample information and parameters
		sca_params.write_summary(figdir=''.join([figdir,'extracted/',label,'/']), extracted=extracted)

		## Save analysis information and relevant AnnData objects to the disk using the Pickle module
		if new_save:
			pickle.dump(sca_params,open(''.join([figdir,'extracted/',label,'/',new_save]),"wb"),protocol=4)

		print("\nAll done!\n")

		return self

	def pipe_cell_rank(self, sca_params, adata_filtered=None, figdir='./figures/', load_save=None):
		import scvelo as scv 
		import cellrank as cr 

		scv.settings.verbosity = 3
		scv.settings.set_figure_params("scvelo")
		cr.settings.verbosity = 2
		# scv.settings.figdir = figdir

		import warnings

		warnings.simplefilter("ignore", category=UserWarning)
		warnings.simplefilter("ignore", category=FutureWarning)
		warnings.simplefilter("ignore", category=DeprecationWarning)

		if adata_filtered:
			adata=adata_filtered.copy()
		else:
			if load_save: # need to read loom files for spliced/unspliced genes - make sure loaded data has that information
				run_save = pickle.load(open(''.join([figdir,load_save]),"rb"))

				adata = run_save.adata.copy()
				# sca_params.vmin_list = run_save.vmin_list
				# sca_params.vmax_list = run_save.vmax_list
				sca_params.adata_preQC = run_save.adata_preQC
				sca_params.adata_unscaled = run_save.adata_unscaled
				sca_params.adata_postQC = run_save.adata_postQC
				sca_params.annotation_dict = run_save.annotation_dict
				sca_params.initial_cell_count = run_save.initial_cell_count
				sca_params.initial_gene_count= run_save.initial_gene_count
			else:
				ld = load_data.load_data(storage_mount_point = sca_params.storage_mount_point,
										 sample_list = sca_params.sample_list)
				adata = ld.load(file_type = '.loom').copy()

			# scv.pl.proportions(adata)
			print(adata)

			qc = quality_control.quality_control(sca_params.qc_params, sca_params.species)
			adata = qc.run_qc(adata).copy()
			sca_params.doublet_clf = qc.doublet_clf
			sca_params.adata_preQC = qc.adata_preQC.copy()

		# ## Normalize the expression matrix to median reads per cell, so that counts become comparable among cells.
		# # This corrects for differences in sequencing depth between cells and samples
		# sc.pp.normalize_total(adata)#,target_sum=10000)

		# ## Log transform the data.
		# sc.pp.log1p(adata)

		# ## Set the .raw attribute of AnnData object to the logarithmized raw gene expression for later use in differential testing and visualizations of gene expression.
		# # We need to do this because the expression matrix will be rescaled and centered which flattens expression too much for some purposes
		# adata.raw = adata.copy()

		# ## Identify highly-variable genes based on dispersion relative to expression level
		# sc.pp.highly_variable_genes(adata, min_mean=sca_params.pp_params.min_mean, max_mean=sca_params.pp_params.max_mean, min_disp=sca_params.pp_params.min_disp)

		# ## Filter the genes to remove non-variable genes since they are uninformative
		# adata = adata[:, adata.var['highly_variable']].copy()

		# scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)

		# ## Run PCA to compute the default number of components
		# sc.tl.pca(adata, svd_solver='arpack')

		# ## Remove batch effects
		# # Note that doing this may override previous sc.pp.neighbors()
		# if sca_params.analysis_params.do_bbknn:
		# 	import bbknn
		# 	bbknn.bbknn(adata, batch_key='sampleName', copy=False)#, 
		# 				# n_pcs=self.n_pcs, neighbors_within_batch=self.n_neighbors)
		# 	#sc.pp.external.mnn_correct(adata,batch_key='sampleName') # Testing another algorithm
		# else:
		# 	## Compute nearest-neighbors
		# 	sc.pp.neighbors(adata, n_neighbors=sca_params.analysis_params.n_neighbors, n_pcs=sca_params.analysis_params.n_pcs)

		# ## Calculate cell clusters via the chosen clustering algorithm
		# getattr(sc.tl, sca_params.analysis_params.clustering_choice)(adata, resolution=sca_params.analysis_params.resolution)

		# sc.tl.umap(
		# 	adata, spread=sca_params.analysis_params.spread, min_dist=sca_params.analysis_params.min_dist, 
		# 	init_pos=sca_params.analysis_params.umap_init_pos)#, n_components=50) # Min_dist needs to be between 0.01 to 0.5

		# scv.pp.moments(adata, n_pcs=None, n_neighbors=None)

		# scv.tl.recover_dynamics(adata, n_jobs=1)

		import pickle 
		# pickle.dump(adata,open(''.join([figdir,'dynamics_adata.p']),"wb"),protocol=4)

		adata = pickle.load(open(''.join([figdir,'dynamics_adata.p']),"rb"))

		scv.tl.velocity(adata, mode="dynamical")
		scv.tl.velocity_graph(adata)
		# scv.pl.velocity_embedding_stream(
		#     adata, basis="umap", legend_fontsize=12, title="", smooth=0.8, min_mass=4,
		#     save='velocity_embedding_stream.png')

		cr.tl.terminal_states(adata, 
			cluster_key=sca_params.analysis_params.clustering_choice,  #weight_connectivities=0.2, 
			softmax_scale=4)
		# cr.pl.terminal_states(adata, save='terminal_states.png')

		cr.tl.initial_states(adata, cluster_key=sca_params.analysis_params.clustering_choice,
			n_states=1)
		# cr.pl.initial_states(adata, discrete=True, save='initial_states.png')

		cr.tl.lineages(adata)
		# cr.pl.lineages(adata, same_plot=False, save='lineages.png')

		k = cr.tl.transition_matrix(
		    adata, weight_connectivities=0.2, softmax_scale=4, show_progress_bar=False
		)
		g = cr.tl.estimators.GPCCA(k)

		g.compute_schur(n_components=4)
		g.compute_macrostates(cluster_key=sca_params.analysis_params.clustering_choice)
		g.set_terminal_states_from_macrostates(["0"])
		g.compute_absorption_probabilities()
		print(g.absorption_probabilities)

		g.compute_lineage_drivers(lineages="0")
		g.lineage_drivers.sort_values("0 corr", ascending=False)

		g.plot_lineage_drivers("0", n_genes=4)

		print(adata)

		# scv.tl.recover_latent_time(adata, root_key="initial_states_probs", end_key="terminal_states_probs")

		# scv.tl.paga(
		#     adata,
		#     groups=sca_params.analysis_params.clustering_choice,
		#     root_key="initial_states_probs",
		#     end_key="terminal_states_probs",
		#     use_time_prior="velocity_pseudotime")

		# cr.pl.cluster_fates(
		#     adata,
		#     mode="paga_pie",
		#     cluster_key="clusters",
		#     basis="umap",
		#     legend_kwargs={"loc": "top right out"},
		#     legend_loc="top left out",
		#     node_size_scale=5,
		#     edge_width_scale=1,
		#     max_edge_width=4,
		#     title="directed PAGA")

		# cr.tl.lineage_drivers(adata)
		print(adata)
		# cr.pl.lineage_drivers(adata, save='lineage_drivers.png')

		return self



	def pipe_phate(self, sca_params, figdir='./figures/',load_save=None):
		import phate
		import scprep
		# Load meta-data table in xls format located at in the storage_mount_point - Change if located elsewhere
		annotation_df = pd.read_excel(''.join([sca_params.storage_mount_point,'01_RNAseq_RAW_Data/single_cell_meta_data_table_excel.xls']),
									  header = None)
		## Creates a dictionary with sample id key, data file location, and relevant metadata
		annotation_dict = dict()
		annotation_dict = annotation_df.set_index(0).T.dropna().to_dict('list')
		print(annotation_dict)

		scprep_array_list = []
		sample_name_list = []
		for sampleID in sca_params.sample_list:
			file_path = ''.join([sca_params.storage_mount_point, annotation_dict[sampleID][0]]).replace("_h5.h5","/hg19").replace(".h5","")
			sample_data = scprep.io.load_10X(file_path, sparse=True, gene_labels='both')
			sample_data = scprep.filter.filter_library_size(sample_data, percentile=20, keep_cells='above')
			sample_data = scprep.filter.filter_library_size(sample_data, percentile=75, keep_cells='below')
			print(sample_data)
			scprep_array_list.append(sample_data)
			sample_name_list.append(str.split(annotation_dict[sampleID][6],':')[1])

		EBT_counts, sample_labels = scprep.utils.combine_batches(
		    scprep_array_list, 
		    sample_name_list,
		    append_to_cell_names=True
		)
		print(EBT_counts)
		print(sample_labels)
		del scprep_array_list # removes objects from memory
		del sample_name_list
		EBT_counts.head()

		EBT_counts = scprep.filter.filter_rare_genes(EBT_counts, min_cells=0)
		EBT_counts = scprep.normalize.library_size_normalize(EBT_counts)

		mito_genes = scprep.select.get_gene_set(EBT_counts, starts_with="MT-") # Get all mitochondrial genes. There are 14, FYI.
		scprep.plot.plot_gene_set_expression(EBT_counts, genes=mito_genes, percentile=90)

		EBT_counts, sample_labels = scprep.filter.filter_gene_set_expression(
		    EBT_counts, sample_labels, genes=mito_genes, 
		    percentile=90, keep_cells='below')

		EBT_counts = scprep.transform.sqrt(EBT_counts)

		phate_operator = phate.PHATE(n_jobs=-2)

		Y_phate = phate_operator.fit_transform(EBT_counts)

		print(Y_phate)

		scprep.plot.scatter2d(Y_phate, c=sample_labels, figsize=(12,8), cmap="Spectral",
                  ticks=False, label_prefix="PHATE")

		return self

	## Pipeline for analysis in which you map labels and embeddings from reference adata to new adata.
	# Extracts clusters to an filtered but unprocessed AnnData object, then reprocesses and reclusters
	def pipe_ingest(self, sca_params, adata, adata_ref, obs = 'leiden', embedding_method = 'umap', figdir='./figures/',
					load_save = None, new_save='ingest_adata_save.p', label=''):
		'''
		sca_params: Class that handles all relevant parameters for setting up a SCARunner session
		adata: The annotated data matrix of shape n_obs × n_vars without labels and embeddings
		adata_ref: The annotated data matrix of shape n_obs × n_vars with labels and embeddings which need to be mapped to adata
		obs: The label of key in adata_ref.obs which need to be mapped to adata.obs, e.g., leiden
		embedding_method: Embeddings in adata_ref which need to be mapped to adata, e.g., umap or pca
		load_save: AnnData saved for the analysis to duplicate, which contains adata and adata_ref
		new_save: Pickle module that contains analysis information and relevant AnnData objects from ingesting
		label: Name of the file folder that contains the output files generated during ingesting
		'''

		if load_save: # See if there is already a save file for the analysis to duplicate
			run_save = pickle.load(open(''.join([figdir,load_save]),"rb"))

			adata = run_save.adata.copy()
			adata_ref = run_save.adata_ref

		## performing ingestion
		sc.tl.ingest(adata, adata_ref, obs=obs, embedding_method = embedding_method)

		## Plot figures
		scplt = plotting.plotting(sca_params.plot_params)
		adata_ingest = scplt.plot_sca(adata,sca_params,figdir = ''.join([figdir,'ingest/',label,'/']))

		sca_params.adata = adata.copy()
		## Write a summary of the analysis to a text file including sample information and parameters
		sca_params.write_summary(figdir=''.join([figdir,'ingest/',label,'/']))

		## Save analysis information and relevant AnnData objects to the disk using the Pickle module
		if new_save:
			pickle.dump(sca_params,open(''.join([figdir,'ingest/',label,'/',new_save]),"wb"),protocol=4)

		print("\nAll done!\n")
