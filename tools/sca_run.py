'''
Applications of single cell data analysis techniques
Written by Josh Wu and Mike Czerwinski
8 September, 2019

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
import matplotlib as mpl
import matplotlib.pyplot as plt
import random
import math
import _pickle as pickle

class sca_run:
	sc.settings.verbosity = 3

	def __init__(self):
		self.storage_mount_point = 'Z:/'
		self.sample_list = []
		self.gene_lists = []
		self.gene_dict = dict()
		
		## Filter Params
		self.min_cells = 0
		self.min_genes = 500
		self.max_genes = 7000
		self.max_counts = 30000
		self.max_mito = 0.1
		
		## Analysis Params
		self.n_neighbors = 30
		self.n_pcs = 15
		self.spread = 1
		self.min_dist = 0.4
		self.resolution = 0.6
		self.remove_batch_effects = False
		self.do_tSNE = False

		## Plot Params
		self.size = 20
		self.umap_obs = ['louvain','sampleName']
		self.exp_grouping = ['louvain']
		self.final_quality = False
		self.umap_color = 'yellow_blue'

		## Summary Params
		self.initial_cell_count = None
		self.initial_gene_count = None
		self.final_cell_count = None
		self.final_gene_count = None
		self.annotation_dict = dict()

		## adata
		self.adata = None
		self.adata_preFiltered = None
		self.adata_preProcessed = None

	def add_gene_list(self, markers=[], label='basic_list', positions=None, groups=None):
		self.gene_lists = self.gene_lists + [label]
		self.gene_dict[label]={'markers':markers,
							  'positions':positions,
							  'groups':groups}

	def set_filter_params(self, min_cells=0, min_genes=500, max_genes=7000, max_counts=30000, max_mito=0.1):
		self.__set_param_attr(locals())

	def set_analysis_params(self, n_neighbors=30, n_pcs=15, spread=1, min_dist=0.4, resolution=0.6):
		self.__set_param_attr(locals())

	def set_plot_params(self, size=20, umap_obs=['louvain','sampleName'], exp_grouping=['louvain'], umap_color='yellow_blue', final_quality=False):
		self.__set_param_attr(locals())

	def __set_param_attr(self,arg_dict):
		for key in arg_dict:
			if key != 'self':
				setattr(self,key,arg_dict[key])

	## Define function to generate a color gradient from a defined starting and ending color
	def __make_cmap(self,colors, position=None, bit=False):
		'''
		make_cmap takes a list of tuples which contain RGB values. The RGB
		values may either be in 8-bit [0 to 255] (in which bit must be set to
		True when called) or arithmetic [0 to 1] (default). make_cmap returns
		a cmap with equally spaced colors.
		Arrange your tuples so that the first color is the lowest value for the
		colorbar and the last is the highest.
		position contains values from 0 to 1 to dictate the location of each color.
		Default sets position 0 of cmap to light gray, and begins proposed gradient
		from >0
		'''
		bit_rgb = np.linspace(0,1,256)
		cdict = {'red':[(0,bit_rgb[220],bit_rgb[220])],'green':[(0,bit_rgb[220],bit_rgb[220])],'blue':[(0,bit_rgb[220],bit_rgb[220])]}
		if position == None:
			position = np.linspace(0.000001,1,len(colors))
		else:
			cdict = {'red':[],'green':[],'blue':[]}
			if len(position) != len(colors):
				sys.exit("position length must be the same as colors")
			elif position[0] != 0 or position[-1] != 1:
				sys.exit("position must start with 0 and end with 1")
		if bit:
			for i in range(len(colors)):
				colors[i] = (bit_rgb[colors[i][0]],
							 bit_rgb[colors[i][1]],
							 bit_rgb[colors[i][2]])
		for pos, color in zip(position, colors):
			cdict['red'].append((pos, color[0], color[0]))
			cdict['green'].append((pos, color[1], color[1]))
			cdict['blue'].append((pos, color[2], color[2]))

		cmap = mpl.colors.LinearSegmentedColormap('my_colormap',cdict,256)
		return cmap

	## Adds metadata from a dictionary to AnnData observations
	def __create_scanpy_anndata(self,storage_mount_point, sampleID, annotation_dict):
		'''
		In:
		storage_mount_point: Data storage mount location
		sampleID: ID numbers of samples from the metadata table (ex: 2235-1)
		annotation_dict: Dictionary of all the sample IDs and metadata

		Out:
		New filled AnnData object
		'''
		metadata_list = annotation_dict[sampleID][1:]
		newAdata = sc.read_10x_h5(''.join([storage_mount_point, annotation_dict[sampleID][0]]))# genome='hg19' or genome='GRCh38'

		## Set gene names to be unique since there seem to be duplicate names from Cellranger
		newAdata.var_names_make_unique()

		## Add metadata for each sample to the observation (cells) annotations in the Anndata objects
		print('\nAdding Metadata to individual samples.\n')
		for field in metadata_list:
			field_list = str.split(field,':')
			meta_name = field_list[0]
			meta_value = field_list[1]
			newAdata.obs[meta_name] = meta_value
		return(newAdata)

	## Swaps the elements at the proposed indices in an applicable data structure
	def __ele_swap(self,structure, index1, index2):
		structure[index1], structure[index2] = structure[index2], structure[index1]
		return structure

	## Writes results of rank genes analysis to multiple csv files, each representing a Louvain cluster
	def __rank_genes(self,adata, groupby, clusters2_compare=None, figdir='./figures/'):
		'''
		groupby: Adata observation metadata categories to compare
		clusters2_compare: Selection of either 2 clusters to compare - if none, then do 1vAll comparison

		Need to add file clearing
		'''

		if clusters2_compare == 'All': # Does 1 to 1 comparison between of all of the clusters
			print("Functionality not developed yet")
			return 0 # Functionality not available yet
		elif clusters2_compare == None: # Do default 1vAll comparison
			print("No clusters selected for comparison. Doing default 1vAll comparison")
			sc.tl.rank_genes_groups(adata,groupby ,method='t-test', rankby_abs=False, n_genes=200)
			self.__write_rank_genes(adata, groupby, figdir)
		else: # Compare 
			adata_temp = adata[adata.obs[groupby].isin(clusters2_compare)]
			sc.tl.rank_genes_groups(adata_temp, groupby, method='t-test', n_genes=200)
			self.__write_rank_genes(adata_temp, groupby, figdir)
		return 0

	## Actually does the writing to csv files of the rank genes analysis
	def __write_rank_genes(self,adata, groupby, figdir='./figures/'):
		rank_genes_data = copy.deepcopy(adata.uns['rank_genes_groups']) # create copy of data to manipulate
		rank_genes_data.pop('params')

		for cluster in adata.obs[groupby].cat.categories:
			csv_fileName = ''.join([figdir,'/csv_files/',cluster,'_',groupby,'_clustercompare.csv'])
			os.makedirs(os.path.dirname(csv_fileName), exist_ok=True) # Make file if it doesn't exist already
			with open(csv_fileName,'w',newline='') as f:
				wr = csv.writer(f)
				wr.writerow(self.__ele_swap(list(rank_genes_data.keys()),0,1))
				wr.writerows(zip(*self.__ele_swap([params[cluster] for params in rank_genes_data.values()],0,1)))

	## Write a summary of the analysis run including sample information, parameters and filtering information
	# Not completely up to date
	def write_summary(self, figdir='./figures/'):

		fileName = ''.join([figdir,'summary.txt'])

		os.makedirs(os.path.dirname(fileName), exist_ok=True) # Create directory if it doesn't exist
		with open(fileName,'w') as f:

			f.write("Summary of single cell sequencing analysis for samples ")
			f.write(''.join([' '.join(self.sample_list),'\n\n']))

			f.write('--------Sample Metadata--------\n')
			for sample in self.sample_list:
				f.write(''.join(['- Sample ',sample,'\n']))
				f.write('\n'.join(self.annotation_dict[sample]))
				f.write('\n\n')

			f.write('--------Basic Run Information--------\n')
			f.write(''.join(['Initial cell count:  ',str(self.initial_cell_count),'\n']))
			f.write(''.join(['Final cell count:  ',str(self.final_cell_count),'\n']))
			f.write(''.join(['Initial gene count:  ',str(self.initial_gene_count),'\n']))
			f.write(''.join(['Final gene count:  ',str(self.final_gene_count),'\n']))

			f.write('\n--------Filter Parameters Used--------\n')
			f.write(''.join(['Min Cells:  ',str(self.min_cells),'\n']))
			f.write(''.join(['Min Genes:  ',str(self.min_genes),'\n']))
			f.write(''.join(['Max Genes:  ',str(self.max_genes),'\n']))
			f.write(''.join(['Max Counts:  ',str(self.max_counts),'\n']))
			f.write(''.join(['Max Mito:  ',str(self.max_mito),'\n']))

			f.write('\n--------Analysis Parameters Used--------\n')
			f.write(''.join(['# Neighbors:  ',str(self.n_neighbors),'\n']))
			f.write(''.join(['# PCs:  ',str(self.n_pcs),'\n']))
			f.write(''.join(['Spread:  ',str(self.spread),'\n']))
			f.write(''.join(['Min Dist:  ',str(self.min_dist),'\n']))
			f.write(''.join(['Resolution:  ',str(self.resolution),'\n']))

		return 0

	## Remove a specified list of genes from AnnData object
	def filter_specific_genes(self,adata, text_file=None, gene_list=[]):
		'''
		List of genes can be in either a line separated text file or a Python list

		Useful for removing unnecessary genes from the analysis such as blood genes
		'''
		if text_file:
			for line in open(text_file,'r'):
				gene_list.append(line.rstrip('\n'))
		a = [(gene not in gene_list) for gene in adata.var_names]
		
		return adata[:, [(gene not in gene_list) for gene in adata.var_names]]

	## Takes a list of genes and determines if they exist within the data set and are variable
	# Appends results to genes_exist or missing_genes list if given
	def __find_genes(self,adata, gene_list, genes_exist=None, missing_genes=None):
		## Check inputs
		if not genes_exist:
			genes_exist = []

		if not missing_genes:
			missing_genes = []

		## Splits given gene list into two based on whether or not they exist or are invariable
		for gene in gene_list:
			gene_zero = (gene in adata.raw.var_names) and np.any(adata.raw[:,gene].X)
			(genes_exist if gene_zero else missing_genes).append(gene)

		missing_genes = list(set(missing_genes))
		if missing_genes:
			print('Sorry, the following genes are not expressed in this dataset or are invariable:',missing_genes,'\n')

		return [genes_exist, missing_genes]

	## Loads data from storage point and creates an AnnData object
	# Adds metadata to adata object 
	def load_data(self, storage_mount_point=None, sample_list=None):
		'''
		Storage structure at the mount point includes 2 folders - processed data and raw
		The raw data folder contains folders of sample runs as well as the meta-data table
		'''

		# If argument is None, set to instance attribute
		param_dict = {k:(v if v else getattr(self,k)) for (k,v) in locals().items()}

		## Location to output the anndata h5ad files
		raw_data_file = ''.join(['./data/Data_','_'.join(param_dict['sample_list']),'.scanpy.raw.h5ad'])  # the file that will store the raw combined data
		results_file = ''.join(['./data/Data_','_'.join(param_dict['sample_list']),'.processed.h5ad'])  # the file that will store the analysis results

		## Creates a dictionary with sample id key, data file location, and relevant metadata
		annotation_dict = dict()

		# Meta-data table located at in the storage_mount_point - Change if located elsewhere
		for line in open(''.join([param_dict['storage_mount_point'],'01_RNAseq_RAW_Data/single_cell_meta_data_table.tsv']),'r'):
			elem = str.split(line.rstrip())
			if elem[0] not in annotation_dict:
				annotation_dict[elem[0]] = elem[1:]

		self.annotation_dict = annotation_dict

		## Read the raw Cellranger filtered data matrices into new Anndata objects
		if Path(raw_data_file).is_file():
			print(''.join(['Data_','_'.join(param_dict['sample_list']),'.scanpy.raw.h5ad']),'found, using this existing raw data file\n')
			adata = sc.read_h5ad(raw_data_file)
		else:
			print('\nNo existing h5ad raw data file found, reading in 10x h5 data for each sample\n')
			adata = 0
			for sample in param_dict['sample_list']:
				if adata:
					adata = adata.concatenate(self.__create_scanpy_anndata(param_dict['storage_mount_point'], sample, annotation_dict))
				else:
					adata = self.__create_scanpy_anndata(param_dict['storage_mount_point'], sample, annotation_dict)
			
			## Make cell names unique by adding _1, _2, _3 sequentially to each duplicated 10x barcode/name
			adata.obs_names_make_unique()
			
			## Write the raw combined dataset to disk so you can skip combining samples next time
			print('\nSaving raw combined sample data to', raw_data_file,'\n')
			adata.write(raw_data_file)

		return adata

	## Filters data based on certain parameters
	# Attempts to remove "bad" data such as dead cells, doublets, etc.
	def filter_data(self, adata, min_cells=None, min_genes=None, max_counts=None, max_genes=None, max_mito=None):
		'''
		Removes cells expressing low to no genes, and genes expressed in few to no cells
		Filters out cells based on mitochondrial genes, UMI and max gene expression
		'''
		# If argument is None, set to instance attribute
		param_dict = {k:(v if v else getattr(self,k)) for (k,v) in locals().items()}

		self.initial_cell_count = len(adata.obs_names)
		self.initial_gene_count=len(adata.var_names)
		## Basic filtering to get rid of useless cells and unexpressed genes
		sc.pp.filter_genes(adata, min_cells=param_dict['min_cells'])
		sc.pp.filter_cells(adata, min_genes=param_dict['min_genes'])

		# Calculate the percent of genes derived from mito vs genome
		# the `.A1` is only necessary as X is sparse (to transform to a dense array after summing)
		mito_genes = adata.var_names.str.startswith('MT-')
		adata.obs['percent_mito'] = np.sum(adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1

		# add the total counts per cell as observations-annotation to adata
		adata.obs['n_counts'] = adata.X.sum(axis=1).A1

		self.adata_preFiltered = adata.copy() # Saving pre-Filtered AnnData
		## Actually do the filtering.
		adata = adata[((adata.obs['n_genes'] < param_dict['max_genes'])   # Keep cells with less than __ genes to remove most doublets
					& (adata.obs['n_counts'] < param_dict['max_counts'])   # Keep cells with less than __ UMIs to catch a few remaining doublets
					& (adata.obs['percent_mito'] < param_dict['max_mito'])), :]   # Keep cells with less than __ mito/genomic gene ratio

		return adata

	## Standardize and normalize the data set
	# Includes additional processing such as removing biased and uninformative data
	def preprocess_data(self,adata):
		## Normalize the expression matrix to 10,000 reads per cell, so that counts become comparable among cells.
		# This corrects for differences in sequencing depth between cells and samples
		sc.pp.normalize_total(adata)

		## Log transform the data.
		sc.pp.log1p(adata)

		## Set the .raw attribute of AnnData object to the logarithmized raw gene expression for later use in differential testing and visualizations of gene expression.
		# We need to do this because the expression matrix will be rescaled and centered which flattens expression too much for some purposes
		adata.raw = adata

		## Identify highly-variable genes based on dispersion relative to expression level.
		sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

		## Filter the genes to remove non-variable genes since they are uninformative
		adata = adata[:, adata.var['highly_variable']].copy()

		## Regress out effects of total reads per cell and the percentage of mitochondrial genes expressed.
		#sc.pp.filter_genes(adata, min_counts=1) # 0 columns negatively affect the convergence of regresion
		sc.pp.regress_out(adata, ['n_counts','percent_mito'])

		print('\nDoing, final filtering...\nKeeping', len(adata.obs_names),'cells and', len(adata.var_names),'genes.\n')

		#sc.pp.filter_cells(adata, min_genes=1) # Remove 0 columns that may have occured in sectioning of data after preprocessing

		## Scale each gene to unit variance. Clip values exceeding standard deviation 10 to remove extreme outliers
		sc.pp.scale(adata, max_value=10)
		return adata

	## Run dimensional reduction analysis and clustering using KNN graph
	def run_analysis(self,adata, n_neighbors=None, n_pcs=None, spread=None,
					 min_dist=None, resolution=None, remove_batch_effects=None, do_tSNE=None):
		# If argument is None, set to instance attribute
		param_dict = {k:(v if v else getattr(self,k)) for (k,v) in locals().items()}

		## Run PCA to compute the default number of components
		sc.tl.pca(adata, svd_solver='arpack')

		## Save the existing data to disk for later
		self.adata_postPCA = adata.copy()

		## Testing of batch correction technique
		#sc.pp.combat(adata, key='sampleName', covariates={'HT239-ESO-2D':1,'HT239_ESO':0})

		## Compute nearest-neighbors
		sc.pp.neighbors(adata, n_neighbors=param_dict['n_neighbors'], n_pcs=param_dict['n_pcs'])

		## Remove batch effects
		# Note that doing this may override previous sc.pp.neighbors()
		if param_dict['remove_batch_effects']:
			import bbknn
			bbknn.bbknn(adata, batch_key='sampleName', copy=False)
			#sc.pp.external.mnn_correct(adata,batch_key='sampleName') # Testing another algorithm

		## Run UMAP Dim reduction
		sc.tl.umap(adata, spread=param_dict['spread'], min_dist=param_dict['min_dist']) # Min_dist needs to be between 0.01 to 0.5

		## Run tSNE analysis
		if param_dict['do_tSNE']:
			sc.tl.tsne(adata, n_pcs=param_dict['n_pcs'])

		## Calculate cell clusters via Louvain algorithm
		sc.tl.louvain(adata, resolution=param_dict['resolution'])

		## Run PAGA to predict non-directional cluster-cluster relationships to infer possible developmental progressions
		sc.tl.paga(adata, groups='louvain', model='v1.2')

		## Do Dendrogram analysis based on PCs
		sc.tl.dendrogram(adata,groupby='louvain',n_pcs=param_dict['n_pcs'],linkage_method="median",use_raw=True)
		return adata

	## Variety of plotting and data display functions
	# Will make more robust in the future
	def plot_sca(self, adata, figdir='./figures/', adata_preFiltered=None,
				annotation_dict=None, adata_postPCA=None, final_quality=None):
		'''
		See the Scanpy visualization library for examples
		'''
		print("Plotting")
		self.final_cell_count = len(adata.obs_names)
		self.final_gene_count=len(adata.var_names)
		param_dict = {k:(v if v else getattr(self,k)) for (k,v) in locals().items()}

		## Create my custom palette for FeaturePlots and define a matlplotlib colormap object
		if self.umap_color=='blue_orange':
			feature_colors = [(35,35,142), (255,127,0)]
			my_feature_cmap = self.__make_cmap(feature_colors,bit=True)
		elif self.umap_color=='yellow_blue':
			feature_colors = [(210,210,210), (210,210,210), (245,245,200), (100,200,225), (0,45,125)]
			position=[0, 0.019999, 0.02, 0.55, 1]
			my_feature_cmap = self.__make_cmap(feature_colors,bit=True,position=position)
		else:
			feature_colors = [(210,210,210), (210,210,210), (245,245,200), (100,200,225), (0,45,125)]
			position=[0, 0.019999, 0.02, 0.55, 1]
			my_feature_cmap = self.__make_cmap(feature_colors,bit=True,position=position)

		gray_cmap = self.__make_cmap([(220,220,220),(220,220,220)], bit=True)

		## Custom color palette for cluster plots and observation plots
		colors = [(0.3,0.3,0.3),(1,0,0),(0,1,0),(0,0,0.9),(1,0,1),(0,1,1),(0.9,0.9,0),
				(0.85,0.5,0.5),(0.5,0.85,0.5),(0.5,0.5,0.85),
				(0.15,0.5,0.5),(0.5,0.15,0.5),(0.5,0.5,0.15),
				(1,0.5,0),(0,0.5,1),(0.5,1,0),(0.5,0,1),(1,0,0.5),(0,1,0.5)]

		## General figure parameters and settings
		sc.set_figure_params(dpi_save=300,dpi=300)#,vector_friendly=False)
		sc.settings.figdir = figdir
		sc.set_figure_params(fontsize=12)
		size = self.size

		# Check to see if user wants publication quality figures
		if param_dict['final_quality']:
			rcParams['figure.figsize'] = 4, 4
			rcParams['savefig.dpi']=1200
			file_type = '.pdf'
		else:
			file_type = '.png'

		## Violin plots for filtering parameters pre and post
		sc.pl.violin(self.adata_preFiltered, ['n_genes','n_counts','percent_mito'],
		 			 jitter=0.4, multi_panel=True, save='_preFiltered_plot.png', show=False)
		sc.pl.violin(adata, ['n_genes','n_counts','percent_mito'],
					 jitter=0.4, multi_panel=True, save='_postFiltered_plot.png', show=False)

		## Draw the PCA elbow plot to determine which PCs to use
		sc.pl.pca_variance_ratio(self.adata_postPCA, log=True, n_pcs=100, save='_elbowPlot.png', show=False)
		## Ranks and displays most contributing genes for each principal component
		components = 4
		loadings_components = range(self.n_pcs-components, self.n_pcs+components+1)
		sc.pl.pca_loadings(self.adata_postPCA, components=loadings_components, save='_rank_genes.png', show=False)

		## Plot results of UMAP dimensional reduction and clustering
		for observation in self.umap_obs:
			legend = 'on data' if (observation=='louvain') else 'right margin'
			sc.pl.umap(adata, color=observation, save=''.join(['_',observation,file_type]), show=False,
					   legend_loc=legend, edges=False, size=size, palette=colors, alpha=0.75)

		## Find marker genes via Wilxocon test based on Louvain cluster assignment
		# Create a simple plot to show the top 25 most significant markers for each cluster
		# Write most significant markers to a csv file
		# adata.obs['is_adult'] = ['Adult' if cell=='ND15989_Fresh_WT_Lung_Adult' else 'Fetal' for cell in adata.obs['sampleName']]
		# rank_grouping = 'age'
		# rank_genes(adata,rank_grouping,figdir=figdir)#,clusters2_compare=['1','4'])
		# sc.pl.rank_genes_groups_heatmap(adata, n_genes=100, use_raw=True, show=False, 
		# 		save=''.join(['_rank_heatmap_',rank_grouping,file_type]), cmap=my_feature_cmap)
		# sc.pl.rank_genes_groups_dotplot(adata, n_genes=5, use_raw=True, show=False, 
		# 		save=''.join(['_rank_dotplot_',rank_grouping,file_type]), color_map=my_feature_cmap)
		# sc.pl.rank_genes_groups_stacked_violin(adata, n_genes=5, use_raw=True, 
		# 		show=False, save=''.join(['_rank_violin_',rank_grouping,file_type]))

		rank_grouping = 'louvain'
		n_genes_rank = 25
		self.__rank_genes(adata,rank_grouping,figdir=figdir)#,clusters2_compare=['1','4'])
		sc.pl.rank_genes_groups_heatmap(adata, n_genes=n_genes_rank, use_raw=True, show=False, 
				save=''.join(['_rank_heatmap_',rank_grouping,file_type]), cmap=my_feature_cmap)
		sc.pl.rank_genes_groups_dotplot(adata, n_genes=n_genes_rank, use_raw=True, show=False, 
				save=''.join(['_rank_dotplot_',rank_grouping,file_type]), color_map=my_feature_cmap)
		sc.pl.rank_genes_groups_stacked_violin(adata, n_genes=n_genes_rank, use_raw=True, 
				show=False, save=''.join(['_rank_violin_',rank_grouping,file_type]))

		## Feature plots and dot plot analysis for each specified set of genes
		#sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, save='_markerPlots.png', show=False)
		if self.gene_lists:
			missing_genes = []
			for gene_list in self.gene_lists:
				gene_dict = self.gene_dict[gene_list]
				genes_to_plot = []
				[genes_to_plot,missing_genes] = self.__find_genes(adata,gene_dict['markers'], missing_genes=missing_genes)

				## Do FeaturePlots for select genes
				print('Plotting standard marker genes: ',genes_to_plot,'\n')
				sc.pl.umap(adata, color=genes_to_plot, save= ''.join(['_featureplots_',gene_list,file_type]), show=False, 
						   cmap=my_feature_cmap, size=size, use_raw=True)
		
				if gene_dict['positions'] and gene_dict['groups']:
					group_positions = gene_dict['positions'] # Manually set and determined
					group_labels = gene_dict['groups']
				else:
					group_positions = None
					group_labels = None

				if len(gene_dict['markers'])!=1:
					for grouping in self.exp_grouping:
						## Dotplot analysis
						# Circle color corresponds to expression level, and circle size corresponds to percentage of cells expressing gene

						## Reordering categories for dotplot or heatmap rows
						#adata_temp = adata.copy()
						#adata_temp.obs['louvain'] = adata.obs['louvain'].cat.reorder_categories(['3','5','0','4','2','1'],inplace = False)

						sc.pl.dotplot(adata, genes_to_plot, groupby=grouping, 
								var_group_positions=group_positions, var_group_labels=group_labels,
								save=''.join(['_markers_',gene_list,'_',grouping,file_type]), show=False, 
								color_map=my_feature_cmap, use_raw=True, dot_max=0.4)
						## Heatmaps
						# Each horizontal line represents expression of one cell
						sc.pl.heatmap(adata, genes_to_plot, groupby=grouping, 
								var_group_positions=group_positions, var_group_labels=group_labels,
								save=''.join(['_markers_',gene_list,'_',grouping,file_type]), show=False, 
								cmap=my_feature_cmap, use_raw=True)

		# Genes that are not expressed or are invariable are plotted using a grayscale
		print('Plotting empty genes: ',missing_genes,'\n')
		sc.pl.umap(adata, color=missing_genes, save=''.join(['_featureplots_gray',file_type]), 
				show=False, cmap=gray_cmap, size=size, use_raw=True)

		if self.do_tSNE:
			# tSNE Plots
			sc.pl.tsne(adata, color='louvain', save = '_clusterIdentity.png', show = False, 
						legend_loc = 'right margin', edges = False, size = 20, 
						palette = colors, alpha = 0.75)
			sc.pl.tsne(adata, color='sampleName', save = '_sample.png', show = False, 
						legend_loc = 'right margin', edges = False, size = 20, 
						palette = colors, alpha = 0.75)
			sc.pl.tsne(adata, color=genes_to_plot, save = '_featureplots.png', show = False, cmap = my_feature_cmap, size=size, use_raw = True)
			sc.pl.tsne(adata, color=missing_genes, save='_featureplots_gray.png', show=False, cmap=gray_cmap, size=size, use_raw=True)
		

		# ## Violin plot for comparing gene expression among different groups/clusters
		# # Create observation field labeled using binary information
		# # Will have option integrated in pipeline in the future
		# adata.obs['CDH5_exp'] = ['CDH5+' if (cell!=0) else 'CDH5-' for cell in adata.raw[:,'CDH5'].X] 

		# # Built in scanpy module
		# sc.pl.violin(adata, genes_to_plot+['CDH5'], groupby='CDH5_exp', jitter=True,
		# 	save='_feature.png', show=False, scale='width',use_raw=True) #order = ['CDH5+','CDH5-'],
		
		# Custom violin plot module -- Not complete/in testing
		df = pd.DataFrame()
		# Add Gaussian y-jitter to better visualize zero expression in violin plots
		for gene in genes_to_plot:
			sigma = np.amax(adata.raw[:,gene].X)/40
			gene_df = [cell if (cell!=0) else np.random.normal(loc=0,scale=sigma) for cell in adata.raw[:,gene].X]
			df[gene] = gene_df

		# df['CDH5_exp']=adata.obs['CDH5_exp'].values
		# vplot, axes = plt.subplots(math.ceil(len(genes_to_plot)/4),4, figsize=(18,12))
		# plt.rcParams.update({'font.size':12})
		# plt.subplots_adjust(left=0.125, right=0.9, bottom=0.1, top=0.9, wspace=0.4, hspace=0.4)
		# for i,gene in enumerate(genes_to_plot):
		# 	sns.violinplot(x='CDH5_exp', y=gene, data=df, inner=None, scale='width', ax=axes[math.floor(i/4),i%4])
		# 	sns.stripplot(x='CDH5_exp', y=gene, data=df, jitter = True, color='black', size=0.4, ax=axes[math.floor(i/4),i%4])
		# vplot.savefig(''.join([figdir,'/violin_feature_jitter.png']))

		## Scatter plots to identify clusters that are high in 
		adata.obs['jitter'] = np.random.rand(len(adata.obs_names))*10
		sc.pl.scatter(adata,x='jitter',y='n_genes',color='louvain',save='_n_genes.png',palette=colors,show=False)
		sc.pl.scatter(adata,x='jitter',y='n_counts',color='louvain',save='_n_counts.png',palette=colors,show=False)
		sc.pl.scatter(adata,x='jitter',y='percent_mito',color='louvain',save='_percent_mito.png',palette=colors,show=False)

		# Set the thresholds and scaling factors for drawing the paga map/plot
		node_size_scale=1.25
		node_size_power=0.9
		edge_width_scale=1
		min_edge_width=0.035
		max_edge_width=2
		threshold=0.08
		# Draw the actual plot 
		sc.pl.paga(adata, layout='fr', threshold=threshold, node_size_scale=node_size_scale, 
			node_size_power=node_size_power, edge_width_scale=edge_width_scale,
			min_edge_width=min_edge_width, max_edge_width=max_edge_width, show=False, save = '_pagaPlot.png',
			title='PAGA: Fruchterman Reingold',
			frameon=False)
		
		return adata

	## Most basic pipeline - Input data and output all figures
	def pipe_basic(self,figdir='./figures/', load_save=None, new_save='adata_save.p', remove_genes=None, final_quality=False):
		'''
		sca_dict is a dictionary of miscellaneous analysis information including
		parameters, sample list and gene_lists

		Uses the pickle module to save adata instances for easy access in future analyses
		*Note that pickled adata files are very large - min 4GB (make sure to clean out unused pickle files)
		'''
		if load_save: # See if there is already a save file for the analysis to duplicate
			run_save = pickle.load(open(''.join([figdir,load_save]),"rb"))

			adata = run_save.adata
			self.adata_preProcessed = run_save.adata_preProcessed
			self.adata_preFiltered = run_save.adata_preFiltered
			self.annotation_dict = run_save.annotation_dict
		else:
			adata = self.load_data()

			## Remove specific genes from the analysis (such as experimentally observed contaminations)
			if remove_genes:
				adata = self.filter_specific_genes(adata,text_file = remove_genes)

			## Filter and process data
			adata = self.filter_data(adata)
			self.adata_preProcessed = adata.copy() # Save for later in case of necessary extraction
			adata = self.preprocess_data(adata)

		## Dimensional reduction and clustering - construction of the neighborhood graph
		adata = self.run_analysis(adata)

		## Plot figures
		adata = self.plot_sca(adata,figdir=figdir)

		self.adata = adata.copy()
		## Save analysis information and relevant AnnData objects to the disk using the Pickle module
		if new_save:
			pickle.dump(self,open(''.join([figdir,new_save]),"wb"),protocol=4)

		## Write a summary of the analysis to a text file including sample information and parameters
		self.write_summary()

		print("\nAll done!\n")

		return sca_run

	## Pipeline for analysis in which you extract interesting clusters/observations after an initial run
	# Extracts clusters to an filtered but unprocessed AnnData object, then reprocesses and reclusters
	def pipe_ext(self, analysis_params_ext, figdir='./figures/', extracted=None, load_save=None, new_save='adata_save.p', final_quality=False, label=''):
		'''
		Allows loading of a saved pickle adata file, file must contained adata that has gone through a complete pipeline
		Otherwise will complete a new full analysis, and then extracts clusters
		'''
		## Ask user for input on which clusters to extract
		if extracted is None:
			print("Which clusters would you like to extract? (Separate your answer with commas) ")
			extracted = list(map(str,input().split(',')))

		if not np.any([not cluster.isdigit() for cluster in extracted]):
			## Load pickle file if one is given
			if load_save:
				run_save = pickle.load(open(''.join([figdir,load_save]),"rb"))
			else: # Otherwise complete a new analysis
				run_save = self.pipe_basic(sca_run,figdir)
		
			adata = run_save.adata
			self.adata_preProcessed = run_save.adata_preProcessed
			self.adata_preFiltered = run_save.adata_preFiltered
			self.annotation_dict = run_save.annotation_dict

			## Create an unprocessed AnnData object with the desired clusters
			adata_ext = self.adata_preProcessed[adata.obs['louvain'].isin(extracted)].copy()
			self.adata_preProcessed = adata_ext.copy()
			## Reprocess and recluster extracted cells
			self.preprocess_data(adata_ext)
		else:
			print("Looking to extract based on gene expression of ",extracted)
			self.load_data()
			self.filter_data(adata)
			self.adata_preProcessed = adata.copy()
			self.preprocess_data(adata)
			adata_ext = adata[adata.raw[:,extracted].X>1.5,:]

		# sca_dict_ext = copy.deepcopy(sca_dict)
		# sca_dict_ext.update(analysis_params = analysis_params_ext) # Use new analysis parameters
		self.set_analysis_params(**(analysis_params_ext))
		self.run_analysis(adata_ext)

		## Plot figures
		adata_ext = self.plot_sca(adata_ext,figdir = ''.join([figdir,'extracted/',label]))

		self.adata = adata_ext.copy()
		## Save analysis information and relevant AnnData objects to the disk using the Pickle module
		if new_save:
			pickle.dump(self,open(''.join([figdir,'extracted/',new_save]),"wb"),protocol=4)

		## Write a summary of the analysis to a text file including sample information and parameters
		self.write_summary()

		print("\nAll done!\n")

		return sca_run