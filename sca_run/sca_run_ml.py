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
from typing import Union, Optional, Tuple, Collection
from anndata import AnnData
from scipy.sparse import issparse, isspmatrix_csr, csr_matrix, spmatrix
from sklearn.utils import sparsefuncs
from sklearn.decomposition import IncrementalPCA, PCA, TruncatedSVD
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier
from sklearn.calibration import CalibratedClassifierCV

class sca_run_ml:
	sc.settings.verbosity = 3
	N_PCS = 50  # default number of PCs

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
		self.cell_score_lists = []
		self.transform_PCA = None
		self.combat = False

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
		self.adata_unscaled = None
		self.adata_preProcessed = None

	## Function to set a list of genes and save relevant information in a dictionary of gene lists
	def add_gene_list(self, markers=['EPCAM','CDH5','VIM','TP63'], label='basic_list', feature_positions=None, feature_groups=None, groupby_positions=None):
		self.gene_lists = self.gene_lists + [label]
		self.gene_dict[label]={'markers':markers,
							  'feature_positions':feature_positions,
							  'feature_groups':feature_groups,
							  'groupby_positions':groupby_positions}

	## Function to set multiple filtering parameters in attributes
	def set_filter_params(self, min_cells=0, min_genes=500, max_genes=7000, max_counts=30000, max_mito=0.1):
		self.__set_param_attr(locals())

	## Function to set multiple analysis parameters in attributes
	def set_analysis_params(self, n_neighbors=30, n_pcs=15, spread=1, min_dist=0.4, resolution=0.6, remove_batch_effects=False, do_tSNE=False, combat=False, cell_score_lists=[]):
		self.__set_param_attr(locals())

	## Function to set multiple plot parameters in attributes
	def set_plot_params(self, size=20, umap_obs=['louvain','sampleName'], exp_grouping=['louvain'], umap_color='yellow_blue', final_quality=False):
		self.__set_param_attr(locals())

	## Function that checks if an argument exists in a function and sets instance attributes based on them
	def __set_param_attr(self,arg_dict):
		for key in arg_dict:
			if key != 'self':
				setattr(self,key,arg_dict[key])

	## Define function to generate a color gradient from a defined starting and ending color
	def make_cmap(self,colors, position=None, bit=False):
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

	## Count number of cells in each louvain cluster as well as the sample splits within each cluster
	def __cell_counter(self,adata):
		#print(adata.obs['louvain'])
		self.sample_counts = {}
		self.louvain_counts = {}
		for sample in adata.obs['sampleName'].cat.categories:
			try:
				self.sample_counts[sample] = {'sample_total':adata.obs['sampleName'].isin([sample]).value_counts()[True]}
			except:
				self.sample_counts[sample] = 0
			for cluster in adata.obs['louvain'].cat.categories:
				if not cluster in self.louvain_counts:
					try:
						self.louvain_counts[cluster] = {'louvain_total':adata.obs['louvain'].isin([cluster]).value_counts()[True]}
					except:
						self.louvain_counts[cluster] = {'louvain_total':0}

				try:
					self.sample_counts[sample][cluster] = self.louvain_counts[cluster][sample] = (adata.obs['sampleName'].isin([sample]) & adata.obs['louvain'].isin([cluster])).value_counts()[True]
				except:
					self.sample_counts[sample][cluster] = 0 
					self.louvain_counts[cluster][sample] = 0

		return 0

	## If doing cell scoring analysis, get average score per cluster given a gene scoring list name
	def __count_cluster_scores(self, adata, cell_score_list):
		cluster_scores = {}
		for cluster in adata.obs['louvain'].cat.categories:
			cluster_scores[cluster] = np.mean(adata[adata.obs['louvain'].isin([cluster])].obs[cell_score_list].values)

		return cluster_scores

	## Write a summary of the analysis run including sample information, parameters and filtering information
	# Not completely up to date
	def write_summary(self, figdir='./figures/'):
		self.__cell_counter(self.adata)

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

			f.write('\n--------Sample Cell Counts Used--------\n')
			f.write(json.dumps(str(self.sample_counts)))
			f.write('\n--------Sample Cell Counts Used--------\n')
			f.write(json.dumps(str(self.louvain_counts)))

			if self.cell_score_lists:
				f.write('\n--------Cluster Scores--------\n')
				for cell_score_list in self.cell_score_lists:
					# f.write(''.join(['--------',cell_score_list,'_raw--------\n']))
					# f.write(json.dumps(str(self.__count_cluster_scores(self.adata, ''.join([cell_score_list,'_raw'])))))
					f.write(''.join(['--------',cell_score_list,'--------\n']))
					f.write(json.dumps(str(self.__count_cluster_scores(self.adata,cell_score_list))))
					f.write('\n')
					# f.write(''.join(['\n--------',cell_score_list,'_processed--------\n']))
					# f.write(json.dumps(str(self.__count_cluster_scores(self.adata, ''.join([cell_score_list,'_processed'])))))

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
	def __find_genes(self,adata, gene_list, genes_exist=[], missing_genes=[]):
		## Check inputs
		if not genes_exist:
			genes_exist = []

		if not missing_genes:
			missing_genes = []

		## Splits given gene list into two based on whether or not they exist or are invariable
		for gene in gene_list:
			gene_zero = (gene in adata.raw.var_names) and np.any(adata.raw[:,gene].X.toarray())
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
			try:
				print('\nSaving raw combined sample data to', raw_data_file,'\n')
				adata.write(raw_data_file)
			except:
				print('\nUnable to save raw combined sample data to', raw_data_file,'\n')

		return adata

	## Filters data based on certain parameters
	# Attempts to remove "bad" data such as dead cells, doublets, etc.
	def filter_data(self, adata, min_cells=None, min_genes=None, max_counts=None, max_genes=None, max_mito=None,
					remove_batch_effects=None):
		'''
		Removes cells expressing low to no genes, and genes expressed in few to no cells
		Filters out cells based on mitochondrial genes, UMI and max gene expression
		'''
		# If argument is None, set to instance attribute
		param_dict = {k:(v if v else getattr(self,k)) for (k,v) in locals().items()}
		self.initial_cell_count = len(adata.obs_names)
		self.initial_gene_count = len(adata.var_names)

		## Basic filtering to get rid of useless cells and unexpressed genes
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

	def __separate_adata(self,adata,obs_field):
		adata_list = []
		for meta_value in adata.obs[obs_field].unique():
			adata_list.append(adata[(adata.obs[obs_field] == meta_value),:].copy())

		return adata_list

	def _scale(self,X, zero_center=True,transform=None):
		# - using sklearn.StandardScaler throws an error related to
		#   int to long trafo for very large matrices
		# - using X.multiply is slower
		#   the result differs very slightly, why?

		if transform:
			mean = transform[0]
			scale = transform[1]
		else:
			mean, var = sc.pp._utils._get_mean_var(X)
			scale = np.sqrt(var)

		if issparse(X):
			print("Is sparse in _scale function")
			if zero_center: raise ValueError('Cannot zero-center sparse matrix.')
			sparsefuncs.inplace_column_scale(X, 1/scale)
		else:
			#print(mean)
			#print(scale)
			X -= mean
			scale[scale == 0] = 1e-12
			X /= scale
		return [mean, scale]

	def scale(self,data, zero_center=True, max_value=None, copy=False, transform=None) -> Optional[AnnData]:
	    """Scale data to unit variance and zero mean.
	    .. note::
	        Variables (genes) that do not display any variation (are constant across
	        all observations) are retained and set to 0 during this operation. In
	        the future, they might be set to NaNs.
	    Parameters
	    ----------
	    data : :class:`~anndata.AnnData`, `np.ndarray`, `sp.sparse`
	        The (annotated) data matrix of shape `n_obs` × `n_vars`. Rows correspond
	        to cells and columns to genes.
	    zero_center : `bool`, optional (default: `True`)
	        If `False`, omit zero-centering variables, which allows to handle sparse
	        input efficiently.
	    max_value : `float` or `None`, optional (default: `None`)
	        Clip (truncate) to this value after scaling. If `None`, do not clip.
	    copy : `bool`, optional (default: `False`)
	        If an :class:`~anndata.AnnData` is passed, determines whether a copy
	        is returned.
	    Returns
	    -------
	    Depending on `copy` returns or updates `adata` with a scaled `adata.X`.
	    """
	    if isinstance(data, AnnData):
	        adata = data.copy() if copy else data
	        # sc.utils.view_to_actual(adata)
	        # need to add the following here to make inplace logic work
	        if zero_center and issparse(adata.X):
	            adata.X = adata.X.toarray()
	        result = self.scale(adata.X, zero_center=zero_center, max_value=max_value, copy=False, transform=transform)
	        mean = result[0]
	        scale = result[1]
	        return [mean, scale, adata] if copy else [mean, scale]
	    X = data.copy() if copy else data  # proceed with the data matrix
	    zero_center = zero_center if zero_center is not None else False if issparse(X) else True
	    if zero_center and issparse(X):
	        X = X.toarray()
	        copy = True
	    mean, scale = self._scale(X, zero_center, transform=transform)
	    if max_value is not None: X[X > max_value] = max_value
	    return [X, mean, scale] if copy else [mean, scale]

	## Standardize and normalize the data set
	# Includes additional processing such as removing biased and uninformative data
	def preprocess_data(self,adata=None, highly_variable=None,transform_X=None,combat=False,adata_unscaled=None,regress_out=['n_counts','percent_mito']):
		if not adata_unscaled:
			## Normalize the expression matrix to 10,000 reads per cell, so that counts become comparable among cells.
			# This corrects for differences in sequencing depth between cells and samples
			sc.pp.normalize_total(adata,target_sum=10000)

			## Log transform the data.
			sc.pp.log1p(adata)

			## Set the .raw attribute of AnnData object to the logarithmized raw gene expression for later use in differential testing and visualizations of gene expression.
			# We need to do this because the expression matrix will be rescaled and centered which flattens expression too much for some purposes
			adata.raw = adata

			## Find cell type score for each cell based on a predefined set of gene lists
			if self.cell_score_lists:
				#adata_score_temp = adata.copy()
				#sc.pp.regress_out(adata_score_temp,['n_counts','percent_mito'])
				adata_scaled = sc.pp.scale(adata, max_value=10, copy=True)
				for file in self.cell_score_lists:
					score_list = [line.rstrip('\n') for line in open(''.join([file,'.txt']),'r')]
					# print(score_list)
					adata.obs[file] = adata_scaled.X[:,adata_scaled.var_names.isin(score_list)].mean(1)

					# sc.tl.score_genes(adata_score_temp, score_list, ctrl_size=50, gene_pool=None, n_bins=25, score_name=file, random_state=0, copy=False, use_raw=False)
					# adata.obs[file]=copy.deepcopy(adata_score_temp.obs[file])

			## Identify highly-variable genes based on dispersion relative to expression level.
			if not np.any(highly_variable):
				sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
				highly_variable = copy.deepcopy(adata.var['highly_variable'])

			## Filter the genes to remove non-variable genes since they are uninformative
			adata = adata[:, highly_variable].copy()

			## Testing of batch correction technique
			if combat:
				print("Conducting combat batch correction")
				sc.pp.combat(adata, key='sampleName')
		else:
			adata = adata_unscaled.copy()
		self.adata_unscaled = adata.copy()

		sc.pp.filter_genes(adata, min_cells=self.min_cells) ## Remove genes expressed in a low number of cells

		## Regress out effects of total reads per cell and the percentage of mitochondrial genes expressed.
		#sc.pp.filter_genes(adata, min_counts=1) # 0 columns negatively affect the convergence of regresion
		if regress_out:
			sc.pp.regress_out(adata, regress_out)

		print('\nDoing, final filtering...\nKeeping', len(adata.obs_names),'cells and', len(adata.var_names),'genes.\n')

		#sc.pp.filter_cells(adata, min_genes=1) # Remove 0 columns that may have occured in sectioning of data after preprocessing

		## Scale each gene to unit variance. Clip values exceeding standard deviation 10 to remove extreme outliers
		transform_X = self.scale(adata, max_value=10, transform=transform_X)

		return [adata,highly_variable,transform_X]

	def pca(self,
	    data: Union[AnnData, np.ndarray, spmatrix],
	    n_comps: int = N_PCS,
	    zero_center: Optional[bool] = True,
	    svd_solver: str = 'auto',
	    random_state: int = 0,
	    return_info: bool = False,
	    use_highly_variable: Optional[bool] = None,
	    dtype: str = 'float32',
	    copy: bool = False,
	    chunked: bool = False,
	    chunk_size: Optional[int] = None,
	    transform_PCA: Union[IncrementalPCA,PCA,TruncatedSVD,None] = None,
	) -> Union[AnnData, np.ndarray, spmatrix]:
	    """Principal component analysis [Pedregosa11]_.
	    Computes PCA coordinates, loadings and variance decomposition. Uses the
	    implementation of *scikit-learn* [Pedregosa11]_.
	    Parameters
	    ----------
	    data
	        The (annotated) data matrix of shape ``n_obs`` × ``n_vars``.
	        Rows correspond to cells and columns to genes.
	    n_comps
	        Number of principal components to compute.
	    zero_center
	        If `True`, compute standard PCA from covariance matrix.
	        If ``False``, omit zero-centering variables
	        (uses :class:`~sklearn.decomposition.TruncatedSVD`),
	        which allows to handle sparse input efficiently.
	        Passing ``None`` decides automatically based on sparseness of the data.
	    svd_solver
	        SVD solver to use:
	        ``'arpack'``
	          for the ARPACK wrapper in SciPy (:func:`~scipy.sparse.linalg.svds`)
	        ``'randomized'``
	          for the randomized algorithm due to Halko (2009).
	        ``'auto'`` (the default)
	          chooses automatically depending on the size of the problem.
	    random_state
	        Change to use different initial states for the optimization.
	    return_info
	        Only relevant when not passing an :class:`~anndata.AnnData`:
	        see “**Returns**”.
	    use_highly_variable
	        Whether to use highly variable genes only, stored in
	        ``.var['highly_variable']``.
	        By default uses them if they have been determined beforehand.
	    dtype
	        Numpy data type string to which to convert the result.
	    copy
	        If an :class:`~anndata.AnnData` is passed, determines whether a copy
	        is returned. Is ignored otherwise.
	    chunked
	        If ``True``, perform an incremental PCA on segments of ``chunk_size``.
	        The incremental PCA automatically zero centers and ignores settings of
	        ``random_seed`` and ``svd_solver``. If ``False``, perform a full PCA.
	    chunk_size
	        Number of observations to include in each chunk.
	        Required if ``chunked=True`` was passed.
	    Returns
	    -------
	    X_pca : :class:`scipy.sparse.spmatrix` or :class:`numpy.ndarray`
	        If `data` is array-like and ``return_info=False`` was passed,
	        this function only returns `X_pca`…
	    adata : anndata.AnnData
	        …otherwise if ``copy=True`` it returns or else adds fields to ``adata``:
	        ``.obsm['X_pca']``
	             PCA representation of data.
	        ``.varm['PCs']``
	             The principal components containing the loadings.
	        ``.uns['pca']['variance_ratio']``)
	             Ratio of explained variance.
	        ``.uns['pca']['variance']``
	             Explained variance, equivalent to the eigenvalues of the covariance matrix.
	    """
	    # chunked calculation is not randomized, anyways
	    # if svd_solver in {'auto', 'randomized'} and not chunked:
	    #     logg.info(
	    #         'Note that scikit-learn\'s randomized PCA might not be exactly '
	    #         'reproducible across different computational platforms. For exact '
	    #         'reproducibility, choose `svd_solver=\'arpack\'.` This will likely '
	    #         'become the Scanpy default in the future.'
	    #     )

	    data_is_AnnData = isinstance(data, AnnData)
	    if data_is_AnnData:
	        adata = data.copy() if copy else data
	    else:
	        adata = AnnData(data)

	    # start = logg.info(f'computing PCA with n_comps = {n_comps}')

	    if adata.n_vars < n_comps:
	        n_comps = adata.n_vars - 1
	        # logg.debug(
	        #     f'reducing number of computed PCs to {n_comps} '
	        #     f'as dim of data is only {adata.n_vars}'
	        # )

	    if use_highly_variable is True and 'highly_variable' not in adata.var.keys():
	        raise ValueError('Did not find adata.var[\'highly_variable\']. '
	                         'Either your data already only consists of highly-variable genes '
	                         'or consider running `pp.filter_genes_dispersion` first.')
	    if use_highly_variable is None:
	        use_highly_variable = True if 'highly_variable' in adata.var.keys() else False
	    # if use_highly_variable:
	    #     logg.info('computing PCA on highly variable genes')
	    adata_comp = adata[:, adata.var['highly_variable']] if use_highly_variable else adata

	    ### Not supported yet
	    if chunked:
	        # if not zero_center or random_state or svd_solver != 'auto':
	        #     logg.debug('Ignoring zero_center, random_state, svd_solver')

	        from sklearn.decomposition import IncrementalPCA

	        X_pca = np.zeros((adata_comp.X.shape[0], n_comps), adata_comp.X.dtype)
	        if transform_PCA is None:
	        	pca_ = IncrementalPCA(n_components=n_comps)
	        else:
	        	pca_ = transform_PCA

	        for chunk, _, _ in adata_comp.chunked_X(chunk_size):
	            chunk = chunk.toarray() if issparse(chunk) else chunk
	            pca_.partial_fit(chunk)

	        for chunk, start, end in adata_comp.chunked_X(chunk_size):
	            chunk = chunk.toarray() if issparse(chunk) else chunk
	            X_pca[start:end] = pca_.transform(chunk)
	    else:
	        if zero_center is None:
	            zero_center = not issparse(adata_comp.X)
	        if zero_center:
	            from sklearn.decomposition import PCA
	            if issparse(adata_comp.X):
	                # logg.debug(
	                #     '    as `zero_center=True`, '
	                #     'sparse input is densified and may '
	                #     'lead to huge memory consumption',
	                # )
	                X = adata_comp.X.toarray()  # Copying the whole adata_comp.X here, could cause memory problems
	            else:
	                X = adata_comp.X

	            if transform_PCA is None:
	            	pca_ = PCA(n_components=n_comps, svd_solver=svd_solver, random_state=random_state)
	            else: 
	            	pca_ = transform_PCA
	        else:
	            from sklearn.decomposition import TruncatedSVD
	            # logg.debug(
	            #     '    without zero-centering: \n'
	            #     '    the explained variance does not correspond to the exact statistical defintion\n'
	            #     '    the first component, e.g., might be heavily influenced by different means\n'
	            #     '    the following components often resemble the exact PCA very closely'
	            # )
	            if transform_PCA is None:
	            	pca_ = TruncatedSVD(n_components=n_comps, random_state=random_state)
	            else:
	            	pca_ = transform_PCA
	            X = adata_comp.X

	        if transform_PCA is None:
	        	X_pca = pca_.fit_transform(X)
	        else:
	        	X_pca = pca_.transform(X)

	    if X_pca.dtype.descr != np.dtype(dtype).descr: X_pca = X_pca.astype(dtype)

	    if data_is_AnnData:
	    	print("does this run?")
	    	adata.obsm['X_pca'] = X_pca
	    	adata.uns['pca'] = {}
	    	adata.uns['pca']['params'] = {'zero_center': zero_center,
	    								  'use_highly_variable': use_highly_variable,}
	    	if use_highly_variable:
	    		adata.varm['PCs'] = np.zeros(shape=(adata.n_vars, n_comps))
	    		adata.varm['PCs'][adata.var['highly_variable']] = pca_.components_.T
	    	else:
	    		adata.varm['PCs'] = pca_.components_.T

	    	adata.uns['pca']['variance'] = pca_.explained_variance_
	    	adata.uns['pca']['variance_ratio'] = pca_.explained_variance_ratio_
	    	return adata, pca_ if copy else None
	    else:
	        # logg.info('    finished', time=start)
	        if return_info:
	            return X_pca, pca_.components_, pca_.explained_variance_ratio_, pca_.explained_variance_, pca_
	        else:
	            return X_pca, pca_

	## Run dimensional reduction analysis and clustering using KNN graph
	def run_analysis(self,adata, n_neighbors=None, n_pcs=None, spread=None, umap_pcs=2, transform_PCA=None,
					 min_dist=None, resolution=None, remove_batch_effects=None, do_tSNE=None):
		# If argument is None, set to instance attribute
		param_dict = {k:(v if v else getattr(self,k)) for (k,v) in locals().items()}

		## If doing dimensional reduction for machine learning
		if True: #transform_PCA is None:
			## Find cell type score for each cell based on a predefined set of gene lists
			# if self.cell_score_lists:
			# 	for file in self.cell_score_lists:
			# 		score_list = [line.rstrip('\n') for line in open(''.join([file,'.txt']),'r')]
			# 		#print(adata)

			# 		# temp_score_genes1 = copy.deepcopy(adata.obs[file+'_processed'])
			# 		# temp_score_genes2 = copy.deepcopy(adata.obs[file+'_raw'])

			# 		# sc.tl.score_genes(adata, score_list, ctrl_size=50, gene_pool=None, n_bins=25, score_name=file, random_state=0, copy=False, use_raw=False)
			# 		# print('processed')
			# 		# print(adata.obs[file+'_processed'] - temp_score_genes1)
			# 		sc.tl.score_genes(adata, score_list, ctrl_size=100, gene_pool=None, n_bins=50, score_name=file, random_state=0, copy=False, use_raw=True)
			# 		# print('raw')
			# 		# print(adata.obs[file+'_raw'] - temp_score_genes2)

			## Run PCA to compute the default number of components
			[adata, transform_PCA] = self.pca(adata, svd_solver='arpack', copy=True, transform_PCA=transform_PCA)
			print(adata.uns['pca']['params']['zero_center'])
			## Save the existing data to disk for later
			self.adata_postPCA = adata.copy()


			## Compute nearest-neighbors
			sc.pp.neighbors(adata, n_neighbors=param_dict['n_neighbors'], n_pcs=param_dict['n_pcs'])
			neighbor_graph = copy.deepcopy(adata.uns['neighbors'])

			## Remove batch effects
			# Note that doing this may override previous sc.pp.neighbors()
			if param_dict['remove_batch_effects']:
				import bbknn
				bbknn.bbknn(adata, batch_key='sampleName', copy=False)#, n_pcs=param_dict['n_pcs'], neighbors_within_batch=param_dict['n_neighbors'])
				#sc.pp.external.mnn_correct(adata,batch_key='sampleName') # Testing another algorithm

			## Run UMAP Dim reduction
			sc.tl.umap(adata, spread=param_dict['spread'], min_dist=param_dict['min_dist'], n_components=umap_pcs) # Min_dist needs to be between 0.01 to 0.5

			## Run tSNE analysis
			if param_dict['do_tSNE']:
				sc.tl.tsne(adata, n_pcs=param_dict['n_pcs'])

			## Calculate cell clusters via Louvain algorithm
			sc.tl.louvain(adata, resolution=param_dict['resolution'])

			## Do Dendrogram analysis based on PCs
			sc.tl.dendrogram(adata,groupby='louvain',n_pcs=param_dict['n_pcs'],linkage_method="median",use_raw=True)
		else:
			self.pca(adata,svd_solver='arpack',transform_PCA = transform_PCA)

		return [adata, transform_PCA]

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
			my_feature_cmap = self.make_cmap(feature_colors,bit=True)
		elif self.umap_color=='yellow_blue':
			feature_colors = [(210,210,210), (210,210,210), (245,245,200), (100,200,225), (0,45,125)]
			position=[0, 0.019999, 0.02, 0.55, 1]
			my_feature_cmap = self.make_cmap(feature_colors,bit=True,position=position)
		else:
			feature_colors = [(210,210,210), (210,210,210), (245,245,200), (100,200,225), (0,45,125)]
			position=[0, 0.019999, 0.02, 0.55, 1]
			my_feature_cmap = self.make_cmap(feature_colors,bit=True,position=position)

		gray_cmap = self.make_cmap([(220,220,220),(220,220,220)], bit=True)

		## Custom color palette for cluster plots and observation plots
		colors = [(1,0.5,0),(0.5,0.5,0.85),(0,1,0),(1,0,0),(0,0,0.9),(0,1,1),
				(0.4,0.4,0.4),(0.5,0.85,0.5),(0.5,0.15,0.5),
				(0.15,0.5,0.5),(0.5,0.5,0.15),(0.9,0.9,0),(1,0,1),
				(0,0.5,1),(0.85,0.5,0.5),(0.5,1,0),(0.5,0,1),(1,0,0.5),(0,0.9,0.6),
				(0.3,0.6,0),(0,0.3,0.6),(0.6,0.3,0),(0.3,0,0.6),(0,0.6,0.3),(0.6,0,0.3)]

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

		# ## Violin plots for filtering parameters pre and post
		# sc.pl.violin(self.adata_preFiltered, ['n_genes','n_counts','percent_mito'],
		#  			 jitter=0.4, multi_panel=True, save='_preFiltered_plot.png', show=False)
		# sc.pl.violin(adata, ['n_genes','n_counts','percent_mito'],
		# 			 jitter=0.4, multi_panel=True, save='_postFiltered_plot.png', show=False)

		## Draw the PCA elbow plot to determine which PCs to use
		sc.pl.pca_variance_ratio(self.adata_postPCA, log=True, n_pcs=100, save='_elbowPlot.png', show=False)
		## Ranks and displays most contributing genes for each principal component
		components = 4
		loadings_components = range(1,self.n_pcs+components+1)
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
		n_genes_rank = 5
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
		
				feature_positions = gene_dict['feature_positions'] # Manually set and determined
				feature_groups = gene_dict['feature_groups']
				groupby_positions = gene_dict['groupby_positions']

				if len(gene_dict['markers'])!=1:
					for grouping in self.exp_grouping:
						## Dotplot analysis
						# Circle color corresponds to expression level, and circle size corresponds to percentage of cells expressing gene

						## Reordering categories for dotplot or heatmap rows
						adata_plots = adata.copy()
						if groupby_positions:
							adata_plots.obs['louvain'] = adata.obs['louvain'].cat.reorder_categories(groupby_positions,inplace = False)

						sc.pl.dotplot(adata_plots, genes_to_plot, groupby=grouping, 
								var_group_positions=feature_positions, var_group_labels=feature_groups,
								save=''.join(['_markers_',gene_list,'_',grouping,file_type]), show=False, 
								color_map=my_feature_cmap, use_raw=True)#, dot_max=0.4)#, dendrogram=True)
						## Heatmaps
						# Each horizontal line represents expression of one cell
						sc.pl.heatmap(adata_plots, genes_to_plot, groupby=grouping, 
								var_group_positions=feature_positions, var_group_labels=feature_groups,
								save=''.join(['_markers_',gene_list,'_',grouping,file_type]), show=False, 
								cmap=my_feature_cmap, use_raw=True)

		# Generate a umap feature plot based on cell scoring
		if self.cell_score_lists:
			vmax = adata.obs.loc[:,self.cell_score_lists].values.max()
			vmin = adata.obs.loc[:,self.cell_score_lists].values.min()
			print(vmax)
			print(vmin)
			sc.pl.umap(adata, color=self.cell_score_lists, 
					   save='_cellType_score.png', show=False, edges=False, color_map=my_feature_cmap, 
					   size=size, vmin=vmin, vmax=vmax)
			sc.pl.umap(adata, color=self.cell_score_lists, 
					   save='_cellType_score_0min.png', show=False, edges=False, color_map=my_feature_cmap, 
					   size=size, vmin=0, vmax=vmax)

			sc.pl.violin(adata,self.cell_score_lists, 
						 jitter=0.4, save='_cell_scores.png',show=False,multi_panel=False,rotation=90)

		# Genes that are not expressed or are invariable are plotted using a grayscale
		print('Plotting empty genes: ',missing_genes,'\n')
		sc.pl.umap(adata, color=missing_genes, save=''.join(['_featureplots_gray',file_type]), 
				show=False, cmap=gray_cmap, size=size, use_raw=True)

		if self.do_tSNE:
			# tSNE Plots
			sc.pl.tsne(adata, color='louvain', save = '_clusterIdentity.png', show = False, 
						legend_loc = 'right margin', edges = False, size = size, 
						palette = colors, alpha = 0.75)
			sc.pl.tsne(adata, color='sampleName', save = '_sample.png', show = False, 
						legend_loc = 'right margin', edges = False, size = size, 
						palette = colors, alpha = 0.75)
			sc.pl.tsne(adata, color=genes_to_plot, save = '_featureplots.png', show = False, cmap = my_feature_cmap, size=size, use_raw = True)
			sc.pl.tsne(adata, color=missing_genes, save='_featureplots_gray.png', show=False, cmap=gray_cmap, size=size, use_raw=True)
		
		return adata

	## Most basic pipeline - Input data and output all figures
	def pipe_basic(self, figdir='./figures/', adata_filtered=None, load_save=None, new_save='adata_save.p', 
		remove_genes=None, final_quality=False, highly_variable=None, transform_X=None, transform_PCA=None, 
		save_transform=False):
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
			self.adata_unscaled = run_save.adata_unscaled
			self.adata_preFiltered = run_save.adata_preFiltered
			self.annotation_dict = run_save.annotation_dict
			self.initial_cell_count = run_save.initial_cell_count
			self.initial_gene_count= run_save.initial_gene_count
		else:
			if not adata_filtered:
				adata = self.load_data()

				## Filter and process data
				adata = self.filter_data(adata)
			else:
				adata=adata_filtered.copy()

			## Remove specific genes from the analysis (such as experimentally observed contaminations)
			if remove_genes:
				print('Removing specified list of genes from analysis')
				adata = self.filter_specific_genes(adata,text_file = remove_genes)

			self.adata_preProcessed = adata.copy() # Save for later in case of necessary extraction
			[adata,highly_variable,transform_X] = self.preprocess_data(adata, highly_variable, transform_X,regress_out=False,combat=self.combat)
		
		## Dimensional reduction and clustering - construction of the neighborhood graph
		[adata, transform_PCA] = self.run_analysis(adata, transform_PCA=transform_PCA)

		## Plot figures
		adata = self.plot_sca(adata,figdir=figdir)

		# ## Run UMAP for machine learning
		# adata = self.run_analysis(adata, umap_pcs=50)

		self.adata = adata.copy()
		## Save analysis information and relevant AnnData objects to the disk using the Pickle module
		if new_save:
			pickle.dump(self,open(''.join([figdir,new_save]),"wb"),protocol=4)

		## Write a summary of the analysis to a text file including sample information and parameters
		print("Writing Summary")
		self.write_summary(figdir=figdir)

		print("\nAll done!\n")

		return [self, highly_variable, transform_X, transform_PCA] if save_transform else self

	## Pipeline for analysis in which you extract interesting clusters/observations after an initial run
	# Extracts clusters to an filtered but unprocessed AnnData object, then reprocesses and reclusters
	def pipe_ext(self, analysis_params_ext, figdir='./figures/', extracted=None, highly_variable=None, transform_X=None,
				 transform_PCA=None, load_save=None, new_save='adata_save.p', final_quality=False, 
				 label='', save_transform=False):
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
			self.adata_unscaled = run_save.adata_unscaled
			self.adata_preFiltered = run_save.adata_preFiltered
			self.annotation_dict = run_save.annotation_dict
			self.initial_cell_count = run_save.initial_cell_count
			self.initial_gene_count= run_save.initial_gene_count

			## Create an unprocessed AnnData object with the desired clusters
			adata_ext = self.adata_preProcessed[adata.obs['louvain'].isin(extracted)].copy()
			self.adata_preProcessed = adata_ext.copy()
			## Reprocess and recluster extracted cells
			[adata_ext,highly_variable,transform_X] = self.preprocess_data(adata_ext, highly_variable, 
														transform_X, regress_out=False, combat=self.combat)
		else: ## Not advised to use this
			print("Looking to extract based on gene expression of ",extracted)
			self.load_data()
			self.filter_data(adata)
			self.adata_preProcessed = adata.copy()
			self.preprocess_data(adata)
			adata_ext = adata[adata.raw[:,extracted].X>1.5,:]

		# sca_dict_ext = copy.deepcopy(sca_dict)
		# sca_dict_ext.update(analysis_params = analysis_params_ext) # Use new analysis parameters
		self.set_analysis_params(**(analysis_params_ext))
		[adata_ext, transform_PCA] = self.run_analysis(adata_ext, transform_PCA=transform_PCA)

		## Plot figures
		adata_ext = self.plot_sca(adata_ext,figdir = ''.join([figdir,'extracted/',label,'/']))

		self.adata = adata_ext.copy()
		## Save analysis information and relevant AnnData objects to the disk using the Pickle module
		if new_save:
			pickle.dump(self,open(''.join([figdir,'extracted/',label,'/',new_save]),"wb"),protocol=4)

		## Write a summary of the analysis to a text file including sample information and parameters
		self.write_summary(figdir=''.join([figdir,'extracted/',label,'/']))

		print("\nAll done!\n")

		return [self, highly_variable, transform_X, transform_PCA] if save_transform else self

	def plot_predictions(self,adata,attr,predicted_key='predicted',figdir='./forest_test/figures/', data_label='all'):
		os.makedirs(os.path.dirname(''.join([figdir,data_label,'/'])), exist_ok=True) # Create directory if it doesn't exist
		
		## Custom color palette for cluster plots and observation plots
		colors = [(1,0.5,0),(0.5,0.5,0.85),(1,0,0),(0,1,0),(0,0,0.9),(0,1,1),
				(0.4,0.4,0.4),(0.5,0.85,0.5),(0.5,0.15,0.5),
				(0.15,0.5,0.5),(0.5,0.5,0.15),(0.9,0.9,0),(1,0,1),
				(0,0.5,1),(0.85,0.5,0.5),(0.5,1,0),(0.5,0,1),(1,0,0.5),(0,0.9,0.6),
				(0.3,0.6,0),(0,0.3,0.6),(0.6,0.3,0),(0.3,0,0.6),(0,0.6,0.3),(0.6,0,0.3)]
		samples = adata.obs['sampleName'].unique()
		df = pd.DataFrame()
		df['predicted'] = adata.obs[predicted_key].values
		#df['predicted_prob_class'] = adata.obs['predicted_prob_class'].values
		#df['predicted_prob_cal'] = adata.obs['predicted_prob_cal'].values
		df['sampleName'] = adata.obs['sampleName'].values
		# for sample in samples:
		# 	df[sample] = [predicted if sampleName==sample for predicted,sampleName in [adata.obs['predicted'],adata.obs['sampleName']]]

		#vplot, axes = plt.subplots(math.ceil(len(samples)/4),4, figsize=(18,12))

		# for i,sample in enumerate(samples):
		vplot = plt.figure(figsize=(24,12))
		sns.violinplot(x='sampleName', y='predicted', data=df, inner=None, scale='width')#, ax=axes[math.floor(i/4),i%4])
		sns.stripplot(x='sampleName', y='predicted', data=df, jitter = True, color='black', size=0.4)#, ax=axes[math.floor(i/4),i%4])
		vplot.savefig(''.join([figdir,data_label,'/violin_',attr,'.png']))

		# vplot = plt.figure(figsize=(24,12))
		# sns.violinplot(x='sampleName', y='predicted_prob_class', data=df, inner=None, scale='width')#, ax=axes[math.floor(i/4),i%4])
		# sns.stripplot(x='sampleName', y='predicted_prob_class', data=df, jitter = True, color='black', size=0.4)#, ax=axes[math.floor(i/4),i%4])
		# vplot.savefig(''.join([figdir,'/violin_all_class_',label,'.png']))

		# vplot = plt.figure(figsize=(24,12))
		# sns.violinplot(x='sampleName', y='predicted_prob_cal', data=df, inner=None, scale='width')#, ax=axes[math.floor(i/4),i%4])
		# sns.stripplot(x='sampleName', y='predicted_prob_cal', data=df, jitter = True, color='black', size=0.4)#, ax=axes[math.floor(i/4),i%4])
		# vplot.savefig(''.join([figdir,'/violin_',data_label,'_calibrated_',label,'.png']))

		# adata.obs['jitter'] = np.random.rand(len(adata.obs_names))*10
		# sc.pl.scatter(adata,x='jitter',y='predicted',color='sampleName',save=''.join([label,'.png']),palette=colors,show=False)

	## Writes random forest regressor feature_importances to a csv file
	def write_feature_importances(self,rf_model, train_adata, attr, figdir='./forest_test/figures/'):
		csv_fileName = ''.join([figdir,'/',attr,'_feature_importances.csv'])
		os.makedirs(os.path.dirname(csv_fileName), exist_ok=True)
		if attr is 'X':
			with open(csv_fileName,'w',newline='') as f:
				wr = csv.writer(f)
				wr.writerow(['Gene', 'Feature Importance'])
				wr.writerows(zip(train_adata.var_names,rf_model.feature_importances_))
		elif attr is 'X_pca':
			with open(csv_fileName,'w',newline='') as f:
				wr = csv.writer(f)
				wr.writerow(['PC', 'Feature Importance'])
				wr.writerows(zip(range(1,51),rf_model.feature_importances_))

	## Function that trains a random forest regressor to model fetal versus adult-ness
	# Plots results of training and results of model on test data
	def forest_regress(self, train_adata, test_dict, obs_key='is_adult', types='X', figdir='./forest_test/figures/', data_label='all'):
		try:
			train_data = getattr(train_adata,types)
		except:
			train_data = train_adata.obsm[types]

		train_labels=train_adata.obs[obs_key].values

		from sklearn.utils import shuffle
		shuffle_train_data, shuffle_train_labels = shuffle(train_data,train_labels)

		## Initialize 
		clf = RandomForestRegressor(n_estimators=150,oob_score=True,verbose=1,n_jobs=10)
		clf = clf.fit(shuffle_train_data, shuffle_train_labels)

		## Check fitting of training data
		train_adata.obs['predicted'] = clf.predict(train_data)

		## Plot training data results in violin plot
		self.plot_predictions(train_adata, attr='_'.join(['train',types]), figdir=figdir, data_label=data_label)
		self.write_feature_importances(clf, train_adata, types, figdir)

		print("Trained Regressor")

		## Check test data
		for key in test_dict:
			try:
				test_data = getattr(test_dict[key],types)
			except:
				test_data = test_dict[key].obsm[types]

			## Check fitting of test data
			test_dict[key].obs['predicted'] = clf.predict(test_data)

			## Plot test data results in violin plot
			self.plot_predictions(test_dict[key], attr='_'.join(['test',key,types]), figdir=figdir, data_label=data_label)

			# print(clf.feature_importances_)
			# print(clf.n_features_)
			# print(clf.n_outputs_)
			# print(clf.oob_score_)
			# print(clf.oob_prediction_)

		return clf

	## Function that trains a random forest regressor to model fetal versus adult-ness
	# Plots results of training and results of model on test data
	def forest_class(self, train_adata, valid_adata, test_dict, obs_key='is_adult', figdir='./forest_test/figures/'):
		train_data = train_adata.X
		train_labels = train_adata.obs[obs_key].values

		valid_data = valid_adata.X
		valid_labels = valid_adata.obs[obs_key].values

		# p = np.random.permutation(len(train_labels))
		# shuffle_train_data = train_data[p,:]
		# shuffle_train_labels = train_labels[p]

		# p = np.random.permutation(len(valid_labels))
		# shuffle_valid_data = valid_data[p,:]
		# shuffle_valid_labels = valid_labels[p]

		from sklearn.utils import shuffle
		print(len(train_data))
		shuffle_train_data, shuffle_train_labels = shuffle(train_data,train_labels)
		shuffle_valid_data, shuffle_valid_labels = shuffle(valid_data,valid_labels)


		clf_class = RandomForestClassifier(n_estimators=150,oob_score=True,verbose=1,n_jobs=10)
		clf_class = clf_class.fit(shuffle_train_data,shuffle_train_labels)

		clf_cal = CalibratedClassifierCV(clf_class,cv="prefit")
		clf_cal = clf_cal.fit(shuffle_valid_data, shuffle_valid_labels)

		## Check fitting of training data
		print(clf_class.predict_proba(train_data))
		print(clf_cal.predict_proba(train_data))
		if obs_key:
			train_adata.obs['predicted_prob_class'] = clf_class.predict_proba(train_data)[:,1]
			train_adata.obs['predicted_prob_cal'] = clf_cal.predict_proba(train_data)[:,1]
		else:
			train_adata.obs['predicted_prob_class'] = clf_class.predict_proba(train_data)
			train_adata.obs['predicted_prob_cal'] = clf_cal.predict_proba(train_data)

		## Plot training data results in violin plot
		self.plot_predictions(train_adata,'train', predicted_key='predicted_prob_class', figdir=figdir, data_label='prob_class')
		self.plot_predictions(train_adata,'train', predicted_key='predicted_prob_cal', figdir=figdir, data_label='prob_cal')

		for key in test_dict:
			test_data = test_dict[key].X

			## Check fitting of test data
			if obs_key:
				test_dict[key].obs['predicted_prob_class'] = clf_class.predict_proba(test_data)[:,1]
				test_dict[key].obs['predicted_prob_cal'] = clf_cal.predict_proba(test_data)[:,1]
			else:
				test_dict[key].obs['predicted_prob_class'] = clf_class.predict_proba(test_data)
				test_dict[key].obs['predicted_prob_cal'] = clf_cal.predict_proba(test_data)

			## Plot test data results in violin plot
			self.plot_predictions(test_dict[key],attr='_'.join(['test',key]), predicted_key='predicted_prob_class', figdir=figdir, data_label='prob_class')
			self.plot_predictions(test_dict[key],attr='_'.join(['test',key]), predicted_key='predicted_prob_cal', figdir=figdir, data_label='prob_cal')

		return [clf_class, clf_cal]
	#def basal_classifier(self, train_adata, test_dict, figdir='./figures/'):
		

	## Fetal adult classifier pipeline
	# Create a double classifier to first classify basal cells, and then use those to classify fetal vs adult data
	def double_class_pipe(self, adata, data_label='combined_cal', figdir='./figures/'):
		## Create my custom palette for FeaturePlots and define a matlplotlib colormap object
		if self.umap_color=='blue_orange':
			feature_colors = [(35,35,142), (255,127,0)]
			my_feature_cmap = self.make_cmap(feature_colors,bit=True)
		elif self.umap_color=='yellow_blue':
			feature_colors = [(210,210,210), (210,210,210), (245,245,200), (100,200,225), (0,45,125)]
			position=[0, 0.019999, 0.02, 0.55, 1]
			my_feature_cmap = self.make_cmap(feature_colors,bit=True,position=position)
		else:
			feature_colors = [(210,210,210), (210,210,210), (245,245,200), (100,200,225), (0,45,125)]
			position=[0, 0.019999, 0.02, 0.55, 1]
			my_feature_cmap = self.make_cmap(feature_colors,bit=True,position=position)

		# self = pickle.load(open(''.join([figdir,'all/','adata_save.p']),"rb"))
		## Run through pipeline and calculate scores
		# self.resolution = 1 # Make higher resolution for more clusters
		self.pipe_basic(adata_filtered=adata, figdir=''.join([figdir,'all/']))#,load_save='adata_save.p') ## batch corrected
		adata_mlp = self.adata_unscaled.copy() #normalized, batch corrected, removed NVGs
		adata = self.adata.copy()

		## Identify clusters in which average score is above 0 and label raw data
		cluster_scores = self.__count_cluster_scores(adata,['basal_markers'])
		basal_score_clusters = [key for key in cluster_scores if cluster_scores[key]>=0]
		adata_mlp.obs['basal_labels'] = [1 if cluster in basal_score_clusters else 0 for cluster in adata.obs['louvain']]

		basal_indices = adata_mlp.obs.index[adata_mlp.obs['basal_labels']==1].tolist()
		not_basal_indices = adata_mlp.obs.index[adata_mlp.obs['basal_labels']==0].tolist()

		basal_dim = min(len(basal_indices),len(not_basal_indices))
		basal_indices = random.sample(basal_indices,basal_dim)
		not_basal_indices = random.sample(not_basal_indices,basal_dim)
		adata_normalized = adata_mlp[adata_mlp.obs_names.isin(basal_indices+not_basal_indices)]
		print(adata_normalized)

		## Split data into training and testing data
		adata_normalized.obs['random_sample'] = [random.random() for x in range(len(adata_normalized.obs_names))]
		print(adata_normalized.obs.shape)
		print(adata_normalized.var.shape)
		print(adata_normalized.X.shape)

		adata_train_reg = adata_normalized[adata_normalized.obs['random_sample']<0.70].copy()
		print("Train Reg")
		print(len(adata_train_reg.obs['basal_labels']))
		print(len(adata_train_reg[adata_train_reg.obs['basal_labels']==1].obs_names))
		print(len(adata_train_reg[adata_train_reg.obs['basal_labels']==0].obs_names))
	
		adata_train = adata_normalized[adata_normalized.obs['random_sample']<0.45].copy()
		print("Train Class")
		print(len(adata_train.obs['basal_labels']))
		print(len(adata_train[adata_train.obs['basal_labels']==1].obs_names))
		print(len(adata_train[adata_train.obs['basal_labels']==0].obs_names))

		adata_valid = adata_normalized[((adata_normalized.obs['random_sample']>=0.45) &
								(adata_normalized.obs['random_sample']<0.7))].copy()
		adata_test = adata_normalized[adata_normalized.obs['random_sample']>=0.7].copy()


		test_dict_basal = {'basal_test': adata_test}

		## Train random forest classifier and calibrate to determine "basal-ness"
		[clf_class, clf_cal] = self.forest_class(adata_train,adata_valid,test_dict_basal,obs_key='basal_labels',figdir=figdir)

		## Compare with the random forest regressor results
		clf_regress = self.forest_regress(adata_train_reg,test_dict_basal,types='X',figdir=figdir,data_label='all',obs_key='basal_labels')

		pickle.dump([clf_class, clf_cal, clf_regress],open(''.join([figdir,'all/clf_basal.p']),"wb"),protocol=4)

		# [clf_class, clf_cal, clf_regress] = pickle.load(open(''.join([figdir,'all/clf_basal.p']),"rb"))

		## Go back to original data, and use classifier to predict basal-ness
		adata.obs['basal_class'] = clf_class.predict_proba(adata_mlp.X)[:,1]
		adata.obs['basal_cal'] = clf_cal.predict_proba(adata_mlp.X)[:,1]
		adata.obs['basal_reg'] = clf_regress.predict(adata_mlp.X)
		sc.settings.figdir = ''.join([figdir,'all/'])
		## Plot basal predictions on a UMAP
		sc.pl.umap(adata, color=['basal_markers','basal_class','basal_cal','basal_reg'], save='_basal_prediction.png', show=False, 
			edges=False, color_map=my_feature_cmap, size=self.size)

		sc.pl.violin(adata, keys=['basal_markers','basal_class','basal_cal','basal_reg'],groupby='dataset', 
					save='_basal_prediction.png', show=False, rotation=90)

		## Extract cell if basal probability is greater than 0.7 and plot
		adata_basal_plot = adata_mlp[(adata.obs['basal_cal']>=0.7)].copy()

		[adata_basal_plot,highly_variable,transform_X] = self.preprocess_data(adata_unscaled=adata_basal_plot)
		[adata_basal_plot,transform_PCA] = self.run_analysis(adata_basal_plot)
		adata_basal_plot = self.plot_sca(adata_basal_plot,figdir=''.join([figdir,'basal/']))
		self.write_summary(figdir=''.join([figdir,'basal/']))

		## Split basal data into training and testing data based on knowledge of datasets
		# print(len(((adata_basal.obs['sampleName']=='HT-187-Tracheal-epi') |
		# 			(adata_basal.obs['sampleName']=='HT234-Airway') |
		# 			(adata_basal.obs['sampleName']=='HT-189-Tracheal-Epi') |
		# 			(adata_basal.obs['sampleName']=='ND16737_Adult_Trachea') |
		# 			(adata_basal.obs['sampleName']=='ND17494_Adult_Trachea') |
		# 			(adata_basal.obs['sampleName']=='ND15989_Fresh_WT_Lung_Adult'))))
		adata_mlp.obs['random_sample'] = [random.random() for x in range(len(adata_mlp.obs_names))]
		adata_mlp.obs['is_adult'] = [1 if (cell=='ND16737_Adult_Trachea' or cell=='ND17494_Adult_Trachea' 
										   or cell=='ND15989_Fresh_WT_Lung_Adult') else 0 for cell in adata_mlp.obs['sampleName']]
		adata_train_reg = adata_mlp[(((adata_mlp.obs['sampleName']=='HT-187-Tracheal-epi') |
								(adata_mlp.obs['sampleName']=='HT234-Airway') |
							    (adata_mlp.obs['sampleName']=='HT-189-Tracheal-Epi') |
						     	(adata_mlp.obs['sampleName']=='ND16737_Adult_Trachea') |
								(adata_mlp.obs['sampleName']=='ND17494_Adult_Trachea') |
								(adata_mlp.obs['sampleName']=='ND15989_Fresh_WT_Lung_Adult')) &
								(adata.obs['basal_cal']>=0.7)), :].copy()
		adata_train = adata_mlp[(((adata_mlp.obs['sampleName']=='HT-187-Tracheal-epi') |
								(adata_mlp.obs['sampleName']=='HT234-Airway') |
							    (adata_mlp.obs['sampleName']=='HT-189-Tracheal-Epi') |
						     	(adata_mlp.obs['sampleName']=='ND16737_Adult_Trachea') |
								(adata_mlp.obs['sampleName']=='ND17494_Adult_Trachea') |
								(adata_mlp.obs['sampleName']=='ND15989_Fresh_WT_Lung_Adult')) &
								(adata.obs['basal_cal']>=0.7) &
								(adata_mlp.obs['random_sample']<0.66)), :].copy()
		adata_valid = adata_mlp[(((adata_mlp.obs['sampleName']=='HT-187-Tracheal-epi') |
								(adata_mlp.obs['sampleName']=='HT234-Airway') |
							    (adata_mlp.obs['sampleName']=='HT-189-Tracheal-Epi') |
						     	(adata_mlp.obs['sampleName']=='ND16737_Adult_Trachea') |
								(adata_mlp.obs['sampleName']=='ND17494_Adult_Trachea') |
								(adata_mlp.obs['sampleName']=='ND15989_Fresh_WT_Lung_Adult')) &
								(adata.obs['basal_cal']>=0.7) &
								(adata_mlp.obs['random_sample']>=0.66)), :].copy()

		print(adata_train_reg.obs['sampleName'])
		print(adata_train_reg.obs[adata_train_reg.obs['is_adult']==1].shape[0])
		print(adata_train_reg.obs[adata_train_reg.obs['is_adult']==0].shape[0])

		# adata_train.obs['random_sample'] = [random.random() for x in range(len(adata_train.obs_names))]
		# adata_train = adata_train[adata_train.obs['random_sample']<0.66, :].copy()
		# adata_valid = adata_train[adata_train.obs['random_sample']>=0.66, :].copy()

		adata_test = {}

		adata_test['1'] = adata_mlp[(((adata_mlp.obs['sampleName']=='HT-182-d125-lung-Prox') |
									(adata_mlp.obs['dataset']=='CFFT_Explant') |
									(adata_mlp.obs['dataset']=='FH')) &
									(adata.obs['basal_cal']>=0.7)),:].copy()

		adata_test['2'] = adata_mlp[(((adata_mlp.obs['dataset']=='CFFT_iBasal') |
						   			(adata_mlp.obs['dataset']=='CFFT_iBasal_AMY')) &
									(adata.obs['basal_cal']>=0.7)),:].copy()

		adata_test['3'] = adata_mlp[(((adata_mlp.obs['dataset']=='CFFT_P1') |
						   			(adata_mlp.obs['dataset']=='CFFT_ALI')) &
									(adata.obs['basal_cal']>=0.7)), :].copy()

		adata_train.obs['is_adult'] = [1 if (cell=='ND16737_Adult_Trachea' or cell=='ND17494_Adult_Trachea' 
										   or cell=='ND15989_Fresh_WT_Lung_Adult') else 0 for cell in adata_train.obs['sampleName']]

		## Train random forest classifer on basal data and calibrate to determin "adult-ness"
		[clf_class, clf_cal] = self.forest_class(adata_train,adata_valid,adata_test,obs_key='is_adult',figdir=figdir)

		## Compare with the random forest regressor results
		clf_regress = self.forest_regress(adata_train_reg,test_dict_basal,types='X',figdir=figdir,data_label='basal',obs_key='is_adult')

		pickle.dump([clf_class, clf_cal, clf_regress],open(''.join([figdir,'basal/clf_adult.p']),"wb"),protocol=4)

		# [clf_class, clf_cal, clf_regress] = pickle.load(open(''.join([figdir,'basal/clf_adult.p']),"rb"))
		print("Predicting")
		## Go back to original data, and use classifier to predict basal-ness
		adata_basal_plot.obs['adult_class'] = clf_class.predict_proba(self.adata_unscaled.X)[:,1]
		adata_basal_plot.obs['adult_cal'] = clf_cal.predict_proba(self.adata_unscaled.X)[:,1]
		adata_basal_plot.obs['adult_reg'] = clf_regress.predict(self.adata_unscaled.X)

		sc.settings.figdir = ''.join([figdir,'basal/'])

		## Plot basal predictions on a UMAP
		sc.pl.umap(adata_basal_plot, color=['adult_class','adult_cal','adult_reg'], save='_adult_prediction.png', show=False, 
			edges=False, color_map=my_feature_cmap, size=self.size)

		sc.pl.violin(adata_basal_plot, keys=['adult_class','adult_cal','adult_reg'],groupby='dataset',
					save='_adult_prediction.png', show=False, rotation=90)


	## Fetal adult classifier pipeline
	# Contains data leakage as training data and testing data go through entire basic pipeline together
	# Not recommended to be used
	def data_leaked_original(self, adata, types=['X','X_pca','X_umap'], data_label='all', figdir='./figures/'):
		## Put train data and test data 
		adata = self.pipe_basic(adata_filtered=adata,figdir=figdir).adata.copy()

		adata_train = adata[((adata.obs['sampleName']=='HT-187-Tracheal-epi') |
					(adata.obs['sampleName']=='HT234-Airway') |
					(adata.obs['sampleName']=='HT-189-Tracheal-Epi') |
					(adata.obs['sampleName']=='ND16737_Adult_Trachea') |
					(adata.obs['sampleName']=='ND17494_Adult_Trachea') |
					(adata.obs['sampleName']=='ND15989_Fresh_WT_Lung_Adult')), :].copy()

		adata_test = {}

		adata_test['1'] = adata[((adata.obs['sampleName']=='HT-182-d125-lung-Prox') |
							(adata.obs['dataset']=='CFFT_Explant') |
							(adata.obs['dataset']=='FH')),:].copy()

		adata_test['2'] = adata[((adata.obs['dataset']=='CFFT_iBasal') |
						   (adata.obs['dataset']=='CFFT_iBasal_AMY')),:].copy()

		adata_test['3'] = adata[((adata.obs['dataset']=='CFFT_P1') |
						   (adata.obs['dataset']=='CFFT_ALI')), :].copy()

		adata_train.obs['is_adult'] = [1 if (cell=='ND16737_Adult_Trachea' or cell=='ND17494_Adult_Trachea' 
								   or cell=='ND15989_Fresh_WT_Lung_Adult') else 0 for cell in adata_train.obs['sampleName']]

		classifiers = []
		for attr in types:
			classifiers.append(self.forest_regress(adata_train, adata_test, types=attr, figdir=figdir, data_label=data_label))

		sc.settings.figdir = ''.join([figdir,'data_leaked_original/'])

		## Plot basal predictions on a UMAP
		sc.pl.umap(adata_basal_plot, color=['adult_class','adult_cal','adult_reg'], save='_adult_prediction.png', show=False, 
			edges=False, color_map=my_feature_cmap, size=self.size)

		sc.pl.violin(adata_basal_plot, keys=['adult_class','adult_cal','adult_reg'],groupby='dataset',
		 			 jitter=0.4, multi_panel=True, save='_adult_prediction.png', show=False)

	## Fetal adult classifier pipeline
	# Fixed data leakage, saves transforms when scaling and doing pca to fix data leakage
	def minimized_dl_original(self, adata_train, adata_test, types=['X','X_pca'], data_label='all', figdir='./figures/'):
		##### No Data Leakage
		if 'X_pca' not in types:
			[adata_train,highly_variable,transform_X] = self.preprocess_data(adata_train)
			#adata_train = self.run_analysis()

			for key in adata_test:
				[adata_test[key],highly_variable,transform_X] = self.preprocess_data(adata_test[key],highly_variable=highly_variable, transform_X=transform_X)
		else:

			[an_run, highly_variable, transform_X, transform_PCA] = self.pipe_basic(adata_filtered=adata_train, figdir=''.join([figdir,'train/']), 
																		new_save=False, save_transform=True)
			adata_train = self.adata.copy()

			## Plot test data
			# for key in adata_test:
			# 	self.pipe_basic(adata_filtered=adata_test[key], figdir=''.join([figdir,key,'/']), new_save=False,
			# 			save_transform=False)#, highly_variable=highly_variable, transform_X=transform_X, transform_PCA=transform_PCA).adata

			for key in adata_test:
				[adata_test[key],highly_variable_test,transform_test] = self.preprocess_data(adata_test[key],highly_variable=highly_variable, transform_X=transform_X)

		for attr in types:
			self.forest_regress(adata_train, adata_test, types=attr, figdir=figdir, data_label=data_label)





