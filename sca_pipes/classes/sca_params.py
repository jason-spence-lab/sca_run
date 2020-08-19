'''
sca_params --
Class that handles all relevant parameters for setting up a SCARunner session

Written by Joshua H Wu
5 August, 2020
'''

from dataclasses import *
from typing import List

class sca_params:
	'''
	Class that handles all relevant paramaters for setting up a SCARunner session
	'''

	def __init__(self):
		self.storage_mount_point = 'Z:/'
		self.sample_list = []
		self.gene_lists = [] ## not needed?
		self.gene_dict = dict()

		## Filter, analysis and plot parameter dataclass objects
		self.filter_params = filter_params()
		self.analysis_params = analysis_params()
		self.plot_params = plot_params()

		## Summary Params - move to singleton
		self.initial_cell_count = None
		self.initial_gene_count = None
		self.final_cell_count = None
		self.final_gene_count = None
		self.annotation_dict = dict()
		self.doublet_clf = None

		## adata
		self.adata = None
		self.adata_preFiltered = None
		# self.adata_unscaled = None
		self.adata_postFiltered = None


	## Function to set a list of genes and save relevant information in a dictionary of gene list objects
	def add_gene_list(self, 
					  markers=['EPCAM','CDH5','VIM','TP63','NEUROG3'], 
					  label='basic_list', 
					  feature_positions=None, 
					  feature_groups=None, 
					  groupby_positions=None):
		'''
		Gene List Params --
			markers: List of interesting genes to be plotted
			label: Name of gene list
			feature_positions: Index positions to highlight groups of markers - will draw a bracket around them
			feature_groups: Labels of groups of markers highlighted
			groupby_positions: Ordering of metadata groups on dot plot y-axis
		'''
		# self.gene_lists = self.gene_lists + [label]
		self.gene_dict[label] = gene_list(markers=markers, 
										  label=label,
										  feature_positions=feature_positions,
										  feature_groups=feature_groups,
										  groupby_positions=groupby_positions)
		return

	## Creates object containing all of the filtering parameters and sets as attribute
	def set_filter_params(self,
						  min_cells=0,
						  min_genes=500,
						  max_genes=7000,
						  max_counts=30000,
						  max_mito=0.1):
		'''
		Filter Params --
			min_cells: Filter out genes with a few number of cells
			min_genes: Filter out cells with fewer genes to remove dead cells
			max_genes: Filter out cells with more genes to remove most doublets
			max_counts: Filter out cells with more UMIs to catch a few remaining doublets
			max_mito: Filter out cells with high mitochondrial gene content
		'''
		self.filter_params = filter_params(min_cells=0,
					  					   min_genes=500,
					  					   max_genes=7000,
					  					   max_counts=30000,
					  					   max_mito=0.1)
		return

	## Return dictionary with filtering parameters
	def get_filter_dict(self):
		return asdict(self.filter_params)

	## Creates object containing all of the analysis parameters and sets as attribute
	def set_analysis_params(self,
							n_neighbors=30,
							n_pcs=15,
							spread=1,
							min_dist=0.4,
							resolution=0.6,
							do_bbknn=False,
							do_tSNE=False,
							cell_score_lists=[]):
		'''
		Analysis Params --
			n_neighbors: Size of the local neighborhood used for manifold approximation
			n_pcs: Number of principle components to use in construction of neighborhood graph
			spread: In combination with min_dist determines how clumped embedded points are
			min_dist: Minimum distance between points on the umap graph
			resolution: High resolution attempts to increases # of clusters identified
			do_bbknn: Run batch balanced k-nearest neighbors batch correction algorithm
			cell_score_lists: List of file names to load cell scoring markers
		'''

		self.analysis_params = analysis_params(n_neighbors=30,
											   n_pcs=15,
											   spread=1,
											   min_dist=0.4,
											   resolution=0.6,
											   do_bbknn=False,
											   do_tSNE=False,
											   cell_score_lists=[])
		return

	## Return dictionary with analysis parameters
	def get_analysis_dict(self):
		return asdict(self.analysis_params)

	## Creates object containing all of the plotting parameters and sets as attribute
	def set_plot_params(self,
						size=20,
						umap_obs=['louvain','sampleName'],
						exp_grouping=['louvain'],
						umap_categorical_color='default',
						umap_feature_color='yellow_blue',
						final_quality=False,
						vmin_list=[],
						vmax_list=[],
						rank_grouping=['louvain'],
						clusters2_compare = ['all']):
		'''
		Plot Params --
			size: Size of dots on all UMAP plots
			umap_obs: Metadata observations by which to color UMAP plots
			exp_grouping: Metadata observations by which to group dot plot rows
			umap_categorical_color: List of colors for UMAP groups (HEX codes or RGB tuples)
			umap_feature_color: Color scheme for UMAP gradient plots (yellow_blue or blue_orange)
			vmin_list: List of y-axis minimums for cell scoring feature plots
			vmax_list:List of y-axis maximums for cell scoring feature plots
			rank_grouping: Metadata observations to group differential gene analyses
			clusters2_compare: Clusters to compare for differential gene analysis
		'''

		self.plot_params = plot_params(size=20,
									   umap_obs=['louvain','sampleName'],
									   exp_grouping=['louvain'],
									   umap_categorical_color='default',
									   umap_feature_color='yellow_blue',
									   final_quality=False,
									   vmin_list=[],
									   vmax_list=[],
									   rank_grouping=['louvain'],
									   clusters2_compare = ['all'])
		return

	## ---DEPRECATED---
	# Function that checks if an argument exists in a function and sets instance attributes based on them
	def __set_param_attr(self,arg_dict):
		for key in arg_dict:
			if key != 'self':
				setattr(self,key,arg_dict[key])


@dataclass
class filter_params:
	'''
	Filter Params DataClass --
		min_cells: Filter out genes with a few number of cells
		min_genes: Filter out cells with fewer genes to remove dead cells
		max_genes: Filter out cells with more genes to remove most doublets
		max_counts: Filter out cells with more UMIs to catch a few remaining doublets
		max_mito: Filter out cells with high mitochondrial gene content
	'''
	min_cells: int=0 
	min_genes: int=500
	max_genes: int=7000
	max_counts: int=30000
	max_mito: float=0.1

@dataclass
class analysis_params:
	'''
	Analysis Params DataClass --
		n_neighbors: Size of the local neighborhood used for manifold approximation
		n_pcs: Number of principle components to use in construction of neighborhood graph
		spread: In combination with min_dist determines how clumped embedded points are
		min_dist: Minimum distance between points on the umap graph
		resolution: High resolution attempts to increases # of clusters identified
		do_bbknn: Run batch balanced k-nearest neighbors batch correction algorithm
		cell_score_lists: List of file names to load cell scoring markers
	'''
	n_neighbors: int=30
	n_pcs: int=15
	spread: int=1
	min_dist: float=0.4
	resolution: float=0.6
	do_bbknn: bool=False
	do_tSNE: bool=False
	cell_score_lists: List[str] = field(default_factory=list)

@dataclass
class plot_params:
	'''
	Plot Params DataClass --
		size: Size of dots on all UMAP plots
		umap_obs: Metadata observations by which to color UMAP plots
		exp_grouping: Metadata observations by which to group dot plot rows
		umap_categorical_color: List of colors for UMAP groups (HEX codes or RGB tuples)
		umap_feature_color: Color scheme for UMAP gradient plots (yellow_blue or blue_orange)
		vmin_list: List of y-axis minimums for cell scoring feature plots
		vmax_list:List of y-axis maximums for cell scoring feature plots
		rank_grouping: Metadata observations to group differential gene analyses
		clusters2_compare: Clusters to compare for differential gene analysis
	'''
	size: int=20
	umap_obs: List[str] = field(default_factory=lambda: ['louvain','sampleName'])
	exp_grouping: List[str] = field(default_factory=lambda: ['louvain'])
	final_quality: bool=False
	umap_categorical_color: List[str] = field(default_factory=lambda: ['default'])
	umap_feature_color: str='yellow_blue'
	vmin_list: List[int] = field(default_factory=list)
	vmax_list: List[int] = field(default_factory=list)
	rank_grouping: List[str] = field(default_factory=lambda: ['louvain'])
	clusters2_compare: List[str] = field(default_factory=lambda: ['all'])

@dataclass
class gene_list:
	'''
	Gene List Params DataClass --
		markers: List of interesting genes to be plotted
		label: Name of gene list
		feature_positions: Index positions to highlight groups of markers - will draw a bracket around them
		feature_groups: Labels of groups of markers highlighted
		groupby_positions: Ordering of metadata groups on dot plot y-axis
	'''
	markers: List[str] = field(default_factory=lambda: ['EPCAM','CDH5','VIM','TP63'])
	label: str='basic_list'
	feature_positions: List[str] = field(default_factory=list)
	feature_groups: List[str] = field(default_factory=list)
	groupby_positions: List[str] = field(default_factory=list)