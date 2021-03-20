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
		self.species = 'human'
		self.sample_list = []
		self.gene_lists = [] 
		self.gene_dict = dict()
		self.cell_score_lists = []

		## Quality control, analysis, and plot parameter dataclass objects
		self.qc_params = qc_params()
		self.analysis_params = analysis_params()
		self.plot_params = plot_params()

		## Summary Params
		self.initial_cell_count = None
		self.initial_gene_count = None
		self.final_cell_count = None
		self.final_gene_count = None
		self.annotation_dict = dict()

		## adata
		self.adata = None
		self.adata_preQC = None
		self.adata_unscaled = None
		self.adata_postQC = None


	## Function to set a list of genes and save relevant information in a dictionary of gene list objects
	def add_gene_list(self, 
					  markers=['EPCAM','CDH5','VIM','TP63','NEUROG3'], 
					  label='basic_list', 
					  feature_positions=None, 
					  feature_groups=None, 
					  groupby_positions=None,
					  cell_score_list=None,
					  load_file=None,):
		'''
		Gene List Params --
			markers: List of interesting genes to be plotted
			label: Name of gene list
			feature_positions: Index positions to highlight groups of markers - will draw a bracket around them
			feature_groups: Labels of groups of markers highlighted
			groupby_positions: Ordering of metadata groups on dot plot y-axis
			cell_score_list: 'None', 'Only', or 'Both'
			load_file: Allows loading a gene list from a line separated text file
		'''
		if load_file:
			markers = [line.rstrip('\n') for line in open(''.join([load_file,'.txt']),'r')]
		if not cell_score_list is 'Only':
			self.gene_lists.append(label)
		self.gene_dict[label] = gene_list(markers=markers, 
										  label=label,
										  feature_positions=feature_positions,
										  feature_groups=feature_groups,
										  groupby_positions=groupby_positions,
										  cell_score_list=cell_score_list)

		return

	## Creates object containing all of the filtering parameters and sets as attribute
	def set_qc_params(self,
					  min_cells=0,
					  min_genes=500,
					  max_genes=7000,
					  max_counts=30000,
					  max_mito=0.1,
					  doublet_detection = False):
		'''
		Quality Control Params --
			min_cells: Filter out genes with a few number of cells
			min_genes: Filter out cells with fewer genes to remove dead cells
			max_genes: Filter out cells with more genes to remove most doublets
			max_counts: Filter out cells with more UMIs to catch a few remaining doublets
			max_mito: Filter out cells with high mitochondrial gene content
			doublet_detection: Run DoubletDetection by Jonathan Shor
		'''
		self.qc_params = qc_params(min_cells=min_cells,
			  					   min_genes=min_genes,
			  					   max_genes=max_genes,
			  					   max_counts=max_counts,
			  					   max_mito=max_mito,
			  					   doublet_detection=doublet_detection)
		return

	## Return dictionary with filtering parameters
	def get_qc_dict(self):
		return asdict(self.qc_params)

	## Creates object containing all of the analysis parameters and sets as attribute
	def set_analysis_params(self,
							n_neighbors=30,
							n_pcs=15,
							spread=1,
							min_dist=0.4,
							resolution=0.6,
							do_bbknn=False,
							do_tSNE=False):
		'''
		Analysis Params --
			n_neighbors: Size of the local neighborhood used for manifold approximation
			n_pcs: Number of principle components to use in construction of neighborhood graph
			spread: In combination with min_dist determines how clumped embedded points are
			min_dist: Minimum distance between points on the umap graph
			resolution: High resolution attempts to increases # of clusters identified
			do_bbknn: Run batch balanced k-nearest neighbors batch correction algorithm
			do_tSNE: Run tSNE dimensional reduction analysis
		'''

		self.analysis_params = analysis_params(n_neighbors=n_neighbors,
											   n_pcs=n_pcs,
											   spread=spread,
											   min_dist=min_dist,
											   resolution=resolution,
											   do_bbknn=do_bbknn,
											   do_tSNE=do_tSNE)
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
						vmin_list=[],
						vmax_list=[],
						rank_grouping=['louvain'],
						clusters2_compare = ['all'],
						final_quality = False):
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
			final_quality: Makes high resolution figures in pdf
		'''

		self.plot_params = plot_params(size=size,
									   umap_obs=umap_obs,
									   exp_grouping=exp_grouping,
									   umap_categorical_color=umap_categorical_color,
									   umap_feature_color=umap_feature_color,
									   vmin_list=vmin_list,
									   vmax_list=vmax_list,
									   rank_grouping=rank_grouping,
									   clusters2_compare = clusters2_compare,
									   final_quality = final_quality)
		return

	## Return dictionary with analysis parameters
	def get_plot_dict(self):
		return asdict(self.plot_params)

	## ---DEPRECATED---
	# Function that checks if an argument exists in a function and sets instance attributes based on them
	def __set_param_attr(self,arg_dict):
		for key in arg_dict:
			if key != 'self':
				setattr(self,key,arg_dict[key])

@dataclass
class qc_params:
	'''
	Quality Control Params DataClass --
		min_cells: Filter out genes with a few number of cells
		min_genes: Filter out cells with fewer genes to remove dead cells
		max_genes: Filter out cells with more genes to remove most doublets
		max_counts: Filter out cells with more UMIs to catch a few remaining doublets
		max_mito: Filter out cells with high mitochondrial gene content
		doublet_detection: Run DoubletDetection by Jonathan Shor
	'''
	min_cells: int=0 
	min_genes: int=500
	max_genes: int=7000
	max_counts: int=30000
	max_mito: float=0.1
	doublet_detection: bool=False

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
	'''
	n_neighbors: int=30
	n_pcs: int=15
	spread: int=1
	min_dist: float=0.4
	resolution: float=0.6
	do_bbknn: bool=False
	do_tSNE: bool=False

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
		final_quality: Makes high resolution figures in pdf
	'''
	size: int=20
	umap_obs: List[str] = field(default_factory=lambda: ['louvain','sampleName'])
	exp_grouping: List[str] = field(default_factory=lambda: ['louvain'])
	umap_categorical_color: List[str] = field(default_factory=lambda: ['default'])
	umap_feature_color: str='yellow_blue'
	vmin_list: List[int] = field(default_factory=list)
	vmax_list: List[int] = field(default_factory=list)
	rank_grouping: List[str] = field(default_factory=lambda: ['louvain'])
	clusters2_compare: List[str] = field(default_factory=lambda: ['all'])
	final_quality: bool=False

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
	markers: List[str] = field(default_factory=lambda: ['EPCAM','CDH5','VIM','TP63','NEUROG3'])
	label: str='basic_list'
	feature_positions: List[str] = field(default_factory=list)
	feature_groups: List[str] = field(default_factory=list)
	groupby_positions: List[str] = field(default_factory=list)
	cell_score_list: str='None'