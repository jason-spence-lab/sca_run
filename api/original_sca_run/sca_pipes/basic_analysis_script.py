'''
BASIC SINGLE CELL ANALYSIS SCRIPT
by Josh Wu
4 June, 2019

Relies heavily on the Scanpy Python module developed by the Theis Lab
Read more about Scanpy at https://scanpy.readthedocs.io/en/latest/index.html

Contains analysis of kidney samples obtained by Emily Holloway
Template for a Basic Analysis 

In progress ---
Moving to encapsulate parameters and relevant functions using class sca_set()
'''

from classes.sca_params import *
import SCARunner as SCARunner

figdir = './figures/'
an_params = sca_params()
#############################################################################
## Change this to point toward your mount location for our MiStorage share ##
#############################################################################
an_params.storage_mount_point = 'Z:/'
an_params.species = 'human'

## IDs of samples as represented in the metadata table
an_params.sample_list = ['2757-1','2761-1']

## List of interesting genes
an_params.add_gene_list(markers = ['CDH5','KDR','FLT1','NOS3','VWF','EMCN','CDH1','KRT8','EPCAM','ITGAM'],
					    label = 'basic_list_1')

an_params.add_gene_list(markers = ['PTPRC','COL1A1','COL1A2','PDGFRA','S100B','STMN2','TUBB3'],
						label = 'basic_list_1',
						cell_score_list = 'Both')

## Parameters used to filter the data - Mainly used to get rid of bad cells
an_params.set_qc_params(min_cells = 0, # Filter out genes with a few number of cells
					    min_genes = 500, # Filter out cells with fewer genes to remove dead cells
					    max_genes = 7000, # Filter out cells with more genes to remove most doublets
					    max_counts = 30000, # Filter out cells with more UMIs to catch a few remaining doublets
					    max_mito = 0.15, # Filter out cells with high mitochondrial gene content
					    doublet_detection = False) # Run DoubletDetection by Jonathan Shor

## Parameters used for initial clustering analysis
an_params.set_analysis_params(n_neighbors = 15, # Size of the local neighborhood used for manifold approximation
						      n_pcs = 11, # Number of principle components to use in construction of neighborhood graph
						      spread = 1, # In combination with min_dist determines how clumped embedded points are
						      min_dist = 0.4, # Minimum distance between points on the umap graph
						      resolution = 0.5, # High resolution attempts to increases # of clusters identified
						      do_bbknn = False, # Run batch balanced k-nearest neighbors batch correction algorithm
						      do_tSNE = False) # Run tSNE dimensional reduction analysis
						      
an_params.set_plot_params(size=20,
						  umap_obs=['louvain','sampleName','age'],
						  exp_grouping=['louvain'],
						  umap_categorical_color='default',
						  umap_feature_color='yellow_blue',
						  final_quality=False,
						  vmin_list=[],
						  vmax_list=[],
						  rank_grouping=['louvain'],
						  clusters2_compare = ['all'])

## Basic pipeline for analysis - will filter data, process, cluster, etc. and output relevant figures
an_run = SCARunner.SCARunner()
an_run.pipe_basic(an_params,figdir=figdir)

## If you find some interesting clusters that you want to "zoom in" on and recluster, you can use the following code

# New analysis parameters for the subset of parameters
analysis_params_ext = dict(n_neighbors = 9,
						n_pcs = 10,
						spread = 1,
						min_dist = 0.4,
						resolution = 0.4)

# an_run.pipe_ext(analysis_params_ext, figdir=figdir, extracted=['2'], load_save='adata_save.p')

