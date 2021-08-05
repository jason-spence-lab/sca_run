'''
Testing basic functionality of sca_params.py

Written by Joshua Wu
5 August, 2020
'''

from sca_params import *

figdir = './figures/'
params_test = sca_params()

sca_params.storage_mount_point = 'X:/'
sca_params.species = 'mouse'

## List of interesting genes
params_test.add_gene_list(markers= ['CDH5','VIM','FLT1','TP63','VWF','EMCN','CDH1','KRT8','EPCAM'],
					 	  label='basic_list',
					 	  feature_positions=[(0,2),(3,5),(6,8)],
					 	  feature_groups=['Group1','Group2','Group3'],
					 	  groupby_positions=['1','2','0','4','5','6'])

print("Gene List")
print(params_test.gene_dict)
print("\n")

## Parameters used to filter the data - Mainly used to get rid of bad cells
params_test.set_filter_params(min_cells = 0, 
							  min_genes = 500, 
							  max_genes = 7000, 
							  max_counts = 30000, 
							  max_mito = 0.15) 

print("Filter Params")
print(params_test.filter_params)
print(params_test.get_filter_dict())
print("\n")

## Parameters used for initial clustering analysis
params_test.set_analysis_params(n_neighbors = 15, 
						   n_pcs = 11, 
						   spread = 1, 
						   min_dist = 0.4, 
						   resolution = 0.5, 
						   do_tSNE=True,
						   cell_score_lists=True)

print("Analysis Params")
print(params_test.analysis_params)
print(params_test.get_analysis_dict())
print("\n")

params_test.set_plot_params(size=5,
					   umap_obs=['age'],
					   exp_grouping=['louvain'],
					   umap_categorical_color='default',
					   umap_feature_color='yellow_blue',
					   final_quality=False,
					   vmin_list=[0,0.1],
					   vmax_list=[1,1.1],
					   rank_grouping=['louvain'],
					   clusters2_compare = ['all'])

print("Plot Params")
print(params_test.plot_params)
print(params_test.get_plot_dict())