import scanpy as sc
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import copy
import os
import csv
import pandas as pd
from matplotlib import gridspec
import services.plotting as plotting

class cellphonedb:

	def __init__(self, cellphonedb_params):
        '''
        Cellphonedb Params --
            max_cbar_height: Maximum colorbar height
            plot_type: 'means' for regular heatmap of significant means per interaction, 'counts' for matrix plot or total significant means for each grouping
            y_labels: Add y axis labels to heatmap
        '''
        self.max_cbar_height = cellphonedb_params.max_cbar_height
        self.plot_type = cellphonedb_params.plot_type
        self.y_labels = cellphonedb_params.y_labels

## Writing metadata and counts csv files to be used as inputs into CellphoneDB analysis
	# Uses log-normalized counts
	def write_cpdb_data(self, adata, sca_params, figdir='/figures/'):

		df_meta = pd.DataFrame(data=[])
		os.makedirs(os.path.dirname(''.join([figdir,'csv_files/'])), exist_ok=True) 
		adata.obs.loc[:,[sca_params.analysis_params.clustering_choice]].to_csv(''.join([figdir,'csv_files/metadata.csv'])) #clustering_choice

		## Export raw counts file
		adata_postQC = sca_params.adata_postQC.copy()
		# sc.pp.normalize_total(adata_postfiltered)#,target_sum=10000)
		adata_raw = adata.raw.copy()

		df = pd.DataFrame(adata_raw.X.T.toarray())
		df.columns = adata_postQC.obs.index
		df.set_index(adata_postQC.var.index, inplace=True)
		print(df)

		df.to_csv(''.join([figdir,'/csv_files/counts.csv']))
		return 0

	def plot_colorbar(self, mappable, fig, subplot_spec, max_cbar_height: float = 4.0):
		"""
		Plots a vertical color bar based on mappable.
		The height of the colorbar is min(figure-height, max_cmap_height)
		Parameters
		----------
		mappable
			The image to which the colorbar applies.
		fig
			The figure object
		subplot_spec
			The gridspec subplot. Eg. axs[1,2]
		max_cbar_height
			The maximum colorbar height
		Returns
		-------
		color bar ax
		"""
		width, height = fig.get_size_inches()
		if height > max_cbar_height:
			# to make the colorbar shorter, the
			# ax is split and the lower portion is used.
			axs2 = gridspec.GridSpecFromSubplotSpec(
				2,
				1,
				subplot_spec=subplot_spec,
				height_ratios=[height - max_cbar_height, max_cbar_height],
			)
			heatmap_cbar_ax = fig.add_subplot(axs2[1])
		else:
			heatmap_cbar_ax = fig.add_subplot(subplot_spec)
		plt.colorbar(mappable, cax=heatmap_cbar_ax)
		return heatmap_cbar_ax

	#### Function for plotting cellphonedb data
	# plot_type is 'means' for regular heatmap of significant means per interaction
	# plot_type is 'counts' for matrix plot or total significant means for each grouping
	def plot_cpdb_heatmap(self, sca_params, plotted_vals, plot_type='means', figsave='cellphone_heatmap.png', y_labels=[]):
		if plot_type=='means':
			feature_colors = [(210,210,210), (210,210,210), (245,245,200), (100,200,225), (0,45,125)]
			position=[0, 0.019999, 0.02, 0.55, 1]
		elif plot_type=='counts':
			feature_colors = [(220,220,220), (25,25,25)]
			position=[0,1]

		scplt = plotting.plotting(sca_params.plot_params)
		my_feature_cmap = scplt.make_cmap(feature_colors,position=position,bit=True)

		colorbar_width = 0.3

		height = 5
		heatmap_width = 5
		width = heatmap_width

		height_ratios = [0, height]

		width_ratios = [
			heatmap_width,
			colorbar_width,
		]
		fig = plt.figure(figsize=(width+2, height))

		axs = gridspec.GridSpec(
			nrows=2,
			ncols=2,
			width_ratios=width_ratios,
			wspace=0.15 / width,
			hspace=0.13 / height,
			height_ratios=height_ratios,
		)

		heatmap_ax = fig.add_subplot(axs[1, 0])
		im = heatmap_ax.imshow(
			plotted_vals.values, aspect='auto', interpolation="nearest",
			cmap=my_feature_cmap
		)
		heatmap_ax.set_ylim(plotted_vals.shape[0] - 0.5, -0.5)
		heatmap_ax.set_xlim(-0.5, plotted_vals.shape[1] - 0.5)
		heatmap_ax.grid(False)
		heatmap_ax.tick_params(axis='x', labelsize='small', labelrotation=90)
		heatmap_ax.set_xticks(np.arange(len(plotted_vals.columns)))
		heatmap_ax.set_xticklabels(plotted_vals.columns)

		if plot_type=='means':
			if y_labels:
				heatmap_ax.tick_params(axis='y', labelsize='small')
				y_label_indices = [plotted_vals.index.get_loc(label) for label in y_labels]
				# print(plotted_vals.index.get_loc('NOTCH4_DLL4','PDGFB-PDGFRB','CCL12_CXCR4','EGFR_HBEGF'))
				heatmap_ax.set_yticks(y_label_indices)
				heatmap_ax.set_yticklabels(y_labels)
			else:
				# heatmap_ax.tick_params(axis='y', left=False, labelleft=False)
				heatmap_ax.tick_params(axis='y', labelsize='small')
				# heatmap_ax.set_ylabel('')
				heatmap_ax.set_yticks(np.arange(len(plotted_vals.index)))
				heatmap_ax.set_yticklabels(plotted_vals.index)
		elif plot_type=='counts':
			heatmap_ax.tick_params(axis='y', labelsize='small')
			heatmap_ax.set_yticks(np.arange(len(plotted_vals.index)))
			heatmap_ax.set_yticklabels(plotted_vals.index)
		else:
			heatmap_ax.tick_params(axis='y', left=False, labelleft=False)
			heatmap_ax.set_ylabel('')
		
		self.plot_colorbar(im, fig, axs[1, 1])
		plt.savefig(figsave,bbox_inches='tight')