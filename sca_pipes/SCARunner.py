'''
Applications of single cell data analysis techniques
Written by Josh Wu and Mike Czerwinski
8 May, 2020

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
	sc.settings.verbosity = 3
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
				   only_plot=False):
		'''
		sca_dict is a dictionary of miscellaneous analysis information including
		parameters, sample list and gene_lists

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
					adata = ld.load().copy()
					sca_params.initial_cell_count = ld.initial_cell_count
					sca_params.initial_gene_count = ld.initial_gene_count
					sca_params.annotation_dict = ld.annotation_dict

				## Filter and process data
				qc = quality_control.quality_control(sca_params.qc_params, sca_params.species)
				adata = qc.run_qc(adata).copy()
				sca_params.doublet_clf = qc.doublet_clf
				sca_params.adata_preQC = qc.adata_preQC.copy()
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
	def pipe_ext(self, sca_params, extracted, figdir='./figures/', load_save=None, new_save='adata_save.p', label='',
				 no_preprocess=False):
		'''
		Allows loading of a saved pickle adata file, file must contained adata that has gone through a complete pipeline
		Otherwise will complete a new full analysis, and then extracts clusters
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
				adata_ext = sca_params.adata_postQC[adata.obs['louvain'].isin(extracted)].copy()
				sca_params.adata_postQC = adata_ext.copy()

				## Reprocess and recluster extracted cells
				pp = preprocess.preprocess(sca_params.gene_dict, sca_params.pp_params)
				adata_ext = pp.run_preprocess(adata_ext)
				sca_params.adata_unscaled = pp.adata_unscaled.copy()
			else:
				sca_params.adata_postQC = sca_params.adata_postQC[adata.obs['louvain'].isin(extracted)].copy()
				# Need to add save sca_params.adata_unscaled
				adata_ext = sca_params.adata[adata.obs['louvain'].isin(extracted)].copy()
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