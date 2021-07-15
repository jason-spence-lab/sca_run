'''
load_data -- 
Loads single cell RNA sequencing data into AnnData object for downstream analyses

Written by Joshua H Wu
09 March, 2021
'''
from pathlib import Path
import scanpy as sc
import pandas as pd
import os

class load_data:
	def __init__(self,
				 storage_mount_point=None,
				 sample_list=None,
				 remove_genes=False):
		self.storage_mount_point=storage_mount_point
		self.sample_list=sample_list
		self.remove_genes=remove_genes

		self.initial_cell_count=None
		self.initial_gene_count=None
		self.annotation_dict={}
	'''
	Load Data Params --
		storage_mount_point: Data storage mount location to turbo or MIstorage
		sample_list: List of interested samples 
		remove_genes: List of unnecessary genes that should not be used in analysis
		initial_cell_count: Inital count of cells in the samples before analysis
		initial_gene_count: Inital count of genes in the samples before analysis
	'''

	## Loads data from storage point and creates an AnnData object 
	# Adds metadata to adata object 
	def load(self, file_type=None):
		'''
		Storage structure at the mount point includes 2 folders - processed data and raw
		The raw data folder contains folders of sample runs as well as the meta-data table
		'''
		print("Loading data into AnnData object")
		
		## Add raw data file directory if it doesn't exist
		os.makedirs(os.path.dirname('./data/'), exist_ok=True)
		## Location to output the anndata h5ad files
		raw_data_file = ''.join(['./data/Data_','_'.join(self.sample_list),
								 '.scanpy.raw.h5ad'])  # the file that will store the raw combined data
		results_file = ''.join(['./data/Data_','_'.join(self.sample_list),
								'.processed.h5ad'])  # the file that will store the analysis results

		if file_type == '.h5':
			# Load meta-data table in xls format located at in the storage_mount_point - Change if located elsewhere
			annotation_df = pd.read_excel(''.join([self.storage_mount_point,'01_RNAseq_RAW_Data/single_cell_meta_data_table_excel.xls']),
										  header = None)
		elif file_type == '.loom':
			annotation_df = pd.read_excel(''.join([self.storage_mount_point,'01_RNAseq_RAW_Data/single_cell_meta_data_table_loom.xls']),
										  header = None)
		else:
			print("No file type selected - assuming .h5")
			annotation_df = pd.read_excel(''.join([self.storage_mount_point,'01_RNAseq_RAW_Data/single_cell_meta_data_table_excel.xls']),
										  header = None)

		## Creates a dictionary with sample id key, data file location, and relevant metadata
		annotation_dict = dict()

		annotation_dict = annotation_df.set_index(0).T.dropna().to_dict('list')

		# Below are the old codes for tsv files. Deprecated
		# for line in open(''.join([self.storage_mount_point,'01_RNAseq_RAW_Data/single_cell_meta_data_table.tsv']),'r'):
		# 	elem = str.split(line.rstrip())
		# 	if elem[0] not in annotation_dict:
		# 		annotation_dict[elem[0]] = elem[1:]

		self.annotation_dict = annotation_dict

		## Read the raw Cellranger filtered data matrices into new Anndata objects
		if Path(raw_data_file).is_file():
			print(''.join(['Data_','_'.join(self.sample_list),'.scanpy.raw.h5ad']),'found, using this existing raw data file\n')
			adata = sc.read_h5ad(raw_data_file)
		else:
			print('\nNo existing h5ad raw data file found, reading in 10x h5 data for each sample\n')
			adata = 0
			for sample in self.sample_list:
				if adata:
					adata = adata.concatenate(self.__create_scanpy_anndata(self.storage_mount_point, sample, annotation_dict, file_type=file_type))
				else:
					adata = self.__create_scanpy_anndata(self.storage_mount_point, sample, annotation_dict, file_type=file_type)
			
			## Make cell names unique by adding _1, _2, _3 sequentially to each duplicated 10x barcode/name
			adata.obs_names_make_unique()
			
			## Write the raw combined dataset to disk so you can skip combining samples next time
			try:
				print('\nSaving raw combined sample data to', raw_data_file,'\n')
				adata.write(raw_data_file)
			except:
				print('\nUnable to save raw combined sample data to', raw_data_file,'\n')

		self.initial_cell_count = len(adata.obs_names)
		self.initial_gene_count = len(adata.var_names)

		## Remove specific genes from the analysis (such as experimentally observed contaminations)
		if self.remove_genes:
			print('Removing specified list of genes from analysis')
			adata = self.filter_specific_genes(adata, text_file = remove_genes)

		return adata

		## Adds metadata from a dictionary to AnnData observations
	def __create_scanpy_anndata(self,
								storage_mount_point,
								sampleID,
								annotation_dict,
								file_type=None):
		'''
		In:
		storage_mount_point: Data storage mount location
		sampleID: ID numbers of samples from the metadata table (ex: 2235-1)
		annotation_dict: Dictionary of all the sample IDs and metadata

		Out:
		New filled AnnData object
		'''
		metadata_list = annotation_dict[sampleID][1:]
		file_path = ''.join([storage_mount_point, annotation_dict[sampleID][0]])
		if file_type == '.h5':
			newAdata = sc.read_10x_h5(file_path)#,genome='GRCh38')# genome='hg19' or genome='GRCh38'
		elif file_type == '.loom':
			newAdata = sc.read_loom(file_path)
		else:
			print("No file type selected")
			newAdata = sc.read(file_path)#,genome='GRCh38')# genome='hg19' or genome='GRCh38'


		## Set gene names to be unique since there seem to be duplicate names from Cellranger
		newAdata.var_names_make_unique()

		## Add metadata for each sample to the observation (cells) annotations in the Anndata objects
		print('\nAdding Metadata to individual samples.\n')
		for field in metadata_list:
			field_list = str.split(field,':')
			meta_name = field_list[0]
			meta_value = field_list[1]
			newAdata.obs[meta_name] = meta_value
		return newAdata.copy()


	## Remove a specified list of genes from AnnData object
	def filter_specific_genes(self, 
							  adata,
							  text_file=None):
		'''
		In:
		adata: AnnData object
		text_file: List of undesired genes that can be in either a line separated text file or a Python list
		
		Out: 
		New filtered AnnData object

		Useful for removing unnecessary genes from the analysis such as blood genes
		'''
		gene_list=[]
		if text_file:
			for line in open(text_file,'r'):
				gene_list.append(line.rstrip('\n'))
		
		return adata[:, [(gene not in gene_list) for gene in adata.var_names]].copy()