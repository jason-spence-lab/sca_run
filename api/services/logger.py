
class logger(object):
	class __logger:
		def __init__(self):
			self.initial_cell_count = 0
			self.initial_gene_count = 0
			self.final_cell_count = 0
			self.final_gene_count = 0
			self.test = 0

			self.sca_params = None

		# def __str__(self):
		# 	return `self + self.initial_cell_count + self.test

		def print_test(self):
			print(self.test)

	instance = None
	def __new__(cls): # __new__ always a classmethod
		if not logger.instance:
			logger.instance = logger.__logger()
		return logger.instance

	def __getattr__(self, name):
		return getattr(self.instance, name)

	def __setattr__(self, name):
		return setattr(self.instance, name)

	# ## Count number of cells in each louvain cluster as well as the sample splits within each cluster
	# def __cell_counter(self,adata):
	# 	#print(adata.obs['louvain'])
	# 	self.sample_counts = {}
	# 	self.louvain_counts = {}
	# 	for sample in adata.obs['sampleName'].cat.categories:
	# 		try:
	# 			self.sample_counts[sample] = {'sample_total':adata.obs['sampleName'].isin([sample]).value_counts()[True]}
	# 		except:
	# 			self.sample_counts[sample] = 0
	# 		for cluster in adata.obs['louvain'].cat.categories:
	# 			if not cluster in self.louvain_counts:
	# 				try:
	# 					self.louvain_counts[cluster] = {'louvain_total':adata.obs['louvain'].isin([cluster]).value_counts()[True]}
	# 				except:
	# 					self.louvain_counts[cluster] = 0

	# 			try:
	# 				self.sample_counts[sample][cluster] = self.louvain_counts[cluster][sample] = (adata.obs['sampleName'].isin([sample]) & adata.obs['louvain'].isin([cluster])).value_counts()[True]
	# 			except:
	# 				self.sample_counts[sample][cluster] = self.louvain_counts[cluster][sample] = 0

	# 	return 0

	# ## If doing cell scoring analysis, get average score per cluster given a gene scoring list name
	# def __count_cluster_scores(self, adata, cell_score_list):
	# 	cluster_scores = {}
	# 	for cluster in adata.obs['louvain'].cat.categories:
	# 		cluster_scores[cluster] = adata[adata.obs['louvain'].isin([cluster])].obs[cell_score_list].mean()

	# 	return cluster_scores

	# ## Write a summary of the analysis run including sample information, parameters and filtering information
	# # Not completely up to date
	# def write_summary(self, figdir='./figures/'):
	# 	self.__cell_counter(self.adata)
	# 	self.final_cell_count = len(self.adata.obs_names)
	# 	self.final_gene_count=len(self.adata.var_names)

	# 	fileName = ''.join([figdir,'summary.txt'])

	# 	os.makedirs(os.path.dirname(fileName), exist_ok=True) # Create directory if it doesn't exist
	# 	with open(fileName,'w') as f:

	# 		f.write("Summary of single cell sequencing analysis for samples ")
	# 		f.write(''.join([' '.join(self.sample_list),'\n\n']))

	# 		f.write('--------Sample Metadata--------\n')
	# 		for sample in self.sample_list:
	# 			f.write(''.join(['- Sample ',sample,'\n']))
	# 			f.write('\n'.join(self.annotation_dict[sample]))
	# 			f.write('\n\n')

	# 		f.write('--------Basic Run Information--------\n')
	# 		f.write(''.join(['Initial cell count:  ',str(self.initial_cell_count),'\n']))
	# 		f.write(''.join(['Final cell count:  ',str(self.final_cell_count),'\n']))
	# 		f.write(''.join(['Initial gene count:  ',str(self.initial_gene_count),'\n']))
	# 		f.write(''.join(['Final gene count:  ',str(self.final_gene_count),'\n']))

	# 		f.write('\n--------Filter Parameters Used--------\n')
	# 		f.write(''.join(['Min Cells:  ',str(self.min_cells),'\n']))
	# 		f.write(''.join(['Min Genes:  ',str(self.min_genes),'\n']))
	# 		f.write(''.join(['Max Genes:  ',str(self.max_genes),'\n']))
	# 		f.write(''.join(['Max Counts:  ',str(self.max_counts),'\n']))
	# 		f.write(''.join(['Max Mito:  ',str(self.max_mito),'\n']))

	# 		f.write('\n--------Analysis Parameters Used--------\n')
	# 		f.write(''.join(['# Neighbors:  ',str(self.n_neighbors),'\n']))
	# 		f.write(''.join(['# PCs:  ',str(self.n_pcs),'\n']))
	# 		f.write(''.join(['Spread:  ',str(self.spread),'\n']))
	# 		f.write(''.join(['Min Dist:  ',str(self.min_dist),'\n']))
	# 		f.write(''.join(['Resolution:  ',str(self.resolution),'\n']))

	# 		f.write('\n--------Sample Cell Counts Used--------\n')
	# 		f.write(json.dumps(str(self.sample_counts)))
	# 		f.write('\n--------Louvain Cell Counts Used--------\n')
	# 		f.write(json.dumps(str(self.louvain_counts)))

	# 		if self.cell_score_lists:
	# 			f.write('\n--------Cluster Scores--------\n')
	# 			for cell_score_list in self.cell_score_lists:
	# 				# f.write(''.join(['--------',cell_score_list,'_raw--------\n']))
	# 				# f.write(json.dumps(str(self.__count_cluster_scores(self.adata, ''.join([cell_score_list,'_raw'])))))
	# 				f.write(''.join(['--------',cell_score_list,'--------\n']))
	# 				f.write(json.dumps(str(self.__count_cluster_scores(self.adata, cell_score_list))))
	# 				f.write('\n')
	# 				# f.write(''.join(['\n--------',cell_score_list,'_processed--------\n']))
	# 				# f.write(json.dumps(str(self.__count_cluster_scores(self.adata, ''.join([cell_score_list,'_processed'])))))
	# 	return 0