'''
helper -- 
Miscellaneous helper functions

Written by Joshua H Wu
09 March, 2021
'''

class helper:

	## Separates an adata object into a list of adata objects split using an obs field
	def separate_adata(self, adata, obs_field):
		adata_list = []
		for meta_value in adata.obs[obs_field].unique():
			adata_list.append(adata[(adata.obs[obs_field] == meta_value),:].copy())

		return adata_list