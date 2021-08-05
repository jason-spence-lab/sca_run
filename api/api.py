import time
from flask import Flask
from flask_cors import CORS
from classes.sca_params import *
import SCARunner as SCARunner
from flask import request
app = Flask(__name__)
CORS(app)
from datetime import date

@app.route('/submit', methods=['POST'])
def analysis_script():
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

	figdir = './figures/'
	an_params = sca_params()
	#############################################################################
	## Change this to point toward your mount location for our MiStorage share ##
	#############################################################################
	an_params.storage_mount_point = 'Z:/'
	an_params.species = request.json.get('species',None)
	if an_params.species is not ('human' or 'mouse'):
		an_params.species = 'human'

	user_name = request.json.get('user_name', None)
	if not user_name:
		print("No username - who are you?")

	## List of interesting genes
	gene_filestr = request.json.get('genes',None)
	print(gene_filestr)
	an_params.read_gene_str(gene_filestr)

	## Parameters used to filter the data - Mainly used to get rid of bad cells
	an_params.set_qc_params(min_cells = 0,
						    min_genes = int(request.json.get('min_genes',None)),
						    max_genes = int(request.json.get('max_genes',None)),
						    max_counts = int(request.json.get('max_counts',None)),
						    max_mito = int(request.json.get('max_mito',None)),
						    doublet_detection = request.json.get('doublet_detection',None))

	regress_out = request.json.get('regress_out', None).replace(" ","").split(sep=",")
	batch_check = request.json.get('batch_check', None)
	if batch_check:
		batch_algo = request.json.get('batch_algo', None)

		combat_boolean = batch_algo.casefold() == 'combat'
		do_bbknn_boolean = batch_algo.casefold() == 'bbknn'
		batch_cat = request.json.get('batch_cat', None)

		if combat_boolean:
			combat = batch_cat
			do_bbknn = None
		elif do_bbknn_boolean:
			combat = None
			do_bbknn = batch_cat
	else:
		combat = None
		do_bbknn = None

	highly_variable = request.json.get('highly_variable', None)
	if highly_variable.casefold() == 'default':
		# print('highly_variable default running')
		an_params.set_pp_params(combat = combat,
								min_mean=0.0125,
								max_mean=3,
								min_disp=0.5,
								n_top_genes=None,
								regress_out=regress_out)
	else:
		highly_variable = highly_variable.replace(" ","").split(sep=",")
		if len(highly_variable)==1:
			an_params.set_pp_params(combat=combat,
									n_top_genes=int(highly_variable[0]),
									regress_out=regress_out)
		elif len(highly_variable)==3:
			an_params.set_pp_params(combat=combat,
									min_mean=float(highly_variable[0]),
									max_mean=float(highly_variable[1]),
									min_disp=float(highly_variable[2]),
									n_top_genes=None,
									regress_out=regress_out)

	dpt=[]
	if request.json.get('dpt_check'):
		dpt = [request.json.get('root_cat'),request.json.get('root_group', None).replace(" ","").split(sep=",")]
		print(dpt)


	## Parameters used for initial clustering analysis
	an_params.set_analysis_params(n_neighbors=int(request.json.get('n_neighbors',None)),
							      n_pcs=int(request.json.get('n_pcs',None)),
							      umap_check=request.json.get('umap_check',None),
							      spread=float(request.json.get('spread',None)),
							      min_dist=float(request.json.get('min_dist',None)),
							      resolution=float(request.json.get('resolution',None)),
							      do_bbknn=do_bbknn,
							      do_tSNE=request.json.get('tsne_check',None),
							      clustering_choice='leiden',
							      dpt=dpt,
							      draw_force_atlas=request.json.get('graph_fa_check',None),
							      umap_init_pos=request.json.get('umap_init_pos',None),
							      phate=request.json.get('phate_check',None))
							      
	umap_feature_color = request.json.get('color_map', None)
	if not umap_feature_color == ('yellow_blue' or 'blue_orange'):
		umap_feature_color = eval(umap_feature_color)

	umap_categorical_color = request.json.get('cat_color',None)
	if not umap_categorical_color.casefold() == 'default':
		umap_categorical_color = umap_categorical_color.replace(" ","").split(sep=",")

	if not request.json.get('diff_exp_check',None):
		rank_grouping=request.json.get('comparison_cats', None).replace(" ","").split(sep=",")
		clusters2_compare=request.json.get('groups_compare', None).replace(" ","").split(sep=",")
	else:
		rank_grouping = None
		clusters2_compare = None

	an_params.set_plot_params(size=int(request.json.get('dot_size',None)),
							  umap_obs=request.json.get('graph_cat', None).replace(" ","").split(sep=","),
							  dot_check=request.json.get('dot_check',None),
							  heatmap_check=request.json.get('heat_check',None),
							  dot_grouping= request.json.get('gene_exp_cat', None).replace(" ","").split(sep=","),
							  umap_categorical_color=umap_categorical_color,
							  umap_feature_color=umap_feature_color,
							  final_quality=request.json.get('figure_quality',None),
							  vmin_list=[],
							  vmax_list=[],
							  rank_grouping=rank_grouping,
							  clusters2_compare=clusters2_compare)

	## Basic pipeline for analysis - will filter data, process, cluster, etc. and output relevant figures
	an_run = SCARunner.SCARunner()

	analysis_label = request.json.get('analysis_label',None)
	if not analysis_label or analysis_label.casefold()==('current date' or 'date'):
		analysis_label = date.today().strftime("%m%d%Y")

	extracted = request.json.get('clusters_extract',None).replace(" ","").split(sep=",")
	sample = request.json.get('sample_list',None)
	figdir = ''.join([an_params.storage_mount_point,'single_cell_analyses/',user_name,'/',analysis_label,'/'])

	print(extracted)
	if ('/' in sample) and extracted:
		an_run.pipe_ext(an_params,figdir=figdir,load_save=sample,extracted=extracted)
	elif ('/' in sample) and not extracted:
		an_run.pipe_basic(an_params,figdir=figdir,load_save=sample)
	else:
		an_params.sample_list = sample.replace(" ","").split(sep=",")
		an_run.pipe_basic(an_params,figdir=figdir)

@app.route('/upload', methods=['POST'])
def upload():
	data_test = request.json.get('data_test', None)
	print(data_test)
	return

