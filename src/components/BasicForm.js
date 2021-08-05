import React from 'react';
import Grid from '@material-ui/core/Grid';
import { withStyles } from '@material-ui/core/styles';
import Container from '@material-ui/core/Container';
import Paper from '@material-ui/core/Paper';
import TextField from '@material-ui/core/TextField';
import Typography from '@material-ui/core/Typography';
import ExpansionPanel from '@material-ui/core/ExpansionPanel';
import Card from '@material-ui/core/Card';
import ExpansionPanelSummary from '@material-ui/core/ExpansionPanelSummary';
import ExpansionPanelDetails from '@material-ui/core/ExpansionPanelDetails';
import FormControlLabel from '@material-ui/core/FormControlLabel';
import Checkbox from '@material-ui/core/Checkbox';
import Input from '@material-ui/core/Input';
import OutlinedInput from '@material-ui/core/OutlinedInput';
import FilledInput from '@material-ui/core/FilledInput';
import InputLabel from '@material-ui/core/InputLabel';
import MenuItem from '@material-ui/core/MenuItem';
import FormHelperText from '@material-ui/core/FormHelperText';
import FormControl from '@material-ui/core/FormControl';
import Select from '@material-ui/core/Select';
import Info from '@material-ui/icons/Info';
import IconButton from '@material-ui/core/IconButton';
import Dialog from '@material-ui/core/Dialog';
import DialogTitle from '@material-ui/core/DialogTitle';
import DialogContent from '@material-ui/core/DialogContent';
import DialogActions from '@material-ui/core/DialogActions';
import Button from '@material-ui/core/Button';
import ButtonBase from '@material-ui/core/ButtonBase';
import DropZone from './DropZone';

const styles = theme => ({
    paper: {
	    margin: theme.spacing(0.5, 3),
	    display: 'stretch',
	    flexDirection: 'column',
	    alignItems: 'center',
  	},
  	root: {
    	height: '80vh',
    	background: theme.palette.common.white,
  	},
    form: {
	    width: '100%', // Fix IE 11 issue.
	    // marginTop: theme.spacing(0),
  	},
  	leftPane: {
  		borderRight: '1px solid lightgray',
  		marginBottom: '5px',
  	},
  	rightPane:{
  		borderLeft: '1px solid lightgray',
  		marginBottom: '5px',
  		overflowY: 'scroll',
  	},
  	expansionGrid: {
  	},
 //  	formControl: {
	//     // margin: theme.spacing(0),
	//     minWidth: 120,
	// },
	submit: {
        margin: theme.spacing(0,0,0),
        background: 'linear-gradient(45deg, #00274C 30%, #00274C 90%)',
        color: "#FFCB05",
        textTransform: "None",
        '&:hover': {
	      	background: '#FFCB05',
	      	color: "#00274C",
	      	// borderColor: '#0062cc',
	      	// boxShadow: 'none',
	    },
    },
    infoTitle: {
    	padding: '8px 24px',
    	marginTop: '4px'
    },
    infoContent: {
    	marginTop: '-6px'
    },
    input: {
    	// display: 'none',
  	},
  	uploadButton: {
  		backgroundColor: '#00274C',
  		textTransform: "None",
  		color: "#FFCB05",
  		'&:hover': {
	      	backgroundColor: '#FFCB05',
	      	color: "#00274C",
	      	// borderColor: '#0062cc',
	      	// boxShadow: 'none',
	    },
  	},
  	button: {
  		color: '#00274C',
  		marginTop: "-5px",
  		fill: '#FFCB05',
  		padding: '9px',
  	},
  	formField: {
  		margin: "0",
  		marginTop:"0px",
  		marginBottom: "10px",
  	},
  	previewChip: {
  		minWidth:160,
  		maxWidth:210,
  	},
});

const basicFields = [
	{
		id:'name',
		label:'Michigan Unique ID',
		value:'',
		placeholder:'wujos',
		width:4,
	},
	{
		id:'label',
		label:'Analysis Label',
		value:'',
		placeholder:'Date',
		width:4,
	},
	{
		id:'species',
		label:'Sample Species',
		value:'human',
		placeholder:'human',
		width:4,
	},
	{
		id:'sampleIDs',
		label:'Sample IDs',
		value:'',
		placeholder:'Ex: 2444-1 or intmed-spence-lab/single_cell_analyses/Josh/figures_07072021/adata_save.p',
		width:12,
	},
];

const qcFields = [
	{
		id:'minGenes',
		label:'Min Genes',
		defaultValue:'500',
		helperText:'Remove cells expressing fewer genes'
	},
	{
		id:'maxGenes',
		label:'Max Genes',
		defaultValue:'7500',
		helperText:'Remove cells expressing more genes'
	},
	{
		id:'maxCounts',
		label:'Max Counts',
		defaultValue:'30000',
		helperText:'Remove cells with more unique molecular identifiers'
	},
	{
		id:'maxMito',
		label:'Max Mito',
		defaultValue:'10',
		helperText:'Remove cells with higher mitochondrail gene expression percentage'
	}
]

const preprocessFields = [
	{
		id:'highly_variable',
		label:'Highly Variable Genes',
		defaultValue:'default',
		helperText:'Number of highly variable genes to extract, if default will automatically calculate top genes'
	},{
		id:'regress_out',
		label:'Partial Variable Regression',
		defaultValue:'total_counts, pct_counts_mito',
		helperText:'Adjust for variables selected using partial linear regression'
	}
]

const analysisFields = [
	{
		id:'n_neighbors',
		label:'N-Neighbors',
		defaultValue:'15',
		helperText:'Size of local neighborhood'
	},
	{
		id:'n_pcs',
		label:'PCs',
		defaultValue:'11',
		helperText:'# of principle components used in construction of neighborhood graph'
	},{
		id:'resolution',
		label:'Resolution',
		defaultValue:'0.4',
		helperText:'Number of clusters identified'
	}
]

const plotFields = [
	{
		id: 'graph_dot_size',
		label:'Graph Dot Size',
		defaultValue:'15',
		helperText:'Size of the dots on the UMAP, t-SNE, Diffusion Map, etc. plots'
	},{
		id: 'categorical_color',
		label:'Categorical Dot Color',
		defaultValue:'default',
		helperText:'Color palette for categorical data to be plotted on graph embeddings'
	},{
		id: 'color_map',
		label:'Color Map',
		defaultValue: 'yellow_blue',
		helperText:'Color map for continuous data to be plotted on graph embeddings, dot plots, etc.'
	},{
		id: 'figure_quality',
		label:'Figure Quality',
		defaultValue:'low',
		helperText:'Resolution of plotted figures (higher quality will result in slower plotting)'
	},{
		id: 'graph_categories',
		label:'Graph Categories',
		defaultValue:'sampleName, leiden',
		helperText:'Categories to be plotted on graph embeddings'
	}
]

const umapFields = [
	{
		id:'spread',
		label:'Spread',
		defaultValue:'1',
		helperText:'How clumped embedded points are'
	},
	{
		id:'min_dist',
		label:'Min Distance',
		defaultValue:'0.4',
		helperText:'Minimum distance between points on the cluster graph'
	},
]

const diffExpFields = [
	{
		id:'diffExpCat',
		label:'Comparison Categories',
		defaultValue:'leiden',
		helperText:'Categories for differential gene expression comparison',
	},
	{
		id:'diffExpGroups2Compare',
		label:'Groups to Compare',
		defaultValue:'All',
		helperText:'Groups within categories to compare'
	},
]

const dptFields = [
	{
		id:'dptCatRoot',
		label:'Category for Root Cell',
		defaultValue:'',
		helperText:'Category for identifying root cell type',
	},
	{
		id:'dptGroupRoot',
		label:'Root Cell Group',
		defaultValue:'',
		helperText:'Group within category containing root cell type'
	},
]

class BasicForm extends React.Component {
	constructor(props) {
		super(props);
		this.handleBasicFieldChange = this.handleBasicFieldChange.bind(this);
		this.handleClickOpen = this.handleClickOpen.bind(this);
		this.handleClose = this.handleClose.bind(this);
		this.handleBasicFieldClickOpen = this.handleBasicFieldClickOpen.bind(this)
		this.handleBasicFieldClose = this.handleBasicFieldClose.bind(this)

		this.handleQCClickOpen = this.handleQCClickOpen.bind(this)
		this.handleQCClose = this.handleQCClose.bind(this)

		this.handlePreprocessFieldChange = this.handlePreprocessFieldChange.bind(this)
		this.handlePreprocessFieldClickOpen = this.handlePreprocessFieldClickOpen.bind(this)
		this.handlePreprocessFieldClose = this.handlePreprocessFieldClose.bind(this)

		this.handleAnalysisFieldClickOpen = this.handleAnalysisFieldClickOpen.bind(this)
		this.handleAnalysisFieldClose = this.handleAnalysisFieldClose.bind(this)

		this.handleMiscFieldClickOpen = this.handleMiscFieldClickOpen.bind(this)
		this.handleMiscFieldClose = this.handleMiscFieldClose.bind(this)

		this.handlePlotFieldChange = this.handlePlotFieldChange.bind(this)
		this.handlePlotFieldClickOpen = this.handlePlotFieldClickOpen.bind(this)
		this.handlePlotFieldClose = this.handlePlotFieldClose.bind(this)

		this.handlePlotFieldChange = this.handlePlotFieldChange.bind(this)
		this.handlePlotFieldClickOpen = this.handlePlotFieldClickOpen.bind(this)
		this.handlePlotFieldClose = this.handlePlotFieldClose.bind(this)

		this.handleBatchCorrClickOpen = this.handleBatchCorrClickOpen.bind(this)
		this.handleBatchCorrClose = this.handleBatchCorrClose.bind(this)

		this.handleUmapFieldChange = this.handleUmapFieldChange.bind(this)
		this.handleUmapFieldClickOpen = this.handleUmapFieldClickOpen.bind(this)
		this.handleUmapFieldClose = this.handleUmapFieldClose.bind(this)

		this.handleDiffExpFieldChange = this.handleDiffExpFieldChange.bind(this)
		this.handleDiffExpFieldClickOpen = this.handleDiffExpFieldClickOpen.bind(this)
		this.handleDiffExpFieldClose = this.handleDiffExpFieldClose.bind(this)

		this.handleGeneExpFieldClickOpen = this.handleGeneExpFieldClickOpen.bind(this)
		this.handleGeneExpFieldClose = this.handleGeneExpFieldClose.bind(this)

		this.handleDPTFieldChange = this.handleDPTFieldChange.bind(this)
		this.handleDPTFieldClickOpen = this.handleDPTFieldClickOpen.bind(this)
		this.handleDPTFieldClose = this.handleDPTFieldClose.bind(this)

		this.handleExtractClickOpen = this.handleExtractClickOpen.bind(this)
		this.handleExtractClose = this.handleExtractClose.bind(this)

		this.handleSubmit = this.handleSubmit.bind(this)
		this.handleUpload = this.handleUpload.bind(this)

		this.handleSetFiles = this.handleSetFiles.bind(this)

		this.handleNoUserNameClose = this.handleNoUserNameClose.bind(this)

		this.state = {
			basicFieldValues:['','','human',''],
			qcFieldValues:['500','7500','30000','10'],
			preprocessFieldValues:['default','total_counts, pct_counts_mito'],
			analysisFieldValues:['15','11','0.4'],
			plotFieldValues:['15','default','yellow_blue','low','sampleName,leiden'],
			umapFieldValues:['1','0.4'],
			diffExpFieldValues:['leiden','all'],
			dptFieldValues:['',''],
			extractClusters:'',
			umapPos:'',
			batchAlgo:'',
			batchCat:'sampleName',
			umapCheck:false,
			tsneCheck:false,
			diffExpCheck:false,
			doubletCheck:false,
			batchCheck:false,
			dptCheck:false,
			phateCheck:false,
			graphFACheck:false,
			extractCheck:false,
			geneExpCheck:false,
			dotCheck:false,
			heatCheck:false,
			violinCheck:false,
			geneExpCat:'leiden, sampleName',
			open:false,
			basicFieldOpen:false,
			qcopen:false,
			ppopen:false,
			anopen:false,
			miscopen:false,
			plotopen:false,
			batchopen:false,
			geneexpopen:false,
			extractopen:false,
			diffexpopen:false,
			submit:false,
			file:"",
			noUserName:false,
		}
	}

	handleBasicFieldChange(value,i) {
		this.setState(prevState => {
			const tempVals = prevState.basicFieldValues;
			tempVals[i] = value;
			console.log(value)
			return{
				basicFieldValues: tempVals
			}
		})
	};

	handleQCFieldChange(value,i) {
		this.setState(prevState => {
			const tempVals = prevState.qcFieldValues;
			tempVals[i] = value;
			console.log(value)
			return{
				qcFieldValues: tempVals
			}
		})
	};

	handleSetFiles(fileString) {
		this.setState({file: fileString})
		console.log(this.state.file)
	}

	handlePreprocessFieldChange(value,i) {
		this.setState(prevState => {
			const tempVals = prevState.preprocessFieldValues;
			tempVals[i] = value;
			console.log(value)
			return{
				preprocessFieldValues: tempVals
			}
		})
	};

	handleAnalysisFieldChange(value,i) {
		this.setState(prevState => {
			const tempVals = prevState.analysisFieldValues;
			tempVals[i] = value;
			console.log(value)
			return{
				analysisFieldValues: tempVals
			}
		})
	};

	handlePlotFieldChange(value,i) {
		this.setState(prevState => {
			const tempVals = prevState.plotFieldValues;
			tempVals[i] = value;
			console.log(value)
			return{
				plotFieldValues: tempVals
			}
		})
	};

	handleUmapFieldChange(value,i) {
		this.setState(prevState => {
			const tempVals = prevState.umapFieldValues;
			tempVals[i] = value;
			console.log(value)
			return{
				umapFieldValues: tempVals
			}
		})
	};

	handleDiffExpFieldChange(value,i) {
		this.setState(prevState => {
			const tempVals = prevState.diffExpFieldValues;
			tempVals[i] = value;
			console.log(value)
			return{
				diffExpFieldValues: tempVals
			}
		})
	};

	handleDPTFieldChange(value,i) {
		this.setState(prevState => {
			const tempVals = prevState.dptFieldValues;
			tempVals[i] = value;
			console.log(value)
			return{
				dptFieldValues: tempVals
			}
		})
	};

	handleChange(e,name) {
		this.setState({[name]:e.target.checked});
	};

	handleClickOpen(e) {
	    this.setState({open:true});
	  }

	handleClose(e) {
	    this.setState({open:false});
	  };

	handleMiscFieldClickOpen(e) {
	    this.setState({miscopen:true});
	  }

	handleMiscFieldClose(e) {
	    this.setState({miscopen:false});
	  };

	handleBasicFieldClickOpen(e) {
	    this.setState({basicFieldOpen:true});
	 }

	handleBasicFieldClose(e) {
	    this.setState({basicFieldOpen:false});
	  };

	handleQCClickOpen(e) {
	    this.setState({qcopen:true});
	}

	handleQCClose(e) {
	    this.setState({qcopen:false});
	};

	handlePreprocessFieldClickOpen(e) {
	    this.setState({ppopen:true});
	}

	handlePreprocessFieldClose(e) {
	    this.setState({ppopen:false});
	};

	handleAnalysisFieldClickOpen(e) {
	    this.setState({anopen:true});
	}

	handleAnalysisFieldClose(e) {
	    this.setState({anopen:false});
	};

	handlePlotFieldClickOpen(e) {
	    this.setState({plotopen:true});
	}

	handlePlotFieldClose(e) {
	    this.setState({plotopen:false});
	};

	handleBatchCorrClickOpen(e) {
	    this.setState({batchopen:true});
	}

	handleBatchCorrClose(e) {
	    this.setState({batchopen:false});
	};

	handleUmapFieldClickOpen(e) {
	    this.setState({umapopen:true});
	}

	handleUmapFieldClose(e) {
	    this.setState({umapopen:false});
	};

	handleGeneExpFieldClickOpen(e) {
	    this.setState({geneexpopen:true});
	}

	handleGeneExpFieldClose(e) {
	    this.setState({geneexpopen:false});
	};

	handleExtractClickOpen(e) {
	    this.setState({extractopen:true});
	}

	handleExtractClose(e) {
	    this.setState({extractopen:false});
	};

	handleDiffExpFieldClickOpen(e) {
	    this.setState({diffexpopen:true});
	}

	handleDiffExpFieldClose(e) {
	    this.setState({diffexpopen:false});
	};

	handleDPTFieldClickOpen(e) {
	    this.setState({dptopen:true});
	}

	handleDPTFieldClose(e) {
	    this.setState({dptopen:false});
	};

	handleNoUserNameClose(e) {
		this.setState({noUserName:false})
	}

	handleSubmit(e) {
		if (this.state.basicFieldValues[0]=='') {
			this.setState({noUserName:true})
		} else {
			const data = {
				user_name: this.state.basicFieldValues[0],
				analysis_label: this.state.basicFieldValues[1],
				species: this.state.basicFieldValues[2],
				sample_list: this.state.basicFieldValues[3],
				min_genes: this.state.qcFieldValues[0],
				max_genes: this.state.qcFieldValues[1],
				max_counts: this.state.qcFieldValues[2],
				max_mito: this.state.qcFieldValues[3],
				highly_variable: this.state.preprocessFieldValues[0],
				regress_out: this.state.preprocessFieldValues[1],
				n_neighbors: this.state.analysisFieldValues[0],
				n_pcs: this.state.analysisFieldValues[1],
				resolution: this.state.analysisFieldValues[2],
				doublet_detection: this.state.doubletCheck,
				dot_size: this.state.plotFieldValues[0],
				cat_color: this.state.plotFieldValues[1],
				color_map: this.state.plotFieldValues[2],
				figure_quality: this.state.plotFieldValues[3],
				umap_check:this.state.umapCheck,
				graph_cat: this.state.plotFieldValues[4],
				umap_init_pos:this.state.umapPos,
				spread: this.state.umapFieldValues[0],
				min_dist: this.state.umapFieldValues[1],
				tsne_check: this.state.tsneCheck,
				// gene_exp_check: this.state.geneExpCheck,
				gene_exp_cat: this.state.geneExpCat,
				dot_check: this.state.dotCheck,
				heat_check: this.state.heatCheck,
				violin_check: this.state.violinCheck,
				diff_exp_check: this.state.diffExpCheck,
				comparison_cats: this.state.diffExpFieldValues[0],
				groups_compare: this.state.diffExpFieldValues[1],
				batch_check: this.state.batchCheck,
				batch_algo: this.state.batchAlgo,
				batch_cat: this.state.batchCat,
				dpt_check: this.state.dptCheck,
				root_cat:this.state.dptFieldValues[0],
				root_group:this.state.dptFieldValues[1],
				phate_check: this.state.phateCheck,
				graph_fa_check: this.state.graphFACheck,
				clusters_extract: this.state.extractClusters,
				genes: this.state.file,
			}

			// call to api
			this.setState({submit:true});

			fetch('http://localhost:5000/submit', {
			  method: 'POST', // or 'PUT'
			  headers: {
			    'Content-Type': 'application/json',
			  },
			  body: JSON.stringify(data), 
			})
		}
	}

	handleUpload(e) {
		const data =  {data_test: 0}
		fetch('http://localhost:5000/upload', {
			method: 'POST',
			headers: {
				'Content-Type': 'application/json',
			},
			body: JSON.stringify(data)
		})

	}

	render() {
		const {classes} = this.props;

		return(
			<Grid container component="main" alignItems='stretch' justify="center" spacing={1}>
			  	<Grid item xs={5} sm={5} md={5}  component={Paper} className={classes.leftPane}>
			  		<Container>
			  			{this.state.noUserName && <Dialog onClose={this.handleNoUserNameClose} aria-labelledby="simple-dialog-title" open={this.state.noUserName}>
					        <DialogTitle id="simple-dialog-title" className={classes.infoTitle}>Please provide your UMich unique ID</DialogTitle>
					    </Dialog>}
			  			<Typography variant="h7" color="inherit" noWrap>
			      			Basic Information
			    		</Typography>
			    		<IconButton className={classes.button} aria-label="info" onClick={this.handleBasicFieldClickOpen}>
			    			<Info fontSize="small" className={classes.iconInfo}/>
			    		</IconButton>
			    		<Dialog onClose={this.handleBasicFieldClose} aria-labelledby="simple-dialog-title" open={this.state.basicFieldOpen}>
					        <DialogTitle id="simple-dialog-title" className={classes.infoTitle}>Basic Fields</DialogTitle>
					        <DialogContent className={classes.infoContent}>
						        <b>Sample IDs</b> - Metadata file run number of the samples to be analyzed
						        <br />
						        <b>Gene List</b> - List of genes to plot in UMAP feature plots, dot plots, etc.
						        <br />
					        </DialogContent>
					    </Dialog>
					    <Grid container spacing={1}>
							{basicFields.map((field,i) => (
								<Grid item xs={field.width}>
									<TextField
								        id={field.id}
								        label={field.label}
								        type="search"
								        value={this.state.basicFieldValues[i]}
								        onChange={(e)=>{
								        		this.handleBasicFieldChange(e.target.value,i);
								        		//console.log(e.target.value)
								        	}
								        }
								        className={classes.formField}
								        placeholder={field.placeholder}
								        fullWidth
								        margin="normal"
								        InputLabelProps={{
								        	shrink: true,
								        }}
									/>
								</Grid>))}
						</Grid>
					    <DropZone handleSetFiles={this.handleSetFiles}/>
			  			<Typography variant="h7" color="inherit">
			      			Quality Control Parameters
			    		</Typography>
			    		<IconButton className={classes.button} aria-label="info" onClick={this.handleQCClickOpen}>
			    			<Info fontSize="small"/>
			    		</IconButton>
			    		<Dialog onClose={this.handleQCClose} aria-labelledby="simple-dialog-title" open={this.state.qcopen}>
					        <DialogTitle id="simple-dialog-title" className={classes.infoTitle}>Quality Control Parameters</DialogTitle>
					        <DialogContent className={classes.infoContent}>
								<b>Min genes</b> - Filter out cells with fewer genes to remove dead cells
								<br />
								<b>Max genes</b> - Filter out cells with more genes to remove most doublets
								<br />
								<b>Max counts</b> - Filter out cells with more UMIs to catch a few remaining doublets
								<br />
								<b>Max mito</b> - Filter out cells with high mitochondrial gene transcript fraction
					        </DialogContent>
					    </Dialog>

						<Grid container spacing={1}>
			    			{qcFields.map((field,i) =>(
			    				<Grid item xs={6} sm={3}>
				    				<TextField
								        id={field.id}
								        label={field.label}
								        className={classes.formField}
								        value={this.state.qcFieldValues[i]}
								        onChange={(e)=>{this.handleQCFieldChange(e.target.value,i)}}
								        defaultValue={field.defaultValue}
								        //helperText={field.helperText}
								        fullWidth
								        margin="normal"
								        InputLabelProps={{
								        	shrink: true,
								        }}
								    />
							    </Grid>
							))}
						</Grid>

						<Typography variant="h7" color="inherit">
			      			Preprocess Parameters
			    		</Typography>
			    		<IconButton className={classes.button} aria-label="info" onClick={this.handlePreprocessFieldClickOpen}>
			    			<Info fontSize="small"/>
			    		</IconButton>
			    		<Dialog onClose={this.handlePreprocessFieldClose} aria-labelledby="simple-dialog-title" open={this.state.ppopen}>
					        <DialogTitle id="simple-dialog-title" className={classes.infoTitle}>Preprocessing Parameters</DialogTitle>
					        <DialogContent className={classes.infoContent}>
						        <b>Highly Variable Genes </b> - Extract highly variable genes for downstream analysis (defaults to self-calculation)
						        <br/>
								<b>Partial Variable Regression </b> - Regress out effects of certain variables such as read depth and mitochondrial transcript fraction
								<br />
					        </DialogContent>
					    </Dialog>
			    		<Grid container spacing={1}>
			    			{preprocessFields.map((field,i) =>(
			    				<Grid item xs={6} sm={6}>
				    				<TextField
								        id={field.id}
								        label={field.label}
								        className={classes.formField}
								        value={this.state.preprocessFieldValues[i]}
								        onChange={(e)=>{this.handlePreprocessFieldChange(e.target.value,i)}}
								        defaultValue={field.defaultValue}
								        //helperText={field.helperText}
								        fullWidth
								        margin="normal"
								        InputLabelProps={{
								        	shrink: true,
								        }}
								    />
							    </Grid>
							))}
						</Grid>

						<Typography variant="h7" color="inherit">
			      			Analysis Parameters
			    		</Typography>
			    		<IconButton className={classes.button} aria-label="info" onClick={this.handleAnalysisFieldClickOpen}>
			    			<Info fontSize="small"/>
			    		</IconButton>
			    		<Dialog onClose={this.handleAnalysisFieldClose} aria-labelledby="simple-dialog-title" open={this.state.anopen}>
					        <DialogTitle id="simple-dialog-title" className={classes.infoTitle}>Analysis Parameters</DialogTitle>
					        <DialogContent className={classes.infoContent}>
						        <b>Nearest Neighbors </b> - Size of the local neighborhood used for manifold approximation
						        <br/>
								<b>Principle Components </b> - Number of principle components to use in construction of neighborhood graph
								<br />
								<b>Resolution </b> - High resolution attempts to increases # of clusters identified
					        </DialogContent>
					    </Dialog>
			    		<Grid container spacing={1}>
			    			{analysisFields.map((field,i) =>(
			    				<Grid item xs={6} sm={4}>
				    				<TextField
								        id={field.id}
								        label={field.label}
								        className={classes.formField}
								        value={this.state.analysisFieldValues[i]}
								        onChange={(e)=>{this.handleAnalysisFieldChange(e.target.value,i)}}
								        defaultValue={field.defaultValue}
								        //helperText={field.helperText}
								        fullWidth
								        margin="normal"
								        InputLabelProps={{
								        	shrink: true,
								        }}
								    />
							    </Grid>
							))}
						</Grid>

						<Typography variant="h7" color="inherit" noWrap>
				      		Extra Analysis Options
				    	</Typography>
				    	<IconButton className={classes.button} aria-label="info" onClick={this.handleMiscFieldClickOpen}>
			    			<Info fontSize="small"/>
			    		</IconButton>
			    		<Dialog onClose={this.handleMiscFieldClose} aria-labelledby="simple-dialog-title" open={this.state.miscopen}>
					        <DialogTitle id="simple-dialog-title" className={classes.infoTitle}>Extra Analysis Options</DialogTitle>
					        <DialogContent className={classes.infoContent}>
					        	Choose the types of analyses and tools you would like to be run
					        	<br/>
						        <b>Doublet Detection </b> - Runs Doublet Detection by Jonathan Shor to predict and remove doublets in the quality control phase
						        <br/>
								<b>UMAP </b> - Uniform Manifold Approximation and Projection to create a graph embedding
								<br />
								<b>t-SNE </b> - t-Distributed Stochastic Neighbor Embedding to create a graph embedding
								<br />
								<b>Differential Gene Expression Analysis </b> - Calculate differential gene expression for selected categories
								<br />
								<b>Batch Correction </b> - Correct for variance among batches using combat or BBKNN
								<br />
								<b>Diffusion Pseudotime </b> - Calculate diffusion pseudotime and diffusion maps to estimate trajectory inference/lineage tracing
								<br />
								<b>Graph Force Atlas </b> - Calculate force atlas graph embedding
								<br />
								<b>PHATE </b> - Run Potential of Heat-diffusion for Affinity-based Trajectory Embedding (PHATE) by the Krishnaswamy Lab
								<br />
								<b>Cluster Extract </b> - Extract clusters from a previous analysis. Will ignore filtering parameters. Must have indicated 
								a file path. All indicated analyses will be applied on the extracted data.
					        </DialogContent>
					    </Dialog>
				    	<Grid container spacing={1}>
				    		<Grid item xs={6} sm={4}>
					    		<FormControlLabel
									control={<Checkbox 
												color="secondary" 
												checked={this.state.doubletCheck} 
												onChange={e=>this.setState({doubletCheck:e.target.checked})}
											/>}
									label={<Typography variant="subtitle2" color="inherit">Doublet Detection</Typography>}
								/>
				    		</Grid>
				    		<Grid item xs={6} sm={4}>
					    		<FormControlLabel
									control={<Checkbox 
												color="secondary" 
												checked={this.state.umapCheck}
												onChange={e=>this.setState({umapCheck:e.target.checked})}
											/>}
									label={<Typography variant="subtitle2" color="inherit">UMAP</Typography>}
								/>
				    		</Grid>
				    		<Grid item xs={6} sm={4}>
					    		<FormControlLabel
									control={<Checkbox 
												color="secondary" 
												checked={this.state.tsneCheck} 
												onChange={e=>this.setState({tsneCheck:e.target.checked})}
											/>}
									label={<Typography variant="subtitle2" color="inherit">t-SNE</Typography>}
								/>
				    		</Grid>
				    		<Grid item xs={6} sm={4}>
					    		<FormControlLabel
									control={<Checkbox 
												color="secondary" 
												checked={this.state.extract}
												onChange={e=>this.setState({geneExpCheck:e.target.checked})}
											/>}
									label={<Typography variant="subtitle2" color="inherit">Gene Expression Plots</Typography>}
								/>
				    		</Grid>
				    		<Grid item xs={6} sm={4}>
					    		<FormControlLabel
									control={<Checkbox 
												color="secondary" 
												checked={this.state.diffExpCheck} 
												onChange={e=>this.setState({diffExpCheck:e.target.checked})}
											/>}
									label={<Typography variant="subtitle2" color="inherit">Differential Expression</Typography>}
								/>
				    		</Grid>
				    		<Grid item xs={6} sm={4}>
					    		<FormControlLabel
									control={<Checkbox 
												color="secondary" 
												checked={this.state.batchCheck} 
												onChange={e=>this.setState({batchCheck:e.target.checked})}
											/>}
									label={<Typography variant="subtitle2" color="inherit">Batch Correct</Typography>}
								/>
				    		</Grid>
				    		<Grid item xs={6} sm={4}>
					    		<FormControlLabel
									control={<Checkbox 
												color="secondary" 
												checked={this.state.dptCheck} 
												onChange={e=>this.setState({dptCheck:e.target.checked})}
											/>}
									label={<Typography variant="subtitle2" color="inherit">Diffusion Pseudotime</Typography>}
								/>
				    		</Grid>
				    		<Grid item xs={6} sm={4}>
					    		<FormControlLabel
									control={<Checkbox 
												color="secondary" 
												checked={this.state.phateCheck} 
												onChange={e=>this.setState({phateCheck:e.target.checked})}
											/>}
									label={<Typography variant="subtitle2" color="inherit"> PHATE</Typography>}
								/>
				    		</Grid>
				    		<Grid item xs={6} sm={4}>
					    		<FormControlLabel
									control={<Checkbox 
												color="secondary" 
												checked={this.state.graphFA} 
												onChange={e=>this.setState({graphFACheck:e.target.checked})}
											/>}
									label={<Typography variant="subtitle2" color="inherit">Graph Force Atlas</Typography>}
								/>
				    		</Grid>
				    		<Grid item xs={6} sm={4}>
					    		<FormControlLabel
									control={<Checkbox 
												color="secondary" 
												checked={this.state.extractCheck}
												onChange={e=>this.setState({extractCheck:e.target.checked})}
											/>}
									label={<Typography variant="subtitle2" color="inherit">Cluster Extract</Typography>}
								/>
				    		</Grid>
				    	</Grid>

			  		</Container>
			  		<br/>
			  	</Grid>

			  	<Grid item xs={7} sm={7} md={7} component={Paper} className={classes.rightPane}>
				  	<Container>
					  	<Grid container spacing={1} className={classes.expansionGrid}>
					  		{(this.state.umapCheck || this.state.tsneCheck || this.state.dptCheck 
					  			|| this.state.phateCheck || this.state.graphFACheck) && 
					  			<Grid item sm={12} m={12} lg={12}>
							  	<Card>
							        <ExpansionPanel square>
								        <ExpansionPanelSummary aria-controls="panel1a-content" id="panel1a-header">
								          	<Typography variant="h7" color="inherit" noWrap>
									      		Graph Embedding Plot Options
									    	</Typography>
									    	<IconButton className={classes.button} aria-label="info" onClick={this.handlePlotFieldClickOpen}>
								    			<Info fontSize="small"/>
								    		</IconButton>
								    		<Dialog onClose={this.handlePlotFieldClose} aria-labelledby="simple-dialog-title" open={this.state.plotopen}>
										        <DialogTitle id="simple-dialog-title" className={classes.infoTitle}>Graph Embedding Plot Options</DialogTitle>
										        <DialogContent className={classes.infoContent}>
										        	<br/>
											        <b>Graph Dot Size </b> - Size of the dots on the UMAP, t-SNE, Diffusion Map, etc. plots
											        <br/>
											        <b>Graph Categories </b> - Categories to be plotted on graph embeddings (UMAP, t-SNE, etc.)
											        <br/>
													<b>Categorical Color </b> - Color palette for categorical 
													<br />
													<b>Color Map </b> - Gradient color map to visualize continuous data for feature plots, dot plots, heatmaps, etc.
													<br />
													<b>Figure Quality </b> - Low will produce regular resolution .png figures, whereas high will produce higher quality .pdf files
										        </DialogContent>
										    </Dialog>

								        </ExpansionPanelSummary>
								        <ExpansionPanelDetails>
									        <Grid container spacing={1}>
								    			{plotFields.map((field,i) =>(
								    				<Grid item xs={4} sm={3}>
									    				<TextField
													        id={field.id}
													        label={field.label}
													        className={classes.formField}
													        value={this.state.plotFieldValues[i]}
													        onChange={(e)=>{this.handlePlotFieldChange(e.target.value,i)}}
													        defaultValue={field.defaultValue}
													        fullWidth
													        margin="normal"
													        InputLabelProps={{
													        	shrink: true,
													        }}
													    />
												    </Grid>
												))}
										    </Grid>
										 
								        </ExpansionPanelDetails>
								    </ExpansionPanel>
						   		</Card>
						   	</Grid>}
						   	{(this.state.umapCheck) && <Grid item sm={12} m={12} lg={12}>
						   		<Card>
							        <ExpansionPanel square>
								        <ExpansionPanelSummary aria-controls="panel1a-content" id="panel1a-header">
								          	<Typography variant="h7" color="inherit" noWrap>
									      		UMAP Options
									    	</Typography>
									    	<IconButton className={classes.button} aria-label="info" onClick={this.handleUmapFieldClickOpen}>
								    			<Info fontSize="small"/>
								    		</IconButton>
								    		<Dialog onClose={this.handleUmapFieldClose} aria-labelledby="simple-dialog-title" open={this.state.umapopen}>
										        <DialogTitle id="simple-dialog-title" className={classes.infoTitle}>UMAP Options</DialogTitle>
										        <DialogContent className={classes.infoContent}>
										        	<br/>
											        <b>Initial Position </b> - Basis from which to intialize the UMAP ("spectral", "random", or "paga")
											        <br/>
													<b>Spread </b> - In combination with min_dist determines how clumped embedded points are
													<br />
													<b>Min Distance </b> - Minimum distance between points on the UMAP figure
										        </DialogContent>
										    </Dialog>

								        </ExpansionPanelSummary>
								        <ExpansionPanelDetails>
									        <Grid container spacing={1}>
									        	<Grid item xs={4} sm={4}>
									    			<FormControl className={classes.formControl} fullWidth>
												        <InputLabel htmlFor="age-simple">UMAP Initial Position</InputLabel>
												        <Select
												          value={this.state.umapPos}
												          onChange={e=>this.setState({umapPos:e.target.value})}
												          inputProps={{
												            name: 'umapPos',
												            id: 'umapPos',
												          }}
												        >
												          <MenuItem value={'spectral'}>Spectral</MenuItem>
												          <MenuItem value={'random'}>Random</MenuItem>
												          <MenuItem value={'paga'}>PAGA</MenuItem>
												        </Select>
												    </FormControl>
												</Grid>
								    			{umapFields.map((field,i) =>(
								    				<Grid item xs={4} sm={4}>
									    				<TextField
													        id={field.id}
													        label={field.label}
													        className={classes.formField}
													        value={this.state.umapFieldValues[i]}
													        onChange={(e)=>{this.handleUmapFieldChange(e.target.value,i)}}
													        defaultValue={field.defaultValue}
													        fullWidth
													        margin="normal"
													        InputLabelProps={{
													        	shrink: true,
													        }}
													    />
												    </Grid>
												))}
										    </Grid>
								        </ExpansionPanelDetails>
								    </ExpansionPanel>
						   		</Card>
						   	</Grid>}
						   	{(this.state.geneExpCheck) && <Grid item sm={12} m={12} lg={12}>
						   		<Card>
							        <ExpansionPanel square>
								        <ExpansionPanelSummary aria-controls="panel1a-content" id="panel1a-header">
								          	<Typography variant="h7" color="inherit" noWrap>
									      		Gene Expression Plot Options
									    	</Typography>
									    	<IconButton className={classes.button} aria-label="info" onClick={this.handleGeneExpFieldClickOpen}>
								    			<Info fontSize="small"/>
								    		</IconButton>
								    		<Dialog onClose={this.handleGeneExpFieldClose} aria-labelledby="simple-dialog-title" open={this.state.geneexpopen}>
										        <DialogTitle id="simple-dialog-title" className={classes.infoTitle}>Gene Expression Plot Options</DialogTitle>
										        <DialogContent className={classes.infoContent}>
										        	<br/>
											        <b>Gene Expression Categories </b> - Categories to plot on vertical axes of gene expression plots
											        <br/>
													<b>Dot Plot </b> - Dot color represents average expression level in a group while dot size represents percentage of cells expressing given gene
													<br />
													<b>Heat Map Plot </b> - Consecutive colored lines representing gene expression of individual cells 
													{/*<br />
													<b>Violin Plot </b> - Violin plot of gene expression levels*/}
										        </DialogContent>
										    </Dialog>

								        </ExpansionPanelSummary>
								        <ExpansionPanelDetails>
									        <Grid container spacing={1}>
									        	<Grid item xs={12} sm={12}>
								    				<TextField
												        id='geneExpCat'
												        label='Gene Expression Categories'
												        className={classes.formField}
												        value={this.state.geneExpCat}
												        onChange={e=>this.setState({geneExpCat:e.target.value})}
												        defaultValue="leiden, sampleName"
												        fullWidth
												        margin="normal"
												        InputLabelProps={{
												        	shrink: true,
												        }}
												    />
											    </Grid>
												<Grid item xs={4} sm={2}>
										    		<FormControlLabel
														control={<Checkbox 
																	color="secondary" 
																	checked={this.state.dotCheck}
																	onChange={e=>this.setState({dotCheck:e.target.checked})}
																/>}
														label={<Typography variant="subtitle2" color="inherit">Dot Plot</Typography>}
													/>
									    		</Grid>
									    		<Grid item xs={4} sm={2}>
										    		<FormControlLabel
														control={<Checkbox 
																	color="secondary" 
																	checked={this.state.heatCheck}
																	onChange={e=>this.setState({heatCheck:e.target.checked})}
																/>}
														label={<Typography variant="subtitle2" color="inherit">Heat Map Plot</Typography>}
													/>
									    		</Grid>
									    		{/*<Grid item xs={4} sm={2}>
										    		<FormControlLabel
														control={<Checkbox 
																	color="secondary" 
																	checked={this.state.violinCheck}
																	onChange={e=>this.setState({violinCheck:e.target.checked})}
																/>}
														label={<Typography variant="subtitle2" color="inherit">Violin Plot</Typography>}
													/>
									    		</Grid>*/}
										    </Grid>
								        </ExpansionPanelDetails>
								    </ExpansionPanel>
						   		</Card>
						   	</Grid>}
							{(this.state.diffExpCheck) && <Grid item sm={12} m={12} lg={12}>
							  	<Card>
							        <ExpansionPanel square>
								        <ExpansionPanelSummary aria-controls="panel1a-content" id="panel1a-header">
								          	<Typography variant="h7" color="inherit" noWrap>
									      		Differential Expression Analysis Options
									    	</Typography>
									    	<IconButton className={classes.button} aria-label="info" onClick={this.handleDiffExpFieldClickOpen}>
								    			<Info fontSize="small"/>
								    		</IconButton>
								    		<Dialog onClose={this.handleDiffExpFieldClose} aria-labelledby="simple-dialog-title" open={this.state.diffexpopen}>
										        <DialogTitle id="simple-dialog-title" className={classes.infoTitle}>Differential Expression Analysis Options</DialogTitle>
										        <DialogContent className={classes.infoContent}>
										        	<br/>
											        <b>Comparison Categories </b> - Category (e.g. leiden, sampleName, age) to calculate differential expression analysis
											        <br/>
											        <b>Groups to Compare </b> - Specific groups within the categories to compare
											        <br/>
										        </DialogContent>
										    </Dialog>

								        </ExpansionPanelSummary>
								        <ExpansionPanelDetails>
									        <Grid container spacing={1}>
								    			{diffExpFields.map((field,i) =>(
								    				<Grid item xs={6} sm={6}>
									    				<TextField
													        id={field.id}
													        label={field.label}
													        className={classes.formField}
													        value={this.state.diffExpFieldValues[i]}
													        onChange={(e)=>{this.handleDiffExpFieldChange(e.target.value,i)}}
													        defaultValue={field.defaultValue}
													        fullWidth
													        margin="normal"
													        InputLabelProps={{
													        	shrink: true,
													        }}
													    />
												    </Grid>
												))}
										    </Grid>
										 
								        </ExpansionPanelDetails>
								    </ExpansionPanel>
						   		</Card>
						   	</Grid>}
						   	{(this.state.batchCheck) && <Grid item sm={12} m={12} lg={12}>
						   		<Card>
							        <ExpansionPanel square>
								        <ExpansionPanelSummary aria-controls="panel1a-content" id="panel1a-header">
								          	<Typography variant="h7" color="inherit" noWrap>
									      		Batch Correction Options
									    	</Typography>
									    	<IconButton className={classes.button} aria-label="info" onClick={this.handleBatchCorrClickOpen}>
								    			<Info fontSize="small"/>
								    		</IconButton>
								    		<Dialog onClose={this.handleBatchCorrClose} aria-labelledby="simple-dialog-title" open={this.state.batchopen}>
										        <DialogTitle id="simple-dialog-title" className={classes.infoTitle}>Batch Correction Options</DialogTitle>
										        <DialogContent className={classes.infoContent}>
										        	<br/>
											        <b>BBKNN </b> - Batch Balanced K-Nearest Neighbors Correction algorithm
											        <br/>
													<b>ComBat </b> - ComBat batch correction
													<br />
													<b>Category </b> - Category across which to batch correct
													<br />
										        </DialogContent>
										    </Dialog>

								        </ExpansionPanelSummary>
								        <ExpansionPanelDetails>
									        <Grid container spacing={1}>
										        <Grid item xs={3}>
									    			<FormControl className={classes.formControl} fullWidth>
												        <InputLabel htmlFor="age-simple">Batch Correction Algorithm</InputLabel>
												        <Select
												          value={this.state.batchAlgo}
												          onChange={e=>this.setState({batchAlgo:e.target.value})}
												          inputProps={{
												            name: 'batch',
												            id: 'batch',
												          }}
												        >
												          <MenuItem value={'bbknn'}>BBKNN</MenuItem>
												          <MenuItem value={'combat'}>ComBat</MenuItem>
												        </Select>
												    </FormControl>
												</Grid>
												<Grid item xs={9}>
													<TextField
													        id='batch_id'
													        label='Category'
													        className={classes.formField}
													        value={this.state.batchCat}
													        onChange={e=>this.setState({batchCat:e.target.value})}
													        defaultValue='sampleName'
													        fullWidth
													        margin="normal"
													        InputLabelProps={{
													        	shrink: true,
													        }}
													    />
												</Grid>
										    </Grid>
								        </ExpansionPanelDetails>
								    </ExpansionPanel>
						   		</Card>
					   		</Grid>}
					   		{(this.state.dptCheck) && <Grid item sm={12} m={12} lg={12}>
							  	<Card>
							        <ExpansionPanel square>
								        <ExpansionPanelSummary aria-controls="panel1a-content" id="panel1a-header">
								          	<Typography variant="h7" color="inherit" noWrap>
									      		Diffusion Pseudotime Options
									    	</Typography>
									    	<IconButton className={classes.button} aria-label="info" onClick={this.handleDPTFieldClickOpen}>
								    			<Info fontSize="small"/>
								    		</IconButton>
								    		<Dialog onClose={this.handleDPTFieldClose} aria-labelledby="simple-dialog-title" open={this.state.dptopen}>
										        <DialogTitle id="simple-dialog-title" className={classes.infoTitle}>Diffusion Pseudotime Options</DialogTitle>
										        <DialogContent className={classes.infoContent}>
										        	<br/>
											        <b>Category for Root Cell </b> - Category (e.g. leiden, sampleName, age, etc.) for identifying root cell type
											        <br/>
											        <b>Root Cell Group </b> - Group (e.g. Cluster 1, Sample_HT239, etc.) within category containing root cell type
											        <br/>
										        </DialogContent>
										    </Dialog>

								        </ExpansionPanelSummary>
								        <ExpansionPanelDetails>
									        <Grid container spacing={1}>
								    			{dptFields.map((field,i) =>(
								    				<Grid item xs={6} sm={6}>
									    				<TextField
													        id={field.id}
													        label={field.label}
													        className={classes.formField}
													        value={this.state.dptFieldValues[i]}
													        onChange={(e)=>{this.handleDPTFieldChange(e.target.value,i)}}
													        defaultValue={field.defaultValue}
													        fullWidth
													        margin="normal"
													        InputLabelProps={{
													        	shrink: true,
													        }}
													    />
												    </Grid>
												))}
										    </Grid>
										 
								        </ExpansionPanelDetails>
								    </ExpansionPanel>
						   		</Card>
						   	</Grid>}
					   		{(this.state.extractCheck) && <Grid item sm={12} m={12} lg={12}>
						   		<Card>
							        <ExpansionPanel square>
								        <ExpansionPanelSummary aria-controls="panel1a-content" id="panel1a-header">
								          	<Typography variant="h7" color="inherit" noWrap>
									      		Extraction Options
									    	</Typography>
									    	<IconButton className={classes.button} aria-label="info" onClick={this.handleExtractClickOpen}>
								    			<Info fontSize="small"/>
								    		</IconButton>
								    		<Dialog onClose={this.handleExtractClose} aria-labelledby="simple-dialog-title" open={this.state.extractopen}>
										        <DialogTitle id="simple-dialog-title" className={classes.infoTitle}>Extraction Options</DialogTitle>
										        <DialogContent className={classes.infoContent}>
											        <b>Clusters to Extract </b> - Clusters to extract from previous analysis
										        </DialogContent>
										    </Dialog>

								        </ExpansionPanelSummary>
								        <ExpansionPanelDetails>
									        <Grid container spacing={1}>
												<Grid item xs={12}>
													<TextField
													        id='extractClusters'
													        label='Clusters to Extract'
													        className={classes.formField}
													        value={this.state.extractClusterst}
													        onChange={e=>this.setState({extractClusters:e.target.value})}
													        defaultValue=''
													        fullWidth
													        margin="normal"
													        InputLabelProps={{
													        	shrink: true,
													        }}
													    />
												</Grid>
										    </Grid>
								        </ExpansionPanelDetails>
								    </ExpansionPanel>
						   		</Card>
					   		</Grid>}
				    	</Grid>
				    </Container>
			  		
			  	</Grid>
			  	<Button
	          	   	fullWidth
	      			type="submit"
	      			variant="contained"
	      			color="secondary"
	      			className={classes.submit}
	                onClick={this.handleSubmit}
	            >
				    Run Analysis
		        </Button>
			</Grid>
		);
	}
}

export default withStyles(styles)(BasicForm);

{/*<Grid item sm={12} m={12} lg={12}>
					  	{this.state.tsneCheck && <Card>
					        <ExpansionPanel square>
						        <ExpansionPanelSummary
						          aria-controls="panel1a-content"
						          id="panel1a-header"
						        >
						          <Typography className={classes.heading}>Additional t-SNE Plot Options</Typography>
						        </ExpansionPanelSummary>
						        <ExpansionPanelDetails>
						          	<TextField
								        id="maxMito"
								        label="Max Mito"
								        style={{ margin: 8 }}
								        defaultValue="0.1"
								        helperText="Remove cells with higher mitochondrial gene expression"
								        fullWidth
								        margin="normal"
								        InputLabelProps={{
								        	shrink: true,
								        }}
								    />
						        </ExpansionPanelDetails>
						    </ExpansionPanel>
				   		</Card>}
				   		</Grid>
				   		<Grid item sm={12} m={12} lg={12}>
				   		{this.state.umapCheck && <Card>
					        <ExpansionPanel square>
						        <ExpansionPanelSummary
						          aria-controls="panel1a-content"
						          id="panel1a-header"
						        >
						          <Typography className={classes.heading}>Additional UMAP Plot Options</Typography>
						        </ExpansionPanelSummary>
						        <ExpansionPanelDetails>
						          	<TextField
								        id="maxMito"
								        label="Max Mito"
								        style={{ margin: 8 }}
								        defaultValue="0.1"
								        helperText="Remove cells with higher mitochondrial gene expression"
								        fullWidth
								        margin="normal"
								        InputLabelProps={{
								        	shrink: true,
								        }}
								    />
						        </ExpansionPanelDetails>
						    </ExpansionPanel>
				   		</Card>}
				   		</Grid>
				   		<Grid item sm={12} m={12} lg={12}>
					  	{this.state.batchCheck && <Card>
					        <ExpansionPanel square>
						        <ExpansionPanelSummary
						          aria-controls="panel1a-content"
						          id="panel1a-header"
						        >
						          	<Typography className={classes.heading}>
						          		Additional Violin Plot Options
						          	</Typography>
						        </ExpansionPanelSummary>
						        <ExpansionPanelDetails>
						        	<Grid container>

						          	<FormControl className={classes.formControl}>
								        <InputLabel htmlFor="age-simple">Age</InputLabel>
								        <Select
								          value={this.state.dotColor}
								          onChange={e=>this.setState({dotColor:e.target.value})}
								          inputProps={{
								            name: 'age',
								            id: 'age-simple',
								          }}
								        >
								          <MenuItem value={10}>Ten</MenuItem>
								          <MenuItem value={20}>Twenty</MenuItem>
								          <MenuItem value={30}>Thirty</MenuItem>
								        </Select>
								    </FormControl>
								    </Grid>
						        </ExpansionPanelDetails>
						    </ExpansionPanel> 
				   		</Card>}
				   		</Grid>*/}
