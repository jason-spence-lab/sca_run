'''
Plotting -- 
Different ways to display data and facilitate data interpretation

Written by Joshua H Wu
09 March, 2021
'''
import scanpy as sc
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import copy
import os
import csv
import pandas as pd

class plotting:

    def __init__(self, plot_params):
        '''
        Plot Params --
            size: Size of dots on all UMAP plots
            umap_obs: Metadata observations by which to color UMAP plots
            dot_check: Plot dot plot figures based off of gene lists
            heatmap_check: Plot heatmap figures based off of gene lists
            dot_grouping: Metadata observations by which to group dot plot rows
            umap_categorical_color: List of colors for UMAP groups (HEX codes or RGB tuples)
            umap_feature_color: Color scheme for UMAP gradient plots (yellow_blue or blue_orange)
            vmin_list: List of y-axis minimums for cell scoring feature plots
            vmax_list:List of y-axis maximums for cell scoring feature plots
            rank_grouping: Metadata observations to group differential gene analyses
            clusters2_compare: Clusters to compare for differential gene analysis
            final_quality: Makes high resolution figures in pdf
        '''
        self.size = plot_params.size
        self.umap_obs = plot_params.umap_obs
        self.dot_check = plot_params.dot_check
        self.heatmap_check = plot_params.heatmap_check
        self.dot_grouping = plot_params.dot_grouping
        self.umap_categorical_color = plot_params.umap_categorical_color
        self.umap_feature_color = plot_params.umap_feature_color
        self.vmin_list = plot_params.vmin_list
        self.vmax_list = plot_params.vmax_list
        self.rank_grouping = plot_params.rank_grouping
        self.clusters2_compare = plot_params.clusters2_compare
        self.final_quality = plot_params.final_quality

    ## Variety of plotting and data display functions
    # Will make more robust in the future
    def plot_sca(self, 
                 adata,
                 sca_params,
                 figdir='./figures/'):
        '''
        See the Scanpy visualization library for examples
        '''
        print("Plotting")

        ## Create my custom palette for FeaturePlots and define a matlplotlib colormap object
        if self.umap_feature_color=='blue_orange':
            feature_colors = [(35,35,142), (255,127,0)]
            my_feature_cmap = self.make_cmap(feature_colors,bit=True)
        elif self.umap_feature_color=='yellow_blue':
            feature_colors = [(210,210,210), (210,210,210), (245,245,200), (100,200,225), (0,45,125)]
            position=[0, 0.019999, 0.02, 0.55, 1]
            my_feature_cmap = self.make_cmap(feature_colors,bit=True,position=position)
        else:
            feature_colors = self.umap_feature_color[0]
            position=self.umap_feature_color[1]
            my_feature_cmap = self.make_cmap(feature_colors,bit=True,position=position)

        gray_cmap = self.make_cmap([(220,220,220),(220,220,220)], bit=True)


        ## Check to see if user specified a color palette for categorical umap plots, ie. leiden, obs_fields
        if self.umap_categorical_color=='default':
            ## Custom color palette for cluster plots and observation plots
            colors = [(1,0.5,0),(0.5,0.5,0.85),(0,1,0),(1,0,0),(0,0,0.9),(0,1,1),
                    (0.4,0.4,0.4),(0.5,0.85,0.5),(0.5,0.15,0.5),
                    (0.15,0.5,0.5),(0.5,0.5,0.15),(0.9,0.9,0),(1,0,1),
                    (0,0.5,1),(0.85,0.5,0.5),(0.5,1,0),(0.5,0,1),(1,0,0.5),(0,0.9,0.6),
                    (0.3,0.6,0),(0,0.3,0.6),(0.6,0.3,0),(0.3,0,0.6),(0,0.6,0.3),(0.6,0,0.3)]
        else:
            colors = self.umap_categorical_color

        ## General figure parameters and settings
        sc.set_figure_params(dpi_save=300,dpi=300)#,vector_friendly=False)
        sc.settings.figdir = figdir
        sc.set_figure_params(fontsize=12)
        size = self.size

        # Check to see if user wants publication quality figures
        if self.final_quality=='high':
            # rcParams['figure.figsize'] = 4, 4
            rcParams['savefig.dpi']=1200
            file_type = '.pdf'
        else:
            file_type = '.png'

        ## Violin plots for filtering parameters pre and post
        sc.pl.violin(adata, ['n_genes','n_counts','percent_mito'],
                     jitter=0.4, multi_panel=True, save='_postFiltered_plot.png', show=False)
        if sca_params.adata_preQC:
            sc.pl.violin(sca_params.adata_preQC, ['n_genes','n_counts','percent_mito'],
                        jitter=0.4, multi_panel=True, save='_preFiltered_plot.png', show=False)

        ## Draw the PCA elbow plot to determine which PCs to use
        sc.pl.pca_variance_ratio(adata, log=True, n_pcs=50, save='_elbowPlot.png', show=False)
        ## Ranks and displays most contributing genes for each principal component
        components = 4
        loadings_components = range(sca_params.analysis_params.n_pcs-components, sca_params.analysis_params.n_pcs+components+1)
        sc.pl.pca_loadings(adata, components=loadings_components, save='_rank_genes.png', show=False)

        ## Plot results of UMAP (and t-SNE) dimensional reduction and clustering
        for observation in self.umap_obs:
            legend = 'on data' if (observation==sca_params.analysis_params.clustering_choice) else 'right margin'

            if sca_params.analysis_params.umap_check:
                sc.pl.umap(adata, color=observation, save=''.join(['_',observation,file_type]), show=False,
                           legend_loc=legend, edges=False, size=size, palette=colors, alpha=0.75)

            if sca_params.analysis_params.do_tSNE:
                sc.pl.tsne(adata, color=observation, save=''.join(['_',observation,file_type]), show=False, 
                            legend_loc=legend, edges=False, size=size, palette=colors, alpha=0.75)

            if sca_params.analysis_params.draw_force_atlas:
                sc.pl.draw_graph(adata, color=observation, save=''.join(['_',observation,file_type]), show=False, 
                            legend_loc=legend, edges=False, size=size, palette=colors, alpha=0.75)

            if sca_params.analysis_params.dpt:
                sc.pl.diffmap(adata, color=observation, save=''.join(['_',observation,file_type]), show=False, 
                            legend_loc=legend, edges=False, size=size, palette=colors, alpha=0.75)

            if sca_params.analysis_params.phate:
                sc.external.pl.phate(adata, color=observation, save=''.join(['_',observation,file_type]), show=False, 
                            legend_loc=legend, edges=False, size=size, palette=colors, alpha=0.75)

        ## Find marker genes via Wilxocon test based on cluster assignment
        # Create a simple plot to show the top 25 most significant markers for each cluster
        # Write most significant markers to a csv file
        for rank_grouping in self.rank_grouping:        
            n_genes_rank = 5
            for comparison in self.clusters2_compare:
                if comparison == 'all':
                    comparison=None
                self.__rank_genes(adata,rank_grouping,figdir=figdir,clusters2_compare=comparison)

            if 'all' in self.clusters2_compare:
                sc.pl.rank_genes_groups_heatmap(adata, n_genes=n_genes_rank, use_raw=True, show=False, 
                        save=''.join(['_rank_heatmap_',rank_grouping,file_type]), cmap=my_feature_cmap)
                sc.pl.rank_genes_groups_dotplot(adata, n_genes=n_genes_rank, use_raw=True, show=False, 
                        save=''.join(['_rank_dotplot_',rank_grouping,file_type]), color_map=my_feature_cmap)
                sc.pl.rank_genes_groups_stacked_violin(adata, n_genes=n_genes_rank, use_raw=True, 
                        show=False, save=''.join(['_rank_violin_',rank_grouping,file_type]))

        if sca_params.qc_params.doublet_detection:
            import doubletdetection
            if sca_params.analysis_params.umap_check:
                sc.pl.umap(sca_params.adata_doublet, color=['doublet_label','doublet_score','n_genes','n_counts','percent_mito'], save='_doublet_test.pdf', show=False, edges=False, size=size)
            f = doubletdetection.plot.convergence(sca_params.doublet_clf, save=''.join([figdir,'convergence_test.pdf']), show=False, p_thresh=1e-16, voter_thresh=0.5)
            f3 = doubletdetection.plot.threshold(sca_params.doublet_clf, save=''.join([figdir,'threshold_test.pdf']), show=False, p_step=6)

        ## Feature plots and dot plot analysis for each specified set of genes
        #sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, save='_markerPlots.png', show=False)
        if sca_params.gene_lists:
            missing_genes = []
            for gene_list in sca_params.gene_lists:
                gene_obj = sca_params.gene_dict[gene_list]
                genes_to_plot = []
                [genes_to_plot,missing_genes] = self.__find_genes(adata,gene_obj.markers, missing_genes=missing_genes)

                ## Do FeaturePlots for select genes
                print('Plotting standard marker genes: ',genes_to_plot,'\n')
                if sca_params.analysis_params.umap_check:
                    sc.pl.umap(adata, color=genes_to_plot, save= ''.join(['_featureplots_',gene_list,file_type]), show=False, 
                               cmap=my_feature_cmap, size=size, use_raw=True, vmin=0)

                if sca_params.analysis_params.do_tSNE:
                    sc.pl.tsne(adata, color=genes_to_plot, save= ''.join(['_featureplots_',gene_list,file_type]), show=False, 
                           cmap=my_feature_cmap, size=size, use_raw=True, vmin=0)

                if sca_params.analysis_params.draw_force_atlas:
                    sc.pl.draw_graph(adata, color=genes_to_plot, save= ''.join(['_featureplots_',gene_list,file_type]), show=False, 
                           cmap=my_feature_cmap, size=size, use_raw=True, vmin=0)

                if sca_params.analysis_params.dpt:
                    sc.pl.diffmap(adata, color=genes_to_plot, save= ''.join(['_featureplots_',gene_list,file_type]), show=False, 
                           cmap=my_feature_cmap, size=size, use_raw=True, vmin=0)

                if sca_params.analysis_params.phate:
                    sc.external.pl.phate(adata, color=genes_to_plot, save= ''.join(['_featureplots_',gene_list,file_type]), show=False, 
                           cmap=my_feature_cmap, size=size, use_raw=True, vmin=0)
        
                feature_positions = gene_obj.feature_positions # Manually set and determined
                feature_groups = gene_obj.feature_groups
                groupby_positions = gene_obj.groupby_positions

                if len(gene_obj.markers)!=0:
                    for grouping in self.dot_grouping:
                        ## Dotplot analysis
                        # Circle color corresponds to expression level, and circle size corresponds to percentage of cells expressing gene

                        ## Reordering categories for dotplot or heatmap rows
                        adata_plots = adata.copy()
                        dendrogram=False
                        if groupby_positions:
                            dendrogram = False
                            clustering_chosen = sca_params.analysis_params.clustering_choice
                            adata_plots.obs[clustering_chosen]=adata.obs[clustering_chosen].cat.reorder_categories(groupby_positions,inplace = False)
                        if self.dot_check:
                            sc.pl.dotplot(adata_plots, genes_to_plot, groupby=grouping, 
                                    var_group_positions=feature_positions, var_group_labels=feature_groups,
                                    save=''.join(['_markers_',gene_list,'_',grouping,file_type]), show=False, 
                                    color_map=my_feature_cmap, use_raw=True, dendrogram=dendrogram)#, figsize=(4,6))#, dot_max=0.4)#, dendrogram=True)
                        ## Heatmaps
                        # Each horizontal line represents expression of one cell
                        if self.heatmap_check:
                            sc.pl.heatmap(adata_plots, genes_to_plot, groupby=grouping, 
                                    var_group_positions=feature_positions, var_group_labels=feature_groups,
                                    save=''.join(['_markers_',gene_list,'_',grouping,file_type]), show=False, 
                                    cmap=my_feature_cmap, use_raw=True)

        # Genes that are not expressed or are invariable are plotted using a grayscale
        sca_params.missing_genes = missing_genes
        print('Plotting empty genes: ',missing_genes,'\n')
        empty_genes = [gene for gene in missing_genes if (gene in adata.raw.var_names)]
        genes_noseq = [gene for gene in missing_genes if (gene not in empty_genes)]
        print('Zero genes: ', empty_genes, '\n')
        print('Gene not in dataset: ', genes_noseq, '\n')
        if empty_genes:
            sc.pl.umap(adata, color=empty_genes, save=''.join(['_featureplots_gray',file_type]), 
                    show=False, cmap=gray_cmap, size=size, use_raw=True)

        # tSNE Plots - should move to integrate in umap code
        if sca_params.analysis_params.do_tSNE:
            sc.pl.tsne(adata, color=missing_genes, save=''.join(['_featureplots_gray',file_type]), 
                    show=False, cmap=gray_cmap, size=size, use_raw=True)

        if sca_params.analysis_params.draw_force_atlas:
            sc.pl.draw_graph(adata, color=missing_genes, save=''.join(['_featureplots_gray',file_type]), 
                    show=False, cmap=gray_cmap, size=size, use_raw=True)

        if sca_params.analysis_params.dpt:
            sc.pl.diffmap(adata, color=missing_genes, save=''.join(['_featureplots_gray',file_type]), 
                    show=False, cmap=gray_cmap, size=size, use_raw=True)

        if sca_params.analysis_params.phate:
            sc.external.pl.phate(adata, color=missing_genes, save=''.join(['_featureplots_gray',file_type]), 
                    show=False, cmap=gray_cmap, size=size, use_raw=True)

        # Generate a umap feature plot based on cell scoring
        if sca_params.cell_score_lists:
            max_list_len = max([len(self.vmax_list),len(self.vmin_list),len(sca_params.cell_score_lists)])
            if not self.vmax_list:
                vmax = [adata.obs.loc[:,sca_params.cell_score_lists].values.max()]*max_list_len
            else:
                vmax = self.vmax_list

            if not self.vmin_list:
                vmin = [adata.obs.loc[:,sca_params.cell_score_lists].values.min()]*max_list_len
            else:
                vmin= self.vmin_list

            for i, score_name in enumerate(sca_params.cell_score_lists):
                sc.pl.umap(adata, color=score_name, 
                           save=''.join(['_',score_name,'_cellType_score.png']), show=False, edges=False, color_map=my_feature_cmap, 
                           size=size, vmin=vmin[i], vmax=vmax[i])
                sc.pl.umap(adata, color=score_name, 
                           save=''.join(['_',score_name,'_cellType_score_0min.png']), show=False, edges=False, color_map=my_feature_cmap, 
                           size=size, vmin=0, vmax=vmax[i])

            sc.pl.violin(adata,sca_params.cell_score_lists, groupby='sampleName',
                         jitter=0.4, save='_cell_scores.png',show=False,multi_panel=True,rotation=90)

        if sca_params.analysis_params.dpt:
            sc.pl.diffmap(adata, color=['dpt_pseudotime', sca_params.analysis_params.clustering_choice], size=self.size, show=False,
                          save=''.join([sca_params.analysis_params.dpt[0],'.png']))
            if sca_params.analysis_params.umap_check:
                sc.pl.umap(adata, color='dpt_pseudotime', size=self.size, show=False,
                           save=''.join(['_','dpt','_',sca_params.analysis_params.dpt[0],'.png']))
            # sc.pl.dpt_groups_pseudotime(adata, color_map=my_feature_cmap, 
            #                           save=''.join([sca_params.analysis_params.dpt[0],sca_params.analysis_params.dpt[1],'.png']))
            # sc.pl.dpt_timeseries(adata, color_map=my_feature_cmap, show=False,
            #                    save=''.join([sca_params.analysis_params.dpt[0],sca_params.analysis_params.dpt[1],'.png']))

        
        # Custom violin plot module -- Not complete/in testing
        df = pd.DataFrame()
        # Add Gaussian y-jitter to better visualize zero expression in violin plots
        for gene in genes_to_plot:
            sigma = np.amax(adata.raw[:,gene].X)/40
            gene_df = [cell if (cell!=0) else np.random.normal(loc=0,scale=sigma) for cell in adata.raw[:,gene].X]
            df[gene] = gene_df

        ## Scatter plots to identify clusters that are high in number of genes, UMI counts, and mito transcript fraction
        # adata.obs['jitter'] = np.random.rand(len(adata.obs_names))*10
        # sc.pl.scatter(adata,x='jitter',y='n_genes',color=sca_params.analysis_params.clustering_choice,save='_n_genes.png',palette=colors,show=False)
        # sc.pl.scatter(adata,x='jitter',y='n_counts',color=sca_params.analysis_params.clustering_choice,save='_n_counts.png',palette=colors,show=False)
        # sc.pl.scatter(adata,x='jitter',y='percent_mito',color=sca_params.analysis_params.clustering_choice,save='_percent_mito.png',palette=colors,show=False)
    

        sc.pl.umap(adata,color=['n_genes','n_counts','percent_mito'],color_map=my_feature_cmap,save='_counts_check.png',show=False)

        # Set the thresholds and scaling factors for drawing the paga map/plot
        node_size_scale=1.25
        node_size_power=0.9
        edge_width_scale=1
        min_edge_width=0.035
        max_edge_width=2
        threshold=0.05
        sc.pl.paga(adata, layout='fr', threshold=threshold, node_size_scale=node_size_scale, 
            node_size_power=node_size_power, edge_width_scale=edge_width_scale,
            min_edge_width=min_edge_width, max_edge_width=max_edge_width, show=False, save = '_pagaPlot.png',
            title='PAGA: Fruchterman Reingold',frameon=False)
        
        return adata

    ## Define function to generate a color gradient from a defined starting and ending color
    def make_cmap(self,colors, position=None, bit=False):
        '''
        make_cmap takes a list of tuples which contain RGB values. The RGB
        values may either be in 8-bit [0 to 255] (in which bit must be set to
        True when called) or arithmetic [0 to 1] (default). make_cmap returns
        a cmap with equally spaced colors.
        Arrange your tuples so that the first color is the lowest value for the
        colorbar and the last is the highest.
        position contains values from 0 to 1 to dictate the location of each color.
        Default sets position 0 of cmap to light gray, and begins proposed gradient
        from >0
        '''
        bit_rgb = np.linspace(0,1,256)
        cdict = {'red':[(0,bit_rgb[220],bit_rgb[220])],'green':[(0,bit_rgb[220],bit_rgb[220])],'blue':[(0,bit_rgb[220],bit_rgb[220])]}
        if position == None:
            position = np.linspace(0.000001,1,len(colors))
        else:
            cdict = {'red':[],'green':[],'blue':[]}
            if len(position) != len(colors):
                sys.exit("position length must be the same as colors")
            elif position[0] != 0 or position[-1] != 1:
                sys.exit("position must start with 0 and end with 1")
        if bit:
            for i in range(len(colors)):
                colors[i] = (bit_rgb[colors[i][0]],
                             bit_rgb[colors[i][1]],
                             bit_rgb[colors[i][2]])
        for pos, color in zip(position, colors):
            cdict['red'].append((pos, color[0], color[0]))
            cdict['green'].append((pos, color[1], color[1]))
            cdict['blue'].append((pos, color[2], color[2]))

        cmap = mpl.colors.LinearSegmentedColormap('my_colormap',cdict,256)
        return cmap

    ## Writes results of rank genes analysis to multiple csv files, each representing a cluster
    def __rank_genes(self, adata, rank_grouping, clusters2_compare=None, figdir='./figures/'):
        '''
        rank_grouping: Adata observation metadata categories to compare
        clusters2_compare: Selection of either 2 clusters to compare - if none, then do 1vAll comparison

        Need to add file clearing
        '''
        if clusters2_compare == None: # Do default 1vAll comparison
            print("No clusters selected for comparison. Doing default 1vAll comparison")
            sc.tl.rank_genes_groups(adata,rank_grouping ,method='t-test', rankby_abs=False, n_genes=200)
            self.__write_rank_genes(adata, rank_grouping, clusters2_compare, figdir)
        else: # Compare 
            adata_temp = adata[adata.obs[rank_grouping].isin(clusters2_compare)]
            sc.tl.rank_genes_groups(adata_temp, rank_grouping, method='t-test', n_genes=200)
            self.__write_rank_genes(adata_temp, rank_grouping, clusters2_compare, figdir)
        return 0

    ## Actually does the writing to csv files of the rank genes analysis
    def __write_rank_genes(self,adata, rank_grouping, clusters2_compare, figdir='./figures/'):
        rank_genes_data = copy.deepcopy(adata.uns['rank_genes_groups']) # create copy of data to manipulate
        rank_genes_data.pop('params')
        if clusters2_compare == None:
            clusters2_compare=['all']

        for cluster in adata.obs[rank_grouping].cat.categories:
            csv_fileName = '/'.join([figdir,'csv_files','_'.join([rank_grouping]+clusters2_compare),
                '_'.join([cluster,'compare.csv'])])
            os.makedirs(os.path.dirname(csv_fileName), exist_ok=True) # Make file if it doesn't exist already
            with open(csv_fileName,'w',newline='') as f:
                wr = csv.writer(f)
                wr.writerow(self.ele_swap(list(rank_genes_data.keys()),0,1))
                wr.writerows(zip(*self.ele_swap([params[cluster] for params in rank_genes_data.values()],0,1)))

    ## Takes a list of genes and determines if they exist within the data set and are variable
    # Appends results to genes_exist or missing_genes list if given
    def __find_genes(self,adata, gene_list, missing_genes=None):

        genes_missing_in_this_list = []
        genes_exist = []

        ## Splits given gene list into two based on whether or not they exist or are invariable
        for gene in gene_list:
            gene_zero = (gene in adata.raw.var_names) and np.any(adata.raw[:,gene].X.toarray())
            (genes_exist if gene_zero else missing_genes).append(gene)
            if not gene_zero: genes_missing_in_this_list.append(gene)

        missing_genes = list(set(missing_genes))
        genes_missing_in_this_list = list(set(genes_missing_in_this_list))
        if genes_missing_in_this_list:
            print('Sorry, the following genes are not expressed in this dataset or are invariable:',genes_missing_in_this_list,'\n')

        return [genes_exist, missing_genes]

    ## Swaps the elements at the proposed indices in an applicable data structure
    def ele_swap(self, structure, index1, index2):
        structure[index1], structure[index2] = structure[index2], structure[index1]
        return structure