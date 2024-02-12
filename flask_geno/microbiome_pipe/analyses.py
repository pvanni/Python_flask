import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms
import os
import matplotlib as mpl
import math
from microbiome_pipe.util import ShellCommand, results_manager
from scipy import stats
from pathlib import Path


@results_manager(True)
def aldex2(project_path, prevalence_threshold, target_column=None):

    def aldex2_visualization(differentials_path, output_dir):
        """Create visualizations according to https:
        //github.com/ggloor/q2-aldex2/blob/master/q2_aldex2/_visualizer.py"""

        effect_statistic_function = 'we.eBH'
        threshold = 0.05

        # Change title sizes
        mpl.rcParams['axes.titlesize'] = 15

        table = pd.read_csv(differentials_path, index_col=0, sep='\t')
        table.drop(index='#q2:types', inplace=True)
        table = table.astype('float')

        fig = plt.figure(figsize=(12.0, 12.0))

        ax1 = fig.add_subplot(221)
        # base effect plot to build on
        ax1.scatter(x="diff.win", y="diff.btw", data=table, color="grey", s=6)

        # get maximum value of diff.btw to dynamically
        # lengthen the effect size line
        btw_max = table['diff.win'].max()

        # use the test that was input as a parameter by
        # taking column from the dict
        # subset features by cutoff
        called = table[effect_statistic_function] <= threshold

        # plot positive and negative effect size line
        ax1.plot([0, btw_max], [0, btw_max], color="grey",
                 linestyle='dashed', linewidth=1)
        ax1.plot([0, btw_max], [0, -btw_max], color="grey",
                 linestyle='dashed', linewidth=1)

        # colour for significant points
        ax1.scatter(x="diff.win", y="diff.btw", data=table[called],
                    color="red", s=6)

        # change titles and labels
        ax1.set_title('A)', loc='left')
        ax1.set_title('Effect plot')
        ax1.set_xlabel('Dispersion')
        ax1.set_ylabel('Difference')

        ########################################################################
        ax2 = fig.add_subplot(222)

        # base effect plot to build on
        ax2.scatter(x="rab.all", y="diff.btw", data=table, color="grey", s=6)

        # colour for significant points
        ax2.scatter(x="rab.all", y="diff.btw", data=table[called],
                    color="red", s=6)

        # change titles and labels
        ax2.set_title('B)', loc='left')
        ax2.set_title('MA plot')
        ax2.set_xlabel('Relative Abundance')
        ax2.set_ylabel('Difference')

        ########################################################################
        ax3 = fig.add_subplot(223)

        # base effect plot to build on
        ax3.scatter(x="diff.btw", y="we.eBH", data=table, color="grey", s=6)

        # colour for significant points
        ax3.scatter(x="diff.btw", y="we.eBH", data=table[called],
                    color="red", s=6)

        # change titles and labels
        ax3.set_title('C)', loc='left')
        ax3.set_title('Volcano plot')
        ax3.set_xlabel('Difference', )
        ax3.set_ylabel('Q score')

        # change p values to log scale
        ax3.set_yscale('log')
        # get min and max to plot for p values
        minimum = table["we.eBH"].min()
        maximum = table["we.eBH"].max()
        ax3.set_ylim([minimum, maximum])

        # plot line where cutoff is located
        ax3.plot([table["diff.btw"].min(), table["diff.btw"].max()],
                 [threshold, threshold], color="grey", linestyle='dashed',
                 linewidth=1)

        ########################################################################
        ax4 = fig.add_subplot(224)

        # base effect plot to build on
        ax4.scatter(x="effect", y="we.eBH", data=table,
                    color="grey", s=6)

        # colour for significant points
        ax4.scatter(x="effect", y="we.eBH", data=table[called],
                    color="red", s=6)

        # change titles and labels
        ax4.set_title('D)', loc='left')
        ax4.set_title('Effect size vs q score')
        ax4.set_xlabel('Effect size')
        ax4.set_ylabel('Q score')

        # change p values to log scale
        ax4.set_yscale('log')
        # get min and max to plot for p values
        minimum = table["we.eBH"].min()
        maximum = table["we.eBH"].max()
        ax4.set_ylim([minimum, maximum])

        # plot line where cutoff is located
        ax4.plot([table["effect"].min(), table["effect"].max()],
                 [threshold, threshold], color="grey",
                 linestyle='dashed', linewidth=1)

        plt.tight_layout()

        if not os.path.isdir(f'{output_dir}'):
            Path(f'{output_dir}').mkdir(parents=True, exist_ok=True)
        for i in ['png']:
            fig.savefig(f'{output_dir}{target_column}_aldex2_output_plots.{i}', dpi=600, bbox_inches='tight')

    def prevalence_filter_tables(project_path, prevalence_threshold=0.1):
        """According to https://docs.qiime2.org/2019.10/tutorials/moving-pictures/
        needs table collapsing to happen first"""

        num_samples = None

        for level in [2, 3, 4, 5, 6, 7]:

            df = pd.read_csv(f'{project_path}/feature_table/level{level}/feature_table.tsv', sep='\t',
                             header=1, index_col=0)
            num_samples = math.floor(df.shape[1] * prevalence_threshold)

            if os.path.isdir(f'{project_path}/feature_table/picrust/pathways/{num_samples}_pathway_table.qza'):
                return num_samples

            # Prevalence filter features that are in less than X% of samples out before ANCOM

            ShellCommand('qiime feature-table filter-features '
                         f'--i-table {project_path}/feature_table/level{level}/collapsed_table.qza '
                         f'--p-min-samples {num_samples} '
                         f'--o-filtered-table {project_path}/feature_table/level{level}/{num_samples}_collapsed_table.qza')

        #Picrust data
        ShellCommand('qiime feature-table filter-features '
                     f'--i-table {project_path}/feature_table/picrust/pathways/pathway_abundance.qza '
                     f'--p-min-samples {num_samples} '
                     f'--o-filtered-table  {project_path}/feature_table/picrust/pathways/{num_samples}_pathway_table.qza')

        return num_samples

    paragraph_list = []
    if not os.path.isdir(f'{project_path}/aldex2'):
        os.mkdir(f'{project_path}/aldex2')
    num_samples = prevalence_filter_tables(prevalence_threshold=prevalence_threshold, project_path=project_path)

    paragraph_list.append("Differentially abundant features in both taxa and pathway data (PICRUSt2 output) "
                          "were analyzed with q2-ALDEx2 -plugin. "
                          f"Feature table was filtered to exclude features "
                          f"not found in {int(prevalence_threshold*100)}% samples")

    # If theres more than 2 groups, return the paragraph
    meta = pd.read_csv(f'{project_path}/start_files/metadata.tsv', sep='\t')
    if len(np.unique(meta[target_column].values)) > 2:
        return paragraph_list

    """https://library.qiime2.org/plugins/q2-aldex2/24/"""
    # export the picrust data and round down the float values to ints
    # Export feature table as biom

    # Dont repeat this part for other target variables
    if not os.path.isfile(f'{project_path}/feature_table/picrust/pathways/rounding_to_int/rounded_feature_table.qza'):
        ShellCommand('qiime tools export '
                     f'--input-path {project_path}/feature_table/picrust/pathways/{num_samples}_pathway_table.qza '
                     f'--output-path {project_path}/feature_table/picrust/pathways/rounding_to_int')
        # Convert biom to tsv (feature table)
        ShellCommand(f'biom convert -i {project_path}/feature_table/picrust/pathways/rounding_to_int/feature-table.biom '
                     f'-o {project_path}/feature_table/picrust/pathways/rounding_to_int/feature_table.tsv --to-tsv')

        # Round the values to ints
        df = pd.read_csv(f'{project_path}/feature_table/picrust/pathways/rounding_to_int/feature_table.tsv', index_col=0, header=1, sep='\t')
        df = df.astype(int)
        df.to_csv(f'{project_path}/feature_table/picrust/pathways/rounding_to_int/rounded_table.tsv', index=True, sep='\t')

        ShellCommand(f'biom convert -i {project_path}/feature_table/picrust/pathways/rounding_to_int/rounded_table.tsv '
                     f'-o {project_path}/feature_table/picrust/pathways/rounding_to_int/rounded_table.biom --to-hdf5')

        # Import feature_table.qza
        ShellCommand("qiime tools import "
                     f"--input-path {project_path}/feature_table/picrust/pathways/rounding_to_int/rounded_table.biom "
                     "--type 'FeatureTable[Frequency]' "
                     "--input-format BIOMV210Format "
                     f"--output-path {project_path}/feature_table/picrust/pathways/rounding_to_int/rounded_feature_table.qza")

    #run aldex on both taxa and picrust data
    #picrust
    if not os.path.isdir(f'{project_path}/aldex2/picrust'):
        Path(f'{project_path}/aldex2/picrust').mkdir(parents=True, exist_ok=True)

    ShellCommand('qiime aldex2 aldex2 '
                 f'--i-table {project_path}/feature_table/picrust/pathways/rounding_to_int/rounded_feature_table.qza '
                 f'--m-metadata-file {project_path}/start_files/metadata.tsv '
                 f'--m-metadata-column {target_column} '
                 f'--o-differentials {project_path}/aldex2/picrust/{target_column}_differentials.qza')

    ShellCommand('qiime aldex2 effect-plot '
                 f'--i-table {project_path}/aldex2/picrust/{target_column}_differentials.qza '
                 f'--o-visualization {project_path}/aldex2/picrust/{target_column}_effect_plot.qzv')

    ShellCommand('qiime tools export '
                 f'--input-path {project_path}/aldex2/picrust/{target_column}_differentials.qza '
                 f'--output-path {project_path}/aldex2/picrust/{target_column}_diff')

    aldex2_visualization(f'{project_path}/aldex2/picrust/{target_column}_diff/differentials.tsv',
                         f'{project_path}/results/picrust/{target_column}/aldex2/')

    data_type_list = [f'level{i}' for i in range(2, 8)]

    for data_type in data_type_list:
        if not os.path.isdir(f'{project_path}/aldex2/{data_type}'):
            Path(f'{project_path}/aldex2/{data_type}').mkdir(parents=True, exist_ok=True)

        ShellCommand('qiime aldex2 aldex2 '
                     f'--i-table {project_path}/feature_table/{data_type}/{num_samples}_collapsed_table.qza '
                     f'--m-metadata-file {project_path}/start_files/metadata.tsv '
                     f'--m-metadata-column {target_column} '
                     f'--o-differentials {project_path}/aldex2/{data_type}/{target_column}_differentials.qza')

        ShellCommand('qiime aldex2 effect-plot '
                     f'--i-table {project_path}/aldex2/{data_type}/{target_column}_differentials.qza '
                     f'--o-visualization {project_path}/aldex2/{data_type}/{target_column}_effect_plot.qzv')

        #Export both differentials files
        ShellCommand('qiime tools export '
                     f'--input-path {project_path}/aldex2/{data_type}/{target_column}_differentials.qza '
                     f'--output-path {project_path}/aldex2/{data_type}/{target_column}_diff')

        aldex2_visualization(f'{project_path}/aldex2/{data_type}/{target_column}_diff/differentials.tsv',
                             f'{project_path}/results/{data_type}/{target_column}/aldex2/')


    return paragraph_list


@results_manager(True)
def beta_diversity_results(project_path, sample_column, target_column=None):
    """Take a specific beta diversity metric and create a collab picture with significance testing results"""

    def confidence_ellipse(x, y, ax, n_std=3.0, facecolor='none', **kwargs):
        """
        Create a plot of the covariance confidence ellipse of *x* and *y*.

        Parameters
        ----------
        x, y : array-like, shape (n, )
            Input data.

        ax : matplotlib.axes.Axes
            The axes object to draw the ellipse into.

        n_std : float
            The number of standard deviations to determine the ellipse's radiuses.

        Returns
        -------
        matplotlib.patches.Ellipse

        Other parameters
        ----------------
        kwargs : `~matplotlib.patches.Patch` properties
        """
        if x.size != y.size:
            raise ValueError("x and y must be the same size")

        cov = np.cov(x, y)
        pearson = cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1])
        # Using a special case to obtain the eigenvalues of this
        # two-dimensionl dataset.
        ell_radius_x = np.sqrt(1 + pearson)
        ell_radius_y = np.sqrt(1 - pearson)
        ellipse = Ellipse((0, 0), width=ell_radius_x * 2, height=ell_radius_y * 2,
                          facecolor=facecolor, **kwargs)

        # Calculating the stdandard deviation of x from
        # the squareroot of the variance and multiplying
        # with the given number of standard deviations.
        scale_x = np.sqrt(cov[0, 0]) * n_std
        mean_x = np.mean(x)

        # calculating the stdandard deviation of y ...
        scale_y = np.sqrt(cov[1, 1]) * n_std
        mean_y = np.mean(y)

        transf = transforms.Affine2D() \
            .rotate_deg(45) \
            .scale(scale_x, scale_y) \
            .translate(mean_x, mean_y)

        ellipse.set_transform(transf + ax.transData)
        return ax.add_patch(ellipse)

    def name_parser(project_path, target_column, sample_column):
        """Opens the metadata file and gets the correct group labels and names of the samples.
        Assumes that there is no categorical/numerical row at the start, for now"""

        metadata = pd.read_csv(f'{project_path}/start_files/metadata.tsv', sep='\t')

        # Contains the sample ids different lists based on groups
        sample_ids = []
        group_names = np.unique(metadata.loc[:, '{}'.format(target_column)].values)

        # Splices the pandas data frame for each unique group found in target column.
        for name in group_names:
            sample_ids.append(metadata.loc[metadata['{}'.format(target_column)] == name]
                              .loc[:, sample_column].values)

        return [sample_ids, group_names]

    def parse_PERMANOVA_p_value(file_path):

        info_dict = {'data_type': '',
                     'target variable': '',
                     'method name': '',
                     'test statistic name': '',
                     'sample size': '',
                     'number of groups': '',
                     'test statistic': '',
                     'p-value': '',
                     'number of permutations': ''}
        # Gather results and BH correct p-values for each data type
        with open(f'{file_path}') as perma_file:
            next_line = False
            current_key = ''
            for line in perma_file:
                if next_line:
                    if current_key == 'p-value':
                        info_dict[current_key] = float(line[line.index('>') + 1:][:line[line.index('>') + 1:].index('<')])
                    else:
                        info_dict[current_key] = line[line.index('>') + 1:][:line[line.index('>') + 1:].index('<')]
                    next_line = False

                for key in info_dict:
                    if key in line:
                        if key == 'test statistic' and 'test statistic name' in line:
                            continue
                        next_line = True
                        current_key = key

        return info_dict

    paragraph_list = []
    fig_list = []
    names = name_parser(project_path, target_column, sample_column)

    for metric, title in zip(['bray_curtis', 'jaccard', 'unweighted_unifrac', 'weighted_unifrac'],
                             ['Bray-Curtis Dissimilarity', 'Jaccard Distance', 'Unweighted UniFrac', 'Weighted UniFrac']):

        if not os.path.isdir(f'{project_path}/diversity/{metric}_results'):
            ShellCommand('qiime tools export '
                         f'--input-path {project_path}/diversity/core-metrics-results/{metric}_pcoa_results.qza '
                         f"--output-path '{project_path}/diversity/{metric}_results'")

            ShellCommand('qiime diversity beta-group-significance '
                         f'--i-distance-matrix {project_path}/diversity/core-metrics-results/{metric}_distance_matrix.qza '
                         f'--m-metadata-file {project_path}/start_files/metadata.tsv '
                         f'--m-metadata-column {target_column} '
                         f"--o-visualization '{project_path}/diversity/{metric}_results/permanova.qzv'")

            ShellCommand('qiime tools export '
                         f"--input-path '{project_path}/diversity/{metric}_results/permanova.qzv' "
                         f"--output-path '{project_path}/diversity/{metric}_results/permanova_export'")

        # Parse the variance explained and pcoa matrix to the results folder
        buffer = []
        with open(f'{project_path}/diversity/{metric}_results/ordination.txt', 'r') as file:
            for line_number, line in enumerate(file):

                # When to stop
                if 'Biplot' in line:
                    buffer = buffer[:len(buffer) - 1]
                    break

                # Write the one line that have all the proportions explained to prop_explained.txt
                if line_number == 4:
                    with open(f'{project_path}/diversity/{metric}_results/prop_explained.txt',
                              'w+') as prop_explained:
                        prop_explained.write(line)

                # Write the values of PC's to different file
                if line_number > 8:
                    buffer.append(line)

        # Write the pcoa matrix as pcoa_matrix.txt file
        with open(f'{project_path}/diversity/{metric}_results/pcoa_matrix.txt',
                  'w+') as file:
            file.write(''.join(buffer))

        fig = plt.figure(figsize=(6.0, 6.0))

        # Change title sizes

        mpl.rcParams['axes.titlesize'] = 15

        # bray curtis
        # Parse the top 2 pcoa variance explained values for each metric
        var_df = pd.read_csv(f'{project_path}/diversity/{metric}_results/prop_explained.txt',
                             header=None, sep='\t')

        bc_df = pd.read_csv(f'{project_path}/diversity/{metric}_results/pcoa_matrix.txt',
                            header=None, sep='\t')
        var_explained = [var_df.iloc[0, :].values[0], var_df.iloc[0, :].values[1]]

        ax1 = fig.add_subplot(111)
        ax1.set_title('D)', loc='left')
        ax1.set_title(f'{title}')

        # list of colors to use
        color_iter = iter(['r', 'b', 'g', 'm', 'c'])



        # Plot each group
        for index, name_group in enumerate(names[0]):
            xy_values = [[], []]
            for i, name in enumerate(bc_df.iloc[:, 0].values):
                if name in name_group:
                    xy_values[0].append(bc_df.iloc[i, 1])
                    xy_values[1].append(bc_df.iloc[i, 2])
            current_color = next(color_iter)
            ax1.scatter(np.array(xy_values[0]), np.array(xy_values[1]), color=current_color, label=names[1][index])
            confidence_ellipse(np.array(xy_values[0]), np.array(xy_values[1]), ax1,
                               n_std=2, facecolor=current_color, alpha=0.2)

        ax1.legend(loc=2)
        ax1.set_xlabel("PC1 (Variance explained {:.1f}%)".format(var_explained[0] * 100))
        ax1.set_ylabel("PC2 (Variance explained {:.1f}%)".format(var_explained[1] * 100))
        ax1.axhline(y=0, color='k', linewidth=0.5)
        ax1.axvline(x=0, color='k', linewidth=0.5)

        # get p-value
        info_dict = parse_PERMANOVA_p_value(f'{project_path}/diversity/{metric}_results/permanova_export/index.html')
        # plot the p-value results
        textstr = f"PERMANOVA test\nadjusted p-value: {round(info_dict['p-value'], 4)}"

        # these are matplotlib.patch.Patch properties
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

        # place a text box in upper left in axes coords
        ax1.text(0.55, 0.10, textstr, transform=ax1.transAxes, fontsize=10,
                verticalalignment='top', bbox=props)

        plt.tight_layout()

        # Create the folder before saving the results
        if not os.path.isdir(f'{project_path}/results/diversity/{target_column}'):
            Path(f'{project_path}/results/diversity/{target_column}').mkdir(parents=True, exist_ok=True)
        for i in ['png']:
            fig.savefig(f'{project_path}/results/diversity/{target_column}/{target_column}_{metric}_beta_scatter.{i}', dpi=600, bbox_inches='tight')
        plt.close('all')

    return []


@results_manager(False)
def diversity_methods(project_path, sampling_depth=1000):

    paragraph_list = []
    # Create tree for phylogenetic metrics

    if not os.path.isdir(f'{project_path}/diversity'):
        os.mkdir(f'{project_path}/diversity')

    ShellCommand('qiime phylogeny align-to-tree-mafft-fasttree '
                 f'--i-sequences {project_path}/feature_table/representative_sequences.qza '
                 f'--o-alignment {project_path}/diversity/aligned_rep_seqs.qza '
                 f'--o-masked-alignment {project_path}/diversity/masked_aligned_rep_seqs.qza '
                 f'--o-tree {project_path}/diversity/unrooted_tree.qza '
                 f'--o-rooted-tree {project_path}/diversity/rooted_tree.qza')

    paragraph_list.append('Rooted tree was created using the q2-phylogeny plugin and FastTree. ')

    # Run the core metrics command on qiime2
    ShellCommand('qiime diversity core-metrics-phylogenetic '
                 f'--i-phylogeny {project_path}/diversity/rooted_tree.qza '
                 f'--i-table {project_path}/feature_table/feature_table.qza '
                 f'--p-sampling-depth {sampling_depth} '
                 f'--m-metadata-file {project_path}/start_files/metadata.tsv '
                 f'--output-dir {project_path}/diversity/core-metrics-results')

    paragraph_list.append(f'Alpha and beta diversity metrics were calculated using the q2-diversity plugin with a'
                          f' sampling rate of: {sampling_depth}. ')

    # Create the files for different metrics
    os.mkdir(f'{project_path}/diversity/shannon_index')
    os.mkdir(f'{project_path}/diversity/faith_pd')
    os.mkdir(f'{project_path}/diversity/observed_OTU')
    os.mkdir(f'{project_path}/diversity/pielou_evenness')

    # Export the tsv files from each metric
    ShellCommand('qiime tools export '
                 f'--input-path {project_path}/diversity/core-metrics-results/shannon_vector.qza '
                 f'--output-path {project_path}/diversity/shannon_index')

    ShellCommand('qiime tools export '
                 f'--input-path {project_path}/diversity/core-metrics-results/faith_pd_vector.qza '
                 f'--output-path {project_path}/diversity/faith_pd')

    ShellCommand('qiime tools export '
                 f'--input-path {project_path}/diversity/core-metrics-results/observed_otus_vector.qza '
                 f'--output-path {project_path}/diversity/observed_OTU')

    ShellCommand('qiime tools export '
                 f'--input-path {project_path}/diversity/core-metrics-results/evenness_vector.qza '
                 f'--output-path {project_path}/diversity/pielou_evenness')

    return paragraph_list


@results_manager(True)
def alpha_boxplots(project_path, sample_column, target_column=None):

    def name_parser(project_path, target_column, sample_column):
        """Opens the metadata file and gets the correct group labels and names of the samples.
        Assumes that there is no categorical/numerical row at the start, for now"""

        metadata = pd.read_csv(f'{project_path}/start_files/metadata.tsv', sep='\t')

        # Contains the sample ids different lists based on groups
        sample_ids = []
        group_names = np.unique(metadata.loc[:, '{}'.format(target_column)].values)

        # Splices the pandas data frame for each unique group found in target column.
        for name in group_names:
            sample_ids.append(metadata.loc[metadata['{}'.format(target_column)] == name]
                              .loc[:, sample_column].values)

        return [sample_ids, group_names]

    paragraph_list = []
    names = name_parser(project_path, target_column, sample_column)

    boxplot_objects = []

    fig = plt.figure(figsize=(12.5, 7.5))

    # Set all title fonts larger
    mpl.rcParams['axes.titlesize'] = 15

    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(223)
    ax4 = fig.add_subplot(224)

    # Create the diversity file in the results
    if not os.path.isdir(f'{project_path}/results/diversity'):
        os.mkdir(f'{project_path}/results/diversity')

    for metric, i_ax, title, letter in zip(['shannon_index', 'faith_pd', 'pielou_evenness', 'observed_OTU'],
                            [ax1, ax2, ax3, ax4],
                            ['Shannon Index', "Faith's Phylogenetic Distance", "Pielou's Evenness", "Observed OTU's"],
                            ['A)', 'B)', 'C)', 'D)']):

        metric_data = pd.read_csv(f'{project_path}/diversity/{metric}/alpha-diversity.tsv', sep='\t')

        data_points = []
        for i in range(len(names[1])):
            # Create boolean array to take correct group samples from data
            boolean_array = np.in1d(metric_data.iloc[:, 0].values, names[0][i])

            # Based on boolean array, pick correct group values
            data_points.append(metric_data.iloc[boolean_array, 1].values)

        # Calculate Kruskal Wallis and export as tsv
        significances = stats.kruskal(*data_points)

        significances_df = pd.DataFrame(data=[['Kruskal-Wallis H-test', significances[0], significances[1]]],
                                        columns=['Test', 'Statistic', 'p-value'])

        significances_df.to_csv(f'{project_path}/results/diversity/{target_column}_{metric}.tsv', index=False, sep='\t')

        # Insert n-counts to the labels
        group_labels = []
        for i, name in enumerate(names[1]):
            group_labels.append('{} (n={})'.format(name, len(names[0][i])))

        values_grouped = [[] for j in names[0]]
        labels = ['' for j in names[1]]
        widths = [0.3 for j in names[1]]
        for index, name_group in enumerate(names[0]):
            for i, value in enumerate(metric_data.iloc[:, 1].values):
                if metric_data.iloc[i, 0] in name_group:
                    values_grouped[index].append(value)

        i_ax.set_title(letter, loc='left')
        i_ax.set_title(title)
        # Box with statistical test results
        kw_text = f'KW H-test\np-value: {round(significances[1], 4)}'
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        i_ax.text(0.05, 0.15, kw_text, transform=i_ax.transAxes, fontsize=8,
                verticalalignment='top', bbox=props)


        if letter in ['C)', 'D)']:
            boxplot_objects.append(i_ax.boxplot(values_grouped, labels=group_labels, widths=widths))
        else:
            boxplot_objects.append(i_ax.boxplot(values_grouped, labels=labels, widths=widths))

    fig.savefig(f'{project_path}/results/{target_column}_Alpha_boxplots.png', dpi=600, bbox_inches='tight')

    return paragraph_list

@results_manager(True)
def ancom(project_path, prevalence_threshold, target_column=None):

    def prevalence_filter_tables(project_path, prevalence_threshold=0.1):
        """According to https://docs.qiime2.org/2019.10/tutorials/moving-pictures/
        needs table collapsing to happen first"""

        num_samples = None

        for level in [2, 3, 4, 5, 6, 7]:

            df = pd.read_csv(f'{project_path}/feature_table/level{level}/feature_table.tsv', sep='\t',
                             header=1, index_col=0)
            num_samples = math.floor(df.shape[1] * prevalence_threshold)

            if os.path.isdir(f'{project_path}/feature_table/picrust/pathways/{num_samples}_pathway_table.qza'):
                return num_samples

            # Prevalence filter features that are in less than X% of samples out before ANCOM

            ShellCommand('qiime feature-table filter-features '
                         f'--i-table {project_path}/feature_table/level{level}/collapsed_table.qza '
                         f'--p-min-samples {num_samples} '
                         f'--o-filtered-table {project_path}/feature_table/level{level}/{num_samples}_collapsed_table.qza')

        #Picrust data
        ShellCommand('qiime feature-table filter-features '
                     f'--i-table {project_path}/feature_table/picrust/pathways/pathway_abundance.qza '
                     f'--p-min-samples {num_samples} '
                     f'--o-filtered-table  {project_path}/feature_table/picrust/pathways/{num_samples}_pathway_table.qza')

        return num_samples

    paragraph_list = []
    if not os.path.isdir(f'{project_path}/ancom'):
        os.mkdir(f'{project_path}/ancom')
    num_samples = prevalence_filter_tables(prevalence_threshold=prevalence_threshold, project_path=project_path)

    # export the picrust data and round down the float values to ints
    # Export feature table as biom

    # Dont repeat this part for other target variables
    if not os.path.isfile(f'{project_path}/feature_table/picrust/pathways/rounding_to_int/rounded_feature_table.qza'):
        ShellCommand('qiime tools export '
                     f'--input-path {project_path}/feature_table/picrust/pathways/{num_samples}_pathway_table.qza '
                     f'--output-path {project_path}/feature_table/picrust/pathways/rounding_to_int')
        # Convert biom to tsv (feature table)
        ShellCommand(
            f'biom convert -i {project_path}/feature_table/picrust/pathways/rounding_to_int/feature-table.biom '
            f'-o {project_path}/feature_table/picrust/pathways/rounding_to_int/feature_table.tsv --to-tsv')

        # Round the values to ints
        df = pd.read_csv(f'{project_path}/feature_table/picrust/pathways/rounding_to_int/feature_table.tsv',
                         index_col=0, header=1, sep='\t')
        df = df.astype(int)
        df.to_csv(f'{project_path}/feature_table/picrust/pathways/rounding_to_int/rounded_table.tsv', index=True,
                  sep='\t')

    paragraph_list.append("Differentially abundant features in both taxa and pathway data (PICRUSt2 output) "
                          "were analyzed with ANCOM2. "
                          f"Feature table was filtered to exclude features "
                          f"not found in {int(prevalence_threshold * 100)}% samples")

    # If theres more than 2 groups, return the paragraph
    meta = pd.read_csv(f'{project_path}/start_files/metadata.tsv', sep='\t')
    if len(np.unique(meta[target_column].values)) > 2:
        return paragraph_list

    # Run the ANCOM rscript to produce results
    # Picrust
    if not os.path.isdir(f'{project_path}/ancom/picrust'):
        Path(f'{project_path}/ancom/picrust').mkdir(parents=True, exist_ok=True)

        # input table: {project_path}/feature_table/picrust/pathways/rounding_to_int/rounded_table.tsv

    # Taxonomic data
    data_type_list = [f'level{i}' for i in range(2, 8)]

    for data_type in data_type_list:
        if not os.path.isdir(f'{project_path}/ancom/{data_type}'):
            Path(f'{project_path}/ancom/{data_type}').mkdir(parents=True, exist_ok=True)