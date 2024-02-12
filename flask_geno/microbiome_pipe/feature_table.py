import pandas as pd
import numpy as np
import os
from microbiome_pipe.util import ShellCommand, results_manager


@results_manager(False)
def filtering_table(project_path, table_path):
    """Take the feature table and apply taxonomic filtering and chimera removal. Replace the original repseq and table
    in the feature_table/ directory with the filtered ones"""

    paragraph_list = []

    os.mkdir(f'{project_path}/feature_table')
    os.mkdir(f'{project_path}/taxonomy')

    # Filter out chimeras based with qiime2 uchime plugin

    ShellCommand('qiime vsearch uchime-denovo '
                 f'--i-table {project_path}/{table_path}/table.qza '
                 f'--i-sequences {project_path}/{table_path}/representative_sequences.qza '
                 f'--output-dir {project_path}/chimera_removal')

    # Remove chimeric and borderline chimeric from table and rep_seqs
    ShellCommand('qiime feature-table filter-features '
                 f'--i-table {project_path}/{table_path}/table.qza '
                 f'--m-metadata-file {project_path}/chimera_removal/nonchimeras.qza '
                 f'--o-filtered-table {project_path}/feature_table/feature_table.qza')

    ShellCommand('qiime feature-table filter-seqs '
                 f'--i-data {project_path}/{table_path}/representative_sequences.qza '
                 f'--m-metadata-file {project_path}/chimera_removal/nonchimeras.qza '
                 f'--o-filtered-data {project_path}/feature_table/representative_sequences.qza')

    ShellCommand('qiime feature-table filter-features '
                 f'--i-table {project_path}/feature_table/feature_table.qza '
                 '--p-min-samples 2 '
                 '--p-min-frequency 100 '
                 f'--o-filtered-table {project_path}/feature_table/feature_table.qza')

    paragraph_list.append("Feature found only in one sample and those with less than 100 reads across all samples "
                          "were removed. ")

    ShellCommand('qiime feature-classifier classify-sklearn '
                 f'--i-classifier {project_path}/taxonomic_classifier/classifier.qza '
                 f'--i-reads {project_path}/feature_table/representative_sequences.qza '
                 f'--o-classification {project_path}/taxonomy/taxonomy.qza')

    paragraph_list.append("Taxonomy was assigned using the naive-bayes classifier. ")

    # Filter out mitochondria and chloroplast sequences
    ShellCommand('qiime taxa filter-table '
                 f'--i-table {project_path}/feature_table/feature_table.qza '
                 f'--i-taxonomy {project_path}/taxonomy/taxonomy.qza '
                 f'--p-exclude mitochondria,chloroplast '
                 f'--o-filtered-table {project_path}/feature_table/feature_table.qza ')

    ShellCommand('qiime taxa filter-table '
                 f'--i-table {project_path}/feature_table/feature_table.qza '
                 f'--i-taxonomy {project_path}/taxonomy/taxonomy.qza '
                 f'--p-include Bacteria '
                 f'--o-filtered-table {project_path}/feature_table/feature_table.qza ')

    # Representative sequences aswell

    ShellCommand('qiime taxa filter-seqs '
                 f'--i-sequences {project_path}/feature_table/representative_sequences.qza '
                 f'--i-taxonomy {project_path}/taxonomy/taxonomy.qza '
                 f'--p-exclude mitochondria,chloroplast '
                 f'--o-filtered-sequences {project_path}/feature_table/representative_sequences.qza ')

    ShellCommand('qiime taxa filter-seqs '
                 f'--i-sequences {project_path}/feature_table/representative_sequences.qza '
                 f'--i-taxonomy {project_path}/taxonomy/taxonomy.qza '
                 f'--p-include Bacteria '
                 f'--o-filtered-sequences {project_path}/feature_table/representative_sequences.qza ')

    paragraph_list.append("Non-bacterial ASV's, mitochondria, and chloroplast "
                          "sequences were removed using the q2-taxa -plugin. ")

    paragraph_list.append("Chimeric sequences were detected and removed using the q2-vsearch -plugin. ")

    return paragraph_list, []


@results_manager(False)
def picrust2_table(project_path):

    paragraph_list = []
    # Remove samples that have lower than 1000 frequency in total

    os.mkdir(f'{project_path}/feature_table/picrust')

    ShellCommand('qiime feature-table filter-samples '
                 f'--i-table {project_path}/feature_table/feature_table.qza '
                 '--p-min-frequency 1000 '
                 f'--o-filtered-table {project_path}/feature_table/picrust/feature_table.qza')

    paragraph_list.append("Samples that had lower than 1000 total ASV's were removed. "
                          "Metabolic pathway composition was predicted using the q2-picrust2 -plugin "
                          "with mp -method and max-nsti 2. ")

    # Save all the file names
    table, repseq = f'{project_path}/feature_table/picrust/feature_table.qza', f'{project_path}/feature_table/representative_sequences.qza'

    ShellCommand('qiime picrust2 full-pipeline '
                 f'--i-table {table} '
                 f'--i-seq {repseq} '
                 f'--output-dir {project_path}/feature_table/picrust/pathways '
                 '--p-threads 4 '
                 '--p-hsp-method mp '
                 '--p-max-nsti 2')

    # Turn to frequency and output the tsv file
    ShellCommand('qiime feature-table relative-frequency '
                 f'--i-table {project_path}/feature_table/picrust/pathways/pathway_abundance.qza '
                 f'--o-relative-frequency-table {project_path}/feature_table/picrust/pathways/rel_table.qza')

    ShellCommand('qiime tools export '
                 f'--input-path {project_path}/feature_table/picrust/pathways/rel_table.qza '
                 f'--output-path {project_path}/feature_table/picrust/export')

    ShellCommand(f'biom convert -i {project_path}/feature_table/picrust/export/feature-table.biom -o {project_path}/feature_table/picrust/feature_table.tsv --to-tsv')

    return paragraph_list, []


@results_manager(False)
def collapse_table(project_path, table_path):

    paragraph_list = []

    for level in [2, 3, 4, 5, 6, 7]:
        ShellCommand('qiime taxa collapse '
                     f'--i-table {project_path}/{table_path} '
                     f'--i-taxonomy {project_path}/taxonomy/taxonomy.qza '
                     f'--p-level {level} '
                     f'--output-dir {project_path}/feature_table/level{level}')

        # Remove samples that have lower than 1000 frequency in total

        ShellCommand('qiime feature-table filter-samples '
                     f'--i-table {project_path}/feature_table/level{level}/collapsed_table.qza '
                     '--p-min-frequency 1000 '
                     f'--o-filtered-table {project_path}/feature_table/level{level}/collapsed_table.qza')

        ShellCommand('qiime feature-table relative-frequency '
                     f'--i-table {project_path}/feature_table/level{level}/collapsed_table.qza '
                     f'--o-relative-frequency-table {project_path}/feature_table/level{level}/relative_table.qza')

        ShellCommand('qiime tools export '
                     f'--input-path {project_path}/feature_table/level{level}/relative_table.qza '
                     f'--output-path {project_path}/feature_table/level{level}/exported')

        ShellCommand('biom convert '
                     f'-i {project_path}/feature_table/level{level}/exported/feature-table.biom '
                     f'-o {project_path}/feature_table/level{level}/feature_table.tsv '
                     '--to-tsv')

    paragraph_list.append(f"Taxa feature table was collapsed to all taxonomic levels from phylum to species. ")

    return paragraph_list, []


@results_manager(False)
def convert_to_ml(project_path, target_col, pos, neg, sample_column='#SampleID'):
    """Convert both relative abundance picrust2 and taxa feature tables"""

    if not os.path.isdir(f'{project_path}/ML_tables/'):
        os.mkdir(f'{project_path}/ML_tables')

    data_type_list = [f'level{i}' for i in range(2, 8)]
    data_type_list.append('picrust')
    for data_type in data_type_list:
        df = pd.read_csv(f'{project_path}/feature_table/{data_type}/feature_table.tsv', sep='\t', header=1)
        meta = pd.read_csv(f'{project_path}/start_files/metadata.tsv', sep='\t')

        # Make it so the header is feature names and rows are samples in df
        df = df.T
        df.columns = df.iloc[0, :].values
        df = df.drop(['#OTU ID'], axis=0)

        # Convert the index col names to strings so they match with the qiime2 export
        index_names = []
        for i in meta[sample_column].values:
            index_names.append(str(i))

        meta.index = index_names
        meta = meta.loc[df.index.values]

        classes = []
        boolean_list = []
        # Take only the target samples
        for i in meta[target_col].values:
            if str(i) in [str(pos), str(neg)]:
                boolean_list.append(True)
            else:
                boolean_list.append(False)

        meta = meta.loc[boolean_list]
        df = df.loc[boolean_list]

        for i in meta[target_col].values:
            if str(i) == str(pos):
                classes.append(1)
            else:
                classes.append(0)

        df['Classes'] = classes
        print(np.unique(classes, return_counts=True))

        df.to_csv(f'{project_path}/ML_tables/{data_type}_ML.tsv', index=False, sep='\t')

    return []
