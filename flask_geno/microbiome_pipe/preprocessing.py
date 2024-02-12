import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from microbiome_pipe.util import ShellCommand, results_manager


@results_manager(False)
def import_function(project_path, path_to_sequences, manifest_type=True, paired=False, error_rate=0, sample_column='#SampleID'):

    """Takes in either multiplexed or demultiplexed formats

    Requirements in manifest_type=True:
    There needs to be a metadata file where #SampleID is the column name for identifiers in start_files called
    "metadata.tsv".
    Each sample must have an {id}.fastq or {id}.fastq.gz file in the .../sequences directory
    If theres no fastq.gz or fastq file for any given sample found in the metadata, then the sample is dropped

    If manifest_type=False, then its assumed for sequences to be Single End demultiplexed format in a single file
    with a barcode sequence. Barcodes are assumed to be found in BarcodeSequence column"""

    paragraph_list = []
    picture_list = []

    # Some of the qiime commands require that there is a specifically called "sample id" column as the first column
    meta = pd.read_csv(f'{project_path}/start_files/start_metadata.tsv', sep='\t')
    columns = list(meta.columns.values)
    columns.remove('#SampleID')
    if sample_column == "#SampleID":
        meta = meta[['#SampleID'] + columns]
    else:
        meta['#SampleID'] = meta[sample_column].copy()
        meta = meta[['#SampleID'] + columns]

    if manifest_type:

        sample_names = list(meta[sample_column].values)

        file_path = path_to_sequences

        # Create for parsing the right samples
        file_list = os.listdir(path_to_sequences)

        # Create the absolute filepath variable
        no_file_list = []
        file_path_list = []
        reverse_path_list = []
        if paired:
            for i in sample_names:
                for orientation in [1, 2]:
                    if os.path.isfile(f'{file_path}/{i}_{orientation}.fastq.gz'):
                        ShellCommand(f'gunzip {file_path}/{i}_{orientation}.fastq.gz')
                    if os.path.isfile(f'{file_path}/{i}_{orientation}.fastq'):
                        if orientation == 1:
                            file_path_list.append(f'{file_path}/{i}_{orientation}.fastq')
                        if orientation == 2:
                            reverse_path_list.append(f'{file_path}/{i}_{orientation}.fastq')
                    else:
                        no_file_list.append(i)
        else:
            for i in sample_names:

                found_file = False
                for file in file_list:
                    if i in file:
                        if not file[file.index(i)+len(i)].isdigit():
                            if '.fastq.gz' in i:
                                ShellCommand(f'gunzip {file_path}/{file}')
                            if os.path.isfile(f'{file_path}/{file}'):
                                file_path_list.append(f'{file_path}/{file}')
                                found_file = True

                if not found_file:
                    no_file_list.append(i)





        meta.index = meta[sample_column].values

        # Report which samples were not found as fastq files and how many there were at the start
        paragraph_list.append(f'There were {meta.shape[0]} samples '
                              f'in the original metadata and {meta.shape[1]} columns. '
                              f'{len(no_file_list)} samples didnt have a fastq file in the sequences/ directory. '
                              f'The following sample IDs were removed from metadata.tsv:')

        paragraph_list.append(f'{no_file_list}')

        meta.drop(index=no_file_list, inplace=True)

        if paired:
            manifest = pd.DataFrame(data={'sample-id': meta[sample_column].values,
                                          'forward-absolute-filepath': file_path_list,
                                          'reverse-absolute-filepath': reverse_path_list},
                                    columns=['sample-id', 'forward-absolute-filepath', 'reverse-absolute-filepath'])
        else:
            manifest = pd.DataFrame(data={'sample-id': meta[sample_column].values, 'absolute-filepath': file_path_list},
                                    columns=['sample-id', 'absolute-filepath'])

        manifest.to_csv(f'{project_path}/sequences/manifest.tsv', index=False, header=True, sep='\t')
        meta.to_csv(f'{project_path}/start_files/metadata.tsv', index=False, header=True, sep='\t')

        if paired:
            ShellCommand('qiime tools import '
                         "--type 'SampleData[PairedEndSequencesWithQuality]' "
                         f"--input-path {project_path}/sequences/manifest.tsv "
                         f"--output-path {project_path}/sequences/demultiplexed_sequences.qza "
                         "--input-format PairedEndFastqManifestPhred33V2")
        else:
            ShellCommand('qiime tools import '
                         "--type 'SampleData[SequencesWithQuality]' "
                         f"--input-path {project_path}/sequences/manifest.tsv "
                         f"--output-path {project_path}/sequences/demultiplexed_sequences.qza "
                         "--input-format SingleEndFastqManifestPhred33V2")

        paragraph_list.append("Sequences were imported into Qiime 2 ")

    else:
        # Compress
        ShellCommand(f'gzip {path_to_sequences}/sequences.fastq')

        # Import
        ShellCommand('qiime tools import '
                     '--type MultiplexedSingleEndBarcodeInSequence '
                     f'--input-path {project_path}/sequences.fastq.gz '
                     f'--output-path {project_path}/sequences/sequences.qza')

        paragraph_list.append("Multiplexed sequences were imported into Qiime2. ")

        ShellCommand('qiime cutadapt demux-single '
                     '--i-seqs sequences/sequences.qza '
                     f'--m-barcodes-file {project_path}/start_files/metadata.tsv '
                     '--m-barcodes-column BarcodeSequence '
                     f'--p-error-rate {error_rate} '
                     f'--output-dir {project_path}/sequences/demultiplexed_sequences.qza')

        paragraph_list.append("Barcode sequences were removed using the q2-cutadapt -plugin. "
                              f"Error rate parameter was set to {error_rate}, while "
                              "other parameters were set to default values. ")

    return paragraph_list, picture_list


@results_manager(False)
def primer_removal(forward, project_path, reverse=None, error_rate=0):

    """Remove primers and create a report + statistics showing how many sequences got shorter"""
    paragraph_list = []
    picture_list = []

    if forward:
        if not os.path.isfile(f'sequences/trimmed_demultiplexed_sequences.qza'):
            if reverse:
                ShellCommand('qiime cutadapt trim-paired '
                             f'--i-demultiplexed-sequences {project_path}/sequences/demultiplexed_sequences.qza '
                             '--p-cores 4 '
                             f'--p-front-f {forward} '
                             f'--p-front-r {reverse} '
                             f'--p-error-rate {error_rate} '
                             f'--o-trimmed-sequences {project_path}/sequences/trimmed_demultiplexed_sequences.qza')
            else:
                ShellCommand('qiime cutadapt trim-single '
                             f'--i-demultiplexed-sequences {project_path}/sequences/demultiplexed_sequences.qza '
                             '--p-cores 4 '
                             f'--p-front {forward} '
                             f'--p-error-rate {error_rate} '
                             f'--o-trimmed-sequences {project_path}/sequences/trimmed_demultiplexed_sequences.qza')

            paragraph_list.append(f"Primer sequences (f: {forward}, r: {reverse}) were removed using the "
                                  f"q2-cutadapt -plugin with "
                                  f"error rate parameter set to {error_rate}. Other parameters were set to default values. ")

    else:
        ShellCommand(f'cp {project_path}/sequences/demultiplexed_sequences.qza '
                     f'{project_path}/sequences/trimmed_demultiplexed_sequences.qza')

        paragraph_list.append(f"No primers were removed. ")

    # Export the sequences before primer trimming
    if not os.path.isdir(f'{project_path}/sequences/demultiplexed_sequences_fastq/'):
        ShellCommand('qiime tools export '
                     f'--input-path {project_path}/sequences/demultiplexed_sequences.qza '
                     f'--output-path {project_path}/sequences/demultiplexed_sequences_fastq')

    if not os.path.isdir(f'{project_path}/sequences/trimmed_demultiplexed_sequences_fastq/'):
        # Export the sequences after primer trimming
        ShellCommand('qiime tools export '
                     f'--input-path {project_path}/sequences/trimmed_demultiplexed_sequences.qza '
                     f'--output-path {project_path}/sequences/trimmed_demultiplexed_sequences_fastq')

    def export_for_quality(before_after):
        # List of files on the "before" directory
        file_before = os.listdir(f'{project_path}/sequences/{before_after}demultiplexed_sequences_fastq/')
        for file in file_before:
            if 'fastq.gz' in file:
                ShellCommand(f'gunzip {project_path}/sequences/{before_after}demultiplexed_sequences_fastq/{file}')

        # Create a list of fastq files in the directory
        file_before = os.listdir(f'{project_path}/sequences/{before_after}demultiplexed_sequences_fastq/')
        before_fastq_list = []
        for file in file_before:
            if 'fastq' in file:
                before_fastq_list.append(f'{file}')

        return before_fastq_list

    # generate a list of fastq files for both before and after trimming
    after_list = export_for_quality('trimmed_')

    def sequence_parser(file_list, path):
        """Parse all the sequences to gather information for quality control"""
        fastq_data = []
        for file_path in file_list:
            name_list = []
            sequence_list = []
            quality_list = []
            length_list = []
            with open(f'{path}{file_path}') as file:
                counter = 1
                for new_line in file:
                    line = new_line.replace("\n", "")
                    if counter == 1:
                        name_list.append(line)
                    if counter == 2:
                        sequence_list.append(line)
                        length_list.append(len(line))
                    if counter == 4:
                        quality_list.append(line)
                        counter = 1
                        continue
                    counter += 1
            fastq_data.append([name_list, sequence_list, quality_list, length_list])

        return fastq_data

    # A list of fastq sequences with quality in the same order as the input list of fastq file paths were
    if len(after_list) > 50:
        after_data = sequence_parser(after_list[:10], f'{project_path}/sequences/trimmed_demultiplexed_sequences_fastq/')
    else:
        after_data = sequence_parser(after_list, f'{project_path}/sequences/trimmed_demultiplexed_sequences_fastq/')

    all_quality_scores = []
    for sample in after_data:
        for quality in sample[2]:
            all_quality_scores.append(quality)

    quality_sample = np.random.choice(all_quality_scores, size=10000, replace=False)
    quality_sample_scores = []
    max_len = 0
    # Transform to quality scores
    for sample in quality_sample:
        q_scores = []
        # Transfer the ascii symbols to PHRED33 quality scores
        for letter in sample:
            q_scores.append(ord(letter) - 33)
        # Record the max lenght needed for plotting
        if len(q_scores) > max_len:
            max_len = len(q_scores)
        quality_sample_scores.append(q_scores)

    # Bin the data into each boxplot
    box_data = [[] for i in range(max_len)]

    for sample in quality_sample_scores:
        for i, value in enumerate(sample):
            box_data[i].append(value)

    # plot the data. Find out in what format the picture needs to be returned for the report
    fig, ax = plt.subplots(figsize=(20.0, 10.0))
    ax.boxplot(box_data, showfliers=False, showcaps=False)
    ax.set_xticks([i for i in np.arange(0, max_len, 20)])
    ax.set_xticklabels([i for i in np.arange(0, max_len, 20)])

    fig.savefig(f"{project_path}/sequence_quality.png")
    plt.close('all')
    picture_list.append(f"{project_path}/sequence_quality.png")

    return paragraph_list, picture_list, 0


@results_manager(False)
def denoising(project_path, paired=False):
    """Denoise the data based on the trunc len found in primer_removal function.
    If more than 10% of the samples are lost during denoising, reduce trunc_len size by 5 and repeat the process"""
    paragraph_list = []
    picture_list = []

    def dada2_denoising(paired=False):
        """Run the dada2"""

        if not paired:
            ShellCommand("qiime dada2 denoise-single "
                         f"--i-demultiplexed-seqs {project_path}/sequences/trimmed_demultiplexed_sequences.qza "
                         f"--p-trim-left 0 "
                         f"--p-trunc-len 0 "
                         f"--p-n-threads 1 "
                         f"--output-dir {project_path}/denoising/dada2_output")
        else:
            if not os.path.isdir(f'{project_path}/denoising/dada2_output/'):
                ShellCommand("qiime dada2 denoise-paired "
                             f"--i-demultiplexed-seqs {project_path}/sequences/trimmed_demultiplexed_sequences.qza "
                             f"--p-trim-left-f 0 "
                             f"--p-trunc-len-f 0 "
                             f"--p-trim-left-r 0 "
                             f"--p-trunc-len-r 0 "
                             f"--p-n-threads 1 "
                             f"--output-dir {project_path}/denoising/dada2_output")

    paragraph_list.append("Sequences were denoised to ASV's using the q2-dada2 -plugin. ")

    dada2_denoising(paired=paired)

    return paragraph_list, []
