from flask_geno import db
from flask_geno.models import RunJob, SequencesDataset
import microbiome_pipe as pipe
import sys
from pathlib import Path


"""
# Arguments
job_id = sys.argv[1]
user_id = sys.argv[2]
dataset_id = sys.argv[3]
"""
job_id = 11
user_id = 4
dataset_id = 10



run_job = RunJob.query.filter_by(id=job_id).first()
dataset = SequencesDataset.query.filter_by(id=job_id).first()

run_path = str(Path(f'/home/pvanni/flask_geno/flask_geno/static/job_files/user_{user_id}/job_{job_id}'))
sequences_path = str(Path(f'/home/pvanni/flask_geno/flask_geno/static/sequences/user_{user_id}/dataset_{dataset_id}'))
target_col = "diagnosis"
"""
pipe.preprocessing.import_function(project_path=run_path,
                                   path_to_sequences=sequences_path,
                                   paired=False,
                                   sample_column=run_job.name_column,
                                   manifest_type=True)

pipe.preprocessing.primer_removal(project_path=run_path, forward=None, reverse=None)

print('starting denoising')

pipe.preprocessing.denoising(project_path=run_path,
                            paired=False)

pipe.taxonomy.train_taxonomic_classifier(project_path=run_path,
                                         f_primer=None,
                                           r_primer="CCGTCAATTCMTTTRAGT",
                                           trunc_len=0,
                                           trim_left=0,
                                           rep_file="/home/pvanni/flask_geno/gg_13_8_otus/rep_set/99_otus.fasta",
                                           tax_file="/home/pvanni/flask_geno/gg_13_8_otus/taxonomy/99_otu_taxonomy.txt",
                                           database_name="GreenGenes")

pipe.feature_table.filtering_table(project_path=run_path,
                                   table_path=f"denoising/dada2_output")


pipe.feature_table.collapse_table(project_path=run_path,
                                  table_path="feature_table/feature_table.qza")

pipe.feature_table.picrust2_table(project_path=run_path)


pipe.feature_table.convert_to_ml(project_path=run_path,
                                 target_col='diagnosis',
                                 pos='Crohn’s disease',
                                 neg='Control')

pipe.analyses.diversity_methods(project_path=run_path)


pipe.analyses.alpha_boxplots(target_column=target_col,
                             project_path=run_path,
                             sample_column=run_job.name_column)
"""

pipe.analyses.beta_diversity_results(target_column=target_col,
                                     project_path=run_path,
                                     sample_column=run_job.name_column)

