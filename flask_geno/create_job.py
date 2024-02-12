from flask_geno import bcrypt, db
from flask_geno.models import User, RunJob
import sys


user = User.query.filter_by(email='admin@genobiomics.com').first()

for title, status in zip(['Mycobiome run #1', 'Saliva Microbiome run // Thesis', 'Thanatomicrobiome in Oulu',
                          'Plant Microbiome analyses // Gold'], ['Complete', 'Queue', 'Running', 'Complete']):
    run_job = RunJob(author=user,
                     job_title=title,
                     sequences_format='Paired-end Manifest',
                     job_description='Description: An unformatted summary '
                                      'describing the purpose, nature, and scope of the data collection, special '
                                      'characteristics of its contents, major subject areas covered, and what questions '
                                      'the PIs attempted to answer when they conducted the study. A listing of major '
                                      'variables in the study is important here.',
                     metadata_file='test_metadata.xlsx',
                     sequences_file='test_sequences.gz',
                     job_status=status)

    db.session.add(run_job)
    db.session.commit()