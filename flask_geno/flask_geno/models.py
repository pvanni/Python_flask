from flask_geno import db, login_manager
from datetime import datetime
from flask_login import UserMixin


@login_manager.user_loader
def load_user(user_id):
    return User.query.get(int(user_id))


class User(db.Model, UserMixin):
    id = db.Column(db.Integer, primary_key=True)
    email = db.Column(db.String(120), unique=True, nullable=False)
    runs_submitted = db.Column(db.Integer, nullable=False, default=0)
    date_joined = db.Column(db.DateTime, nullable=False, default=datetime.utcnow)
    password = db.Column(db.String(60), nullable=False)
    runs = db.relationship('RunJob', backref='author', lazy=True)
    datasets = db.relationship('SequencesDataset', backref='author', lazy=True)

    def __repr__(self):
        return f"User('Email: {self.email}, ID: {self.id}, Runs submitted: {self.runs_submitted}, " \
               f"Date joined: {self.date_joined}')"


class RunJob(db.Model):
    # Meta variables
    id = db.Column(db.Integer, primary_key=True)
    user_id = db.Column(db.Integer, db.ForeignKey('user.id'), nullable=False)
    date_created = db.Column(db.DateTime, nullable=False, default=datetime.utcnow)
    job_title = db.Column(db.String(100), nullable=False, default='New Run')
    job_status = db.Column(db.String(100), nullable=False, default='Unfinished')

    # Job files and other related information
    sequences_format = db.Column(db.String(50), nullable=False, default='No sequences format')
    job_description = db.Column(db.Text, nullable=False, default='No description')
    metadata_file = db.Column(db.String(100), nullable=False, default='No metadata file')
    run_dataset_id = db.Column(db.Integer)
    f_primer = db.Column(db.String(100))
    r_primer = db.Column(db.String(100))
    name_column = db.Column(db.String(100))


    # Preprocessing and methods texts for results view
    importing_text = db.Column(db.Text)
    primer_removal_text = db.Column(db.Text)
    denoising_text = db.Column(db.Text)
    taxonomic_classifier_text = db.Column(db.Text)
    filtering_feature_table = db.Column(db.Text)
    collapsing_text = db.Column(db.Text)
    picrust2_text = db.Column(db.Text)


class SequencesDataset(db.Model):
    # Meta variables
    id = db.Column(db.Integer, primary_key=True)
    user_id = db.Column(db.Integer, db.ForeignKey('user.id'), nullable=False)
    name = db.Column(db.String(250), nullable=False)
    description = db.Column(db.Text)
    sequences_format = db.Column(db.String(100))
    paired_end = db.Column(db.Boolean, nullable=False, default=False)
    data_uploaded = db.Column(db.Boolean, nullable=False, default=False)
