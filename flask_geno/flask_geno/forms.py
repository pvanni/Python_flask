from flask_wtf import FlaskForm
from flask_wtf.file import FileField, FileAllowed, FileRequired
from flask_login import current_user
from wtforms import StringField, PasswordField, SubmitField, BooleanField, SelectField
from wtforms.validators import DataRequired, Email, ValidationError
from flask_geno.models import User


class LoginForm(FlaskForm):

    email = StringField('Email', validators=[DataRequired(), Email()])
    password = PasswordField('Password', validators=[DataRequired()])
    remember = BooleanField('Remember Me')
    submit = SubmitField('Login')


class UpdateAccountForm(FlaskForm):

    email = StringField('Email', validators=[DataRequired(), Email()])
    submit = SubmitField('Update')

    def validate_email(self, email):
        if email.data != current_user.email:
            user = User.query.filter_by(email=email.data).first()
            if user:
                raise ValidationError('That email is taken. Please choose a different one.')


class MetadataForm(FlaskForm):
    """For creating a new run job"""
    metadata_file = FileField('Metadata file', validators=[FileRequired(), FileAllowed(['xlsx', 'tsv', 'csv'])])
    submit = SubmitField('Preview metadata')


class RunDataForm(FlaskForm):
    """For creating a new run job"""
    name_column = StringField('Sample identifier column', validators=[DataRequired()])
    job_name = StringField('Run name', validators=[DataRequired()])
    job_description = StringField('Run description', validators=[DataRequired()])
    f_primer = StringField('Sequencing Forward primer', validators=[DataRequired()])
    r_primer = StringField('Sequencing Reverse primer', validators=[DataRequired()])
    target_column = StringField('Target variable column', validators=[DataRequired()])
    dataset_choice = SelectField('Dataset')
    submit = SubmitField('Submit run')


class NewDatasetForm(FlaskForm):
    """For the first part of creating a new dataset entry. Data upload is handled by dropzone.js"""
    name = StringField('Dataset name', validators=[DataRequired()])
    description = StringField('Description')
    sequences_format = BooleanField('Paired-end sequences')
    submit = SubmitField('To data upload')


class SubmitForm(FlaskForm):
    """General submit button for everything"""
    submit = SubmitField('Submit')

