from flask import render_template, url_for, flash, redirect, request, make_response
from flask_geno.models import User, RunJob, SequencesDataset
from flask_geno.forms import LoginForm, UpdateAccountForm, RunDataForm, NewDatasetForm, SubmitForm, MetadataForm
from flask_geno import app, db, bcrypt
from flask_login import login_user, current_user, logout_user, login_required
from flask_geno.util import create_dataset_directory, clear_dataset_temp, \
                            create_job_directory, clear_job_temp, parse_dataset_dir, parse_documentation_folder
from werkzeug.utils import secure_filename
from pathlib import Path
import os
import pandas as pd
import shutil

@app.context_processor
def inject_documentation():
    doc_path = Path(f'/home/pvanni/flask_geno/flask_geno/static/documentation')
    title_files = os.listdir(doc_path)

    title_list = []
    subtitle_list = []
    # Loop through all titles
    for i in range(1, len(title_files) + 1):
        for file in title_files:
            title_dict = {}
            title_dict['subtitle_list'] = []
            if f'Title {i} - ' in file:
                title = file.replace(f'Title {i} - ', '')
                #title = title.replace('.txt', '')
                title_dict['title'] = title
                title_dict['title_file'] = file

                title_list.append(title_dict)
                # Loop through all subtitles
                current_subtitles = []
                subtitle_files = os.listdir(Path(f'/home/pvanni/flask_geno/flask_geno/static/documentation/{file}'))
                for j in range(1, len(title_files) + 1):
                    for subfile in subtitle_files:
                        subtitle_dict = {}
                        if f'Subtitle {i} - ' in subfile:
                            subtitle = subfile.replace(f'Subtitle {i} - ', '')
                            #subtitle = subtitle.replace('.txt', '')
                            subtitle_dict['subtitle'] = subtitle
                            subtitle_dict['sub_file'] = subfile

                            current_subtitles.append(subtitle_dict)

                subtitle_list.append(current_subtitles)

    return dict(title_list=title_list, subtitle_list=subtitle_list)

@app.context_processor
def inject_enumerate():
    return dict(enumerate=enumerate)

@app.route('/info')
@app.route('/')
@login_required
def info():
    return render_template('info.html', title=info)


@app.route('/new_run', methods=['GET', 'POST'])
@login_required
def new_run():
    metadata_form = MetadataForm()

    if metadata_form.validate_on_submit():
        metadata_file = metadata_form.metadata_file.data
        new_name = os.path.splitext(metadata_file.filename)
        file_path = Path(f'/home/pvanni/flask_geno/flask_geno/static/job_files/user_{current_user.id}/temp/' + 'start_metadata' + new_name[1])
        metadata_file.save(file_path)
        run_job = RunJob(author=current_user,
                         metadata_file=str(file_path))
        db.session.add(run_job)
        db.session.commit()
        return redirect(url_for('new_run_fill', runjob_id=run_job.id))

    return render_template('new_run.html', metadata_form=metadata_form)


@app.route('/new_run_fill/<int:runjob_id>', methods=['GET', 'POST'])
@login_required
def new_run_fill(runjob_id):

    run_job = RunJob.query.get_or_404(runjob_id)
    form = RunDataForm()
    user_datasets = SequencesDataset.query.filter_by(author=current_user).order_by(SequencesDataset.id.desc())
    form.dataset_choice.choices = [(dataset.id, dataset.name) for dataset in user_datasets]

    df = pd.read_csv(run_job.metadata_file, sep='\t')

    if form.validate_on_submit():
        # Move the metadata file from temp to job files
        create_job_directory(current_user.id, run_job.id)
        new_path = run_job.metadata_file.replace('temp', f'job_{run_job.id}/start_files')
        shutil.move(run_job.metadata_file, new_path)

        # Update the run_job database entry
        run_job.metadata_file = new_path
        run_job.name_column = form.name_column.data
        run_job.job_title = form.job_name.data
        run_job.job_description = form.job_description.data
        run_job.f_primer = form.f_primer.data
        run_job.r_primer = form.r_primer.data
        run_job.run_dataset_id = form.dataset_choice.data
        run_job.job_status = 'Running'
        db.session.add(run_job)
        db.session.commit()

        # Clear temp dir
        clear_job_temp(current_user.id)

        return redirect(url_for('user_home'))

    return render_template('new_run_fill.html', metadata_form=form, df=df)


@app.route('/home')
@login_required
def user_home():
    """This template will have a list of all the run jobs and info about them"""
    run_jobs = RunJob.query.filter_by(author=current_user).order_by(RunJob.id.desc())
    return render_template('user_home.html', run_jobs=run_jobs)

@app.route('/docs/<title>/<subtitle>')
@login_required
def doc_page(title, subtitle):
    """Goes into the documentations folder and finds the correct materials and passes it to the template"""
    # List all topics documentation titles
    doc_path = Path(f'/home/pvanni/flask_geno/flask_geno/static/documentation/{title}/{subtitle}')
    subject_dict = parse_documentation_folder(doc_path)
    return render_template('documentation.html', subject_dict=parse_documentation_folder(doc_path),
                           file_path=str(Path(f'documentation/{title}/{subtitle}/')) + '/')


@app.route('/login', methods=['GET', 'POST'])
def login():
    if current_user.is_authenticated:
        return redirect(url_for('new_run'))
    form = LoginForm()
    if form.validate_on_submit():
        user = User.query.filter_by(email=form.email.data).first()
        if user and bcrypt.check_password_hash(user.password, form.password.data):
            login_user(user, remember=form.remember.data)
            flash('You have succesfully logged in!', 'success')
            return redirect(url_for('user_home'))
        else:
            return render_template('login.html', title='Login', form=form, error_message='Login failed, try again')
    return render_template('login.html', title='Login', form=form)


@app.route('/logout')
def logout():
    logout_user()
    flash('You have logged out', 'success')
    return redirect(url_for('info'))


@app.route('/account', methods=['GET', 'POST'])
@login_required
def account():
    form = UpdateAccountForm()
    if form.validate_on_submit():
        if form.email.data not in User.query.all():
            current_user.email = form.email.data
            db.session.commit()
            return redirect(url_for('new_run'))
        else:
            return render_template('account.html', form=form, error_message='Invalid email')

    return render_template('account.html', form=form)


@app.route("/run/<int:runjob_id>")
@login_required
def run_results(runjob_id):
    """
    Splice the results page into sections. Sections contain paragraphs of text, pictures and datatables

        Input the picture paths and texts into the template in a dictionary

        [
            example = {
                'data_type': str,
                'target_variable': str'
                'heading': str,
                'paragraphs': list of str's,
                'pictures': list of paths in str,
                'datatables': dictionary of 'dataframe': pandas.DataFrames and 'df_id': id for jquery

            }
        ]
    """
    run_job = RunJob.query.get_or_404(runjob_id)


    example_results = [{
        'data_type': 'level6',
        'target_variable': 'diagnosis',
        'heading': 'Differential abundance results - ALDEx2',
        'paragraphs': ['The learning phase and the subsequent prediction of machine learning algorithms can be affected by the problem of imbalanced data set. '
                       'The balancing issue corresponds to the difference of the number of samples in '
                       'the different classes. We illustrate the effect of training a linear SVM classifier with different level of class balancing.'],
        'pictures': ['job_files/user_1/job_4/results/diagnosis_Alpha_boxplots.png'],
        'datatables': [{'dataframe': pd.read_csv('/home/pvanni/flask_geno/flask_geno/static/job_files/user_1/job_4/aldex2/level6/diagnosis_diff/differentials.tsv', sep='\t'),
                       'df_id': 1}]

    }]

    if run_job.author == current_user:
        return render_template('run_results.html', results_data=example_results)
    else:
        return redirect(url_for('user_home'))


@app.route('/sequence_uploader/<int:sequencesdataset_id>', methods=['GET', 'POST'])
@login_required
def sequence_uploader(sequencesdataset_id):

    dataset = SequencesDataset.query.get_or_404(sequencesdataset_id)
    if dataset.data_uploaded:
        flash('Data has already been added for this dataset', 'danger')
        return redirect(url_for('user_home'))
    form = SubmitForm()
    if form.validate_on_submit():
        dataset.data_uploaded = True
        db.session.commit()

        create_dataset_directory(current_user.id, sequencesdataset_id)
        shutil.move(str(Path(f'/home/pvanni/flask_geno/flask_geno/static/sequences/user_{current_user.id}/temp')),
                    str(Path(f'/home/pvanni/flask_geno/flask_geno/static/sequences/user_{current_user.id}/dataset_{sequencesdataset_id}')))
        clear_dataset_temp(current_user.id)

        # Parse the dataset directory after uploading
        parse_dataset_dir(user_id=current_user.id,
                          dataset_id=dataset.id)

        flash('Dataset created!', 'success')
        return redirect(url_for('user_home'))

    if dataset.author == current_user:
        return render_template('sequence_uploader.html', dataset=dataset, form=form)
    else:
        flash("You dont have access to this dataset", 'danger')
        return redirect(url_for('user_home'))


@app.route('/upload_seq/<int:sequencesdataset_id>', methods=['POST'])
@login_required
def upload_seq(sequencesdataset_id):

    file = request.files['file']

    save_path = Path(f'/home/pvanni/flask_geno/flask_geno/static/sequences/user_{current_user.id}/temp/' + secure_filename(file.filename))
    current_chunk = int(request.form['dzchunkindex'])

    if os.path.exists(save_path) and current_chunk == 0:
        # 400 and 500s will tell dropzone that an error occurred and show an error
        return make_response(('File already exists', 400))

    try:
        with open(save_path, 'ab+') as f:
            f.seek(int(request.form['dzchunkbyteoffset']))
            f.write(file.stream.read())
    except OSError:
        return make_response(("Not sure why,"
                              " but we couldn't write the file to disk", 500))

    total_chunks = int(request.form['dztotalchunkcount'])

    if current_chunk + 1 == total_chunks:
        # This was the last chunk, the file should be complete and the size we expect
        if os.path.getsize(save_path) != int(request.form['dztotalfilesize']):
            return make_response(('Size mismatch', 500))

    return make_response(("Chunk upload successful", 200))


@app.route("/new_dataset", methods=['GET', 'POST'])
@login_required
def new_dataset():
    form = NewDatasetForm()
    if form.validate_on_submit():
        dataset = SequencesDataset(name=form.name.data,
                                   description=form.description.data,
                                   author=current_user,
                                   sequences_format=form.sequences_format.data)
        db.session.add(dataset)
        db.session.commit()
        clear_job_temp(current_user.id)
        create_dataset_directory(current_user.id, dataset.id)
        return redirect(url_for('sequence_uploader', sequencesdataset_id=dataset.id))
    return render_template('new_dataset.html', form=form)
