from pathlib import Path
import os
import shutil


def create_job_directory(user_id, job_id):
    """Create job directory files"""
    if os.path.isdir(Path(f'/home/pvanni/flask_geno/flask_geno/static/job_files/user_{user_id}/job_{job_id}')):
        shutil.rmtree(Path(f'/home/pvanni/flask_geno/flask_geno/static/job_files/user_{user_id}/job_{job_id}'))

    os.mkdir(str(Path(f'/home/pvanni/flask_geno/flask_geno/static/job_files/user_{user_id}/job_{job_id}')))
    # start_files
    os.mkdir(str(Path(f'/home/pvanni/flask_geno/flask_geno/static/job_files/user_{user_id}/job_{job_id}/start_files')))
    # sequences
    os.mkdir(str(Path(f'/home/pvanni/flask_geno/flask_geno/static/job_files/user_{user_id}/job_{job_id}/sequences')))
    # methods
    os.mkdir(str(Path(f'/home/pvanni/flask_geno/flask_geno/static/job_files/user_{user_id}/job_{job_id}/methods')))
    # results
    os.mkdir(str(Path(f'/home/pvanni/flask_geno/flask_geno/static/job_files/user_{user_id}/job_{job_id}/results')))


def create_dataset_directory(user_id, dataset_id):
    """Create job directory files"""
    if os.path.isdir(Path(f'/home/pvanni/flask_geno/flask_geno/static/sequences/user_{user_id}/dataset_{dataset_id}')):
        shutil.rmtree(Path(f'/home/pvanni/flask_geno/flask_geno/static/sequences/user_{user_id}/dataset_{dataset_id}'))

    os.mkdir(str(Path(f'/home/pvanni/flask_geno/flask_geno/static/sequences/user_{user_id}/dataset_{dataset_id}')))


def clear_dataset_temp(user_id):
    if os.path.isdir(Path(f'/home/pvanni/flask_geno/flask_geno/static/sequences/user_{user_id}/temp')):
        shutil.rmtree(Path(f'/home/pvanni/flask_geno/flask_geno/static/sequences/user_{user_id}/temp'))

    os.mkdir(str(Path(f'/home/pvanni/flask_geno/flask_geno/static/sequences/user_{user_id}/temp')))


def clear_job_temp(user_id):
    if os.path.isdir(Path(f'/home/pvanni/flask_geno/flask_geno/static/job_files/user_{user_id}/temp')):
        shutil.rmtree(Path(f'/home/pvanni/flask_geno/flask_geno/static/job_files/user_{user_id}/temp'))

    os.mkdir(str(Path(f'/home/pvanni/flask_geno/flask_geno/static/job_files/user_{user_id}/temp')))

def submit_run_to_batch(user, job_id):

    return 0
    # Submit the run


def parse_dataset_dir(user_id, dataset_id, current_dir=None, target_extensions=['.fastq', '.fasta', '.qual']):
    """Iterative directory crawler for uploaded datasets. Will move every .fastq file to the dataset root directory"""
    dataset_root = Path(f'/home/pvanni/flask_geno/flask_geno/static/sequences/user_{user_id}/dataset_{dataset_id}')

    root_level = False
    if current_dir == None:
        root_level = True
        current_dir = dataset_root

    # unpack archives first
    current_dir_files = os.listdir(current_dir)
    for dir_or_file in current_dir_files:
        for extension in ['.zip', '.tar']:
            if extension in dir_or_file:
                shutil.unpack_archive(os.path.join(current_dir, dir_or_file), current_dir)
                os.remove(os.path.join(current_dir, dir_or_file))
                break

    # Get all the files
    current_dir_files = os.listdir(current_dir)

    for dir_or_file in current_dir_files:
        if os.path.isdir(os.path.join(current_dir, dir_or_file)):
            next_current_dir = os.path.join(current_dir, dir_or_file)
            parse_dataset_dir(user_id, dataset_id, current_dir=next_current_dir)
        else:
            remove = True
            for extension in target_extensions:
                if extension in dir_or_file:
                    if dataset_root == current_dir:
                        remove = False
                    else:
                        remove = False
                        shutil.move(os.path.join(current_dir, dir_or_file), os.path.join(dataset_root, dir_or_file))
                        break
            if remove:
                os.remove(os.path.join(current_dir, dir_or_file))

    if root_level:
        current_dir_files = os.listdir(current_dir)
        for dir_or_file in current_dir_files:
            if os.path.isdir(os.path.join(current_dir, dir_or_file)):
                shutil.rmtree(os.path.join(current_dir, dir_or_file))

def get_text_file_content(file_path):
    file_paragraphs = []
    with open(file_path) as file:
        paragraph = ''
        for line in file:

            if '\n' in line:
                paragraph += line.replace('\n', '')
                if paragraph != '':
                    file_paragraphs.append(paragraph)
                paragraph = ''
            else:
                paragraph += line

    return file_paragraphs

def parse_documentation_folder(folder_path):
    information_list = []

    files = os.listdir(folder_path)
    #print(files)
    for i in range(1, len(files) + 1):
        found_portion = False
        for file in files:
            if f'Portion {i} -' in file:
                found_portion = True
                break

        if not found_portion:
            break

        # gather information into this dictionary about the portion
        portion_dict = {}

        for file in files:
            if f'Portion {i} - Main -' in file:
                # Title
                title = file.replace(f'Portion {i} - Main - ', '')
                title = title.replace('.txt', '')
                portion_dict['title'] = title
                # Main text that includes picture names and legends
                portion_dict['main_text'] = get_text_file_content(file_path=f'{folder_path}/{file}')

        information_list.append(portion_dict)

    return information_list