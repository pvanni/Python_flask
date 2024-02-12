from flask_geno import bcrypt, db
from flask_geno.models import User
import sys, os
from pathlib import Path


email = sys.argv[1]
password = sys.argv[2]
reset = sys.argv[3]

if reset == 'Drop':
    db.drop_all()
    db.create_all()

hashed_pw = bcrypt.generate_password_hash(password).decode('utf-8')

new_user = User(email=email, password=hashed_pw)

db.session.add(new_user)
db.session.commit()

os.mkdir(Path(f'/home/pvanni/flask_geno/flask_geno/static/job_files/user_{new_user.id}'))
os.mkdir(Path(f'/home/pvanni/flask_geno/flask_geno/static/sequences/user_{new_user.id}'))
os.mkdir(Path(f'/home/pvanni/flask_geno/flask_geno/static/job_files/user_{new_user.id}/temp'))
os.mkdir(Path(f'/home/pvanni/flask_geno/flask_geno/static/sequences/user_{new_user.id}/temp'))