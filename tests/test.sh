#!/usr/bin/env bash

#pip install "${dir_repo}"; amalgkit metadata -h

dir_repo='/Users/kef74yk/Dropbox_w/repos/amalgkit'
dir_work="${dir_repo}/tests/test"
mkdir -p ${dir_work}
cd ${dir_work}

pip install "${dir_repo}"; \
amalgkit metadata \
--config_dir "${dir_repo}/config/test" \
--work_dir ${dir_work} \
--entrez_email 'kfuku52@gmail.com' \
--overwrite 'no'




dir_repo='/Users/kef74yk/Dropbox_w/repos/amalgkit'
dir_work="${dir_repo}/tests/vertebrate"
mkdir -p ${dir_work}

pip install "${dir_repo}"; \
amalgkit metadata \
--config_dir "${dir_repo}/config/vertebrate" \
--work_dir ${dir_work} \
--entrez_email 'kfuku52@gmail.com' \
--publication_date 1900/01/01:2019/03/12 \
--overwrite 'no'



pip install "${dir_repo}"; \
amalgkit getfastq \
--entrez_email 'kfuku52@gmail.com' \
--id 'PRJDB4514' \
--threads 8 \
--work_dir ${dir_repo}/tests/getfastq \
--save_metadata 'yes' \
--pfd 'yes' \
--max_bp '75000' \
--fastp 'yes' \
--layout 'auto' \
--remove_tmp 'yes' \
--remove_sra 'yes' \
--read_name 'trinity' \
--concat 'yes'

# a small-scale sample
pip install "${dir_repo}"; \
amalgkit getfastq \
--entrez_email 'kfuku52@gmail.com' \
--id 'SRR7764276' \
--threads 8 \
--work_dir ${dir_repo}/tests/getfastq \
--save_metadata 'yes' \
--pfd 'yes' \
--max_bp '75000' \
--fastp 'yes' \
--layout 'auto' \
--remove_tmp 'yes' \
--remove_sra 'yes' \
--concat 'no' \
--ascp 'yes' \
--ascp_key "${HOME}/Applications/Aspera Connect.app/Contents/Resources/asperaweb_id_dsa.openssh"

