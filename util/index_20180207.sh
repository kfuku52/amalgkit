#!/bin/sh
#$ -S /bin/bash
#$ -cwd
#$ -l s_vmem=16G
#$ -l mem_req=16G
#$ -pe def_slot 1
#$ -t 1-95

# 1-95

source /home/kfuku/.bashrc
ulimit -s unlimited

echo running on `hostname`
echo starting at `date`

dir_mask_fasta='/lustre1/home/kfuku/my_db/Ensembl/release-91/mask/fasta'
dir_cdna='/lustre1/home/kfuku/my_db/Ensembl/release-91/cdna'
dir_cdna_mask='/lustre1/home/kfuku/my_db/Ensembl/release-91/cdna.mask'
dir_genome='/lustre1/home/kfuku/my_db/Ensembl/release-91/genome'
#wget -r -nd -np -A .cdna.all.fa.gz ftp://ftp.ensembl.org/pub/release-91/fasta/
#wget -r -nd -np -A .cds.all.fa.gz ftp://ftp.ensembl.org/pub/release-91/fasta/
#wget -r -nd -np -A .gtf.gz ftp://ftp.ensembl.org/pub/release-91/gtf/

files=( `ls ${dir_cdna}/*.fa.gz` )
sci_name=`basename ${files[$[${SGE_TASK_ID}-1]]} | sed -e "s/\..*//" | sed -e "s/_/PLACEHOLDER/" | sed -e "s/_.*//g" | sed -e "s/PLACEHOLDER/_/"`
echo ${sci_name}

cd ${dir_cdna}
fasta=`ls ${sci_name}*.fa.gz`
if [ -s ${fasta}.index ]; then
	echo ${fasta} is already indexed, skipped.
else	
	echo indexing ${fasta}
	kallisto index --index=${fasta}.index --kmer-size=31 ${fasta}
fi

cd ${dir_cdna_mask}
fasta=`ls ${sci_name}*.fasta`
if [ -s ${fasta}.index ]; then
	echo ${fasta} is already indexed, skipped.
else	
	echo indexing ${fasta}
	kallisto index --index=${fasta}.index --kmer-size=31 ${fasta}
fi

#cd ${dir_genome}
#fasta=`ls ${sci_name}*.fa.gz`
#if [ -e ${fasta}.index ]; then
#	echo ${fasta} is already indexed, skipped.
#else	
#	echo indexing $fasta
#	star \
#	--runThreadN ${NSLOTS} \
#	--runMode genomeGenerate \
#	--genomeDir
#fi

###################
echo ending at `date`

