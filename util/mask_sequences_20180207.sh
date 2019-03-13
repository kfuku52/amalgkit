#!/bin/sh
#$ -S /bin/bash
#$ -cwd
#$ -l s_vmem=8G
#$ -l mem_req=8G
#$ -pe def_slot 1
#$ -t 1-95

# 1-95

source /lustre1/home/kfuku/.bashrc
ulimit -s unlimited

echo running on `hostname`
echo starting at `date`

#wget -r -c -nd -np -A .gtf.gz ftp://ftp.ensembl.org/pub/release-91/gtf/
#wget -r -c -nd -np -A .dna.toplevel.fa.gz ftp://ftp.ensembl.org/pub/release-91/fasta/
dir_work=/lustre1/home/kfuku/my_project/convergence_duplication/20180207_kallisto
dir_gtf=/lustre1/home/kfuku/my_db/Ensembl/release-91/gtf
dir_cdna=/lustre1/home/kfuku/my_db/Ensembl/release-91/cdna
dir_cds=/lustre1/home/kfuku/my_db/Ensembl/release-91/cds
dir_genome=/lustre1/home/kfuku/my_db/Ensembl/release-91/genome
dir_mask=/lustre1/home/kfuku/my_db/Ensembl/release-91/mask

cd ${dir_work}

ext=".gtf.gz"
files=( `ls ${dir_gtf}/*${ext}` )
file=${files[$[${SGE_TASK_ID}-1]]}
bn=`basename ${file} .gtf.gz`
sci_name=`echo ${bn} | sed -e "s/\..*//"`
cdna=`ls ${dir_cdna}/${sci_name}*.fa.gz`
genome=`ls ${dir_genome}/${sci_name}*.fa.gz`
echo file: ${file}
echo basename: ${bn}
echo sci_name: ${sci_name}
echo cdna: ${cdna}
echo genome: ${genome}

if [ ! -e ${dir_work}/${bn} ]; then
	mkdir ${dir_work}/${bn}
fi
cd ${dir_work}/${bn}

if [ -e ${dir_mask}/fasta/${bn}.mask.fasta.1.bt2 ]; then
	echo ${file} is already processed, skipped.
else	
	python /lustre1/home/kfuku/my_script/ensembl_gtf2geneid.py \
	--gtf ${file} \
	--search_attrs 'gene_biotype|transcript_name' \
	--search_values 'lincRNA|macro_lncRNA|miRNA|misc_RNA|Mt_rRNA|Mt_tRNA|rRNA|scaRNA|scRNA|snoRNA|snRNA|sRNA|tRNA|mt-.*' \
	--out_attr 'transcript_id'

	gunzip -c ${cdna} | sed -e "s/\..*//" -e "s/[[:space:]].*//" > tmp.cdna.fasta
	fatt extract --file ensembl_gtf2geneid.out_attr.tsv tmp.cdna.fasta > ${dir_mask}/fasta.cdna_overlap/${bn}.mask.cdna_overlap.fasta

	gunzip -c ${genome} > tmp.genome.fasta
	bedtools getfasta \
	-s \
	-fi tmp.genome.fasta \
	-bed ensembl_gtf2geneid.gtf \
	-fo ${dir_mask}/fasta/${bn}.mask.fasta

	bowtie2-build -f ${dir_mask}/fasta/${bn}.mask.fasta ${dir_mask}/fasta/${bn}.mask.fasta

	cp ensembl_gtf2geneid.gtf ${dir_mask}/gtf/${bn}.mask.gtf
	cp ensembl_gtf2geneid.out_attr.tsv ${dir_mask}/transcript_id/${bn}.mask.txt
fi

if [ ! -e ${dir_cdna}.mask/${sci_name}.cdna.mask.fasta ]; then
	transcript_id=`ls ${dir_mask}/transcript_id/${sci_name}.*`
	cdna=`ls ${dir_cdna}/${sci_name}.*.fa.gz`
	zcat ${cdna} | sed -e "s/[[:space:]].*//" > ${sci_name}.cdna.id_formatted.fasta
	echo cDNA: num line before mask: `wc -l ${sci_name}.cdna.id_formatted.fasta`
	python -c "import sys,re; keys = open(sys.argv[1]).read().split('\n'); entries = open(sys.argv[2]).read().split('>'); [ sys.stdout.write('>'+e) for e in entries if not re.sub('\\\..*','',e,flags=re.DOTALL) in keys ]" \
	${transcript_id} ${sci_name}.cdna.id_formatted.fasta > ${sci_name}.cdna.mask.fasta
	echo cDNA: num line after mask: `wc -l ${sci_name}.cdna.mask.fasta`
	cp ${sci_name}.cdna.mask.fasta ${dir_cdna}.mask/${sci_name}.cdna.mask.fasta
fi

if [ ! -e ${dir_cds}.mask/${sci_name}.cds.mask.fasta ]; then
	transcript_id=`ls ${dir_mask}/transcript_id/${sci_name}.*`
	cds=`ls ${dir_cds}/${sci_name}.*.fa.gz`
	zcat ${cds} | sed -e "s/[[:space:]].*//" > ${sci_name}.cds.id_formatted.fasta
	echo CDS: num line before mask: `wc -l ${sci_name}.cds.id_formatted.fasta`
	python -c "import sys,re; keys = open(sys.argv[1]).read().split('\n'); entries = open(sys.argv[2]).read().split('>'); [ sys.stdout.write('>'+e) for e in entries if not re.sub('\\\..*','',e,flags=re.DOTALL) in keys ]" \
	${transcript_id} ${sci_name}.cds.id_formatted.fasta > ${sci_name}.cds.mask.fasta
	echo CDS: num line after mask: `wc -l ${sci_name}.cds.mask.fasta`
	cp ${sci_name}.cds.mask.fasta ${dir_cds}.mask/${sci_name}.cds.mask.fasta
fi

rm -rf ${dir_work}/${bn}

###################
echo ending at `date`

