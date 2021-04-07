#!/bin/sh
#$ -S /bin/bash
#$ -cwd
##$ -l month
##$ -l d_rt=1400:00:00
##$ -l s_rt=1400:00:00
#$ -l s_vmem=4G
#$ -l mem_req=4G
#$ -pe def_slot 4
#$ -t 1-2695

#1-2695

source /home/kfuku/.bashrc
ulimit -s unlimited

echo running on `hostname`
echo `date`: start

sra_date=2018_5_1

PATH=/lustre1/home/kfuku/klab/build/:${PATH}
wd="/lustre1/home/kfuku/my_project/convergence_duplication/20180207_kallisto"
dir_out="/lustre1/home/kfuku/my_db/Ensembl/release-91/kallisto"
dir_sra="/lustre1/home/kfuku/ncbi/public/sra"
dir_cdna_mask="/lustre1/home/kfuku/my_db/Ensembl/release-91/cdna.mask"
dir_mask_fasta='/lustre1/home/kfuku/my_db/Ensembl/release-91/mask/fasta'
sra_table="/lustre1/home/kfuku/my_db/Ensembl/release-91/sra/${sra_date}/sra_table_qualified_${sra_date}.tsv"
key_ascp="/lustre1/home/kfuku/.aspera/connect/etc/asperaweb_id_dsa.openssh"

overwrite=0 # 0 = False, 1 = True
row=$[$SGE_TASK_ID+1]
sci_name=`cat ${sra_table} | awk -F "\t" "NR==${row}" | awk -F "\t" '{print $1}'` ; echo sci_name: ${sci_name}
sci_name_ub=`echo ${sci_name} | sed -e "s/[[:space:]]/_/"` ; echo sci_name_ub: ${sci_name_ub}
index_file=`ls ${dir_cdna_mask}/${sci_name_ub}*.fasta.index` ; echo index_file: ${index_file}
sra_primary=`cat ${sra_table} | awk -F "\t" "NR==${row}" | awk -F "\t" '{print $18}'` ; echo sra_primary: ${sra_primary}
sra_exp=`cat ${sra_table} | awk -F "\t" "NR==${row}" | awk -F "\t" '{print $16}'` ; echo sra_exp: ${sra_exp}
sra_run=`cat ${sra_table} | awk -F "\t" "NR==${row}" | awk -F "\t" '{print $17}'` ; echo sra_run: ${sra_run}
lib_layout=`cat ${sra_table} | awk -F "\t" "NR==${row}" | awk -F "\t" '{print $27}'` ; echo lib_layout: ${lib_layout}
nominal_length=`cat ${sra_table} | awk -F "\t" "NR==${row}" | awk -F "\t" '{print $35}' | sed -e "s/\..*//"` ; echo nominal_length: ${nominal_length}
if [ -z "${nominal_length}" ] || [ ${nominal_length} -lt 200 ]; then
	nominal_length=200
fi
echo adjusted nominal_length: ${nominal_length}
nominal_sd=$[${nominal_length}/10]; echo nominal_sd: ${nominal_sd}

# set directory
if [ ! -e ${dir_out}/${sci_name_ub} ]; then
	mkdir ${dir_out}/${sci_name_ub}
fi
if [ ! -e ${dir_out}/${sci_name_ub}/${sra_run} ]; then
	mkdir ${dir_out}/${sci_name_ub}/${sra_run}
fi
cd ${dir_out}/${sci_name_ub}/${sra_run}
echo "working in: `pwd`"
echo "`date`: start kallisto quantification"
	if [ ${fastq_file_n} -eq 1 ]; then
		kallisto_infiles_fastp="--single --fragment-length=${nominal_length} --sd=${nominal_sd} ${sra_run}.fastp.fastq"
	elif [ ${fastq_file_n} -gt 1 ]; then
		kallisto_infiles_fastp="${sra_run}_1.fastp.fastq ${sra_run}_2.fastp.fastq"
	fi
 amalgkit quant -
	#kallisto quant \
	#--index=${index_file} \
	#--output-dir=. \
	#--plaintext \
	#--threads=${NSLOTS} \
	#--bias \
	#${kallisto_infiles_unfiltered}
	#mv abundance.tsv ${sra_run}.unfiltered.kallisto.tsv
	#mv run_info.json ${sra_run}.unfiltered.kallisto.json

#	echo "`date`: start deleting all fastq files"
#	rm ${sra_run}*.fastq
else
	echo "Output files exist already. Exiting."
fi

echo `date`: end





: <<'#_______________CO_______________'

# single end example /usr/local/ftp/public/ddbj_database/dra/fastq/SRA142/SRA142856/SRX477531
# paired end example /usr/local/ftp/public/ddbj_database/dra/fastq/SRA059/SRA059960/SRX196362

# debug
grep -nw ./* -e "ERROR" | sed -e "s/:.*//" | uniq
grep -nw ./* -e "integer expression expected" | sed -e "s/:.*//" | uniq
grep -nw ./* -e "Warning:" | sed -e "s/:.*//" | uniq
grep -nw ./* -e "not found while parsing within network system module" | sed -e "s/:.*//" | uniq
grep -nw ./* -e "timeout exhausted while reading file within network system module" | sed -e "s/:.*//" | uniq
grep -nw ./* -e "path not found while resolving tree within virtual file system module" | sed -e "s/:.*//" | uniq
grep -nw ./* -e "item not found while constructing within virtual database module" | sed -e "s/:.*//" | uniq

	echo "`date`: start fastq subsampling"
	for fastq in ${fastq_files[*]}
	do
		fatt_out="${fastq}_fatt_count_before_trim.tsv"
		fatt count ${fastq} > ${fatt_out}
		num_reads=`cat ${fatt_out} | awk -F "\t" "NR==2" | awk -F "\t" '{print $2}'`
		avg_read_len=`cat ${fatt_out} | awk -F "\t" "NR==2" | awk -F "\t" '{print $4}'`
		min_read_len=`cat ${fatt_out} | awk -F "\t" "NR==2" | awk -F "\t" '{print $5}'`
		max_read_len=`cat ${fatt_out} | awk -F "\t" "NR==2" | awk -F "\t" '{print $6}'`
		if [ ${max_read_len} -ne ${min_read_len} ]; then
			echo "Warning: maximum read length was different from minimum read length"
		fi
		if [ ${num_reads} -gt ${subsample_size} ]; then
			num_row=$[${subsample_size}*4]
			head -${num_row} ${fastq} > ${fastq}.subsample_tmp
			rm ${fastq}
			mv ${fastq}.subsample_tmp ${fastq}
		fi
		if [ ${max_read_len} -gt ${trim_length} ]; then
			num_trim=$[${max_read_len}-${trim_length}]
			/lustre1/home/kfuku/bin/seqtk trimfq -e ${num_trim} ${fastq} > ${fastq}.trim_tmp
			rm ${fastq}
			mv ${fastq}.trim_tmp ${fastq}
		fi
		fatt count ${fastq} > ${fastq}_fatt_count_after_trim.tsv
		/lustre1/home/kfuku/FastQC/fastqc --nogroup ${fastq}
	done


#_______________CO_______________
