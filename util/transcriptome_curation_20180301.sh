#!/bin/sh
#$ -S /bin/bash
#$ -cwd
#$ -l month
#$ -l d_rt=1400:00:00
#$ -l s_rt=1400:00:00
#$ -l s_vmem=8G
#$ -l mem_req=8G
#$ -pe def_slot 1
#$ -t 1-1

#1-66

echo `date`: start

dist_method=pearson
mapping_rate_cutoff=0.2

mode=imac
sra_date=2018_5_1
if [ ${mode} = 'nig' ]; then
	sra_table="/lustre1/home/kfuku/my_db/Ensembl/release-91/sra/${sra_date}/sra_table_mapped_${sra_date}.tsv"
	dir_work="/lustre1/home/kfuku/my_project/convergence_duplication/20180207_kallisto"
	dir_out="/lustre1/home/kfuku/my_db/Ensembl/release-91/curated_transcriptome/${sra_date}"
	dir_exp="/lustre1/home/kfuku/my_db/Ensembl/release-91/kallisto_summary/${sra_date}/tpm.masked.kallisto.tsv"
	dir_exp_curated="/lustre1/home/kfuku/my_db/Ensembl/release-91/kallisto_summary/${sra_date}/tpm.masked.kallisto.gene.log.tsv"
	dir_idm="/lustre1/home/kfuku/my_db/Ensembl/release-91/id_mapping"
	dir_script="/lustre1/home/kfuku/my_script"
	id=$[${SGE_TASK_ID}-1]
elif [ ${mode} = 'imac' ]; then
	sra_table="/Users/kf/Dropbox/kfdata/02_Data/my_db/Ensembl/release-91/sra/${sra_date}/sra_table_mapped_${sra_date}.tsv"
	dir_work="/Users/kf/Dropbox/kfdata/02_Data/04_Convergence_Duplication/20180207_kallisto"
	dir_out="/Users/kf/Dropbox/kfdata/02_Data/my_db/Ensembl/release-91/curated_transcriptome/${sra_date}"
	dir_exp="/Users/kf/Dropbox/kfdata/02_Data/my_db/Ensembl/release-91/kallisto_summary/${sra_date}/tpm.masked.kallisto.tsv"
	dir_exp_curated="/Users/kf/Dropbox/kfdata/02_Data/my_db/Ensembl/release-91/kallisto_summary/${sra_date}/tpm.masked.kallisto.gene.log.tsv"
	dir_idm="/Users/kf/Dropbox/kfdata/02_Data/my_db/Ensembl/release-91/id_mapping"
	dir_script="/Users/kf/Dropbox/kfdata/02_Data/my_script"
fi

dir_out_tau=${dir_out}/tau
dir_out_tc=${dir_out}/tc
dir_out_r2=${dir_out}/r2
dir_out_sra=${dir_out}/sra
dir_out_uncorrected_tissue_mean=${dir_out}/uncorrected_tissue_mean
dir_out_tissue_mean=${dir_out}/tissue_mean
dir_out_plot=${dir_out}/plot
dir_out_sva=${dir_out}/sva
infiles=( `ls ${dir_exp}` )
echo "infiles: ${infiles[@]}"

for infile in ${infiles[@]}; do
	sci_name=`echo ${infile} | sed -e "s/\..*//"`
	dir_tmp=${dir_work}/${sci_name}
	file_exp=${infile}
	sci_name_idm=`echo ${sci_name:0:1} | awk '{print tolower($0)}'``echo ${sci_name} | sed -e "s/\..*//" -e "s/.*_//"`
	file_idm=`ls ${dir_idm} | grep -e "^${sci_name_idm}"`

	# prepare and set directories
	directories=( `set | grep "^dir_" | sed -e "s/=.*//"` )
	for d in ${directories[@]}; do
		if [ ! -e `eval echo '$'${d}` ]; then
			echo creating: `eval echo '$'${d}`
			mkdir `eval echo '$'${d}`
		fi
	done
	cd ${dir_tmp}

	if [ ! -s ${dir_exp_curated}/${sci_name}.gene.log.tsv ]; then
		echo `date`: Start: transcriptome_data_preparation.r

		Rscript ${dir_script}/transcriptome_data_preparation.r \
		${dir_exp}/${file_exp} \
		${sci_name}.gene.log.tsv \
		${dir_idm}/${file_idm} \
		1 \
		1

		cp ${sci_name}.gene.log.tsv ${dir_exp_curated}/${sci_name}.gene.log.tsv
	else
		echo `date`: Skipped: transcriptome_data_preparation.r
	fi

	if [ ! -s ${dir_out_tc}/${sci_name}.tc.tsv ]; then
		echo `date`: Start: transcriptome_curation.r
		
		Rscript ${dir_script}/transcriptome_curation.r \
		${dir_exp_curated}/${sci_name}.gene.log.tsv \
		${sra_table} \
		${dir_tmp} \
		${dist_method} \
		${mapping_rate_cutoff} \
		0 \
		1 \
		'brain|heart|kidney|liver|ovary|testis'

		cp *.pdf ${dir_out_plot}
		cp *.sva.*.RData ${dir_out_sva}
		cp *.r2.tsv ${dir_out_r2}
		cp ${sci_name}.tau.tsv ${dir_out_tau}
		cp ${sci_name}.uncorrected.tissue.mean.tsv ${dir_out_uncorrected_tissue_mean}
		cp ${sci_name}.tissue.mean.tsv ${dir_out_tissue_mean}
		cp ${sci_name}.sra.tsv ${dir_out_sra}
		cp ${sci_name}.tc.tsv ${dir_out_tc}
	else
		echo `date`: Skipped: transcriptome_curation.r
	fi

	rm -r ${dir_tmp}
done

echo `date`: end





: <<'#_______________CO_______________'



#_______________CO_______________



