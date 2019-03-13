#!/bin/sh
#$ -S /bin/bash
#$ -cwd
#$ -l month
#$ -l d_rt=1400:00:00
#$ -l s_rt=1400:00:00
#$ -l s_vmem=6G
#$ -l mem_req=6G
#$ -pe def_slot 1
#$ -t 1-66

#1-66

echo `date`: start

date_sra='2018_5_1'
dir_kallisto="/lustre1/home/kfuku/my_db/Ensembl/release-91/kallisto"
dir_summary="/lustre1/home/kfuku/my_db/Ensembl/release-91/kallisto_summary/${date_sra}/"
sci_names=( `ls ${dir_kallisto}` )
id=$[${SGE_TASK_ID}-1]
sci_name=${sci_names[$id]}
cd ${dir_kallisto}/${sci_name}
sras=`ls ${dir_kallisto}/${sci_name}`
echo working on `pwd`
echo SRA IDs: ${sras[@]}

directories=( `set | grep "^dir_" | sed -e "s/=.*//"` )
for d in ${directories[@]}; do
	if [ ! -e `eval echo '$'${d}` ]; then
		echo creating: `eval echo '$'${d}`
		mkdir `eval echo '$'${d}`
	fi
done


ext_kallisto=( '.unfiltered.kallisto.tsv' '.fastp.kallisto.tsv' '.masked.kallisto.tsv' )
echo `date`: checking output directories
for ext in ${ext_kallisto[@]}; do
	for ext2 in eff_length est_counts tpm; do
		if [ ! -e ${dir_summary}/${ext2}${ext} ]; then
			mkdir ${dir_summary}/${ext2}${ext}
		fi
	done
done
if [ ! -e ${dir_summary}/fatt.tsv ]; then
	mkdir ${dir_summary}/fatt.tsv
fi

for ext in ${ext_kallisto[@]}; do
	echo `date`: start making table: ${ext}
	table_esize=""
	table_count=""
	table_tpm=""
	for sra in ${sras[@]}; do
		cd ${dir_kallisto}/${sci_name}/${sra}
		if [ -s ${sra}${ext} ]; then
			if [ -z "${table_esize}" ]; then
				table_esize="`cat ${sra}${ext} | cut -f 1,2,3 | sed -e s/eff_length/${sra}/`"
				table_count="`cat ${sra}${ext} | cut -f 1,2,4 | sed -e s/est_counts/${sra}/`"
				table_tpm="`cat ${sra}${ext} | cut -f 1,2,5 | sed -e s/tpm/${sra}/`"
				echo num_line: ${sra}: `echo "${table_count}" | wc -l`
			else
				tmp_size="`cat ${sra}${ext} | cut -f 3 | sed -e s/eff_length/${sra}/`"
				tmp_count="`cat ${sra}${ext} | cut -f 4 | sed -e s/est_counts/${sra}/`"
				tmp_tpm="`cat ${sra}${ext} | cut -f 5 | sed -e s/tpm/${sra}/`"
				table_esize=`paste <(echo "${table_esize}") <(echo "${tmp_size}")`
				table_count=`paste <(echo "${table_count}") <(echo "${tmp_count}")`
				table_tpm=`paste <(echo "${table_tpm}") <(echo "${tmp_tpm}")`
				echo num_line: ${sra}: `echo "${tmp_count}" | wc -l`
			fi
		else
			>&2 echo Warning: no ${ext} file found: ${sci_name}: ${sra}
		fi
		cd ${dir_kallisto}/${sci_name}
	done
	cd ${dir_summary}
	echo -e "${table_esize}" > ${dir_summary}/eff_length${ext}/${sci_name}.eff_length${ext}
	echo -e "${table_count}" > ${dir_summary}/est_counts${ext}/${sci_name}.est_counts${ext}
	echo -e "${table_tpm}" > ${dir_summary}/tpm${ext}/${sci_name}.tpm${ext}
done

echo `date`: start making table: fatt
table_fatt=""
for sra in ${sras[@]}; do
	cd ${dir_kallisto}/${sci_name}/${sra}
	if [ -s ${sra}.fatt.tsv ]; then
		if [ -z "${table_fatt}" ]; then
			table_fatt="`cat ${sra}.fatt.tsv`"
		else
			tmp_fatt="`tail -n +2 ${sra}.fatt.tsv`"
			table_fatt=${table_fatt}"\n"${tmp_fatt}
		fi
	else
		>&2 echo Warning: no fatt file found: ${sci_name}: ${sra}
	fi
	cd ${dir_kallisto}/${sci_name}
done
cd ${dir_summary}
echo -e "${table_fatt}" > ${dir_summary}/fatt.tsv/${sci_name}.fatt.tsv

echo `date`: end

