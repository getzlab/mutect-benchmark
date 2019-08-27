#!/bin/bash

# usage: ./slim_maf.sh xxxx.maf
# get patient ID
tcga="TCGA"
patient=$(grep "## tumor_sample=" $1  | head -n 1 | tail -c 25)
patient_id="$tcga$patient"


# get cohort
bn=$(basename $1)
cohort=$(echo $bn | awk '{print substr($0,0,5)}')


grep -v '^#' $1 | while read -r file ; do
	awk -F $'\t' \
	-v cc="$cohort" -v pid="$patient_id" \
	-v OFS='\t' \
	'{ if($6 == $7) {
		ref=$11
		mut_type=$9

		if (ref==$12) {alt=$13}
		if (ref==$13) {alt=$12}
		if (mut_type=="Missense_Mutation") {mut_type="Missense"}
		if (mut_type=="Silence") {mut_type="Synonymous"}
		print $5,$6,$1,ref,alt,mut_type,$10,pid,cc
	} }' <<< "$file" 
	#echo $filenew
	#new=cut -f2 $file
	#echo $new
done





