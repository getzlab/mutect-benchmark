#!/bin/bash
# export PID="ACC-OR-A5J1"
# export mem=10
# export interval_list="/demo-mount/refs/interval1_22.list"
# export ref="/demo-mount/refs/Homo_sapiens_assembly19.fasta"
# export ref_dict="/demo-mount/refs/Homo_sapiens_assembly19.dict"
# export data_source_folder="/demo-mount/funcotator_dataSources.v1.6.20190124s"
set -eo pipefail

outdir="/demo-mount/M2full/${pid}"
vcfs="${outdir}/vcfs/*.vcf"
stats="${outdir}/vcfs/*.stats"
normal_piles="${outdir}/normal_pile/*.table"
tumor_piles="${outdir}/tumor_pile/*.table"
f1r2s="${outdir}/f1r2/*.tar.gz"

# DOUBLE check if scattered work has been finished
[ `ls -1q $vcfs | wc -l` -eq 22 ] || { echo "not 22 vcfs, exiting .."; exit 1; }
[ `ls -1q $stats | wc -l` -eq 22 ] || { echo "not 22 vcf stats files, exiting .."; exit 1; }
[ `ls -1q $normal_piles | wc -l` -eq 22 ] || { echo "not 22 normal_piles, exiting .."; exit 1; }
[ `ls -1q $tumor_piles | wc -l` -eq 22 ] || { echo "not 22 tumor piles, exiting .."; exit 1; }
[ `ls -1q $f1r2s | wc -l` -eq 22 ] || { echo "not 22 f1r2 tar files, exiting .."; exit 1; }

all_vcfs_input=`for file in $vcfs; do printf -- " -I ${file}"; done`
all_stats_input=`for file in $stats; do printf -- " -stats ${file}"; done`
all_normal_piles_input=`for file in $normal_piles; do printf -- " -I ${file}"; done`
all_tumor_piles_input=`for file in $tumor_piles; do printf -- " -I ${file}"; done`
all_f1r2_input=`for file in $f1r2s; do printf -- " -I ${file}"; done`

# command and output formatter
gatkm="gatk --java-options -DGATK_STACKTRACE_ON_USER_EXCEPTION=true"
merged_unfiltered_vcf="${outdir}/merged_unfiltered.vcf"
merged_stats="${outdir}/merged_unfiltered.stats"
artifact_prior="${outdir}/artifact_prior.tar.gz"
normal_pile_table="${outdir}/normal_pile.tsv"
tumor_pile_table="${outdir}/tumor_pile.tsv"
contamination_table="${outdir}/contamination.table"
segments_table="${outdir}/segments.table"
filtering_stats="${outdir}/filtering.stats"
merged_filtered_vcf="${outdir}/merged_filtered.vcf"
annot_merged_filtered_maf="${outdir}/annot_merged_flitered.vcf" 

echo $annot_merged_filtered_maf
# check if this pair has been finished
[ -f $annot_merged_filtered_maf ] && { echo "${pid} has finished, exiting .."; exit 0; }


echo "----------------------- merge vcfs----------------------------"
$gatkm MergeVcfs \
	$all_vcfs_input \
	-O $merged_unfiltered_vcf

echo "------------------------merge mutect stats-------------------------------"
$gatkm MergeMutectStats \
	$all_stats_input \
	-O $merged_stats

echo "-----------------------learn read orientation model-----------------------"
$gatkm LearnReadOrientationModel \
	$all_f1r2_input \
	-O $artifact_prior

echo "-----------------------gather pile up summaries-----------------------"
$gatkm GatherPileupSummaries \
	$all_normal_piles_input \
	--sequence-dictionary $ref_dict \
	-O $normal_pile_table

echo "-----------------------gather tumor piles-----------------------"
$gatkm GatherPileupSummaries \
	$all_tumor_piles_input \
	--sequence-dictionary $ref_dict \
	-O $tumor_pile_table

echo "-----------------------calc contamination-----------------------"
$gatkm CalculateContamination \
	-I $tumor_pile_table \
	-O $contamination_table \
	--tumor-segmentation $segments_table \
	-matched $normal_pile_table

echo "-----------------------filter calls---------------------------"
$gatkm FilterMutectCalls \
	-V $merged_unfiltered_vcf \
	-R $ref \
	--contamination-table $contamination_table \
	--tumor-segmentation $segments_table \
	--ob-priors $artifact_prior \
	-stats $merged_stats \
	--filtering-stats $filtering_stats \
	-O $merged_filtered_vcf

echo "-----------------------annotate-----------------------"
$gatkm Funcotator \
	--data-sources-path $data_source_folder \
	--ref-version hg19 \
	--output-file-format MAF \
	-R $ref \
	-V $merged_filtered_vcf \
	-O $annot_merged_filtered_maf \
	-L $interval_list \


