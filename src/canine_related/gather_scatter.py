#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  6 08:59:50 2019

The M2 workflow per tumor-normal pair

The staging dir is where the scattered outputs are stored



@author: qingzhang
"""


import os
import yaml
import glob

# global params
expected_chr_counts = 22
filepath = "/demo-mount/refs/mpairs.tsv"
chr_l = list(range(1, 23))
ref = "/demo-mount/refs/Homo_sapiens_assembly19.fasta"
germ = "/demo-mount/refs/af-only-gnomad.raw.sites.b37.vcf"
vfc = "/demo-mount/refs/gnomad_var_for_contamination/variants_for_contamination.vcf"
pon = "/demo-mount/M2pon/makePON/AllinOne/merged_vcfs/fin_concat/AIO_merged_PON.vcf"
heap_mem = 18

# check if the sample pair is runned by canine at all
with open(filepath) as f:

    tumor_l = []
    normal_l = []
    pid_l = []
    cohort_l = []
    for l in f:
        tumor, normal, pid, cohort = l.strip().split("\t")
        scatter_stage = "/demo-mount/canine_m2full/{}".format(pid)
        if os.path.exists(scatter_stage):
            
            # how many files are there
            vcfs = glob.glob("{}/vcfs/*.stats".format(scatter_stage))
            tumor_piles = glob.glob("{}/tumor_pile/*.table".format(scatter_stage))
            normal_piles = glob.glob("{}/normal_pile/*.table".format(scatter_stage))
            f1r2s = glob.glob("{}/f1r2/*.tar.gz".format(scatter_stage))
            
            if len(vcfs) == 22 and len(tumor_piles) == 22 and len(normal_piles) == 22 and len(f1r2s) == 22:    
                tumor_l.append(tumor)
                normal_l.append(normal)
                pid_l.append(pid)


outdir = "/demo-mount/M2full"
stage_dir = "/demo-mount/canine_m2merge"
ostrdir = "/demo-mount/canine_M2full_ostr"

if os.path.exists(stage_dir):
    print("Warning! stage dir exist!")


"""
# what to run
gatk --java-options "-Xmx${command_mem}m" MergeVcfs \
    -I ${sep=' -I ' input_vcfs} \
    -O ${output_vcf}

gatk --java-options "-Xmx${command_mem}m" MergeMutectStats \
    -stats ${sep=" -stats " stats} \
    -O merged.stats

gatk --java-options "-Xmx${command_mem}m" LearnReadOrientationModel \
    -I ${sep=" -I " f1r2_tar_gz} \
    -O "artifact-priors.tar.gz"
    
gatk --java-options "-Xmx${command_mem}m" GatherPileupSummaries \
    --sequence-dictionary ${ref_dict} \
    -I ${sep=' -I ' input_tables} \
    -O ${output_name}.tsv
gatk --java-options "-Xmx${command_mem}m" GatherPileupSummaries \
    --sequence-dictionary ${ref_dict} \
    -I ${sep=' -I ' input_tables} \
    -O ${output_name}.tsv

gatk --java-options "-Xmx${command_mem}m" CalculateContamination \
    -I ${tumor_pileups} \
    -O contamination.table \
    --tumor-segmentation segments.table \
    ${"-matched " + normal_pileups}
    
gatk --java-options "-Xmx${command_mem}m" FilterMutectCalls -V ${unfiltered_vcf} \
    -R ${ref_fasta} \
    -O ${output_vcf} \
    ${"--contamination-table " + contamination_table} \
    ${"--tumor-segmentation " + maf_segments} \
    ${"--ob-priors " + artifact_priors_tar_gz} \
    ${"-stats " + mutect_stats} \
    --filtering-stats filtering.stats \
    ${m2_extra_filtering_args}

         gatk --java-options "-Xmx${command_mem}m" Funcotator \
             --data-sources-path $DATA_SOURCES_FOLDER \
             --ref-version ${reference_version} \
             --output-file-format ${output_format} \
             -R ${ref_fasta} \
             -V ${input_vcf} \
             -O ${output_file} \
             ${interval_list_arg} ${default="" interval_list} \
             --annotation-default normal_barcode:${default="Unknown" control_id} \
             --annotation-default tumor_barcode:${default="Unknown" case_id} \
             --annotation-default Center:${default="Unknown" sequencing_center} \
             --annotation-default source:${default="Unknown" sequence_source} \
             ${"--transcript-selection-mode " + transcript_selection_mode} \
             ${transcript_selection_arg}${default="" sep=" --transcript-list " transcript_selection_list} \
             ${annotation_def_arg}${default="" sep=" --annotation-default " annotation_defaults} \
             ${annotation_over_arg}${default="" sep=" --annotation-override " annotation_overrides} \
             ${excluded_fields_args}${default="" sep=" --exclude-field " funcotator_excluded_fields} \
             ${filter_funcotations_args} \
             ${extra_args_arg}
"""

# input vcfs to merge
tumor_piles = " ".join(["-I {}/{pid}/tumor_pile/tumor-{}-pileups.table".format(outdir,"{pid}",chr) for chr in chr_l])
normal_piles = " ".join(["-I {}/{pid}/normal_pile/normal-{}-pileups.table".format(outdir,"{pid}",chr) for chr in chr_l])
f1r2s = " ".join([" -I {}/{}/f1r2/{}_chr{}-f1r2.tar.gz".format(outdir,"{pid}","{pid}",chr) for chr in chr_l])
vcfs = " ".join([" -I {}/{}/vcfs/{}_chr{}_unfilter.vcf".format(outdir,"{pid}","{pid}",chr) for chr in chr_l])
stats = " ".join([" -stats {}/{}/vcfs/{}_chr{}_unfilter.vcf.stats".format(outdir,"{pid}","{pid}",chr) for chr in chr_l])
# input 
# the dictionary is scattered by pids
yy = dict(
        name = "M2_merge",
        script = ["set -eo pipefail",
                  "$gatk MergeVcfs $vcfs_args -O $merged_vcf",
                  "$gatk MergeMutectStats $stats_args -O $merged_stats",
                  "$gatk LearnReadOrientationModel $f1r2_args -O $merged_model",
                  "$gatk GatherPileupSummaries $tumor_pile_args -O $merged_tumor_pile",
                  "$gatk GatherPileupSummaries $normal_pile_args -O $merged_normal_pile",
                  "$gatk CalculateContamination -I $merged_tumor_pile -matched $merged_normal_pile -O $contamination_table --tumor-segmentation $segments_table",
                  "$gatk FilterMutectCalls -V $merged_vcf -R $ref -O $filtered_vcf --contamination-table $contamination_table --tumor-segmentation $segments_table --ob-priors $merged_model -stats $merged_stats",
                  "$gatk Funcotator -I $filtered_vcf -O $annot_vcf -R $ref --data-sources-path $dataSourcesFolder --ref-version hg19"],
        inputs = {
                "gatk":"gatk --java-options -Xmx{}g".format(heap_mem),
                "vcfs_args":[vcfs.format(pid = pid) for pid in pid_l],
                "merged_vcf":["{out}/{pid}/{pid}-merged_unfiltered.vcf.gz".format(out=outdir, pid=pid) for pid in pid_l],
                "stats_args":[stats.format(pid = pid) for pid in pid_l],
                "merged_stats":["{out}/{pid}/{pid}-merged.stats".format(out=outdir, pid=pid) for pid in pid_l],
                "f1r2_args":[f1r2s.format(pid = pid) for pid in pid_l],
                "merged_model":["{out}/{pid}/{pid}-merged_model.vcf.gz".format(out=outdir, pid=pid) for pid in pid_l],
                "merged_tumor_pile":["{out}/{pid}/{pid}-tumor-pile.table".format(out=outdir, pid=pid) for pid in pid_l],
                "merged_normal_pile":["{out}/{pid}/{pid}-normal-pile.table".format(out=outdir, pid=pid) for pid in pid_l],
                "normal_pile_args":[normal_piles.format(pid = pid) for pid in pid_l],
                "tumor_pile_args":[tumor_piles.format(pid = pid) for pid in pid_l],
                "contamination_table":["{out}/{pid}/{pid}-contamination.table".format(out=outdir, pid=pid) for pid in pid_l],
                "segments_table":["{out}/{pid}/{pid}-segments.table".format(out=outdir, pid=pid) for pid in pid_l],
                "annot_vcf":["{out}/{pid}/{pid}-annot.vcf.gz".format(out=outdir, pid=pid) for pid in pid_l],
                "dataSourcesFolder":"/demo-mount/funcotator_dataSources.v1.6.20190124s",
                "ref":ref,
                "vfc":vfc
                }
        
        )













