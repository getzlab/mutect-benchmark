"""
Input a tsv with three columns - <tumor, normal, pid>
generate a yaml with exact replicate of gatk-best-practice[Tag:2.5.0]
that calls <gatk Mutect2> scattered on a predefined <interval_list>

Qing Zhang 09/25/2019

"""

import argparse
import pandas as pd
import bp_vars  # common best-practice variables
import itertools
from datetime import datetime
import yaml

parser = argparse.ArgumentParser(description="-pfile filepath [-test]")
parser.add_argument("-pfile", dest="pfile", type=str,
                    default="~/Documents/igv_remote/helpers/mpairs.tsv",
                    help="a tsv that contains the paired gs paths")
parser.add_argument("-test", dest="iftest",
                    action="store_true",
                    help="if you want to test on first 3 samples",
                    default=True)
args = parser.parse_args()
if args.iftest:
    print("Under testing mode, only first three pairs will be run...")
    pair_df = pd.read_csv(args.pfile,
                          sep="\t",  nrows=3)
else:
    pair_df = pd.read_csv(args.pfile, sep='\t')

if False:
    # only runs in debug mode
    pair_df = pd.read_csv("~/Documents/igv_remote/helpers/mpairs.tsv",
                          sep="\t", header=None)
    pair_df.columns = ["tumor", "normal", "pid", "cohort"]
    pair_df = pair_df.iloc[:4, :]  # testing with the first three samples

############## make params for yaml config ###########

# input files for each pair
tumor_paths = pair_df["tumor"].tolist()
normal_paths = pair_df["normal"].tolist()

# pid defines where output files go
pids = pair_df["pid"].tolist()
cohort = pair_df["cohort"].tolist()

# intervals input
subints = ['{:0>4}'.format(i) for i in range(M2["scatter_count"])]
intervals = ["{}/{}-scattered.interval_list".format(
    M2["sub_intervals_dir"], int1) for int1 in subints]

# directory structure (somewhere in NFS)
res_dirs = ["{}/M2bp/res/{}/{}".format(
    M2["NFS"],
    cohort1, pid1) for (cohort1, pid1) in zip(cohort, pids)]

# the full permutation of each samples wrt all pieces
file_indices = list(range(len(pids)))
intervals_indices = list(range(50))
comb_indices = list(itertools.product(file_indices, intervals_indices))

# an utlity function to unfold index tuples to list


def unfold(uniq_vec, idx, tups=comb_indices):
    newlist = [uniq_vec[tup[idx]] for tup in tups]
    return(newlist)


input_dict = dict(
    tumor=unfold(tumor_paths, idx=0),
    normal=unfold(normal_paths, idx=0),
    resdir=unfold(res_dirs, idx=0),
    interval=unfold(intervals, idx=1),
    chr=unfold(subints, idx=1)
)
input_dict.update(M2)  # append dict with default params


def gen_dict(input_dict=input_dict, job_name="sca_m2"):
    """
#!/bin/bash
set -eo pipefail


# define the format of output files
out_tpile="${resdir}/${tpile}/${chr}-tpile.table"
out_npile="${resdir}/${npile}/${chr}-npile.table"
out_f1r2="${resdir}/${f1r2}/${chr}-f1r2.tar.gz"
out_vcf="${resdir}/${subvcfs}/${chr}_unfilter.vcf"

# make file structure
mkdir -p $resdir
mkdir -p $resdir/$f1r2
mkdir -p $resdir/$tpile
mkdir -p $resdir/$npile
mkdir -p $resdir/$subvcfs

command_mem=$heap_mem
[ -f $normal_name.txt ] && gatk --java-options "-Xmx${command_mem}g" GetSampleName -R $ref_fasta \
    -I ${tumor_bam} -O tumor_name.txt

tumor_command_line="-I ${tumor_bam} -tumor `cat tumor_name.txt`"
cat tumor_name.txt


[ -f $out_vcf ] && gatk --java-options "-Xmx${command_mem}g" GetSampleName -R $ref_fasta \
    -I ${normal_bam} -O normal_name.txt

normal_command_line="-I ${normal_bam} -normal `cat normal_name.txt`"
cat normal_name.txt


[ -f $out_tpile] && gatk --java-options "-Xmx${command_mem}g" Mutect2 \
    -R $ref_fasta \
    $tumor_command_line \
    $normal_command_line \
    --germline-resource $gnomad \
    -pon $default_pon \
    -L  $interval \
    -O $out_vcf \
    --f1r2-tar-gz $out_f1r2
echo finish Mutect2 ----------



[ -f $out_npile ] && gatk --java-options "-Xmx${command_mem}g" GetPileupSummaries \
    -R $ref_fasta -I ${tumor_bam} \
    --interval-set-rule INTERSECTION \
    -L $interval \
    -V $vfc \
    -L $vfc \
    -O $out_tpile
echo Getpiles tumor ---------

gatk --java-options "-Xmx${command_mem}g" GetPileupSummaries \
    -R $ref_fasta -I ${normal_bam} \
    --interval-set-rule INTERSECTION \
    -L $interval \
    -V $vfc \
    -L $vfc \
    -O $out_npile
echo Getpiles normal --------

echo allfin!

    """
    input_key_list = list(input_dict.keys())
    yy = dict(
        name="M2sca",
        resources={
            "cpus_per_task": 2,
            "mem_per_cpu": "{}G".format(input_dict["heap_mem"])
        },
        inputs=input_dict,
        backend={"type": "Local"},
        localization={
            # staging directory will be a random folder in home
            "staging_dir": "{}/{}/{}_{}".format(input_dict["NFS"], "M2bp", job_name, datetime.now().strftime('%m%d%H%M')),
            "overrides": dict(zip(input_key_list, [None]*len(input_key_list)))
        }
    )
    return(yy)


"""
#!/bin/bash
# export PID="ACC-OR-A5J1"
export mem=20
export interval_list="gs://broad-references/hg19/v0/wgs_calling_regions.v1.interval_list"
export ref="/demo-mount/refs/Homo_sapiens_assembly19.fasta"
export ref_dict="/demo-mount/refs/Homo_sapiens_assembly19.dict"
export data_source_folder="/demo-mount/refs/funco_slim"
export tumor_bam="gs://fc-a6cde48f-6248-4428-a950-a256d46a828b/ichor_patients/7_FC19663997_HFW3TBBXX.7.aligned.duplicates_marked.bam"
set -eo pipefail

outdir="/demo-mount/for_Ziao/res"
vcfs="${outdir}/subvcfs/*.vcf"
stats="${outdir}/subvcfs/*.stats"
normal_piles="${outdir}/normal_piles/*"
tumor_piles="${outdir}/tumor_piles/*"
f1r2s="${outdir}/f1r2/*.tar.gz"

all_vcfs_input=`for file in $vcfs; do printf -- " -I ${file}"; done`
all_stats_input=`for file in $stats; do printf -- " -stats ${file}"; done`
all_normal_piles_input=`for file in $normal_piles; do printf -- " -I ${file}"; done`
all_tumor_piles_input=`for file in $tumor_piles; do printf -- " -I ${file}"; done`
all_f1r2_input=`for file in $f1r2s; do printf -- " -I ${file}"; done`

# command and output formatter
gatkm="gatk --java-options -Xmx20g"
merged_unfiltered_vcf="${outdir}/merged_unfiltered.vcf"
merged_stats="${outdir}/merged_unfiltered.stats"
artifact_prior="${outdir}/artifact_prior.tar.gz"
normal_pile_table="${outdir}/normal_pile.tsv"
tumor_pile_table="${outdir}/tumor_pile.tsv"
contamination_table="${outdir}/contamination.table"
segments_table="${outdir}/segments.table"
filtering_stats="${outdir}/filtering.stats"
merged_filtered_vcf="${outdir}/merged_filtered.vcf"
aligned_merged_filtered_vcf="${outdir}/aligned_merged_filtered.vcf"
annot_merged_filtered_maf="${outdir}/annot_merged_filtered.vcf" 
bwa_img="/demo-mount/refs/Homo_sapiens_assembly19.fasta.img"

# check if this pair has been finished
[ -f $annot_merged_filtered_maf ] && { echo "${pid} has finished, exiting .."; exit 0; }


echo "----------------------- merge vcfs----------------------------"
[ -f $merged_stats ] && $gatkm MergeVcfs \
	$all_vcfs_input \
	-O $merged_unfiltered_vcf

echo "------------------------merge mutect stats-------------------------------"
[ -f $artifact_prior ] && $gatkm MergeMutectStats \
	$all_stats_input \
	-O $merged_stats

echo "-----------------------learn read orientation model-----------------------"
[ -f $normal_pile_table ] && $gatkm LearnReadOrientationModel \
	$all_f1r2_input \
	-O $artifact_prior

echo "-----------------------gather pile up summaries-----------------------"
[ -f $tumor_pile_table ] && $gatkm GatherPileupSummaries \
	$all_normal_piles_input \
	--sequence-dictionary $ref_dict \
	-O $normal_pile_table

echo "-----------------------gather tumor piles-----------------------"
[ -f $contamination_table ] && $gatkm GatherPileupSummaries \
	$all_tumor_piles_input \
	--sequence-dictionary $ref_dict \
	-O $tumor_pile_table

echo "-----------------------calc contamination-----------------------"
[ -f $merged_filtered_vcf ] && $gatkm CalculateContamination \
	-I $tumor_pile_table \
	-O $contamination_table \
	--tumor-segmentation $segments_table \
	-matched $normal_pile_table

echo "-----------------------filter calls---------------------------"
[ -f $aligned_merged_filtered_vcf ] && $gatkm FilterMutectCalls \
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
	-V $aligned_merged_filtered_vcf \
	-O $annot_merged_filtered_maf \
	-L $interval_list \
	--remove-filtered-variants true

echo "-----------all finished-----------------"


"""





if __name__ == "__main__":
    job_name = "sca_m2"
    yaml_name="k9_{}.yaml".format(job_name)
    script_name="k9_{}.sh".format(job_name)
    yy=gen_dict(job_name=job_name)

    # the yaml config (worth to have it as to document!)
    with open(yaml_name, "w") as outfile:
        yaml.dump(yy, outfile, default_flow_style = False)

    # generate the companion script
    with open(script_name, "w") as outfile:
        outfile.write(gen_dict.__doc__)
