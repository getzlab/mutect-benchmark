M2 = dict(
    # default gatk best practices input
    artifact_modes=["G/T", "C/T"],
    ref_fasta="gs://gatk-best-practices/somatic-b37/Homo_sapiens_assembly19.fasta",
    ref_dict="gs://gatk-best-practices/somatic-b37/Homo_sapiens_assembly19.dict",
    default_pon="gs://gatk-best-practices/somatic-b37/Mutect2-exome-panel.vcf",
    gnomad="gs://gatk-best-practices/somatic-b37/af-only-gnomad.raw.sites.vcf",
    vfc="gs://gatk-best-practices/somatic-b37/small_exac_common_3.vcf",



    # these variables are only used once for preparation
    scatter_count=50,
    wes_interval_list="gs://gatk-best-practices/somatic-b37/whole_exome_agilent_1.1_refseq_plus_3_boosters.Homo_sapiens_assembly19.baits.interval_list",

)



# command history of generating the input files
bp_logs = dict(
    gen_interval_list="""
    gatk SplitIntervals -R Homo_sapiens_assembly19.fasta -L whole_exome_agilent_1.1_refseq_plus_3_boosters.Homo_sapiens_assembly19.baits.interval_list --scatter-count 50 -O gatk_bp_intervals/
    """
)


PoN = dict(
    interval_list="gs://gatk-best-practices/somatic-b37/whole_exome_agilent_1.1_refseq_plus_3_boosters.Homo_sapiens_assembly19.baits.interval_list",
    gnomad="gs://gatk-best-practices/somatic-b37/af-only-gnomad.raw.sites.vcf",
    artifact_mode=["G/T", "C/T"]

M2scatter_sh = 
"""#!/bin/bash
set -eo pipefail

# make file structure
subvcfs_dir=${base_dir}/${pid}/${subvcfs_dir}
fir2_dir=${base_dir}/${pid}/${f1r2_dir}
tpile_dir=${base_dir}/${pid}/${tpile_dir}
npile_dir=${base_dir}/${pid}/${npile_dir}

mkdir -p $subvcfs_dir
mkdir -p $tpile_dir
mkdir -p $npile_dir
mkdir -p $f1r2_dir


command_mem=$heap_mem
[ -f normal_name.txt ] && gatk --java-options "-Xmx${heap_mem}g" GetSampleName -R $ref_fasta \
    -I ${tumor_bam} -O tumor_name.txt

tumor_command_line="-I ${tumor_bam} -tumor `cat tumor_name.txt`"
cat tumor_name.txt


[ -f $out_vcf ] && gatk --java-options "-Xmx${heap_mem}g" GetSampleName -R $ref_fasta \
    -I ${normal_bam} -O normal_name.txt

normal_command_line="-I ${normal_bam} -normal `cat normal_name.txt`"
cat normal_name.txt


[ -f $out_tpile] && gatk --java-options "-Xmx${heap_mem}g" Mutect2 \
    -R $ref_fasta \
    $tumor_command_line \
    $normal_command_line \
    --germline-resource $gnomad \
    -pon $default_pon \
    -L  $interval \
    -O ${subvcfs_dir}/${out_vcf} \
    --f1r2-tar-gz ${f1r2_dir}/${out_f1r2}
echo finish Mutect2 ----------



[ -f $out_npile ] && gatk --java-options "-Xmx${heap_mem}g" GetPileupSummaries \
    -R $ref_fasta -I ${tumor_bam} \
    --interval-set-rule INTERSECTION \
    -L $interval \
    -V $vfc \
    -L $vfc \
    -O ${tpile_dir}/${out_tpile}
echo Getpiles tumor ---------

gatk --java-options "-Xmx${heap_mem}g" GetPileupSummaries \
    -R $ref_fasta -I ${normal_bam} \
    --interval-set-rule INTERSECTION \
    -L $interval \
    -V $vfc \
    -L $vfc \
    -O ${npile_dir}/${out_npile}
echo Getpiles normal --------

echo allfin!
"""

M2merge_sh="""
#!/bin/bash
# export PID="ACC-OR-A5J1"
export mem=20
export interval_list="gs://broad-references/hg19/v0/wgs_calling_regions.v1.interval_list"
export ref="/demo-mount/refs/Homo_sapiens_assembly19.fasta"
export ref_dict="/demo-mount/refs/Homo_sapiens_assembly19.dict"
export data_source_folder="/demo-mount/refs/funco_slim"
export tumor_bam="gs://fc-a6cde48f-6248-4428-a950-a256d46a828b/ichor_patients/7_FC19663997_HFW3TBBXX.7.aligned.duplicates_marked.bam"
set -eo pipefail

outdir="${base_dir}/${pid}"
vcfs="${outdir}/${subvcfs_dir}/*.vcf"
stats="${outdir}/${subvcfs_dir}/*.stats"
normal_piles="${outdir}/${npile_dir}/*"
tumor_piles="${outdir}/${tpile_dir}/*"
f1r2s="${outdir}/${f1r2_dir}/*.tar.gz"

all_vcfs_input=`for file in $vcfs; do printf -- " -I ${file}"; done`
all_stats_input=`for file in $stats; do printf -- " -stats ${file}"; done`
all_normal_piles_input=`for file in $normal_piles; do printf -- " -I ${file}"; done`
all_tumor_piles_input=`for file in $tumor_piles; do printf -- " -I ${file}"; done`
all_f1r2_input=`for file in $f1r2s; do printf -- " -I ${file}"; done`

# command and output formatter
gatkm="gatk --java-options -Xmx${heap_mem}g"
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

M2normal_sh="""
#!bin/bash
set -eo pipefail
outdir="${base_dir}/${pid}"

gatk Mutect2 \ # no germline arg here
    -R $ref_fasta \
    -I $normal_bam \
    -max-mnp-distance 0 \ # uses this only when making PON
    -O ${outdir}/normal.vcf
"""

PONcreate_sh="""
#!bin/bash
set -eo pipefail

mkdir -p ${basedir}/sub_pons
input_stats="${base_dir}/*/*.stats"
all_vcfs_input=`for file in $input_stats; do printf -- " -V ${file%.*}"; done`

# --scatter by interval
gatk GenomicsDBImport \
    --genomicsdb-workspace-path gendb://pon_db \
    -R $ref_fasta \
    --batch-size 50 --merge-input-intervals \ # make things faster
    $all_vcfs_input

echo "genomics db import finished"

gatk CreateSomaticPanelOfNormals 
	-R $ref_fasta \
	--germline-resource ${gnomad} \ # germline arg only appears here
    -V gendb://pon_db \
    -O ${base_dir}/sub_pons/${chr}_sub.vcf 

echo "create PON finished!"
"""

