#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 19 07:57:01 2019
Mutect2 on per interval basis
@author: qingzhang

    
"""


gnomad="gs://gatk-best-practices/somatic-b37/af-only-gnomad.raw.sites.vcf"
tumor_bam_index="gs://fc-a6cde48f-6248-4428-a950-a256d46a828b/ichor_patients/7_FC19663997_HFW3TBBXX.7.aligned.duplicates_marked.bai"
tumor_bam="gs://fc-a6cde48f-6248-4428-a950-a256d46a828b/ichor_patients/7_FC19663997_HFW3TBBXX.7.aligned.duplicates_marked.bam"
normal_bam_index="gs://fc-a6cde48f-6248-4428-a950-a256d46a828b/healthy_donor_wgs/FC15884027.markDuplicates.bai"
normal_bam="gs://fc-a6cde48f-6248-4428-a950-a256d46a828b/healthy_donor_wgs/FC15884027.markDuplicates.bam"
ref_fasta="/demo-mount/refs/Homo_sapiens_assembly19.fasta"
pon="gs://fc-a6cde48f-6248-4428-a950-a256d46a828b/pon/wgs_125_plus_targeted_dream.vcf"
vfc="gs://gatk-best-practices/somatic-b37/small_exac_common_3.vcf"
command_mem=4
# varying: inteval, vcf output and f1r2 out
subints = ["/demo-mount/for_Ziao/intevals/sublist/" + '{:0>4}'.format(i) + "-scattered.interval_list" for i in range(1000)]



scripts = """
set -e
mkdir -p /demo-mount/for_Ziao/res/
mkdir -p /demo-mount/for_Ziao/res/f1r2
mkdir -p /demo-mount/for_Ziao/res/tumor_piles
mkdir -p /demo-mount/for_Ziao/res/normal_piles
mkdir -p /demo-mount/for_Ziao/res/subvcfs

command_mem=4
gatk --java-options "-Xmx${command_mem}g" GetSampleName -R $ref \
    -I ${tumor_bam} -O tumor_name.txt \

tumor_command_line="-I ${tumor_bam} -tumor `cat tumor_name.txt`"
cat tumor_name.txt


gatk --java-options "-Xmx${command_mem}g" GetSampleName -R $ref \
    -I ${normal_bam} -O normal_name.txt \
    -encode normal_command_line="-I ${normal_bam} -normal `cat normal_name.txt`"
cat normal_name.txt


gatk --java-options "-Xmx${command_mem}g" Mutect2 \
    -R $ref \
    $tumor_command_line \
    $normal_command_line \
    --germline-resource  $gnomad \
    -pon $pon \
    -L  $subint \
    -O "$output_vcf \
    --f1r2-tar-gz $f1r2.tar.gz 
echo finish Mutect2 ----------



gatk --java-options "-Xmx${command_mem}g" GetPileupSummaries \
    -R $ref -I ${tumor_bam} \
    --interval-set-rule INTERSECTION \
    -L $subint \
    -V $vfc -L $vfc \
    -O $tumor-pileups.table
echo Getpiles tumor ---------

gatk --java-options "-Xmx${command_mem}g" GetPileupSummaries \
    -R $ref -I ${normal_bam} \
    --interval-set-rule INTERSECTION \
    -L $subint \
    -V $vfc -L $vfc \
    -O $normal-pileups.table
echo Getpiles normal --------

echo allfin!

"""




yy = dict(
         name="ziao",
    resources={
            "cpus_per_task":2,
            "mem_per_cpu":"5G"
            },
    inputs={
            "ref":ref_fasta,
            "gnomad":gnomad,
            "tumor_bam":tumor_bam,
            "normal_bam":normal_bam,
            "subint": subints,
                #0499-scattered.interval_list
            "output_vcf":["/demo-mount/for_Ziao/res/subvcfs/"+ int1 + "_out.vcf" for int1 in subints],
            "f1r2.tar.gz":["/demo-mount/for_Ziao/res/f1r2/" + int1 + "_f1r2.tar.gz" for int1 in subints],
            "pon":pon,
            "tumor-pileups.table":["/demo-mount/for_Ziao/res/tumor_piles/" + int1 + "_f1r2.tar.gz" for int1 in subints],
            "normal-pileups.table":["/demo-mount/for_Ziao/res/normal_piles/" + int1 + "_f1r2.tar.gz" for int1 in subints]
            },
    backend={
            "type":"Local"
            },
    localization={
            "staging_dir":"/demo-mount/for_Ziao/staging0/",
            "startegy":"NFS",
            "overrides":{
                    "tumor_bam":None,
                    "normal_bam":None,
                    "ref":None,
                    "gnomad":None,
                    "pon":None,
                    "vfc":None,
                    "subint":None
                    }
            }
        
        
        )


