#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 25 15:47:18 2019

create a set of yaml file to run all 
tumor-only M2 calls to build PON.

# caveats: 
1. do not specify nodelist for preempt machines!
2. must use --max-mnp-distance 0
3. use None to create a quoteless "null" options for each vcf/bam
4. 4g heap size and 6g mem alloc will work for (supposedly) all samples

@author: qingzhang
"""


# for some reason my local pandas does not work, should definately go with pandas!
# the file is generated by a helper script in igv_remove repo under dalmatian helper.
filepath="/home/qingzhang/Documents/igv_remote/helpers/mpairs.tsv"
with open(filepath) as f:

    tumor_l = []
    normal_l = []
    pid_l = []
    cohort_l = []
    for l in f:
        tumor, normal, pid, cohort = l.strip().split("\t")
        tumor_l.append(tumor)
        normal_l.append(normal)
        pid_l.append(pid)
        cohort_l.append(cohort)




import yaml
TCGA_abbr=["ACC", "BLCA", "LGG", "BRCA", 
           "CESC", "CHOL", "COAD","THYM",
           "ESCA", "GBM", "HNSC", "KICH",
           "KIRC", "KIRP", "LIHC", "LUAD", 
           "LUSC", "DLBC", "MESO", "OV",
           "PAAD", "PCPG", "PRAD", "READ",
           "SKCM", "STAD", "TGCT", "SARC",
           "THCA", "UCS", "UCEC", "UVM", "MESO"]

# IMPORTANT! must set --max-mnp-distance 0 for tumor only mode
for cohort in TCGA_abbr:
    normals = [n for n in normal_l if cohort in n]
    pids = [pid for (pid, n) in zip(pid_l, normal_l) if cohort in n]
    pon_dict =dict(
        name=cohort+"-PON", # need 4g heap size
        script=["gatk --java-options -Xmx4g Mutect2 -R $ref -I $normal --germline-resource $germ -O $outfile --max-mnp-distance 0"],
        resources={
            "cpus_per_task":1,
            "mem_per_cpu":"6G" # add some buffer mem size
        },
        inputs={
            "ref":"/demo-mount/refs/Homo_sapiens_assembly19.fasta",
            "germ":"/demo-mount/refs/af-only-gnomad.raw.sites.b37.vcf",
            "normal":normals,
            "outfile":["/demo-mount/M2pon/"+ cohort + "/" + pid + "_normal.vcf.gz" for pid in pids]
            },
        backend={
            "type":"Local",
            },
        localization={
            "staging_dir":"/demo-mount/canine_pon/"+cohort+"/",
            "overrides":{ 
			# set all params to None so that canine will not attempt 
			# to copy files / create simlinks that disrupts the dependency 
			# with its index files
                    "normal":None,
                    "ref":None, 
                    "germ":None,
                    "outfile":None
                    }
            }
        )
    with open(cohort + "_pon.yaml","w") as outfile:
        yaml.dump(pon_dict, outfile, default_flow_style=False)