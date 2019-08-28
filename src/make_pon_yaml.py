#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 25 15:47:18 2019

@author: qingzhang
"""

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
           "CESC", "CHOL", "LCML", "COAD",
           "ESCA", "GBM", "HNSC", "KICH",
           "KIRC", "KIRP", "LIHC", "LUAD", 
           "LUSC", "DLBC", "MESO", "OV",
           "PAAD", "PCPG", "PRAD", "READ",
           "SKCM", "STAD", "TGCT", "SARC",
           "THCA", "UCS", "UCEC","UVM"]

for cohort in TCGA_abbr:
    normals = [n for n in normal_l if cohort in n]
    pids = [pid for (pid, n) in zip(pid_l, normal_l) if cohort in n]
    pon_dict =dict(
        name=cohort+"-PON",
        script=["gatk Mutect2 -R $ref -I $normal --germline-resource $germ -O $outfile --max-mnp-distance 0"],
        resources={
            "cpus_per_task":1,
            "mem_per_cpu":"4G"
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
                    "normal":None,
                    "ref":None,
                    "germ":None,
                    "outfile":None
                    }
            }
        )
    with open(cohort + "_pon.yaml","w") as outfile:
        yaml.dump(pon_dict, outfile, default_flow_style=False)
