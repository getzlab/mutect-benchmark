#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  5 19:24:21 2019

For each scattered pieces, 
    get Mutect2 calls for a pair of tumor-normal
    estimate f1r2
    get pileup summaries for tumor
    get pileup summaries for normal
    
@author: qingzhang
"""
import yaml

filepath = "/demo-mount/refs/mpairs.tsv"

# filepath="/home/qingzhang/Documents/igv_remote/helpers/mpairs.tsv"
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

batch_size=200
n_batch=len(pid_l) // batch_size+ 1
chr_l = list(range(1,23))

# files
vfc="/demo-mount/refs/gnomad_var_for_contamination/variants_for_contamination.vcf"
ref="/demo-mount/refs/Homo_sapiens_assembly19.fasta"
germ="/demo-mount/refs/af-only-gnomad.raw.sites.b37.vcf"
pon="/demo-mount/M2pon/makePON/AllinOne/merged_vcfs/fin_concat/AIO_merged_PON.vcf"


for i in list(range(1, n_batch+1)):
    id_start=(i-1)*batch_size
    id_end=min(id_start+batch_size, len(tumor_l))
    print((id_start, id_end))
    tumor_in_batch=tumor_l[id_start:id_end]
    normal_in_batch=normal_l[id_start:id_end]
    pid_in_batch=pid_l[id_start:id_end]

    chr_n=22
    yy = dict(
            name="m2scatter-{}".format(i),
            resources={
                    "cpus_per_task":1,
                    "mem_per_cpu":"4G"
                    },
                    
            inputs={
                    "pid":[item for item in pid_in_batch for i in range(chr_n)],
                    "tumor":[item for item in tumor_in_batch for i in range(chr_n)],
                    "normal":[item for item in normal_in_batch for i in range(chr_n)],
                    "chr":list(range(1,23))*(id_end-id_start),
                    "ref":ref,
                    "germ":germ,
                    "vfc":vfc,
                    "pon": pon
                    },
            backend={
                    "type":"Local",
                    },
            localization={
                    "staging_dir":"/demo-mount/canine_M2scatter/batch_{}".format(i),
                    "overrides":{
                            "pid":None,
                            "chr":None,
                            "tumor":None,
                            "normal":None,
                            "ref":None,
                            "germ":None,
                            "pon":None,
                            "vfc":None
                            }
                    } 
            )
    with open("/demo-mount/canine_M2scatter/yamls_p{}/batch-{}.yaml".format(batch_size, i),"w") as outfile:
        yaml.dump(yy, outfile, default_flow_style=False)



    
