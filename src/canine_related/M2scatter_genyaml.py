#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  5 19:24:21 2019

@author: qingzhang
"""
import os
import yaml

filepath = "/demo-mount/refs/mpairs.tsv"
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


basedir = "/demo-mount/M2full"
stagedir = "/demo-mount/canine_m2full"
ostrdir = "/demo-mount/canine_M2full_ostr"

if os.path.exists(stagedir):
    raise Exception("staging dir exist!")

for pid in pid_l:
    scatter_dir = "/".join([basedir,pid])
    if not os.path.exists(scatter_dir):
        os.makedirs(scatter_dir)
    
    yy = dict(
            name="M2full"+pid,
            script=["gatk --java-options -Xmx2g Mutect2 -R $ref -I $tumor -I $normal -L $chr --germline-resource $germ -O $outfile -pon $pon"],
            resources={
                    "cpus_per_task":1,
                    "mem_per_cpu":"2G"
                    },
            inputs={
                    "ref":"/demo-mount/refs/Homo_sapiens_assembly19.fasta",
                    "germ":"/demo-mount/refs/af-only-gnomad.raw.sites.b37.vcf",
                    "tumor":tumor_l,
                    "normal":normal_l,
                    "chr": list(range(1,23)),
                    "pon": "/demo-mount/M2pon/makePON/AllinOne/merged_vcfs/fin_concat/AIO_merged_PON.vcf",
                    "outfile":["{}/{}_chr{}_somatic.vcf".format(scatter_dir,pid,chr) for chr in list(range(23))]
            },
    backend={
            "type":"Local",
            },
    localization={
            "staging_dir":stagedir,
            "overrides":{
                    "tumor":None,
                    "normal":None,
                    "ref":None,
                    "germ":None,
                    "pon":None,
                    }
            }
    )
    with open("{}/canine_con_full.yaml".format(ostrdir),"w") as outfile:
        yaml.dump(yy, outfile, default_flow_style=False)
    
