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

chr_l = list(range(1, 23))+["X","Y"]
basedir = "/demo-mount/M2full"
stagedir = "/demo-mount/canine_m2full"
ostrdir = "/demo-mount/canine_M2full_ostr"

if os.path.exists(stagedir):
    raise Exception("staging dir exist!")

for pid, tumor, normal in zip(pid_l, tumor_l, normal_l):
    scatter_dir = "/".join([basedir,pid])
    if not os.path.exists(scatter_dir):
        os.makedirs(scatter_dir)
    
    yy = dict(
            name="M2full"+pid,
            script=["gatk --java-options -Xmx2g Mutect2 -R $ref -I $tumor -I $normal -L $chr --germline-resource $germ -O $outfile -pon $pon --f1r2-tar-gz $f1r2"],
            resources={
                    "cpus_per_task":1,
                    "mem_per_cpu":"4G"
                    },
            inputs={
                    "ref":"/demo-mount/refs/Homo_sapiens_assembly19.fasta",
                    "germ":"/demo-mount/refs/af-only-gnomad.raw.sites.b37.vcf",
                    "tumor":tumor,
                    "normal":normal,
                    "chr": chr_l,
                    "pon": "/demo-mount/M2pon/makePON/AllinOne/merged_vcfs/fin_concat/AIO_merged_PON.vcf",
                    "f1r2": ["{}/{}_chr{}-f1r2.tar.gz".format(scatter_dir,pid,chr) for chr in chr_l],
                    "outfile":["{}/{}_chr{}_unfilter.vcf".format(scatter_dir,pid,chr) for chr in chr_l]
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
    with open("{}/{}-bychr.yaml".format(ostrdir,pid),"w") as outfile:
        yaml.dump(yy, outfile, default_flow_style=False)
    
