#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 26 16:32:58 2019

bulid a yaml file for dispatching genomicsDBimport and 
CreatePON for each subintervals of interval.list

@author: qingzhang
"""
import os
import glob
import yaml
import argparse

parser = argparse.ArgumentParser(description='Example with long option names')

# parse cohort & separation by chromsome or pieces
parser.add_argument('--cohort', action="store", dest="cohort")
parser.add_argument('--bychr', action="store_true", 
                    dest="bychr", default=False)
parser.add_argument("--memG", action="store", dest="mem",
                    default=20)
parser.add_argument("--endwith", action="store", dest="endwith",
                    default="vcf.gz")

results = parser.parse_args()
cohort = results.cohort
bychr = results.bychr
mem = results.mem
endwith = results.endwith

changeDir_script = "cd ~"
RemovePreviousDB_script = "rm -rf $dbname"

loc="/demo-mount/M2pon/{}/*{}".format(cohort, endwith)
if endwith == "vcf":
    allvcfs = glob.glob(loc)
else: 
    allvcfs = [vcf[:-6] for vcf in glob.glob(loc+".stats")]

outpath = "/demo-mount/M2pon/makePON/byCohort/{}".format(cohort)
if bychr:
    allintervals = list(range(1,23))
    outfiles = ["{}/chr{}pon.vcf".format(outpath, chro) for chro in allintervals]
    dbnames = ["chr{}{}_db".format(cohort,chro) for chro in allintervals]
else:
    allintervals = glob.glob("/demo-mount/refs/sub_intervals/*")
    outfiles = ["{}/{}pon.vcf".format(outpath, os.path.basename(sub)[:4]) for sub in allintervals]
    dbnames = ["{}{}_db".format(cohort, os.path.basename(sub)[:4]) for sub in allintervals]
    
vcf_opts = ["-V "+vcf for vcf in allvcfs]
genomicsDBImport_script = "gatk --java-options -Xmx{}g GenomicsDBImport -R $ref -L $chr {} --genomicsdb-workspace-path $dbname --batch-size 50 --merge-input-intervals".format(mem, " ".join(vcf_opts))
createPon_script = "gatk --java-options -Xmx{}g CreateSomaticPanelOfNormals -R $ref --germline-resource $germ -V gendb://$dbname -O $outvcf".format(mem)
dictyaml = dict(
            name="{}-GDBI".format(cohort),
            script=[changeDir_script, RemovePreviousDB_script, genomicsDBImport_script, createPon_script],
            inputs={
                    "ref":"/demo-mount/refs/Homo_sapiens_assembly19.fasta",
                    "germ":"/demo-mount/refs/af-only-gnomad.raw.sites.b37.vcf",
                    "dbname":dbnames,
                    "chr":allintervals,
                    "outvcf":outfiles
                    },
            resources={
                    "cpus-per-task":1,
                    "mem-per-cpu":"30G"
                    },
            backend={
                    "type":"Local"
                    },
            localization={
                    "staging_dir":"/demo-mount/canine_genomicsDBI_PON/{}/".format(cohort),
                    "overrides":{
                                "dbname":None,
                                "ref":None,
                                "germ":None,
                                "outvcf":None,
                                "chr":None
                                }
                    }
        )
with open("/demo-mount/canine_genomicsDBI_PON/gdb_yamls/{}make_pon.yaml".format(cohort),"w") as outfile:
    yaml.dump(dictyaml, outfile, default_flow_style=False)

# make sure the target directory exists
import os.path
if not os.path.exists(outpath):
    os.makedirs(outpath)

