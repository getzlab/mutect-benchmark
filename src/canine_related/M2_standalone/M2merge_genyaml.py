#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  6 08:59:50 2019

@author: qingzhang
"""


import os
import yaml
import glob
import datetime

now=datetime.datetime.now()
hash_time=now.strftime("%f")
# global params
filepath = "/demo-mount/refs/mpairs.tsv"
chr_l = list(range(1, 23))
ref = "/demo-mount/refs/Homo_sapiens_assembly19.fasta"
interval_list="/demo-mount/refs/interval1_22.list"
ref_dict="/demo-mount/refs/Homo_sapiens_assembly19.dict"
data_source_folder="/demo-mount/funcotator_dataSources.v1.6.20190124s"
heap_mem = 8

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
            output_dir = "/demo-mount/M2full/{}".format(pid)
            # how many files are there
            vcfs = glob.glob("{}/vcfs/*.stats".format(output_dir))
            tumor_piles = glob.glob("{}/tumor_pile/*.table".format(output_dir))
            normal_piles = glob.glob("{}/normal_pile/*.table".format(output_dir))
            f1r2s = glob.glob("{}/f1r2/*.tar.gz".format(output_dir))
            
            if len(vcfs) == 22 and len(tumor_piles) == 22 and len(normal_piles) == 22 and len(f1r2s) == 22:
                pid_l.append(pid)

print(pid_l)
print(hash_time)



yy = dict(
        name = "M2{}".format(len(pid_l)),
        inputs = {
                "mem": heap_mem,
                "pid": pid_l,
                "ref":ref,
                "interval_list":interval_list,
                "ref_dict":ref_dict,
                "data_source_folder":data_source_folder
                },
         resources={
                "cpus_per_task":1,
                "mem_per_cpu":"10G"
                },
         backend={
                "type":"Local"
                },
         localization={
                "staging_dir":"/demo-mount/canine_merge_ostr/merged-{}".format(hash_time),
                "overrides":{
                             "pid":None,
                             "ref":None,
                             "ref_dict":None,
                             "interval_list":None,
                             "data_source_folder":None
                             }
                                                                                                                                                                                             }
                                                         
        )

with open("/demo-mount/canine_merge_ostr/merge.yaml","w") as outfile:
    yaml.dump(yy, outfile, default_flow_style=False)
