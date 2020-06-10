"""
Exact replicate of gatk best practice 2.5.0 for somatic indel and snv

1. only requires three vectors & the name for the PoN
    a) the paths for tumor bams
    b) the paths for normal bams
    c) the unique identifiers for each sample

2. bp_var contains the default input for best-practice, can be overwriten
    by user-input dict.

3. each workflow is defined as methods. The docstring for each function
    is the companion .sh script that canine depends.

4. The configuation is designed to be done within python console, which
    could be convenient for those who prefer to call canine within python.

Qing Zhang

"""

import yaml
from bp_vars import *
import pandas as pd
import argparse
from datetime import datetime
import itertools


class k9m2:


"""
Given three list (tumor, normal, id) and the name of the pon
0. Check that normal and pid are unique
    a) One-time script: generate interval lists
1. Generate make-pon workflow that
    a) scatter M2 tumor-only calls by normals
    b) scatter GenomicsDBImport & CreateSomaticPoN
    c) generate a script Merges the individual vcfs
       to get a single PoN named by our input.
2. Run M2 for each tumor-normal pair
    a) scatter M2 paired calls on each pair X interval
    b) merge result for all intervals
"""
    resource_dict = {
        "cpus_per_task": {
            "M2_scatter": 1,
            "M2_merge": 4,  # parallel funcotator
            "PoN_normal": 1,
            "PoN_create": 1
           }
        "mem_per_cpu": {
            "M2_scatter": 4,
            "M2_merge": 4,
            "PoN_normal": 10,  # since not splitted by interval
            "PoN_create": 20  # this requires a lot
           }
    }



   def __init__(self, tumor_list, normal_list, pid_list,
                 name_conf, default_dict = M2, resource_dict=resource_dict,
                 pon_name = "PoN_{}.vcf"):
        self.__dict__.update(default_dict)
        self.__dict__.update(name_conf)
        self.__dict__.update(resource_dict)
        self.timestamp = datetime.now().strftime('%m%d%H%M')
        self.tumor_bam = tumor_list
        self.normal_bam = normal_list
        self.pids = pid_list
        self.pon_name = pon_name.format(self.timestamp)

    @classmethod
    def from_table(cls, table_path, test=True):
        if test:
            df = pd.read_csv(table_path, sep="\t", header=None, nrows=4)
        else:
            df = pd.read_csv(table_path, sep="\t", header=None)

        # we assume the first three columns are
        # "tumor"/"normal"/"identifiers"
        df.columns = ["tumor", "normal", "pid"]
        normal_list = df["normal"].tolist()
        cls.validate(normal_list)
        cls.validate(pid_list)

        return cls(tumor_list, normal_list, pid_list)

    @staticmethod
    def validate(a_list):
        """
        The normals must be unique to be used in CreateSomaticPoN
        The pids must be unique for the file structure.
        """
        seen = set()
        for x in a_list:
            if x in seen:
                raise Exception()
            seen.add(x)
    
    # define jobs
    def make_config(self, job_name):
        # this dict will work for most jobs (may be 
        # redundant but never bad to have some entries)
        general_conf = dict(

            # reference file and all inputs
            ref_fasta = self.ref_fasta
            gnomad = self.gnomad
            ref_dict = self.ref_dict
            vfc = self.vfc
            pon = self.pon_name

            # about file structure
            subvcfs_dir=self.subvcfs_dir
            f1r2_dir=self.f1r2_dir
            tpile_dir=self.tumor_pile
            npile_dir=self.normal_pile   

        )

        # the name of interval files
        subints = ['{:0>4}'.format(i) for i in range(self.scatter_count)]
        intervals = ["{}/{}-scattered.interval_list".format(self.sub_intervals_dir,int1) for int1 in subints]

        # generate map indices for permutated pairs
        file_indices = list(range(len(self.pids)))
        intervals_indices = list(range(int(self.scatter_count)))
        comb_indices = list(itertools.product(file_indices, intervals_indices))

        # a general function that unfolds a vector by list of indices
        def unfold(uniq_vec, idx, tups=comb_indices):
            newlist = [uniq_vec[tup[idx]] for tup in tups]
            return(newlist)


        chr=unfold(subints, idx=1)

        if job_name == "M2scatter":
            base_dir = "/".join(self.NFS, self.res_folder, "paired_calls")
            cfdict = dict(
                # most essential input: itervals and paired data
                tumor_bam=unfold(self.tumor_bam, idx=0),
                normal_bam=unfold(self.normal_bam, idx=0),
                pid = unfold(self.pids, idx=0) # pid must go first!
                interval=unfold(intervals, idx=1),

                # the output structure
                out_vcf = "{}.vcf".format(chr)
                out_f1r2 = "{}.f1r2.tar.gz".format(chr )
                out_tpile = "{}.table".format(chr )
                out_npile =  "{}.table".format(chr )
            )

            cfsh = M2scatter_sh
            
        elif  job_name == "M2merge":
            cfdict = dict(
                tumor_bam = self.tumor_bam
                normal_bam = self.normal_bam
                base_dir = "/".join(self.NFS, self.res_folder, "paired_calls")
                pid = self.pids
            )
            cfsh = M2merge_sh

        elif job_name == "M2normal":
            # scatter by normals
            base_dir = "/".join(self.NFS, self.res_folder, "normal_calls", "${pid}")
            cfdict = dict(
                normal_bam = self.normal_bam
                base_dir = "/".join(self.NFS, self.res_folder, "normal_calls")
                pid = self.pids
            )
            cfsh = M2normal_sh

        elif job_name== 'createPON':
            cfdict  =dict(
                chr = subints
                base_dir = "/".join(self.NFS, self.res_folder, "normal_calls", "PON")
            )
            cfsh = createPON_sh
        else:
            print("job name must be one of [M2scatter / M2merge / M2normal / createPON]")
            raise Exception()

        return cfsh.update(general_conf), cfsh
    
    # ACTION!
    def input_to_k9conf(self, input_conf,
                              job_name,
                              cpus_per_task=1):
        input_conf, script = make_config(job_name = job_name)
        input_conf.update(heap_mem= self.mem_per_cpu[job_name]-2)
        input_key_list = list(input_conf.keys())
        staging_dir = "/".join(self.NFS, self.res_folder, self.staging)
        yy = dict(
            name=job_name,
            resources={
                "cpus_per_task": self.cpus_per_task[job_name],
                "mem_per_cpu": "{}G".format(self.mem_per_cpu[job_name])
            },
            inputs=input_conf,
            backend={"type": "Local"},
            localization={
                # staging directory will be a random folder in home
                "staging_dir": "{}/{}_{}".format(staging_dir, job_name, self.timestamp),
                "overrides": dict(zip(input_key_list, [None]*len(input_key_list)))
            }
        )
        return yy, script



if __name__ == "__main__":

    fs_dict = dict(
        # the NFS disk
        NFS="/demo-mount",
        # used for define file structure
        tpile_dir="tumor_pile",
        npile_dir="normal_pile",
        f1r2_dir="f1r2",
        subvcfs_dir="subvcfs"
        res_folder="M2bp",
        staging = "staging"
        sub_intervals_dir="/demo-mount/refs/gatk_bp_intervals"
    )

    cc = k9m2()
    make_config = cc.make_config

