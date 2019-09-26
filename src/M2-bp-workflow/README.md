# M2-best-practice in Canine



The best-practices suggested variables are stored as python dict in `bp_vars.py` and will be read and merged into input configs (also python dicts, but with user-input e.g. tumor/normal paths) to generate yaml and scripts.

## What is different?

* There should be a `gatk GetSampleName` call before M2 that marks the sample name for tumor and normal as suggested [here](https://github.com/gatk-workflows/gatk4-somatic-snvs-indels/blob/2.4.0/mutect2_nio.wdl#L606). When `-normal` argument is not provided, M2 will run as if it is multiple tumor samples from the same individual and calls mutation as if there was no matched normal. **This is a real *bug* and I apologize for this.** Further readings [here](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.2.0/org_broadinstitute_hellbender_tools_walkers_mutect_Mutect2.php)

  ```bash
  gatk Mutect2 \
       -R reference.fa \
       -I tumor.bam -I normal.bam \
       -normal normal_sample_name \ # this was missing!
       --germline-resource af-only-gnomad.vcf.gz \
       --panel-of-normals pon.vcf.gz \
       -O somatic.vcf.gz
  ```

* When creating the PoN it is suggested to run M2 without `--germline-resource`  under tumor only mode as (since there is only one sample, no need for `-tumor` or `-normal`):

  ```bash
  # -- scatter by sample
  gatk Mutect2 \ # no germline arg here
      -R reference.fa \
      -I sample.bam \
      -max-mnp-distance 0 \ # uses this only when making PON
      -O input.vcf
  
  # --scatter by interval
  gatk GenomicsDBImport \
      --genomicsdb-workspace-path gendb://pon_db \
      -R reference.fa \
      --batch-size 50 --merge-input-intervals \ # make things faster
      -V ${sep=' -V ' input_vcfs}
      
  gatk CreateSomaticPanelOfNormals 
  	-R reference.fa \
  	--germline-resource ${gnomad} \ # germline arg only appears here
      -V gendb://pon_db \
      -O single_sample_pon.vcf 
  ```

* To merge single sample vcfs into a PON by `MergeVcfs` (previously used `GatherVcfsCloud` which cannot regenerate index, not sure if it really matters here).

  ```bash
  gatk MergeVcfs -I ${sep=' -I ' input_vcfs} -O ${output_vcf}
  ```

* an extra argument for `SplitIntervals` for `--subdivision-mode`

  ```bash
  gatk SplitIntervals
  	-R ${ref_fasta} \
      -L  intervals \
      -scatter ${scatter_count} \
      -O interval-files \
      --subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION 
      --min-contig-size 1000000 # %need double check%
  ```

* There is a new tool `FilterAlignmentArtifacts` but to my knowledge this has not been incorporated into best-practice workflow (default not to run this). 
* For `Funcotator`  with output format set to `MAF` , best-practice uses `touch $outpt_maf_index` in the end.

## BP: M2-scatter

```
python3 Mutect2 -pfile mpairs.txt -test
```

which will generate a script named `k9_sca_m2.sh` and a `k9_sca_m2.yaml` config for canine. Then run by

```
canine k9_sca_m2.yaml --script k9_sca_m2.sh
```



## BP: CreatePON



## BP: M2-gather

