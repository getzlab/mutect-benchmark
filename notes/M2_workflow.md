# M2 workflow

This assumes you have got the PON with substantial normals

## Params

```bash
ref="/demo-mount/refs/Homo_sapiens_assembly19.fasta"
# means the only info field is the "alleic freq"
germ="/demo-mount/refs/af-only-gnomad.raw.sites.b37.vcf" 
pon="/demo-mount/M2pon/makePON/AllinOne/merged_vcfs/fin_concat/AIO_merged_PON.vcf"
interval=[1-22,"X","Y"] # a number denoting chromosome
tumor="" # a gs path
normal="" # a gs path
```

## Prep: common contamination

- Get `variants_for_contamination.vcf`

```bash
time srun -N 1 -n 4 --mem 80g \
	gatk --java-options -Xmx70g SelectVariants \
	-V af-only-gnomad.raw.sites.b37.vcf \
	-select-type SNP \
	-restrict-alleles-to BIALLELIC -select "AF > 0.05" \
	-O gnomad_var_for_contamination/variants_for_contamination.vcf
```

## Scatter: M2 for each pair by chr

Run `python3 M2scatter_genyaml.py` and `canine XXX.yaml` for each pair, whose script is essentially - 

```bash
gatk --java-options -Xmx2g \ # so total 
	-L 1 \ # need to scatter!
	-R $ref \
	-I $tumor \
	-I $normal \
	--germline-resources $germ \
	-pon $pon \
	-O xxx.vcf # do not use gz!
	--f1r2-tar-gz $f1r2
```

use `.py` to generate a yaml for all tumor-normal pairs by `python3 M2scatter_genyaml.py` which will also produce the output directory.

*[Note on directory structure]*

* staging dir = `/demo-mount/canine_m2full/{pid}`
* instruction yaml dir = `/demo-mount/canine_M2full_ostr`
* output dir = `/demo-mount/M2full/{pid}`
  * f1r2 dir = `{output_dir}/f1r2/`
  * vcfs/stats_dir = `{output_dir}/vcfs/`
  * normal_pile = `{output_dir}/normal_pile/`
  * tumor_pile=`{output_dir}/tumor_pile/`

## Gather by pair

### Merge outputs

* merge vcfs

  ```
  gatk --java-options "-Xmx${command_mem}m" MergeVcfs \
      -I ${sep=' -I ' input_vcfs} \
      -O ${output_vcf}
  ```

* merge stats

  ```
  gatk --java-options "-Xmx${command_mem}m" MergeMutectStats \
      -stats ${sep=" -stats " stats} \
      -O merged.stat
  ```

* merge pileups (for tumor and normal separately)

  ```
  gatk --java-options "-Xmx${command_mem}m" GatherPileupSummaries \
      --sequence-dictionary ${ref_dict} \
      -I ${sep=' -I ' input_tables} \
      -O ${output_name}.tsv
  gatk --java-options "-Xmx${command_mem}m" GatherPileupSummaries \
      --sequence-dictionary ${ref_dict} \
      -I ${sep=' -I ' input_tables} \
      -O ${output_name}.tsvl
  ```

### learn orientation bias model

```bash
gatk --java-options "-Xmx${command_mem}m" LearnReadOrientationModel \
    -I ${sep=" -I " f1r2_tar_gz} \
    -O "artifact-priors.tar.gz"
```

### Calculate contamination

```bash
gatk --java-options "-Xmx${command_mem}m" CalculateContamination \
    -I ${tumor_pileups} \
    -O contamination.table \
    --tumor-segmentation segments.table \
    ${"-matched " + normal_pileups}
```

### Filter Mutect Calls

```bash
gatk --java-options "-Xmx${command_mem}m" FilterMutectCalls -V ${unfiltered_vcf} \
    -R ${ref_fasta} \
    -O ${output_vcf} \
    ${"--contamination-table " + contamination_table} \
    ${"--tumor-segmentation " + maf_segments} \
    ${"--ob-priors " + artifact_priors_tar_gz} \
    ${"-stats " + mutect_stats} \
    --filtering-stats filtering.stats \
    ${m2_extra_filtering_args}
```

### Funcotate

```bash
gatk --java-options "-Xmx${command_mem}m" Funcotator \
    --data-sources-path $DATA_SOURCES_FOLDER \
    --ref-version ${reference_version} \
    --output-file-format ${output_format} \
    -R ${ref_fasta} \
    -V ${input_vcf} \
    -O ${output_file} \
    ${interval_list_arg} ${default="" interval_list} \
    --annotation-default normal_barcode:${default="Unknown" control_id} \
    --annotation-default tumor_barcode:${default="Unknown" case_id} \
    --annotation-default Center:${default="Unknown" sequencing_center} \
    --annotation-default source:${default="Unknown" sequence_source} \
	${"--transcript-selection-mode " + transcript_selection_mode} \
	${transcript_selection_arg}${default="" sep=" --transcript-list " transcript_selection_list} \
	${annotation_def_arg}${default="" sep=" --annotation-default " annotation_defaults} \
	${annotation_over_arg}${default="" sep=" --annotation-override " annotation_overrides} \
	${excluded_fields_args}${default="" sep=" --exclude-field " funcotator_excluded_fields} \
	${filter_funcotations_args} \
	${extra_args_arg}
```

## Refs

[Mutect2-WDL-resources](https://github.com/broadinstitute/gatk/blob/master/scripts/mutect2_wdl/mutect_resources.wdl) prepares the *variants_for_contamination* file that used for `getPileUpSummary` and `CalculateContanmiantion`

