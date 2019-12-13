# Get bams

## from Google bucket

If the google bucket for the TCGA legacy bams can be accessed (with your ERA commons account), and assume you have a working installation of `dalmation`, please use

```
python3 dalmation_helper.py
```

which will generate a table of paired gs urls in a file named `mpairs_update.tsv`. Since most GATK tools can stream google buckets without localization, direct using such bucket would be optimal practice.


### from NCI data commons

Otherwise, in case the google bucket is not accessible due to any reasons, please use `get_paired.R` which utilizes a R API for GDC `TCGAbiolinks` to get **UUID** for each patient (this is similar to a manifest file). 

This table, composing columns for `cohort`/`patient_id`/`tumor_uuid`/`normal_uuid`, will be the input for [snakemake pipeline with *localization* step](https://github.com/hurrialice/smk-m2). In short, NCI data commons will help generate signed URLs for both `bam` and `bai`. 
