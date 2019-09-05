# M2 w PON

## params

```bash
ref="/demo-mount/refs/human_g1k_v37.fasta"
germ="/demo-mount/refs/af-only-gnomad.raw.sites.b37.vcf"
sublist="/demo-mount/refs/sub_simple/*"
```

## workflow

1. split the chromosome into smaller sub intervals

   ```bash
   seq 1 22 > interval.list
   gatk SplitIntervals \
   	-R $ref \
   	-L interval.list \
   	--scatter-count 10 \
   	-O /demo-mount/refs/sub_simple
   ls /demo-mount/refs/sub_simple > input.list
   ```

2. run M2 on normals with tumor-only mode

   ```
   python3 make_pon_yaml.py
   ```

   will generate canine's yaml config for all cohorts. **NOTE: we need to specify `--max-mnp-distance 0`  for tumor-only mode, and set heapmem ~4g **

   The most script field inside this yaml is 

   ```bash
   gatk --java-options -Xmx4g Mutect2 -R $ref -I $normal --germline-resource $germ -O $outfile --max-mnp-distance 0
   ```

   ```
   canine LUAD_pon.yaml
   ```

   will run all LUAD normals

3. run GenomicsDBImport and CreatePON for all normals scattered by 1K subintervals

   ```
   python3 genomicsdb_dispatched.py --cohort SARC
   canine SARCmake_pon.yaml
   ```

   will generate a set of cohort-specific scattered PONs. Similarly, for a big PON that contains all normals we have, replace `genomicsdb_dispatched.py` with `makeBigPON.py`. Some of the normal M2 results seems corrupted - so the `all_in_onemakepon.yaml` is the current working biggest PON (~3k normals).

   To merge the scattered PONs, use 

   ```
   gatk --java-options Xmx13g GatherVcfsCloud \
   	-I input.list \
   	-O SARC_merged.vcf
   ```

   `GatherVcfsCloud` is a beta tool whose error msg is not informative. When it noted some weird line, try increase memory first (if controller node cannot suffice, try `srun -N 1 -n 2 --mem 50g gatk ...`)! If still not working, double check that output format is `.vcf` instead of `.vcf.gz`!

4. run M2 with matched tumor-normal (to be continued!)