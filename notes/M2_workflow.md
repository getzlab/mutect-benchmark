# M2 workflow

This assumes you have got the PON with substantial normals

## params

```
ref="/demo-mount/refs/Homo_sapiens_assembly19.fasta"
germ="/demo-mount/refs/af-only-gnomad.raw.sites.b37.vcf"
pon="/demo-mount/M2pon/makePON/AllinOne/merged_vcfs/fin_concat/AIO_merged_PON.vcf"
interval=1 # a number denoting chromosome
tumor= # a gs path
normal=# a gs path
```

### M2 for each pair

```bash
gatk --java-options -Xmx2g \
	-L 1 \ # need to scatter!
	-R $ref \
	-I $tumor \
	-I $normal \
	--germline-resources $germ \
	-pon $pon \
	-O xxx.vcf # do not use gz!
```

use `.py` to generate a yaml for all tumor-normal pairs by `python3 M2scatter_genyaml.py` which will also produce the output directory.