#!/bin/bash
# The companion script


gatk="gatk --java-options -Xmx2g"
outdir="/demo-mount/M2full/${pid}"
tpile="${outdir}/tumor_pile/tumor-${chr}-pileups.table"
npile="${outdir}/normal_pile/normal-${chr}-pileups.table"
f1r2="${outdir}/f1r2/${chr}-f1r2.tar.gz"
unfilter="${outdir}/vcfs/chr${chr}_unfilter.vcf"

set -eo pipefail
mkdir -p $outdir/tumor_pile
mkdir -p $outdir/normal_pile
mkdir -p $outdir/vcfs
mkdir -p $outdir/f1r2


echo -------------- mutect2 ----------------------
$gatk Mutect2 \
    -R $ref \
    -I $tumor \
    -I $normal \
    -L $chr \
    --germline-resource $germ \
    -pon $pon \
    --f1r2-tar-gz $f1r2 \
    -O $unfilter

echo -------------------------get tumor piles ------------
$gatk GetPileupSummaries \
    -R $ref \
    -I $tumor \
    -L $chr \
    --interval-set-rule INTERSECTION \
    -V $vfc \
    -L $vfc \
    -O $tpile

echo -----------------------get normal piles--------------------
$gatk GetPileupSummaries \
    -R $ref \
    -I $normal \
    -L $chr \
    --interval-set-rule INTERSECTION \
    -V $vfc \
    -L $vfc \
    -O $npile
