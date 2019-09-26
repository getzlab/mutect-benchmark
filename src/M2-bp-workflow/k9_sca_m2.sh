
#!/bin/bash
set -eo pipefail


# define the format of output files
out_tpile="${resdir}/${tpile}/${chr}-tpile.table"
out_npile="${resdir}/${npile}/${chr}-npile.table"
out_f1r2="${resdir}/${f1r2}/${chr}-f1r2.tar.gz"
out_vcf="${resdir}/${subvcfs}/${chr}_unfilter.vcf"

# make file structure
mkdir -p $resdir
mkdir -p $resdir/$f1r2
mkdir -p $resdir/$tpile
mkdir -p $resdir/$npile
mkdir -p $resdir/$subvcfs

command_mem=$heap_mem
gatk --java-options "-Xmx${command_mem}g" GetSampleName -R $ref_fasta     -I ${tumor_bam} -O tumor_name.txt

tumor_command_line="-I ${tumor_bam} -tumor `cat tumor_name.txt`"
cat tumor_name.txt


gatk --java-options "-Xmx${command_mem}g" GetSampleName -R $ref_fasta     -I ${normal_bam} -O normal_name.txt

normal_command_line="-I ${normal_bam} -normal `cat normal_name.txt`"
cat normal_name.txt


gatk --java-options "-Xmx${command_mem}g" Mutect2     -R $ref_fasta     $tumor_command_line     $normal_command_line     --germline-resource $gnomad     -pon $default_pon     -L  $interval     -O $out_vcf     --f1r2-tar-gz $out_f1r2
echo finish Mutect2 ----------



gatk --java-options "-Xmx${command_mem}g" GetPileupSummaries     -R $ref_fasta -I ${tumor_bam}     --interval-set-rule INTERSECTION     -L $interval     -V $vfc     -L $vfc     -O $out_tpile
echo Getpiles tumor ---------

gatk --java-options "-Xmx${command_mem}g" GetPileupSummaries     -R $ref_fasta -I ${normal_bam}     --interval-set-rule INTERSECTION     -L $interval     -V $vfc     -L $vfc     -O $out_npile
echo Getpiles normal --------

echo allfin!

    