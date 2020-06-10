set -e
mkdir -p /demo-mount/for_Ziao/res/
mkdir -p /demo-mount/for_Ziao/res/f1r2
mkdir -p /demo-mount/for_Ziao/res/tumor_piles
mkdir -p /demo-mount/for_Ziao/res/normal_piles
mkdir -p /demo-mount/for_Ziao/res/subvcfs

command_mem=4
gatk --java-options "-Xmx${command_mem}g" GetSampleName -R $ref \
    -I ${tumor_bam} -O tumor_name.txt \
    -encode tumor_command_line="-I ${tumor_bam} -tumor `cat tumor_name.txt`"
cat tumor_name.txt


gatk --java-options "-Xmx${command_mem}g" GetSampleName -R $ref \
    -I ${normal_bam} -O normal_name.txt \
    -encode normal_command_line="-I ${normal_bam} -normal `cat normal_name.txt`"
cat normal_name.txt


gatk --java-options "-Xmx${command_mem}g" Mutect2 \
    -R $ref \
    $tumor_command_line \
    $normal_command_line \
    --germline-resource  $gnomad \
    -pon $pon \
    -L  $subint \
    -O "$output_vcf \
    --f1r2-tar-gz $f1r2.tar.gz 
echo finish Mutect2 ----------



gatk --java-options "-Xmx${command_mem}g" GetPileupSummaries \
    -R $ref -I ${tumor_bam} \
    --interval-set-rule INTERSECTION \
    -L $subint \
    -V $vfc -L $vfc \
    -O $tumor-pileups.table
echo Getpiles tumor ---------

gatk --java-options "-Xmx${command_mem}g" GetPileupSummaries \
    -R $ref -I ${normal_bam} \
    --interval-set-rule INTERSECTION \
    -L $subint \
    -V $vfc -L $vfc \
    -O $normal-pileups.table
echo Getpiles normal --------

echo allfin!
