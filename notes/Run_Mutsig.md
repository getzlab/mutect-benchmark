# Run Mutsig

### Separate a table by one column

MC3 is the callset from Mutect1. The last column patient_type is used to split the file into sub tumor types.

```
$ head MC3.maf
chr	pos	gene	ref_allele	newbase	type	classification	patientttype
1	13308		T	G	IGR	SNP	TCGA-L5-A88V-01A-11D-A351-09	ESCA
1	13372		G	T	IGR	SNP	TCGA-EB-A3Y7-01A-11D-A23B-08	SKCM
1	13501		G	A	IGR	SNP	TCGA-CR-7386-01A-11D-2012-08	HNSC
1	13515		C	A	IGR	SNP	TCGA-06-0221-01A-01D-1491-08	GBM
1	13516		C	A	IGR	SNP	TCGA-06-0221-01A-01D-1491-08	GBM
1	13516		C	A	IGR	SNP	TCGA-EO-A3AY-01A-12D-A19Y-09	UCEC
1	13521		C	A	IGR	SNP	TCGA-D1-A103-01A-11D-A10M-09	UCEC
1	13557		G	A	IGR	SNP	TCGA-ZN-A9VS-01A-11D-A39R-32	MESO
1	14612		G	A	IGR	SNP	TCGA-OR-A5LJ-01A-11D-A29I-10	ACC
```

To split the files and add header to each small chunks - 

```bash
#!/bin/bash
addname="head_"

# make a header file
rm *maf
head -n 1 ../MC3.maf > header.txt 
awk '{print>$NF ".maf"}' ../MC3.maf

filename="*maf"
for file in $filename
do
        echo "start with $file"
        cat header.txt "$file" >> "$addname$file"
        rm "$file"
        mv "$addname$file" "$file"
done

```



### Run Mutsig on cluster

I am using the Broad cluster "cga2", which uses GridEngine as job scheduler. 

- load the Matlab runtime

`reuse .matlab_2016a_mcr`

Put it under `reuse -q .matlab-2016a` in bashrc!

- make a text file to dispatch on UGER

```bash
#!/bin/bash
filenames="sub_ttypes_mafs/*.maf"
for file in $filenames
do 
        basename="${file##*/}"
        ttype=${basename%.*}
        mkdir "results/$ttype"
        echo "./MutSig2CV_v3_coding $file results2/$ttype params.txt"

done
```

Note that `params.txt` are used to avoid permission issues. *Thank you Julian!*

However there are some problematic nodes  - it can run on some nodes but not on the others. To repeat the submitting jobs that are not finished, we need to figure out what is finished and what is not. A reasonable proxy is the "number of hard link" in `ls -al` of the result folder, when it equals to 793 then it is finished. Under this logic, we can write - 

```bash
#!/bin/bash
ls -l  ../ | awk '{ if ($5 == 171) { print $NF } }' > try.txt
cat try.txt |  
while IFS= read -r line
do
        ttype="${line%/*}"
        echo "../../MutSig2CV_v3_coding ../../sub_ttypes_mafs/$ttype.maf ../../results2/$ttype ../../params.txt"
done
```

However a better proxy is to check if **finished.txt** exists in each output folder. 

```
%TODO%
```

- Dispatch on cluster

in a scheduler-independent way:

```
qsubb.sh <name of runs file> --pre=". broad/software/scripts/useuse;reuse .matlab_2016a_mcr" -V -cwd -N mutsig -j y -o mutsig_log -l h_vmem=4G -l h_rt=2:00:00
```

For some cancers - LUAD etc, we need 20h to run. Most of cancers tend to finish by 5h. `qstat` can be used to track job status.