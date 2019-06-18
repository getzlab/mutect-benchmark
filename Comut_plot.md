# CoMut plot

Based on the Mutsig outputs we want something as this - 

![](<http://gdac.broadinstitute.org/runs/awg_gbm__latest/mutsignozzlereport2.0/GBM-TP_coMut.png>)

*(credit to <http://gdac.broadinstitute.org/runs/awg_gbm__latest/mutsignozzlereport2.0/nozzle.html>)*

I made minor edits to make it compatible with current mutsig outputs. (without allelic fraction plot)

- Load R

```
use R-2.15
```

- Run R executables by - 

```bash
Rscript ./coMut.R -v -o . --sort.by.mutation.status -a <tag> -p <mutsig output directory> -s sig_genes.txt -c patient_counts_and_rates.txt -m final_analysis_set.maf -q 0.1 --firehose.mutsig.mutcategs mutcategs.txt
```

- To generate comut plot for every cancer we have finished

```bash
#!/bin/bash
ls -l ../results2 | awk '{if ($5 == 793) {print $NF}}' > finished.txt
echo "#!/bin/bash"
cat finished.txt |
while IFS= read -r line
do 
        ttype="${line%/*}"
        echo "Rscript ./coMut.R -v -o . --sort.by.mutation.status -a $ttype -p ../results2/$ttype -s sig_genes.txt -c patient_counts_and_rates.txt -m final_analysis_set.maf -q 0.1 --firehose.mutsig.mutcategs mutcategs.txt"
done
```

Note that the name of text files are invariant across all tumor types.

and run the command file you just generated.