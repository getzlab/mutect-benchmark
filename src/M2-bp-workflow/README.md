# M2-best-practice in Canine



The best-practices suggested variables are stored as python dict in `bp_vars.py` and will be read and merged into input configs (also python dicts, but with user-input e.g. tumor/normal paths) to generate yaml and scripts.



## M2-scatter

```
python3 Mutect2 -pfile mpairs.txt -test
```

which will generate a script named `k9_sca_m2.sh` and a `k9_sca_m2.yaml` config for canine. Then run by

```
canine k9_sca_m2.yaml --script k9_sca_m2.sh
```



## CreatePON



## M2-gather

