# Somatic Variant Calling benchmark with Mutect and GATK Mutect2



This repo has visualization code / thesis / slides for defense for my thesis on somatic SNV caller comparison between Mutect1 and GATK Mutect2. The code to reproduce Mutect2 (in GATK 4.1.4.0) in Snakemake can be found at [smk-m2](https://github.com/getzlab/smk-m2). 

Major finding is summarized as below:

- Significance analysis is crucial to find scientifically important patterns that drives the difference in call sets. The seemingly concordant call sets might produce very divergent sets of significant genes.
- M2's postfilters are able to filter some mapping artifacts and base quality artifacts compared to M2, however some of the realigned haplotypes could be false positives
- M1 tends to be conservative at regions with more variations, thus would have missed SNVs near SVs / indels, although the portion affected is negligible.

I would like to thank Gaddy, Julian and Chip for their generous input to this work.  I would also like to thank GATK team (David Benjamin, Soo Hee Lee) for their nice documentation and clear explanation on GATK forum. The computational resource is supported by GDAN grant.





