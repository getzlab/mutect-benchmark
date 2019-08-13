# Fishy genes summarized 

## General

Most of the questionable calls come from low-coverage regions and limited number of supporting reads.  Many lies within low mappability reads - but I am not sure if there is a way that we can ever safely call mutation from such region (where there is lots of repetitive sequences / homologs).

* how do we define "enrichment of truncating events"? I think maybe there will be some models behind - assuming homogeneous base mutation prob.
  * PSIP1 & ATAD5 in ACC
  * MAP2K4 in BRCA
  * COL11A1 in STAD
  * EIF4G3 in SARC
  * NF1 for LUSC
* some gene lacks coverage in tumor portal
  * General question: the definition for transcripts used in tumorportal is for both intron and exon or exon only?
  * SMTNL2 in BRCA
    * ![](/home/qingzhang/.config/Typora/typora-user-images/1564502488262.png)
    * This region (4475311-4475900) is present in both tumor and normal, however i cannot find the corresponding transcript in refseq track.
    * The first exon lies really far - I do not know exactly how first several bases disappear - at least the coverage is extremely low.
  * ADCY1 in HNSC
    * I am suspecting if tumor portal x-axis includes intron region? There is no reads at all in first 100bp and there is no supporting refseq exons.
  * PDGFRA in GBM
    * This one I am not very sure - four isoforms spans many exon-introns. 
* How to faithfully call mutation for regions with low mappability? Do we need a separate algorithm for this?
* When we talk about "alignment artifact" can we safely predict what is the correct way of alignment so that this alternation will not be called? I remembered you showed me an example where gap should prioritize mis alignments - but cannot recall how you did that.
* Refseq track: why there are reads aligned to regions where there is no exon notion? (and some even has no intron notion...) How should we deal with those reads?



## Revisit previous identified fishy genes 

> copied from the update with Gaddy:

This is the summary of fishy genes found by Mutsig - we try to dissect them into either Mutsig False Positives (Mutsig FPs in short, meaning that the sequencing calls are okay but Mutsig falsely regard them as true signal) and sequencing artifacts (perhaps the mutation does not even exist, i.g. due to problems with upstream mutation caller).

- Mutsig FPs
  \- Too few mutations to confidently classify (C12orf65, CYHR1, PTMA [one mutation!])

  \- Indistinct pattern of missense mutations (RBM26, ARHGAP4)
  \- dN/dS ~ 1 (FAT4, ASPM)

- Mutation calling FPs
  \- Mismappings (as inferred by mappability track)
  \- Other possible recurrent artifacts [tissue specificity, strange contexts]

**Slides on mismappings:**

**PNPLA4** — most egregious example; alignment artifact hotspot pervades many tumor types.

> The number of supporting reads are less than 5 in ~400X

**RNF5** — another bad example

> In BRCA , the sequencing depth is limited within 50X, and only within 5 reads support the ALT while two has MAPQ=0. (maybe should not be called?)

**NBPF3**

>  The mutation at 6606 position is poorly aligned, and there are less than three reads supporting it (out of 300 reads) I do not think it should be called!

**PTPN11** — difference between MC3 callset and PanCan14k

> in LGG for the recurrent 8199 alt, the signals seem to be okay -  although there are mapping issue, there are about 50-60% reads supporting the existence of ALT.

**HSPA8** — no hotspots, but gene has generally poor mappability (significant by pCV): 

> Is there a technical writeup for how to get pCV? And the coverage is limited

**NONO**  

> NONO has a consensus region on chr16, but based on the number of supporting reads I think it is okay. 



**Slides on tissue specific hotspots:**



**MED12**: strong hotspot only in prostate; mapping looks fine. Could not be validated with Sanger sequencing:

https://www.ncbi.nlm.nih.gov/pubmed/23661306 (include screenshot of article title/abstract)

> The sequencing coverage is limited to 60X (12 supported)

**EPAS1**: strong hotspots only in PCPG and testicular. Testicular hotspot is in poorly mappable region.

> in testicular there are only 3-4 reads supporting the existence of the alt, the mappabiliy is actually okay here

**AP2B1**: poor mapping hotspot only in esophageal

> The number of supporting reads are limited and it is in a not well mapped region. Is there some possibilities that we first try to separate chromosomes and then do sequencing on each chromosome individually? I feel that this way we can eliminate most of the mapping issue...??

**WDR89**: strong hotspot only in mesothelioma

> There are only 3-4 supporting reads for the ALT (coverage ~90) in every samples. I assume it is a batch effect.

**PPCDC**: strong hotspot only in mesothelioma

> There are only 3-4 supporting reads for the ALT. Meanwhile for all three observed calls, there is a germline variant at chr15:75,336,729. Seems like a batch effect

**SOHLH2**: strong hotspot in odd context (A->G) only in stomach. MSI?

> How to determine MSI in IGV? It seems that the genome is relatively stable (there is not much germline/tumor variation) in the neighborhood.

**Thymoma seems to have a cohort-wide problem: many THYM-specific hotspots in poorly mappable regions**

**GTF2I**

> This is a very poorly aligned region (top 2-3 alignments with score 1000 on the same chromosome) Many reads blanked. Question: how do we study mutations within a homolog-rich region?

**SLC9B1**

> mapping quality issue, the number of supporting reads are limited (with in 7/300)

**SLC22A2**

> The coverage is low, alignment issue. Question - is there a way to color by F2R1 or F1R2  in IGV (are those handled by OXOG filter)?

**USP6**

> 10/250 reads support, but overall  in a low-mappable region.

**ZNF658**

> <10 reads, low mappability

**Possible artifacts:**

**IGLL5**  (lots of synonymous mutations [AID?], strange coverage pattern)

> In BLCA there seems to be a lot of SSNV (some with low phred score) but there seem to be a clear contrast between tumor and normal in terms of called ALT. I still cannot really  understand why there are no mutation in the latter portion of the gene - the chromosome location seems confusing to me - In tumor portal, the first dot in BLCA has genomic location 23230234 but the first red dot in SKCM has genomic location 23237863 - which is 7kb apart.

**LCP1**: hotspot is in atypical context, G(C->T)G. Could be misaligned microsatellite?

> The background mutation rate is high, but mappability is okay. supported by >30 reads. I saw lots of SNPs in this region, but I am not sure if there are microsatelite nearby?

------



**PCDHGA8** (faded)
This is a definite calling artifact.

> My question is why there seems no exsome in 140772815? The coverage is about 70X. Is it because the transcript we defined here is broader than exome?

**PSIP1** (faded, lack of pattern)

enrichment of truncating events

**ATAD5**  (long gene)
This is indeed a long gene, but there definitely seems to be an enrichment of truncating events. We can discuss a method to quantify this.

> Question: what is the distribution of truncating events under H0? 

**TRRAP** (most likely)
Seems to be one of the rare examples in which a driver event occurs in a highly mutable gene. Here's a study that specifically assessed the functional effect of the S722F hotspot: https://stm.sciencemag.org/content/3/83/83ec74

> There are some tymide dimers (double T in tumor)

- BRCA

**MAP2K4** 
Another example of a gene that, although it contains some problematic regions, has an enrichment for truncating events in non-problematic regions.

**SMTNL2**
This looks like a MutSig FP.  However, we should check to see why the first ~150 amino acids of this gene entirely lack coverage — it could simply be that whole exome capture kits do not target it, but it also might indicate other problems.

> The hotspot at 6467 calls are of poor base quality. I am not sure how to look for the first 150 amino acids - it seems to have coverage < 15. On the refseq track, there is no solid bloc indicating exon (does it mean exon?)

**BTG2** (an mutational cluster appear only in BTG2, which overlaps with some syn muts)
I agree, it's weird that this locus is only mutated in DLBCL. Maybe this is the result of AID target hypermutation?

![1564438659708](/home/qingzhang/.config/Typora/typora-user-images/1564438659708.png)

> Yes I can identify some AID signatures!

IGLL5 (too many syn muts, short, faded)
We should definitely survey these mutations in IGV — this looks like alignment artifacts. 

> This region has lots of SSNVs (single nucleotide in both tumor and normal)  - but mapQ okay in DLBC.  I am quite confused about the alignment in tumor-portal, on the same genomic location (same vertical position) we have 23237723 and 23230403 - does it mean there are two transcript defined 7kb apart?

- ESCA

  **AP2B1** (faded) P246T is probably an alignment artifact.

> Yes it has a homolog sequence at chrX and the fraction of ALT is about 2%, not very convincing!

GBM: many genes related to vesicle transport appears, not sure if there is a problem.

**PDGFRA** (sudden clear cut??? is that possible? I want to check)
Could indicate a sequencing artifact, or later exons in the gene might not be targeted for capture. Definitely something to check in IGV.

> I suspect it is the problem of how we define transcripts - The alignments agree with refseq genes but the display on tumor portal is truncated.

**SEMG1** (faded, lack pattern)
Might check IGV.

> most of the mutations seem to be okay

> COL1A2 (collagen, not sure) 
> This one is borderline; there does seem to be a slight enrichment of truncating events

HNSC

**NCOR1** (long, many syn, lack pattern) 
That's true, but there does seem to be an overall pan-cancer enrichment for truncating events. I don't really see that in head+neck.

> I really want to know how do you mathematically identify "pan-cancer enrichment for truncating events".

**LCP1** (plastin, border)
Significance is probably due to R488H hotspot; might check this in IGV.

> There are lots of sporatic passgener events near hotspot, overall seems okay

**ADCY1** (many syn) [coverage]
Also first part of gene lacks coverage.

**CCDC50** (border)
Probably driven by R19G hotspot; might check in IGV. Also recurrent events around AAs 102/103.

> seems okay now. The R19G hotspot has really low coverage.

LUAD
**COL11A1** (collagen, long, syn)
Looks enriched for truncating events.

MYH7 (myosin, long syn faded)
Worth checking in IGV

> There seems some mapping issue (lots of blank reads aligned) - and the fraction of mutation is not that much

> HSPA8 (cancer impact protein folding?? faded) 
> Definitely looks artifactual. Worth checking in IGV

**CPED1** 
Check if cluster of mutations around AA 940-950 are artifactual in IGV

> There are two adjacent mutation occurring simultaneously - 2 adjacent G->T mutation(fraction=45%). In other cases the mutation fraction is less than 10%

**AACS** (border)
Two G645 mutations occur in the same patient (44-8120). Worth looking at.

> ![1564447911201](/home/qingzhang/.config/Typora/typora-user-images/1564447911201.png)

> 44-8120 patient has two adjacent G->C mutation (occur on the same read). Note that there is one homologue on chromosome 5



LUSC

**CPED1** (border)
See previous comment on this gene.

> coverage is slightly low, but overall seems okay with the hotspot

**MYH7** (myosin, long, faded)
Worth looking at in IGV

> The mutations are usually of very little allelic fraction (less than 7%) . Myosin has some repeating domains so there could be some mappbility issue (first alignment scoring 1000 while the second scoring 750 at 400b apart)

**PTPN11** (faded)
Worth looking at in IGV

> contains region with low mappability

- MESO
  **WDR89** (rank 1st but only reported in MESO & not hallmark - I want to check) 
  This is definitely worth looking at. Very strong hotspot only in mesothelioma.

> although there is an aligned hotspot, the area is of low coverage(~70) and the number of supporting reads is usually less than 5. It seems weired that it presents in almost all samples.

**PPCDC** (sparse and syn, very conserved)
R82H hotspot is worth looking at in IGV — only really occurs in MESO (one event in UCEC)

> same. The mutant allele fraction is very limited but the coverage is okay(~200).

**PAAD** 

**NBPF3** (faded, is reported to expand in primates - is it a bias from uniform PON?) 
That F424C hotspot is almost certainly an artifact. Check in IGV.

> Yes this region is really poorly aligned

**PCPG**
**EPAS1** (border, syn)
Large hotspot at AA ~530. Might be a real PCPG-specific driver? https://www.ncbi.nlm.nih.gov/pubmed/23418310
Worth checking in IGV anyway.

> Yes there is a NTAD domain near this hotspot that might be related to PCPG-specific mutagenesis.https://cancer.sanger.ac.uk/cosmic3d/protein/EPAS1?pdb=5KIZ

**PRAD**
**MED12** (long sporatic) 
L1224F hotspot might be an artifact: https://www.ncbi.nlm.nih.gov/pubmed/23661306
Definitely check in IGV

**SARC**
**PNPLA4** (aligned faded spots)
Hotspot is almost definitely an alignment artifact. Check in IGV.

> The number of supporting reads are very limited (<2%) and sometimes are also observed in normal. Mapping is okay here (the second best alignment only has score of 200+ while the dorminant one has score of 1000)



> SKCM: contains many border genes, hard to draw the line.
> SLC9C2 (border)
> Hard to say if those missense events are significantly clustered in certain protein domains.

STAD

**BNC2** (not enough pattern, skin pigment)
Check S575R hotspot in IGV.

> The mapping seems okay but it is in a low coverage region.

**RIMKLB** (low quality calls, no clustering)
Check for artifacts in IGV.

> There is a similar sequence in chr 21

**CNBD1**
Check L396 hotspot in IGV.

> The coverage in this region is very low (~40) but seems no mapping issue - the portion of reads supporting ALT is fairly high

**BTBD1**
Could be a real driver; hotspot is borderline significant. Check IGV.

> seems okay but the number of supporting reads quite low, but there are significantly more mutations in callset than in tumor portal - maybe different set of transcript definition is used.

**MYCT1**
Significance probably driven by hotspot.

> Yes this one seems okay

**SOHLH2** (testis-specific TF, questionable)
Hotspot mutation is in interesting context; check in IGV

**ARMC4** (expressed in respiratory cilia, mutation causes imorbility, questionable)
Significance probably driven by F305; check in IGV (since region adjacent to it is faded)

> This gene uses the negative strand and there is similar sequence on the positive strand.

**TGCT**
**SPIN2A** (faded overall)
Check in IGV

> In a very poorly aligned region

**EPAS1** (border, syn, sporatic)
Check S478C mutation in IGV

**THYM** 

**GTF2I** (rank 1st but all faded dots aligned) 
Check hotspot in IGV

> low mappability, 25% supporting reads

**SLC9B1** (faded overall)
Check hotspot in IGV

> I do not think they are real. There are only ~10 reads in ~400, and many of the reads are of low MAPQ.

**SLC22A2** (many syn muts)
Check hotspot in IGV

> not much supporting reads (4 out of 100)

**USP6** (faded overall)
I think this cohort might have an issue … check hotspot in IGV

> the number of supporting reads is less than 15 (out of 300), and it is on a poorly mapped region.

**ZNF658** (faded overall) 
Check IGV

> Homologs prevails entire gene, many of the hotspots are not well supported

