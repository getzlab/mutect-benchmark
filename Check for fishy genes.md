# Check for fishy genes

This human-intelligence training tasks - loop through each significant genes in the call-set and check for some unusual behavior. The major "dubious" aspects are -

* lie in a poorly-aligned region
  * examples, with in-page references
* artifact in sequencing: the limited calls are from "hypermutation patients"
  * examples, 
* Not enough support

All the potential fishy gene-cancer pairs are flagged for further checking with IGV. potential erroneous patient-wise mutation burden are also marked for further check.

marks in  the document

* **#fishy** potentially problematic
* **#interesting** great pattern, uniqueness
* **#broader** broader cases, the support is not as strong
* **#error** gene symbol cannot be found
* **#specific** cancer-type specific
* **#question** simple thought on something ;) 

## ACC: Adrenocortical carcinoma

![comut-plot](/home/qingzhang/Downloads/comut_pack/comut_plots_pngs/ACC_coMut.png)

- MEN1: overwhelming amount of nonsense mutations (and indels). (alternate chromosome structures and epigenetic gene regulation, encodes fro menin, a potential TSG)
- C12orf65: This gene is very short (120bp), and the mutations seems to be scattered. **#broader** This gene encodes an mitochrondrial matrix protein that contributes to mitochondrial translation machinery, translation release factor activity.

Three patients have simliar albeit different mutation spectrum from others.

## BLCA: Bladder Urothelial Carcinoma

![plot](/home/qingzhang/Downloads/comut_pack/comut_plots_pngs/BLCA_coMut.png)

#### Patients

**#fishy-patient** check for the one with hypermutation.

#### Spectrum

![1561310115659](/home/qingzhang/Downloads/comut_pack/comut_plots_pngs/1561310115659.png)

The LEGO plot indicates a prevalence of C>T and C>G mutation, which agrees with our spectrum.

#### Genes

* KDM6A: hallmark gene, modulates H3K27 acetylation and performs H3K27 demethylation, TSG+oncogene, leukemia-maintenance
  * prevalence of truncation indels Q555* compared to missense! #surprise
* PIK3CA: hallmark, ras pathway
* FBXW7: encodes a F-box protein that constitutes ubiquitin ligase complex, P-dependent ubiquitination
* RHOB: **#interesting** there are a lot of missense in BLCA but not in any other cancers (E172K / P75T <proline>) **#specific** 
* ARID1A: long gene, hallmark (when depleted cannot go through cell-cycle arrest). 

* RB1: 900bp, negative regulator in cell cycle, mostly frameshift and non-sense.
* CREBBP: 2250bp **#broader** transcriptional coactivation
* KMT2C **#error** cannot find the symbol. The previous name MLL3 cannot be displayed (perhaps nothing?) and there are multiple forms of mutations spanning over many individuals **#fishy**
* KMT2B **#error** cannot be found - The previous name MLL2 cannot be displayed (perhaps nothing found?)
* KMT2D **#error** same
* STAG2: 1200, nonsense Q593*, cohesion complex at centromeres, mediates cell-type specific contacts between enhancers and promoters.
* ELF3: evenly scattered mutation, short (350) **#broader**
* ERCC2: **#specific** short gene but very concentrated mutation spots - Y72C and N238S (**#interesting**, because N/S are chemically similar). In transcription-coupled nucleotide excision repair pathway.
* CDKN2A: interact and stabilize P53, cell cycle G1 control
  * ![Simplified diagram of how p53 halts the cell cycle at the G1/S checkpoint. p53 is activated by DNA damage and causes production of a Cdk inhibitor, which binds to the Cdk-G1/S cyclin complex and inactivates it. This halts the cell in G1 and prevents it from entering S phase, allowing time for the DNA damage to be fixed.](/home/qingzhang/Downloads/comut_pack/comut_plots_pngs/portal_plots/G1-checkpoint.png)
* GNA13 **#specific** R200G, Gprotein
* CTNNB1**#overlap** gene product constitute adheren junctions,  which create and maintain epithelial cell layers by regulating cell growth and adhesion between cells.
* ERBB3: mostly missense. Encodes a EGFR family tyrosine kinase, it can bind its ligand but not convey signal to cell through protein phosphorylation (no kinase domain).
* ERBB2: S310Y S310F missense. This gene encodes a member of the epidermal growth factor (EGF) receptor family of receptor tyrosine kinases. This protein has no ligand binding domain of its own and therefore cannot bind growth factors. However, it does bind tightly to other ligand-bound EGF receptor family members to form a heterodimer, stabilizing ligand binding and enhancing kinase-mediated activation of downstream signalling pathways, such as those involving mitogen-activated protein kinase and phosphatidylinositol-3 kinase. (ERBB2 does not have ligand-binding domain while ERBB3 does not contain kinase domain)
* OGDH: **#fishy** no clustering signal, 1kb, function mostly with energy metabolism, not a hallmark gene
* PTPN4: protein tyrosine phosphatase, all missense, no clustering signal, potentially **#fishy** <online note:This RefSeq record was created from transcript and genomic sequence data to make the sequence consistent with the reference genome assembly. The genomic coordinates used for the transcript record were based on transcript alignments.>
* HIST1H3B: **#error** basic nuclear proteins that are responsible for nucleosome structure, no MC3 records found **#fishy**
* TMCO4: **#overlap** Q13* truncation. 
* TFPI2: clustered in R222C/R222H but not a hallmark gene. sparse elsewhere **#interesting**
* RBM26: no clustering , 1kb, **#fishy**
* RIPK4: **#error** nothing displayed in MC3, no signal from other callsets
* ZBTB7B: codes for a zinc-finger containing transcription fatcor that acts as a key regulator of lineage commitment of immature T-cell precursors **#broader**
* PCDHGA8: **#fishy** sequencing artifacts (not saturated), no clustering
* SF1: splicing factor 1: recognizes intron branch point sequence and is required for the early stages of spliceosome assembly
* SF3B1: **#spefic** E902K in bladder (charge changes!) **#interesting** also splicing factors
* PSIP1: **#fishy** sequencing artifacts, no clustering. hallmark: embryonic development, lysosomal statibility
* ATAD5: scattered **#fishy** no hotspots, ATPase
* AHR: aryl hydrocarbon receptor  - a ligand-activated tf to respond to aromatic hydrocarbons. Q383H presents in bladder / kidney / liver and lung. **#interesting**
* TRRAP: 3.6kb, potentially **#fishy** , many missense
  * ![1561319212821](/home/qingzhang/Downloads/comut_pack/comut_plots_pngs/portal_plots/TRRAP.png)

* DIAPH2: no hotspots, perhaps over-representation, potentially **#broader, #fishy** 1kb
  * ![1561319440421](/home/qingzhang/Downloads/comut_pack/comut_plots_pngs/portal_plots/DIAPH2.png)
* CUL1: E485K, ubiqitin protein ligase binding, activated TLR signaling, quite promising..
* RXRA: **#specific** S427Y S427F retinoic acid-mediated gene activation. ~~retinoic acid related to bladder??~~
* EP300: p300 transctiptional coactivator histone acetyltransferase via chromatin remodeling and important for proliferation and differentiation. 2.25kb Q182**#specific** truncation, bladder-specific mutation hotspots in 1100-1500bp

## BRCA: Breast cancer

![brca](/home/qingzhang/Downloads/comut_pack/comut_plots_pngs/BRCA_coMut.png)

:warning: sequencing artifact in at least one individual!

:stop_sign: Individual ID: to be fillled

* PIK3CA: Classic! ras effector, catalytic subunit of kinase, many activating mutations :paperclip:worth explore! Most missense mutation but there are inframe deletion and splice-site as well.
  * ![img](/home/qingzhang/Downloads/comut_pack/comut_plots_pngs/portal_plots/rtk-pathway.jpg)
* TP53
* AKT1: serine protein kinase. catalytically inactivate in serum-starved primary and immortalized fibroblasts. AKT1 and AKT2 are activated by platelet-derived growth factor. Intense hotspot in E17K.
* CDH1(Cadherin1): **#specific** intense clustering of frameshift indel / nonsense in BRCA, but not other tumors! **#interesting** Q23*. cadherin, Ca-dependent cell cell adhesion glycoprotein. Increasing proliferatoin, invasion and metastasis.
* CBFB: forms heterodimer with RUNX transcription factor to enhance their affinity for DNA, negatively regulated ribosomal gene expression during mitosis. Required for transcription of a subset of Runx2-target genes that are sufficient to maintain the invasive phenotype of the cells in breast cancer. (TSG/fusion) There are no specific mutation hotspot, but it seems to be an over-representation of CBFB mutation in BRCA compared to other cancers.
* PTEN: seems the sequencing quality is not very good. classic TSG. **#interesting Q:** ~~why quality so low? special conformation??~~
  * ![The Akt/PTEN pathway. Oncogenic and mitogenic stimuli that activate PI3kinase can lead to Akt activation, either directly, via the actions of phosphatidylinositol 3,4,5-trisphosphate (/home/qingzhang/Downloads/comut_pack/comut_plots_pngs/portal_plots/pten-pathway.jpg ) and phosphati-Â ](https://www.researchgate.net/profile/Iain_Mcneish2/publication/8887434/figure/fig2/AS:277685025689615@1443216677879/The-Akt-PTEN-pathway-Oncogenic-and-mitogenic-stimuli-that-activate-PI3kinase-can-lead-to.png)
* MAP3K1: mostly indel and truncating events in breast. The first 150 bp are error-free ~~could be vitally important??~~! **#broader** but have literature support in breast.... hmmm. It is a mitogen activated protein kinase kinase kinase that activate JNK signalling pathway, stabilizes Myc by promoting N-terminal Phos and enhancing its transcriptional activity (Oncogene + TSG)
  * ![1561324327829](/home/qingzhang/Downloads/comut_pack/comut_plots_pngs/portal_plots/MAP3K1.png)
* MAP2K4: potentially **#fishy** sequencing quality not good and it is not a hallmark gene. mutations are not clustered.
  * ![1561324530318](/home/qingzhang/Downloads/comut_pack/comut_plots_pngs/portal_plots/MAP2K4.png)
  * ![MAP2K4 with MAP3K1](/home/qingzhang/.config/Typora/typora-user-images/1561325119350.png)
* CASP8: encodes cysteine-aspartic acid protease (caspase) family, sequential activation of caspase plays central role in execution-phase of cell apoptosis. Q482H mutation of procaspase-8 abolishes caspase-8-mediated apoptosis by impairing procaspase-8 dimerization in AML.
* FOXA1: oncogene hallmark. The major signal comes from prostate cancer cluster, it is very dense at 250bp! siliencing inhibits EMT and suppresses proliferation. (basically most findings come from gene perturbation? ~~we are bayesian~~) 
* KMT2C **#error**
* KRAS: very clean plot! G12R/C, G13D, E62K, A146T spanning multiple cancers
* SF3B1: mostly clusted in 600-720 region, but not a hallmark gene
* ERBB2: as mentioned in BLCA
* RNF5: potentially **#fishy** sequencing calls are not of good quality, not a hallmark gene
  * ![1561326718061](/home/qingzhang/Downloads/comut_pack/comut_plots_pngs/portal_plots/RNF5.png)
* RB1: classic TSG, chromatin-associated protein that limits the transcription of cell cycle genes, primarily via regulation of the E2F transcription factor
* ARID1A: TSG/fusion, classical, frequent inactivating mutations, fuse with MAST2 in triple-negative breast cancer. Maintains genome statibility, promote differentiation.
* SMTNL2: potentially **#fishy** not enough support
  * ![1561327055800](/home/qingzhang/.config/Typora/typora-user-images/1561327055800.png)

## CESC: Cervical squamous cell carcinoma and endocervical adenocarcinoma

![cesc](/home/qingzhang/Downloads/comut_pack/comut_plots_pngs/CESC_coMut.png)

There are much more variation in mutation spectrum. :warning: check the first patient for sequencing issue!

* PIK3CA
* FBXW7
* EP300: coactivator lf HIF1A and thus plays a role in hypoxia-induced genes such as VEGF (~~cervial cancer is related to hypoxia??~~)
* NFE2L2: nuclear factor like 2, encodes a TF which is a member of basic leucine zipper proteins, it regulates genes which contain antioxidant response elements (ARE) in their promoters. Very sharp peak! (it worths to check cervical cancer and oxygene related responses!)
  * ![1561327820111](/home/qingzhang/.config/Typora/typora-user-images/1561327820111.png)
* KMT2C **#error**
* KRAS
* KMT2D **#error**
* ERBB2
* TP53
* ARHGAP4: Rho GTPase activating protein 4, regulating proteins from Ras family. But the signals are weak - **#broader**
  * ![1561328018499](/home/qingzhang/.config/Typora/typora-user-images/1561328018499.png)

## CHOL: Cholangiocarcinoma

![chol](/home/qingzhang/Downloads/comut_pack/comut_plots_pngs/CHOL_coMut.png)

It seems the mutational burden for each CHOL patients are not as much, could be worth to check the patient with 12 muts(maybe not). The mutation spectrum are consistent.

* IDH1: isocitrate dehydrogenase 1 (NADP+), soluble, significant R132G/H, especially in brain lower grade glioma. Hallmarkgene: cytosolic krebs cycyle enzyme, catalyzing the formation of alpha-ketoglutarate from isocitrate. R132H alters histone marks and induces extensive DNA hypermethylation. Missense R132H causes a shift in the active site resulting in decreased binding affinity for NADPH. ​ :question:  ~~is CHOL pathologically related to energy metabolism?~~
* BAP1: BRCA1 associated protein-1 (ubiquitin carboxy-terminal hydrolase) important in the G1/S transition, regulates the levels of cell cycle proteins by affecting E2F1-dependent S-phase gene expression
* PBRM1: **#interesting** mostly nonsense. Similar scattered nonsense mutation observed in Kidney renal clear cell carcinoma. TSG chromatin remodeling complex. Mutation occurs mostly in bromodomain and BAH domain.

## COAD: Colon adenocarcinoma

![COAD](/home/qingzhang/Downloads/comut_pack/comut_plots_pngs/COAD_coMut.png)

:warning:  There are lot of patients whose mutation burden is very high >150 mutations/mb, they might be sequencing artifacts! The "hyper-mutated" patients have elevated number of C->A in T_T and A->G

* APC: intense signal from clustering of nonsense and inframe indels in colon adenocarcinoma (perhaps suggesting similar etiology). repressor of Wnt signaling pathway part of beta-catenin destruction complex, promotes cell migration, promotes apoptosis, classic hallmark gene, mutation mostly in unstructured protein sequence. 
  * ![1561329323273](/home/qingzhang/.config/Typora/typora-user-images/1561329323273.png)
  * ![Addiction to the Wnt signalling pathway through gene-silencing events.a | In normal colon epithelial cells, secreted frizzled-related proteins (SFRPs) function as antagonists of Wnt signalling by competing with Wnt proteins for binding to their receptor, Frizzled (FRZ). Expression of SFRPs is therefore the epigenetic gatekeeper step. When Wnt signalling is inactive, the adenomatosis polyposis coli (APC) complex phosphorylates -catenin, leading to its degradation. This prevents the accumulation of nuclear -catenin and therefore its ability to engage its transcription factor partners (TGF), which results in the differentiation and homeostasis of colon epithelial cells. Expression of APC is therefore a genetic gatekeeper step. b | When SFRP expression is lost, through epigenetic silencing of the gene that encodes it (loss of the epigenetic gatekeeper), Wnt signalling becomes activated through the receptor FRZ. This Wnt signalling potentially inactivates the APC complex (loss of the genetic gatekeeper), allowing -catenin to accumulate in the cytoplasm and eventually in the nucleus. In the nucleus, -catenin activates transcription of genes such as MYC, cyclin D and other genes whose products promote cell proliferation and survival rather than differentiation. This results in the expansion of colon epithelial stem and progenitor cells and formation of atypical crypt foci (ACF). c | Persistent activation of the Wnt pathway allows mutations to occur in other pathway components, such as those that permanently disable the APC complex and promote nuclear accumulation of -catenin (loss of the genetic gatekeeper, as indicated by the bold cross). These cells are selected for because of their survival and proliferative advantages. This combination of epigentic and genetic events fully activates the Wnt pathway to promote tumour progression. Without the epigenetic events that silence the SFRP genes, mutations that disrupt the APC complex might not be sufficient to promote tumorigenesis or tumour progression.](https://www.researchgate.net/profile/Joyce_Ohm/publication/7286336/figure/fig4/AS:271656811757600@1441779439889/Addiction-to-the-Wnt-signalling-pathway-through-gene-silencing-eventsa-In-normal-colon.png)
* TP53
* KRAS
* PIK3CA
* BRAF
* FBXW7
* SMAD4
* TCF7L2: transcription factor 7-like 2 (T-cell specific, HMG-box), Wnt signaling pathway, blood glucose homeostasis
* PCBP1: **#interesting #specific** L100Q for colon cancer, KH domain, mrna splicing.
* NRAS: signature G12D G13R Q61H/R/L/K
* SMAD2:  This protein mediates the signal of the transforming growth factor (TGF)-beta, and thus regulates multiple cellular processes, such as cell proliferation, apoptosis, and differentiation.**#interesting** it has a consensus truncation at the end of transcript, (what is that for???)
  * ![1561346526779](/home/qingzhang/.config/Typora/typora-user-images/1561346526779.png)
* AMER1: clustered truncation near 500 bp. R358*. binds to the tumour suppressor APC and acts as an inhibitor of Wnt signalling by inducing beta-catenin degradation [[Pubmed\]![img](https://cancer.sanger.ac.uk/core/gfx/blank.gif)](https://www.ncbi.nlm.nih.gov/pubmed/21498506); inhibits the ubiquitination of NRF2; AMER1 and NRF2 compete for binding to KEAP1 [[Pubmed\]![img](https://cancer.sanger.ac.uk/core/gfx/blank.gif)](https://www.ncbi.nlm.nih.gov/pubmed/22215675); enhances p53 acetylation by CBP/p300 [[Pubmed](https://www.ncbi.nlm.nih.gov/pubmed/22285752)]. hallmark TSG.
* PDYN: prodynorphin, The protein encoded by this gene is a preproprotein that is proteolytically processed to form the secreted opioid peptides beta-neoendorphin, dynorphin, leu-enkephalin, rimorphin, and leumorphin. These peptides are ligands for the kappa-type of opioid receptor. Dynorphin is involved in modulating responses to several psychoactive substances, including cocaine. Mutations tends to happen at the end of protein. (~~psychoactive substance with colon cancer???~~) **#broader** 
* DAB2: **#specific** V172E disabled homolog 2, mitogen-responsive phosphoprotein, but it is not a classical hallmark
* B2M: beta-2-microglobulin **#question** (why the frameshift indel is not flagged?) hallmark: recurrent missense mutations of the 1st amino acid, recurrent L15fs, inactivating mutations in 5' region of the gene lead to impaired expression of HLA class I antigen, thus helping tumour to avoid the immune destruction in melanoma FO-1 cell line.
* MAP2K4
* CCDC30: **#broader** V396E, not classic hallmark
* CHST15: glycosaminoglycan, This gene has also been identified as being co-expressed with RAG1 in B-cells and as potentially acting as a B-cell surface signaling receptor. **#broader** not well-supported, no obvious clustering, potentially **#fishy**
  * ![1561345520545](/home/qingzhang/.config/Typora/typora-user-images/1561345520545.png)
* CCT6B: **#fishy** molecular chaperone to achieve optimal fold, not functionally-likely???? **#broader**
  * ![1561345735572](/home/qingzhang/.config/Typora/typora-user-images/1561345735572.png)
* PTEN
* ATM: classical hallmark. phosphorylates histone H2AX in response to DNA double-strand breaks [[Pubmed\]![img](https://cancer.sanger.ac.uk/core/gfx/blank.gif)](https://www.ncbi.nlm.nih.gov/pubmed/11571274); plays a central role in the repair of DNA double-strand breaks, recruits and cooperates with other sensor proteins. Although  overall mutational burden is not very high (nor is it shows sign of clustering) we can believe the signals.
* SNX13: **#specific** L663F Intracellular trafficking, G-protein signaling. **#broader**
* FAT4: FAT tumor suppressor homolog 4 (Drosophila) Hallmark TSG, protein belonging to cadherin (calcium-dependent cell adhesion protein) family [[Pubmed\]![img](https://cancer.sanger.ac.uk/core/gfx/blank.gif)](https://www.ncbi.nlm.nih.gov/pubmed/23076869); suppresses phosphorylation and nuclear accumulation of Yap, a protein associated to the promoted proliferation, migration and cell cycle progression in gastric cancer.However the gene is very long 4.5kb and the mutation does not show any clustering. Therefore do not know much about colon-specific implications.
  * ![1561346244141](/home/qingzhang/.config/Typora/typora-user-images/1561346244141.png)
* ASPM: Human ortholog of Drosophila abnormal spindle gene. mitotic spindle regulation, with a preferential role in regulating neurogenesis. Very long (3.5kb) high burden of syn mutations. **#fishy**
  * ![1561346426785](/home/qingzhang/.config/Typora/typora-user-images/1561346426785.png)
* CYHR1: not enough information, no clustering mutation hotspots
  * ![1561346495312](/home/qingzhang/.config/Typora/typora-user-images/1561346495312.png)

## DLBC: Lymphoid Neoplasm Diffuse Large B-cell Lymphoma

![dlbc](/home/qingzhang/Downloads/comut_pack/comut_plots_pngs/DLBC_coMut.png)

:warning: Per patient mutation burden relatively stable, quite some A->G and C->A mutations

The following genes have a expected tendency for immune function (features CDxxx), however will it be a result of higher power from overly expressed genes? **#question** The mutational plot for many genes here is sparse, but given the fact that they are mostly short protein - it could be expected.

* B2M
* MYD88: This gene encodes a cytosolic adapter protein that plays a central role in the innate and adaptive immune response. This protein functions as an essential signal transducer in the interleukin-1 and Toll-like receptor signaling pathways. These pathways regulate that activation of numerous proinflammatory genes. **#question** why cosmic says the post prevalent substitution is L265* but we have most L252P (why is it labeled non/mis)??
* TP53
* BTG2: **#specific but #fishy**. However it is structually related to proteins that appear to have antiproliferative properties, involves in G1/S transition of cell cycle. **#question** why is this cluster only appear in DLBC, if the mechanism is so fundamental???
  * ![1561346996388](/home/qingzhang/.config/Typora/typora-user-images/1561346996388.png)
* TMSB4X: really short!!! an actin sequestering protein which plays a role in regulation of actin polymerization. The protein is also involved in cell proliferation, migration, and differentiation. This gene escapes X inactivation and has a homolog on chromosome Y. But viewing from the plot, we do not have enough confidence due to lack of support.**#fishy**
  * ![1561347197852](/home/qingzhang/.config/Typora/typora-user-images/1561347197852.png)
* IGLL5: many syn mutation, potentially **#fishy** Not much support from literature
  * ![1561347297861](/home/qingzhang/.config/Typora/typora-user-images/1561347297861.png)
* KLHL6: This gene encodes a member of the kelch-like (KLHL) family of proteins, which is involved in B-lymphocyte antigen receptor signaling and germinal-center B-cell maturation. The encoded protein contains an N-terminal broad-complex, tramtrack and bric a brac (BTB) domain that facilitates protein binding and dimerization, a BTB and C-terminal kelch (BACK) domain, and six C-terminal kelch repeat domains. Naturally occurring mutations in this gene are associated with chronic lymphocytic leukemia. **#broader** there is no clustering and the background syn rate is high
  * ![1561347744236](/home/qingzhang/.config/Typora/typora-user-images/1561347744236.png)
* UBE2A: ubiquitin pathway - The modification of proteins with ubiquitin is an important cellular mechanism for targeting abnormal or short-lived proteins for degradation. Ubiquitination involves at least three classes of enzymes: ubiquitin-activating enzymes, or E1s, ubiquitin-conjugating enzymes, or E2s, and ubiquitin-protein ligases, or E3s. This gene encodes a member of the E2 ubiquitin-conjugating enzyme family. This enzyme is required for post-replicative DNA damage repair. However the mutation is very sparse and lacking clustering - it worth to note that we do not have syn here - only splice and nonsense. **#broader**
  * ![1561347902537](/home/qingzhang/.config/Typora/typora-user-images/1561347902537.png)
* CD70: The protein encoded by this gene is a cytokine that belongs to the tumor necrosis factor (TNF) ligand family. This cytokine is a ligand for TNFRSF27/CD27. It is a surface antigen on activated, but not on resting, T and B lymphocytes. It induces proliferation of costimulated T cells, enhances the generation of cytolytic T cells, and contributes to T cell activation. This cytokine is also reported to play a role in regulating B-cell activation, cytotoxic function of natural killer cells, and immunoglobulin sythesis. **#broader** mutations are so sparse
  * ![1561348160032](/home/qingzhang/.config/Typora/typora-user-images/1561348160032.png)
* CD79B: **#specific** A68N for DLBC (A196N) **#question** different gene model? Ig beta, activator of Syk kinase [[Pubmed\]![img](https://cancer.sanger.ac.uk/core/gfx/blank.gif)](https://www.ncbi.nlm.nih.gov/pubmed/15219998); truncating mutation in Ig beta prevents the assembly of the IgM BCR on the cell surface, resulting in the absence of peripheral B cells and low/absent immunoglobulin serum levels [[Pubmed\]![img](https://cancer.sanger.ac.uk/core/gfx/blank.gif)](https://www.ncbi.nlm.nih.gov/pubmed/17709424)
* PIM1: The protein encoded by this gene belongs to the Ser/Thr protein kinase family, and PIM subfamily. This gene is expressed primarily in B-lymphoid and myeloid cell lines, and is overexpressed in hematopoietic malignancies and in prostate cancer. It plays a role in signal transduction in blood cells, contributing to both cell proliferation and survival, and thus provides a selective advantage in tumorigenesis. Both the human and orthologous mouse genes have been reported to encode two isoforms (with preferential cellular localization) resulting from the use of alternative in-frame translation initiation codons, the upstream non-AUG (CUG) and downstream AUG codons. 
* FAS: TNF receptor superfamily. The interaction of this receptor with its ligand allows the formation of a death-inducing signaling complex that includes Fas-associated death domain protein (FADD), caspase 8, and caspase 10. The autoproteolytic processing of the caspases in the complex triggers a downstream caspase cascade, and leads to apoptosis. This receptor has been also shown to activate NF-kappaB, MAPK3/ERK1, and MAPK8/JNK, and is found to be involved in transducing the proliferating signals in normal diploid fibroblast and T cells. Several alternatively spliced transcript variants have been described, some of which are candidates for nonsense-mediated mRNA decay (NMD). The isoforms lacking the transmembrane domain may negatively regulate the apoptosis mediated by the full length isoform.
* HLA-B **#error** nothing displayed
* CARD11: molecular scaffolds, This protein is also a member of the CARD protein family, which is defined by carrying a characteristic caspase-associated recruitment domain (CARD). This protein has a domain structure similar to that of CARD14 protein. The CARD domains of both proteins have been shown to specifically interact with BCL10, a protein known to function as a positive regulator of cell apoptosis and NF-kappaB activation. When expressed in cells, this protein activated NF-kappaB and induced the phosphorylation of BCL10. **#broader**
* KMT2D **#error**

## ESCA: Esophageal carcinoma



![esca](/home/qingzhang/Downloads/comut_pack/comut_plots_pngs/ESCA_coMut.png)

Two patients are significantly hypermutated. 

* TP53
* NFE2L2: Mostly missense. This gene encodes a transcription factor which regulates genes which contain antioxidant response elements (ARE) in their promoters; many of these genes encode proteins involved in response to injury and inflammation which includes the production of free radicals. (**#question** is there some mutational signatures associated with free radicals?)
* CDKN2A: inhibitor of CDK4 kinase. cell cycle G1 control. 
* SMAD4
* PIK3CA
* AP2B1: protein product links clathrin to receptors in coated vesicles potentially **#fishy**  - could be sequencing artifacts
  * ![1561349321459](/home/qingzhang/.config/Typora/typora-user-images/1561349321459.png)

## GBM

* TP53
* PTEN
* EGFR
* PIK3CA
* RB1
* NF1
* IDH1
* PIK3R1
* BRAF
* PDGFRA
* MAP3K1
* SEMA3C
* SEMG1
* SLC26A3
* EXOC2
* LZTR1
* COL1A2
* KRT15

## HNSC

* TP53
* PIK3CA
* CDKN2A
* FAT1
* NOTCH1
* CASP8
* NSD1
* EP300
* HRAS
* TGFBR2
* KMT2D
* HLA-B
* RAC1
* FBXW7
* RHOA
* ZNF750
* RB1
* KEAP1
* EPHA2
* FGFR3
* MAPK1
* MB21D2
* NFE2L2
* ZNF623
* DOK6
* PTEN
* EEF1A1
* NCOR1
* LCP1
* CCDC50
* ADCY1
* GIGYF2
* CUL3

## KICH

TP53

## KIRC

* MTOR

* BAP1
* VHL
* PBRM1
* STED2
* SETD2
* TP53 (why so less important)
* PIK3CA
* PTEN
* ATM

## KIRP

* MET
* KRAS

## LGG

* IDH1
* TP53
* ATRX
* CIC
* PIK3CA
* EGFR
* IDH2
* PTEN
* FUBP1
* NOTCH1
* ZBTB20
* PIK3R1
* SMARCA4
* NF1
* NRAS
* PTPN11
* DNMT3A

# General Questions

## chromosome translocations

mutation happens on an unstructed protein segments or a well-characterized one? bent is more vulnerable?

if truncating events and frameshift prevails, what is it implies? very obvious selection advantage?

# General pathways

### RAS signaling

