
# Population genetic structure and clonal analyses of _Racomitrium lanuginosum_ in the forest-tundra transition zone.

## 1. Overview
This repository contains the scripts and dataset related to the project “Hide-and-seek in the tundra: Accounting for cryptic species uncovers molecular variation in moss clonal populations”. The main goal of the code is to provide detailed information on the project/article analyses (reproducibility).

> **Abstract**: Plant genetic variation can sometimes be structured in complex patterns where several processes intervene, like cryptic speciation, masking the real diversity across habitats. When not taken into account, cryptic species can over-inflate the genetic diversity estimates by the confluence of two or more distinct likely similar lineages. Such a situation is exacerbated in clonal lineages spread and intermingled in mats, and as a consequence, a multilevel pattern of interspecific and intraspecific genetic variation arises. This research delves into the population genetics and clonal structure of the moss Racomitrium lanuginosum in the Canadian tundra to prove the importance of decoupling genetic variation for an accurate biodiversity assessment in cold environments. We performed genotyping-by-sequencing to produce genomic data. Cryptic molecular lineages within populations were identified using haplotype-based coancestry matrices and phylogenetic analyses of a combined dataset that included samples of the north distribution range. We tested the accuracy of genetic diversity estimates with different clustering scenarios: one species or two cryptic lineages. The clonal structure was characterized using multilocus genotypes at population and fine scales. The correct characterization of cryptic lineages in R. lanuginosum revealed genetic diversity differences between habitats and proved the inflation of genetic estimates when treating differentiated cryptic lineages as a single species. Clonal growth influences the genetic structure of moss populations at the local scale with no apparent sexual reproduction. This plant exemplifies the multilevel genetic variation of cryptic species in the Northern Hemisphere. We highlight the value of determining genetic partition in complex biological entities for biodiversity assessments.

## 2. Software information
The code was run on the bioinformatic platform of the [Institut de Biologie Intégrative et des Systèmes (IBIS)](https://www.ibis.ulaval.ca/en/services-2/bioinformatics/documentation-servers/) at Laval University and on the statistic environment R. Check each script for details.
The following software and R packages were employed for the analyses:

* Software
  - fastqc/0.11.8
  - Trimmomatic-0.36
  - gcc/6.2.0
  - stacks/2.5
  - bwa/0.7.17 
  - samtools/1.8
  - BBTools/36.92
  - python/2.7 
  - fineRADstructure/1.0
  - RAxML/8.2.9
  - vcftools/0.1.16
 
* R/v.4.1.3 packages
    - poppr/v.2.9.3 
    - adegenet/v.2.0.2
    - vcfR/v.1.12.0 
    - ade4
    - rlang
    - ggplot2/v.3.4.2 
    - tidyr
    - MASS/v.7.3-55 
    - dplyr
    - dartR/v.2.0.4 
    - ggpubr


## 3. Data availability
The sequence data is available on the Sequence Read Archive repository (NCBI) under the BioProject PRJNA1023519. The sequence accession numbers for the samples include SRR26265052 to SRR26265159. Previously published sequences added to our dataset can be retrieved from the [BioProject PRJNA735773](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA735773).

Files used for analyses are provided in the dataset folder.

Reference-based alignments were performed using two transcriptomes of _R. elongatum_ Ehrh. ex Frisvoll [(ID: ABCD)](https://www.ncbi.nlm.nih.gov/biosample/SAMEA104170995/) and _R. varium_ (Mitt.) A.Jaeger [(ID: RDOO)](https://www.ncbi.nlm.nih.gov/biosample/?term=Racomitrium+varium) (Carpenter et al., 2019; Leebens-Mack et al., 2019).

## 4. Citation
If you use part or the entire code in your work, please cite it using the following referece: AVAILABLE SOON.

> [!NOTE]
Be aware that the code may contain issues related to updated software or packages and IS NOT intended to serve as an optimized pipeline but as supplementary information for the research.
