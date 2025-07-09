
# Population genetic structure and clonal analyses of _Racomitrium lanuginosum_ in the forest-tundra transition zone.

## 1. Overview
This repository contains the scripts and dataset related to the project _“Decoding cryptic diversity of moss populations in a forest-tundra ecotone”. The main goal of the code is to provide detailed information on the project/article analyses (reproducibility)_.

> **Abstract**: Cryptic speciation is widespread among bryophytes, and it appears common in arctic and subarctic mosses where cryptic lineages may occur sympatrically in the environment. However, cryptic lineages are rarely considered in genetic diversity assessments. This situation poses a challenge, as the complex population structure resulting from cryptic speciation, along with factors like the predominance of clonality, can lead to inaccurate estimates of biodiversity. In this sense, we studied two populations of the subarctic moss _Racomitrium lanuginosum_ in the forest-tundra ecotone to test the impact of accounting for differentiated genetic groups on moss genetic diversity estimates, clonal structure and microbial covariation. We performed genotyping-by-sequencing to infer genetic diversity and structure in the forest tundra and the shrub tundra. Genetic groups were identified using haplotype-based coancestry matrices and phylogenetic analyses. The clonal structure was explored by determining multilocus genotypes at the population and a finer scale (225 cm<sup>2</sup>). The covariation between genetic groups and microbial communities (bacterial and diazotrophic) was explored. The recognition of cryptic lineages in genetic diversity estimations revealed differences between habitats that remained undetected when treating _R. lanuginosum_ as a single species. Clonal growth seemed to predominate and affect the genetic structure at the local scale. Finally, genetic groups did not host specific microbiomes, suggesting that moss microbial associations in the forest-tundra ecotone did strongly respond to a genetic component. This study highlights the importance of accounting for cryptic lineages in genetic diversity estimations for a precise biodiversity assessment in subarctic and arctic ecosystems.

## 2. Software information
The code was run on the bioinformatic platform of the [Institut de Biologie Intégrative et des Systèmes (IBIS)](https://www.ibis.ulaval.ca/en/services-2/bioinformatics/documentation-servers/computer-description/) at Laval University and on the statistic environment R. Check each [script](Scripts/) for details.
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
The sequence data is available on the Sequence Read Archive repository (NCBI) under the [BioProject PRJNA1023519](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA1023519). The sequence accession numbers for the samples include SRR26265052 to SRR26265159. Previously published sequences added to our dataset can be retrieved from the [BioProject PRJNA735773](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA735773).

Reference-based alignments were performed using two transcriptomes: _R. elongatum_ Ehrh. ex Frisvoll [(ID: ABCD)](https://www.ncbi.nlm.nih.gov/biosample/SAMEA104170995/) and _R. varium_ (Mitt.) A.Jaeger [(ID: RDOO)](https://www.ncbi.nlm.nih.gov/biosample/?term=Racomitrium+varium) (Carpenter et al., 2019; Leebens-Mack et al., 2019).

The [files](Data/) used for all the analyses are provided in the repository. The [scripts](Scripts/) for the genomic analyses are stored as .txt files or similar formats. The detailed analyses of the covariation between genetic groups and microbial communities (bacterial and diazotrophic) are presented  in separate workflows:

### [Bacterial workflow: Analyses based on metabarcoding  of the 16S rRNA gene.](https://escolastico-ortizda.github.io/Cryptic_diversity_Bacteria/)
### [Diazotroph workflow: Analyses based on metabarcoding of the nifH gene.](https://escolastico-ortizda.github.io/Cryptic_diversity_Diazotroph/)

## 4. Citation
If you use part or the entire code in your work, please cite it using the following reference: **[Escolastico-Ortiz, D.A., Derome, N., & Villarreal-A., J.C. 2025. Decoding cryptic diversity of moss populations in a forest-tundra ecotone. Preprint.](https://doi.org/10.1101/2025.06.25.660410)**.

> [!NOTE]
Be aware that the code may contain issues related to updated software or packages and IS NOT intended to serve as an optimized pipeline but as supplementary information for the research.
