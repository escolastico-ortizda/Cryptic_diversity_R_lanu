#!/bin/bash

#SBATCH -J fineRAD
#SBATCH -o log_files/fineRAD_-%j.out
#SBATCH -c 1
#SBATCH -p small
#SBATCH --mail-type=ALL
#SBATCH --mail-user=
#SBATCH --time=1-00:00
#SBATCH --mem=5G

module load python/2.7 
module load fineRADstructure/1.0

### --------------------- 0. Reorder according to linkage disequilibrium when you have UNMAPPED LOCI  --------------------- ###

#Reorder according to linkage disequilibrium when you have a UNMAPPED LOCI. 
#Rscript ./fineRADstructure -s 1 -n 500 ./stacks.ref/populations_ref_2020/R20_127_samples/populations.haps.radpainter \
#./Softwares/fineRADstructure-master/R40_127_samples_bwa/R40_127_samples_bwa_reorder.finerad

### --------------------- 1. Run fineRADstructure on the Species dataset (230 samples) --------------------- ###

#RADpainter
#RADpainter paint -p 1 ./6_populations/R20_230_samples_t125/populations.haps.radpainter

#fineRADStructure analysis
#finestructure -x 100000 -y 100000 -z 1000 ./6_populations/R20_230_samples_t125/populations.haps_chunks.out \
#./6_populations/R20_230_samples_t125/populations_haps_chunks.mcmc.xml

finestructure -m T -x 10000 ./6_populations/R20_230_samples_t125/populations.haps_chunks.out  \
./6_populations/R20_230_samples_t125/populations_haps_chunks.mcmc.xml \
./6_populations/R20_230_samples_t125/populations_haps_chunks.mcmc.Tree.xml

# The name of the tree file is important. Always preserve the name : "*.mcmc.Tree.xml. Otherwise fineRAD will not storage the tree.


### --------------------- 2. Run fineRADstructure on the Tundra dataset (125 samples) --------------------- ###

#RADpainter
#RADpainter paint -p 1 ./6_populations/all_r80_125samples/populations.haps.radpainter

#fineRADStructure analysis
#finestructure -x 100000 -y 100000 -z 1000 ./6_populations/all_r80_125samples/populations.haps_chunks.out \
#./6_populations/all_r80_125samples/populations_haps_chunks.mcmc.xml

#finestructure -m T -x 10000 ./6_populations/all_r80_125samples/populations.haps_chunks.out  \
#./6_populations/all_r80_125samples/populations_haps_chunks.mcmc.xml \
#./6_populations/all_r80_125samples/populations_haps_chunks_Tree.xml



