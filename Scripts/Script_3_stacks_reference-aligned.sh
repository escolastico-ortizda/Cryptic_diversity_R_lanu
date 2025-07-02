#!/bin/bash

#SBATCH -J stacks_reference-aligned
#SBATCH -o log_files/populations_%j.out
#SBATCH -c 1
#SBATCH -p small
#SBATCH --mail-type=ALL
#SBATCH --mail-user=
#SBATCH --time=1-00:00
#SBATCH --mem=5G

# Source path
src=~/amsip1/Dennis_Clonality_Chapter_II 

### --------------------- 1. Species dataset (230 samples) -process_radtgas t125 aligning reads per sample, calling variant sites in the patch, genotyping in each individual --------------------- >

#populations -P $src/4_alignment/gstacks_reference_t125_all_235_samples/ -M $src/2_info/3_popmaps/popmap_all_230_samples.txt --max-obs-het 0 -R 20 --min-maf 0.05 --vcf \
#--radpainter --structure \
#-O $src/6_populations/R20_230_samples_t125/ &> populations.oe

### Unlinked SNPs

#populations -P $src/4_alignment/gstacks_reference_t125_all_235_samples/ -M $src/2_info/3_popmaps/popmap_all_230_samples.txt --max-obs-het 0 -R 20 --min-maf 0.05 --vcf \
#--fasta-samples --structure --write-random-snp --phylip --phylip-var-all \
#-O $src/6_populations/R20_230_samples_t125_single/ &> populations.oe


### --------------------- 2. Tundra samples (1254 samples)- Aligning reads per sample, calling variant sites in the patch, genotyping in each individual --------------------- ###

#gstacks -I $src/4_alignment/bam_files/ -M $src/2_info/3_popmaps/popmap_all.txt -O $src/gstacks_reference/ -t 4

#populations -P $src/gstacks_reference/ -M $src/2_info/3_popmaps/popmap_all.txt --max-obs-het 0 -r 80 -H --min-maf 0.05 --vcf --treemix \
#-fstats --fst-correction --radpainter --bootstrap-fst --bootstrap-pifis --bootstrap-phist \
#-O $src/6_populations/all_r80_125samples/ &> populations.oe

### --------------------- 3. Population scale - All groups (66 samples). Even number of samples for habitat type on each genetic group (AB, CD) --------------------- >

#populations -P $src/5_gstacks_reference/ -M $src/2_info/3_popmaps/popmap_population_scale_even.txt --max-obs-het 0 -r 80 -H --min-maf 0.05 --vcf --treemix -fstats --fst-correction --bootstrap-fst --bootstrap-pifis --bootstrap-phist --radpainter -O /$src/6_populations/population_scale_r80_75_samples_even/

### --------------------- 4. Population scale - Genetic groups. Even number of samples for habitat type on each genetic group (AB, CD) --------------------- >

#populations -P $src/4_alignment/gstacks_reference_t125_all_235_samples/ -M $src/2_info/3_popmaps/popmap_population_scale_AB.txt --max-obs-het 0 -r 80 --min-maf 0.05 --vcf \
#-fstats --fst-correction  --bootstrap-fst --bootstrap-pifis --bootstrap-phist \
#-O $src/6_populations/population_scale_even_AB_habitat/

populations -P $src/4_alignment/gstacks_reference_t125_all_235_samples/ -M $src/2_info/3_popmaps/popmap_population_scale_CD.txt --max-obs-het 0 -r 80 --min-maf 0.05 --vcf \
-fstats --fst-correction  --bootstrap-fst --bootstrap-pifis --bootstrap-phist \
-O $src/6_populations/population_scale_even_CD_habitat/
