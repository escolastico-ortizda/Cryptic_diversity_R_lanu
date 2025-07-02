#!/bin/bash

#SBATCH -J stacks_reference-aligned
#SBATCH -o log_files/missing_data_populations_%j.out
#SBATCH -c 1
#SBATCH -p small
#SBATCH --mail-type=ALL
#SBATCH --mail-user=
#SBATCH --time=1-00:00
#SBATCH --mem=1G

### --------------------- Extract missing data per individual in a population --------------------- ###

module load vcftools/0.1.16

# Species dataset (230 samples; -R20)

#vcftools --vcf ./6_populations/R20_230_samples_t125/populations.snps.vcf --missing-indv 
#vcftools --vcf ./6_populations/R20_230_samples_t125/populations.snps.vcf populations.snps.vcf --missing-site 

# Tundra dataset (125 samples; -r80)

vcftools --vcf ./6_populations/all_r80_125samples/populations.snps.vcf --missing-indv 
#vcftools --vcf ./6_populations/all_r80_125samples/populations.snps.vcf populations.snps.vcf --missing-site 
