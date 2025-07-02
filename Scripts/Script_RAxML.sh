#!/bin/bash

#SBATCH -J RAxML
#SBATCH -o RAxML.out
#SBATCH -c 4
#SBATCH -p medium
#SBATCH --mail-type=ALL
#SBATCH --mail-user=
#SBATCH --time=4-00:00
#SBATCH --mem=10G

module load RAxML/8.2.9
src=~/amsip1/Dennis_Clonality_Chapter_II

#mkdir -p ls _single_RAxML
raxmlHPC-PTHREADS-SSE3 -n all_R20_235_samples.output -s $src/6_populations/all_R20_235_samples_single/populations.fixed.phylip \
-m ASC_GTRGAMMA --asc-corr=lewis -f a -p 164655 -x 12345 -# 100 -T 6
