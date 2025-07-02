#!/bin/bash

#SBATCH -J GBS-cleaning
#SBATCH -o log_files/process_radtags-%j.out
#SBATCH -c 2
#SBATCH -p small
#SBATCH --mail-type=ALL
#SBATCH --mail-user=
#SBATCH --time=1-00:00
#SBATCH --mem=3G

# Source path
src=~/amsip1/Dennis_Clonality_Chapter_II

# Process_radtags needs only 4G and 2 cores !!! 

### --------------------- 1. Check read quality using fastqc --------------------- ###

module load fastqc/0.11.8

# Genotyping-by-sequencing data of Racomitrium lanuginosum
# Use no more than 15 GB of memory for this software

#fastqc $src/0_raw_data/Villarreal_Dennis_20211129_P1.fastq.gz -outdir $src/1_Fastqc_output
#fastqc $src/0_raw_data/Villarreal_Dennis_20211129_P2.fastq.gz -outdir $src/1_Fastqc_output/

### --------------------- 2. Trimmomatic for trimming (optional) --------------------- ###

# Trimmomatic: Upload directory to the server. 
# 10 GB memory

#java -jar $src/Trimmomatic-0.36/trimmomatic-0.36.jar SE -phred33 $src/0_raw_data/Villarreal_Dennis_20211129_P1.fastq.gz $src/0_raw_data/Villarreal_Dennis_20211129_P1_trimmed125.fastq.gz \
#ILLUMINACLIP:$src/Trimmomatic-0.36/adapters/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:125

#java -jar $src/Trimmomatic-0.36/trimmomatic-0.36.jar SE -phred33 $src/0_raw_data/Villarreal_Dennis_20211129_P2.fastq.gz $src/0_raw_data/Villarreal_Dennis_20211129_P2_trimmed125.fastq.gz \
#ILLUMINACLIP:$src/Trimmomatic-0.36/adapters/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:125

#fastqc $src/0_raw_data/Villarreal_Dennis_20211129_P1_trimmed125.fastq.gz -outdir $src/1_Fastqc_output
#fastqc $src/0_raw_data/Villarreal_Dennis_20211129_P2_trimmed125.fastq.gz -outdir $src/1_Fastqc_output/


### --------------------- 3. Process_radtags to  cleaning, trimming and demultiplexing --------------------- ###

# Quality check, trimming, demultiplexing using STACKS. Use 20 GB of memory. Check the paths before running each analysis. 

module load gcc/6.2.0
module load stacks/2.5


# Test with --adapter command (optional)

# GBS data plate 1
#process_radtags -f $src/0_raw_data/Villarreal_Dennis_20211129_P1.fastq.gz -b $src/2_info/1_barcodes/Barcodes_2021_Plate1.txt -o $src/3_process_radtags/Plate1_t125 --renz_1 pstI --renz_2 mspI \
#-c -r -q -s 10 -t 125 --adapter_1 CGAGATCGGAAGAGCGGGGAGCTTAAGC --adapter_mm 1

# GBS data plate 2
#process_radtags -f $src/0_raw_data/Villarreal_Dennis_20211129_P2.fastq.gz -b $src/2_info/1_barcodes/Barcodes_2021_Plate2.txt -o $src/3_process_radtags/Plate2_t125 --renz_1 pstI --renz_2 mspI \
#-c -r -q -s 10 -t 125 --adapter_1 CGAGATCGGAAGAGCGGGGAGCTTAAGC --adapter_mm 1

# Results: Demultiplexing with the --adapter command will produce no-adapter contamination but some samples will have very low reads (e.g. 247 reads in JC3091_399-A5).
# The default process_radtags (no --adapter command) produced better results: less than 10% of adapter containt and a good amount of reads per sample.


# Main demultiplexing 

# GBS data 2021 plate 1 
#process_radtags -f $src/0_raw_data/Villarreal_Dennis_20211129_P1.fastq.gz -b $src/2_info/1_barcodes/Barcodes_2021_Plate1.txt -o $src/3_process_radtags/Plate1_t125 --inline_null --renz_1 pstI --renz_2 mspI \
#-c -r -t 125 --len_limit 125 -q -s 10 --barcode_dist_1 1

# GBS data 2021 plate 2
#process_radtags -f $src/0_raw_data/Villarreal_Dennis_20211129_P2.fastq.gz -b $src/2_info/1_barcodes/Barcodes_2021_Plate2.txt -o $src/3_process_radtags/Plate2_t125 --inline_null --renz_1 pstI --renz_2 mspI \
#-c -r -t 125 --len_limit 125 -q -s 10 --barcode_dist_1 1

# GBS data 2018 (120 length)
#process_radtags -f $src/0_raw_data/racomitrium_p01_c01_copy.fastq.gz -b $src/2_info/1_barcodes/R_lanu_2018_barcodes.txt -o $src/3_process_radtags/GBS_2018_t125 --inline_null --renz_1 pstI --renz_2 mspI \
#-c -r -t 125 --len_limit 125 -q -s 10 --barcode_dist_1 1

# GBS data 2020 (120 length)
process_radtags -p $src/0_raw_data/GBS_2020 -b $src/2_info/1_barcodes/R_lanu_2020_barcodes.txt -o $src/3_process_radtags/GBS_2020_t125 --inline_null --renz_1 pstI --renz_2 mspI \
-c -r -t 125 --len_limit 125 -q -s 10 --barcode_dist_1 1



# Check the read quality of some of the resulting demultiplexed files

# Plate 1 with --adapter command
#fastqc $src/3_process_radtags/Plate1_t125/JC3091_399-A5.fq.gz -outdir $src/1_Fastqc_output
#fastqc $src/3_process_radtags/Plate1_t125/JC2950_407-E1.fq.gz -outdir $src/1_Fastqc_output

# Plate 2 with --adapter command
#fastqc $src/3_process_radtags/Plate2_t125/JC3302_393-D5.fq.gz -outdir $src/1_Fastqc_output
#fastqc $src/3_process_radtags/Plate2_t125/JC3256_393-E3.fq.gz -outdir $src/1_Fastqc_output

# Plate 1 processed with Trimmomatic
#fastqc $src/3_process_radtags/Plate1_trimmo_125/JC3303_416-A5.fq.gz -outdir $src/1_Fastqc_output
#fastqc $src/3_process_radtags/Plate1_trimmo_125/JC2950_407-E1.fq.gz -outdir $src/1_Fastqc_output

# Plate 1 
#fastqc $src/3_process_radtags/Plate1_t125_Marta_script/JC3091_399-A5.fq.gz -outdir $src/1_Fastqc_output
#fastqc $src/3_process_radtags/Plate1_t125_Marta_script/JC2950_407-E1.fq.gz -outdir $src/1_Fastqc_output

