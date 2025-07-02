#!/bin/bash

#SBATCH -J stacks_reference
#SBATCH -o log_files/mapping_bwa%j.out
#SBATCH -c 2
#SBATCH -p medium
#SBATCH --mail-type=ALL
#SBATCH --mail-user=
#SBATCH --time=1-00:00
#SBATCH --mem=2G

cd ~/amsip1/Dennis_Clonality_Chapter_II

### --------------------- 1. Create index for transcriptomes of Racomitrium varium and R. elongatum --------------------- ###

# Load the software with module if applicable:
module load bwa/0.7.17 
module load samtools/1.8
#module load bowtie/2.3.4.1
module load BBTools/36.92

# Place the transcriptomes in "info" folder. If you already combine the transcriptomes or genomes into one file, place it inside the folder.
# Make index for alignement

# BOWTIE Aligment 
#bowtie2-build -f ./2_info/2_reference_based/RefTrans_2Racomitrium.fa ./2_info/2_reference_based/RacoTrans

# BWA Alignment
#bwa index 2_info/2_reference_based/RefTrans_2Racomitrium.fa 2_info/2_reference_based/RefTrans_2Racomitrium

# Create folders
#mkdir -p sam_files
#mkdir -p bam_files
#mkdir -p fastq_files
#mkdir -p fastq_reformatted

### --------------------- 2. Mapping reads to the reference file --------------------- ###

# In this case, we performed the mapping per run for memmory reasons.
FILES=~/amsip1/Dennis_Clonality_Chapter_II/3_process_radtags/GBS_2018_t125/*.fq.gz
#echo -e 'Sample\tread_no\tmapped_no\tper_mapped' > bowtie2_mapping.out
echo -e 'Sample\tread_no\tmapped_no\tper_mapped' > bwa_mapping.out
for f in $FILES
do
	name=${f%%.fq.gz}
    	name=${name##*/}
	FqNo=$(zcat $f | echo $((`wc -l`/4)))
#	#bowtie2 -p 6 -x ~/Dennis_GBS_protocol/info/RacoTrans -S ./sam_files/${name}.sam -U $f --al-gz ./fastq_files/${name}.mapped.fastq.gz -D 20 -R 3 -N 0 -L 20 -i S,1,0.50 -N 1
	bwa mem -t 6 -o ./4_alignment/sam_files/${name}.sam ./2_info/2_reference_based/RefTrans_2Racomitrium.fa $f
	samtools view --threads 6 -b -S ./4_alignment/sam_files/${name}.sam | samtools sort --threads 6 -o ./4_alignment/bam_files_t125/${name}.bam
	reformat.sh in=./4_alignment/bam_files_t125/${name}.bam out=./4_alignment/fastq_reformatted/${name}.reformatted.fastq mappedonly &>> bbtools_fastq.out
	gzip ./4_alignment/fastq_reformatted/${name}.reformatted.fastq
	SamNo=`samtools view -c -F 260 ./4_alignment/bam_files_t125/${name}.bam`
	Per=$(bc <<<"scale=2; $SamNo / $FqNo")
#	#echo -e "$name\t$FqNo\t$SamNo\t$Per" >> bowtie2_mapping.out
	echo -e "$name\t$FqNo\t$SamNo\t$Per" >> bwa_mapping.out
done
#echo "bowtie2 mapping is done."
echo "bwa mapping is done."

