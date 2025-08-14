#!/bin/bash
#SBATCH --job-name=WNV_script
#SBATCH --mail-user=jfauver@unmc.edu
#SBATCH --mail-type=ALL
#SBATCH --time=1-00:00:00
#SBATCH --output=jvar.out
#SBATCH --error=jvar.err
#SBATCH --nodes=1
#SBATCH --mem=24g
#SBATCH --partition=guest

module purge
module load fastp
module load samtools
module load bwa
module load fgbio
module load picard
module load viral_consensus
module load seqkit

#Unique name given to each library. This will be provided in a list file for each library.  
file_base=$1

#Creates a log for each processing step for each library 
log=${file_base}.pipeline.log

{

echo "***********************************" 
echo "begin consensus generation for sample: $file_base" 
echo "***********************************" 

#Hard coded paths to the reference genome you are aligning to, the .bed file for primer postions, .fasta file containing all primer sequences, and .tsv file containg primer pair info 
#NOTE must provide this path for your own system
#Easy way to determine path is to enter directory where each file is located, type "pwd". Copy path below and add /your_file_name for each 
reference_genome=/work/fauverlab/shared/WNV_Nebraska/IDT_Files/WestNile_1_50_consensus_40.fa
IDT_trimprimers=/work/fauverlab/shared/WNV_Nebraska/IDT_Files/WNV.220525.fgbio.trimprimers.txt
IDT_intervals=/work/fauverlab/shared/WNV_Nebraska/IDT_Files/WMV.220525.target.interval_list

#Change into directory with symlinks to paired read data for specific sample
cd ${file_base}

#input variables
f1=${file_base}_L004_R1_001.fastq.gz
f2=${file_base}_L004_R2_001.fastq.gz

#run QC with fastp
fastp --in1 $f1 --in2 $f2 --out1 ${file_base}_L004_R1_001_trimmed.fastq --out2 ${file_base}_L004_R2_001_trimmed.fastq -l 50 -h ${file_base}_QC.html

#new read output
f1=${file_base}_L004_R1_001_trimmed.fastq
f2=${file_base}_L004_R2_001_trimmed.fastq

#align reads to reference genome
bwa mem -t 16 $reference_genome $f1 $f2 -I 95.0,200.0[1000,60] | samtools view -Shb -F 0x100 -F 0x800 | samtools sort -o ${file_base}_sorted.bam

#index sorted .bam file
samtools index ${file_base}_sorted.bam

#trim primers from 5' end of aligned reads

fgbio -Xmx8g TrimPrimers -i ${file_base}_sorted.bam -o ${file_base}_sorted_primertrim.bam -p $IDT_trimprimers -r $reference_genome

## generate QC stats
echo ""
echo "***********************************************"
echo "** SUMMARIZING QC STATS FOR" ${file_base}"..."
echo "***********************************************"
echo ""

samtools flagstat ${file_base}_sorted_primertrim.bam | tee ${file_base}.flagstats.txt
samtools coverage ${file_base}_sorted_primertrim.bam | tee ${file_base}.coverage.txt

picard -Xmx10g CollectHsMetrics -I ${file_base}_sorted_primertrim.bam -TARGET_INTERVALS $IDT_intervals -BAIT_INTERVALS $IDT_intervals -O ${file_base}.hsmetrics.txt -PER_BASE_COVERAGE ${file_base}.per.base.coverage.txt -NEAR_DISTANCE 0 -COVERAGE_CAP 10000 -CLIP_OVERLAPPING_READS false

picard -Xmx10g CollectAlignmentSummaryMetrics -I ${file_base}_sorted_primertrim.bam -R $reference_genome -O ${file_base}.alignment.metrics.txt

##final assembly
echo ""
echo "***********************************************"
echo "** FINAL GENOME ASSEMBLY FOR" ${file_base}"..."
echo "***********************************************"
echo ""

samtools view -b -F 4 ${file_base}_sorted_primertrim.bam | samtools sort -o ${file_base}_sorted_primertrim_mapped.bam

viral_consensus --in_reads ${file_base}_sorted_primertrim_mapped.bam --ref_genome $reference_genome --min_depth 10 --min_freq 0.5 --out_consensus ${file_base}_consensus_genome.fa

#rename header of final consensus with seqkit
seqkit replace -p "\s.+" ${file_base}_consensus_genome.fa | seqkit replace -p "viral_consensus" -r "${file_base}" > ${file_base}_consensus_genome_final.fa

##clean up before we are finished
gzip $f1
gzip $f2
rm *.json
rm ${file_base}_sorted.bam
rm ${file_base}_sorted_primertrim.bam

} 2>&1 | tee -a $log
