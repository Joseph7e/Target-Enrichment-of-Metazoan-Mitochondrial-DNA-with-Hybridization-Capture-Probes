#!/usr/bin/sh
#SBATCH --ntasks=1
#SBATCH --job-name="enrich"
#SBATCH --output="enrich"

module purge
module load linuxbrew/colsa

#sample="Sample_199440773"
#
#contig_ref=$sample"contigs.fasta"
#reads_for=/mnt/lustre/hcgs/joseph7e/insert_size_experiment/Project_Fr*/$sample/paired-1.fastq.gz
#reads_rev=/mnt/lustre/hcgs/joseph7e/insert_size_experiment/Project_Fr*/$sample/paired-2.fastq.gz





# Assembly Reference
contig_ref=$1
# Forward Reads
reads_for=$2
# Reverse Reads
reads_rev=$3
# Unique Sample Name
sample=$4

echo mapping $reads_for $reads_rev to $contig_ref

# Make output directory
outdir=bwa_mapping_$sample
mkdir $outdir
echo contig_ref=$contig_ref >> $outdir/mapping_log.txt
echo forward_reads=$reads_for >> $outdir/mapping_log.txt
echo reverse_reads=$reads_rev >> $outdir/mapping_log.txt
echo sample_name=$sample >> $outdir/mapping_log.txt

#index genome assembly
bwa index -a bwtsw $contig_ref

# run bwa
bwa mem -M -t 24  $contig_ref $reads_for $reads_rev > $outdir/raw_mapped.sam

#sort pe sam file and retain only mapped reads
samtools view -@ 24 -Sb -F 4  $outdir/raw_mapped.sam  | samtools sort -@ 24 - -o $outdir/sorted_mapped.bam
samtools index $outdir/sorted_mapped.bam

#extract_mapped_reads
#nohup bam2fastq -o $outdir/mapped#.fastq --aligned $outdir/sorted_mapped.bam &
####nohup bam2fastx -q -M -P -o $outdir/mapped_reads.fastq $outdir/sorted_mapped.bam &

# find mapping matrics
samtools flagstat $outdir/raw_mapped.sam > $outdir/mapping_metrics.txt

#find insert sizes
picard CollectInsertSizeMetrics \
     I=$outdir/sorted_mapped.bam \
     O=$outdir/name_insert_size_metrics.txt \
     H=insert_size_histogram.pdf \
     M=0.5
#rm $outdir/raw_mapped.sam
#gzip $outdir/*.fastq
mv insert_size_histogram.pdf $outdir/insert_size_histogram.pdf

#determine coverage_info
bioawk -c fastx '{ print $name, length($seq) }' < $contig_ref > $outdir/genome.txt
genomeCoverageBed -ibam $outdir/sorted_mapped.bam -g genome.txt > $outdir/coverage.txt



bedtools genomecov -ibam $outdir/sorted_mapped.bam > $outdir/coverage.out

# Generate concoct input table
python3 /mnt/lustre/software/linuxbrew/colsa/Cellar/concoct/0.4.0/libexec/scripts/gen_input_table.py  --isbedfiles $contig_ref $outdir/coverage.out >  $outdir/concoct_inputtable.tsv

#echo nohup bowtie2 $outdir/ --un-conc-gz $outdir/ --al-conc-gz $outdir/ --al-gz $outdir/ --un-gz \
#-x $contig_ref -1 $reads_for -2 $reads_rev -r $unpaired\
#-p 20 -S sample_a8.sam &