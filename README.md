# Target-Enrichment-of-Metazoan-Mitochondrial-DNA-with-Hybridization-Capture-Probes-
code and tutorials for data analyses related to the manuscript submitted to Ecological Indicators

### Downloading the mitochondrial genomes available on NCBI
The complete set of reference mitochondrial genomes is available from https://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/

```bash
## downlaod the data in BASH
wget https://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/mitochondrion.1.1.genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/mitochondrion.2.1.genomic.fna.gz
```

The R script edited from Dunn et al and for construction of Figure1 is called. phylo-availability.r

### FASTA and BED file with probe locations

BED file with probe coordinates and coverage. *mitochondrial-hybridcapture-targets.bed.gz*

The FASTA file used to design the probes was divided into two halfs to adhear to githubs 25 MB per file size limit. Use the  'cat' command in BASH to combine these files together before any analyses.

```bash
cat mitochondrial_genomes_for_probe_design_masked_repeats_half1.fasta.gz mitochondrial_genomes_for_probe_design_masked_repeats_half1.fasta.gz > mitos.fasta
```

### Genome assembly and assessment pipeline.

This script will take a set of paired-end read files and run genome assembly and assessment tools (QUAST, BWA, BUSCO, BLOBTOOLS, BLAST).
It assumes you have a dircetory named 'Sample_SAMPLENAME' that contains two fastq sequences with _R1_ and _R2_ in the names. 

```bash
assembly_pipeline.sh <Sample_directory>
```


### Calculate coverage per contig using BWA and samtools

```bash
bwa_index_and_mapV2.sh <reference_fasta> <forward_reads> <reverse_reads> <sample_name>
```

### Coverage per million reads (CPM) from BED

Coverage data was then normalized based on the starting read counts to produce coverage per million reads (CPM) using an in-house python script called annotate_mitochondria.py. 

```bash
python3 annotate_mitochondria.py <fasta> <coverage_table1> <coverage_table2> <taxonomy_table> <blast_results>
```

### Download mitochondrial annotations
Annotations were completed with the mitos web server. http://mitos.bioinf.uni-leipzig.de/index.py.

The following script will download the results to the server.

```bash
python3 mitos_dowload_files.py <results_url>
```

### Creation of circos plots

Circular figures depicting mitochondrial genome annotations and coverage were constructed using custom scripts in Python v3.6.1 and the Circos software package (Kryzywinski et al. 2009). 

Tested and applied with circos v 0.69.3. A help menu is available for this program, use the --help option.


```bash
python3 circos-mitogenome.py <fasta> <gff> <bed1> <bed2>
```

#### Rscripts for box and scatter plots.
Creation of boxplots, and scatter plots were completed in R Studio v 3.5.1 using ggplot2 (Wickham 2016).

```bash
mito-erichment-box-scatters.r
```
#### Average nucleotide identity with FASTANI

Starting Data: 

1.) Reference sequences used for probe design (seperated into individual fastas)
2.) Mitochondrial genomes from this study


```bash
#step 1: divide mutlifasta into individual fastas
python3 mitochondrial-enrichment-paper/multi-fasta-to-single-fasta.py <reference_multi_fasta>
#move all the output fastas into a directory called 'reference_fastas'

#step 2: create list of reference sequences and wuery sequences
readlink -f reference_fastas/* > probe_list.txt
readlink -f study_fastas/* > query_list.txt

# step 3: run fastaANI on all sequences and identify the ANI for the top match
fastANI --output matrix --rl sequence.fasta --ql sequence.fasta
```
