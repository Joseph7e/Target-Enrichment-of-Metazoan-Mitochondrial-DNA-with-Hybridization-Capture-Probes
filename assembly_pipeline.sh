#!/usr/bin/sh
#SBATCH --ntasks=1
#SBATCH --mem=409600
#SBATCH --exclude=node117,node118
#SBATCH --job-name="41907B"
#SBATCH --output="41907B"

module purge
module load linuxbrew/colsa

# this script is setup for a fresh Sample directory produuced from HCGS Hiseq2500. Don't have anything else in there! okay?
# HCGS changed its demuxing. Use this.
# for f in *_R1_*; do echo $f; s=$(echo $f | awk -F'_' '{print $1"-"$2"-"$3}'); echo $s; r=$(echo $f | sed 's/_R1_/_R2_/g'); echo $r; mkdir Sample_$s; mv $f Sample_$s; mv $r Sample_$s; done


# Use multiple paths if you need to combine read datasets
path_to_round1="" 
path_to_round2="" 
path_to_round3="" 

reverse_read_identifier=_R2_

# CAREFUL MODE IS NOT ACTIVE BELOW

for dir in $@
do
    case "$dir" in
    */)
        dir="${dir::-1}"
        ;;
    esac
    echo $dir
    #sample=$(echo $dir | cut -d'_' -f 2)
    sample=${dir/Sample_/}
    sample=${sample//_/-}
    echo "___________________________________________"
    echo '#'working on: Sample $sample

    ### check all directories for like files
    forward_growing_path=$(echo $dir/*_R1_*)

    if [ -d $path_to_round1/$dir ] ; then
        forward_growing_path=$forward_growing_path" $(echo $path_to_round1/$dir/*_R1_*)"
    fi

    if [ -d $path_to_round2/$dir ] ; then
        forward_growing_path=$forward_growing_path" $( echo $path_to_round2/$dir/*_R1_*)"
    fi

    if [ -d $path_to_round3/$dir ] ; then
        forward_growing_path=$forward_growing_path" $( echo $path_to_round3/$dir/*_R1_*)"
    fi


    echo '#'concatenating the following files:
    echo '#'FORWARD: $forward_growing_path

    reverse_growing_path=$(echo $dir/*$reverse_read_identifier*)

    if [ -d $path_to_round1/$dir ] ; then
        reverse_growing_path=$reverse_growing_path" $(echo $path_to_round1/$dir/*$reverse_read_identifier*)"
    fi
        echo $path_to_round1/$dir
    if [ -d $path_to_round2/$dir ] ; then
        reverse_growing_path=$reverse_growing_path" $( echo $path_to_round2/$dir/*$reverse_read_identifier*)"
    fi

    if [ -d $path_to_round3/$dir ] ; then
        reverse_growing_path=$reverse_growing_path" $( echo $path_to_round3/$dir/*$reverse_read_identifier*)"
    fi

    echo '#'REVERSE: $reverse_growing_path

    #### concatenate files
    raw_forward=$dir/raw_forward_$sample.fastq.gz
    raw_reverse=$dir/raw_reverse_$sample.fastq.gz
    echo cat $forward_growing_path ">" $raw_forward
    echo cat $reverse_growing_path ">" $raw_reverse


    ###############################################################XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX customize reads going into the workflow
    ## Be sure to remove this when not in use
#    forward_growing_path=$dir/combined-1_round3.fastq.gz
#    reverse_growing_path=$dir/combined-2_round3.fastq.gz

   ########################################################### remove this


    cat $forward_growing_path > $raw_forward
    cat $reverse_growing_path > $raw_reverse

    echo '#'running TRIMMOMATIC:
    forward=$dir/trimmed_forward_$sample.fastq.gz
    reverse=$dir/trimmed_reverse_$sample.fastq.gz


    echo trimmomatic PE -threads 24 $raw_forward $raw_reverse $forward $dir/unpaired-1.fastq.gz $reverse $dir/unpaired-2.fastq.gz ILLUMINACLIP:/mnt/lustre/hcgs/joseph7e/program_files/NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:5 MINLEN:36
    trimmomatic PE -threads 24 $raw_forward $raw_reverse $forward $dir/unpaired-1.fastq.gz $reverse $dir/unpaired-2.fastq.gz ILLUMINACLIP:/mnt/lustre/hcgs/joseph7e/program_files/TruSeq2-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:5 MINLEN:36



    echo concatenating unpaired_reads:
    unpaired=$dir/trimmed_unpaired_$sample.fastq.gz
    echo cat $dir/unpaired-1.fastq.gz $dir/unpaired-2.fastq.gz ">" $unpaired
    cat $dir/unpaired-1.fastq.gz $dir/unpaired-2.fastq.gz > $unpaired



    echo running SPADES:
    spades_dir=$dir/spades_assembly
    echo spades.py -1 $forward -2 $reverse -s $unpaired -t 24 -o $spades_dir -m 400 --careful
    # CAREFUL MODE IS ACTIVE
    spades.py -1 $forward -2 $reverse -s $unpaired -t 24 -o $spades_dir -m 400 --careful
    #spades.py -1 $forward -2 $reverse -s $unpaired -t 24 -o $spades_dir -m 400
    #spades.py -t 24 -o $spades_dir -m 600 --continue
    #srun spades.py -o $dir/spades_assembly --continue -t 24 -m 400
    fasta=$spades_dir/contigs.fasta


    echo '#'running QUAST:
    echo mkdir ../QUAST_results
    echo quast.py -o QUAST_results/quast_$sample $fasta
    quast.py -o QUAST_results/quast_$sample $fasta

    ### construct new directpry for downstream analysis
    mkdir initial_spades_contigs
    cp $fasta initial_spades_contigs/$sample.fasta
    cd initial_spades_contigs

    echo '#'running BUSCO:
    echo ~/scripts/quality_check_genome_busco_bbacteria.sh $sample.fasta
    #~/scripts/quality_check_genome_busco_bacteria.sh $sample.fasta
    ~/scripts/quality_check_genome_busco_metazoa.sh $sample.fasta

    echo '#'running MAPPING:
    ~/scripts/bwa_index_and_mapV2.sh $sample.fasta ../$forward ../$reverse $sample

    echo '#'running BLOBTOOLS:
    ~/scripts/blob_plots/blob_blast.sh $sample.fasta
    ~/scripts/blob_plots/blob_create.sh $sample.fasta $sample.fasta.vs.nt.cul5.1e5.megablast.out bwa_mapping_$sample/sorted_mapped.bam $sample
    cd ../


done
#sort pe sam file and retain only mapped reads
