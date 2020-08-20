#!/usr/bin/python3

import os

fasta_dir = "reference_fastas/"
mapping_dir = "mapping_files/"
table_dir = "reference_tables/"
blast_dir = "reference_blasts/"


for file in os.listdir(fasta_dir):
    print (file)
    sample = file.split('.')[0]
    enriched_mapping = mapping_dir + "bwa_mapping_" + sample + "-enriched/concoct_inputtable.tsv"
    original_mapping = mapping_dir + "bwa_mapping_" + sample + "-original/concoct_inputtable.tsv"
    table = table_dir + sample + ".taxonomy." + sample + ".blobDB.table.txt"
    blast = blast_dir + sample + ".fasta.vs.nt.cul5.1e5.megablast.out"
    os.system("~/scripts/mitochondrial_metagenomics/annotate_mitochondria.py " + file + " " + enriched_mapping + " " + original_mapping + " " + table + " " + blast + "> " + sample+".tsv")