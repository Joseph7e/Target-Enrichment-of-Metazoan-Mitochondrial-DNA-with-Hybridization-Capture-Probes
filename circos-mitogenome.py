#!/usr/bin/python3
# Author: Joseph Sevigny
# Date: Feb, 28th 2020
# Purpose: automatic construction of circos plot from gff, fasta, and bed file
# tested for circos v 0.69.3

import os
import shutil
import sys
import argparse
from Bio import SeqIO




# Parse arguments
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
#OPTIONAL ARGUMENTS
parser.add_argument("-v", "--verbose", help="increase output verbosity", action="store_true")
# Circos installation stuff
parser.add_argument("--circos_dir", help="path to circos install", type=str, default="/home/genome/joseph7e/bin/circos-0.69-3/")
parser.add_argument("--circos", help="path to circos program", type=str, default="/home/genome/joseph7e/bin/circos-0.69-3/bin/circos")

# Output options
parser.add_argument("-o", "--outdir", help="directory name for output", type=str, default="circos-mitogenome/")
parser.add_argument("--force", help="force delete and rewrite output dir", action="store_true")

# Required positional arguments
parser.add_argument("fasta", help="FASTA file of mitochondrial genome")
parser.add_argument("gff", help="path to standard gff file associated with FASTA")
parser.add_argument("bed", help="""bed graph file, make with this  
bedtools genomecov -ibam sorted_mapped.bam -g ../../reference_fastas/MDD02-FG02-Achotines-A1.fasta -bga > coverage_histogram.bedg
grep NODE_6_ coverage_histogram.bedg | sed 's/NODE_6_length_15196_cov_22.680403/MDD02-FG02-Achotines-A1/g' mito_contig_coverage.bedg | sed 's/\t/ /g' > mito_contig_coverage.bedg.fixed
""")
parser.add_argument("bed2", help="another bed graph file, make as above")

args = parser.parse_args()

if os.path.isdir(args.outdir):
    if args.force:
        print ("deleting old directory and making new one")
        shutil.rmtree(args.outdir)
    else:
        print("Output directory already exists, please remove, use a new name, or force with --force")
        sys.exit()
os.mkdir(args.outdir)

# construct circos input files

# construct karotype from fasta and gff
output_karyotype = open(args.outdir + 'karyotype.txt','w')
for seq_record in SeqIO.parse(args.fasta, "fasta"):
    length = len(seq_record)
    name = seq_record.id.split('_')[0]
    name2 = seq_record.id
    output_karyotype.writelines("chr - {} 1 0 {} {}\n".format(name, length, name2))


output_genes = open(args.outdir + 'genes_protein.txt', 'w')
output_tRNA = open(args.outdir + 'genes_trna.txt', 'w')
output_rRNA = open(args.outdir + 'genes_rrna.txt', 'w')
output_labels = open(args.outdir + 'gene_labels.txt', 'w')
saved_starts_and_stops = [] # ensure no overlapping bands
for line in open(args.gff):
    # the following is based off of mitos2 gff files.
    elements = line.rstrip().split('\t')
    contig = elements[0].split('_')[0]
    type = elements[2]
    start = elements[3]
    stop = elements[4]
    gene = elements[-1].split('=')[-1]

    for pairs in saved_starts_and_stops:
        if int(start) >= pairs[0] and int(start) <= pairs[1]:
            start = str(pairs[1] + 1)
        if int(stop) >= pairs[0] and int(stop) <= pairs[1]:
            stop = str(pairs[0] - 1)

    saved_starts_and_stops.append([int(start), int(stop)])
    if type == 'gene':
        output_genes.writelines("{} {} {} {}\n".format(contig, start, stop, gene))
    elif type == 'tRNA':
        output_tRNA.writelines("{} {} {} {}\n".format(contig, start, stop, gene))
    elif type == 'rRNA':
        output_rRNA.writelines("{} {} {} {}\n".format(contig, start, stop, gene))
    output_labels.writelines("{} {} {} {}\n".format(contig, start, stop, gene))

# write circos config file
output_config = open(args.outdir + 'config.txt', 'w')
output_config.writelines("""
# circos.conf
karyotype = karyotype.txt

<ideogram>
<spacing>
default = 0.005r
</spacing>

radius    = 0.9r
thickness = 5p
fill      = yes
</ideogram>


<highlights>
# the default value for z-depth and fill_color for all highlights
z = 0
fill_color = green

# the first set will be drawing from 0.6x 1x-25pixels of the ideogram
# radius and will be green (color by default)

<highlight>
file       = genes_protein.txt
r0         = .95r
r1         = 1r
fill_color = dblue
stroke_color = black
stroke_thickness = 3p
</highlight>


<highlight>
file       = genes_trna.txt
r0         = .95r
r1         = 1r
fill_color = black
stroke_color = black
stroke_thickness = 0.2
</highlight>


<highlight>
file       = genes_rrna.txt
r0         = .95r
r1         = 1r
fill_color = red
stroke_color = black
stroke_thickness = 3p
</highlight>
</highlights>
########################################################## gene labels
<plots>
<plot>
type             = text
color            = black
file             = gene_labels.txt

r0 = 1.01r
r1 = 1.01r+200p

link_dims      = 0p,0,50p,0p,10p
link_thickness = 2p
link_color     = red

label_size   = 34p
label_font   = condensed

padding  = 0p
rpadding = 0p

</plot>
################################################################# bedg
<plot>
# construct the histogram based on bedg file.
type      = histogram
file      = {}


r1        = 0.73r
r0        = 0.54r

stroke_type = outline
thickness   = 4
color       = vdgrey
extend_bin  = no

<backgrounds>
<background>
color = vvlgrey
</background>
</backgrounds>

<axes>
<axis>
spacing   = 0.1r
color     = lgrey
thickness = 2
</axis>
</axes>

############################################ bedg 2
</plot>
<plot>
# construct the histogram based on bedg file.
type      = histogram
file      = {}

r1        = 0.94r
r0        = 0.75r

# min = 0
# max = 5000

stroke_type = outline
thickness   = 4
color       = dgreen
extend_bin  = no

<backgrounds>
<background>
color = vvlgrey
</background>
</backgrounds>

<axes>
<axis>
spacing   = 0.1r
color     = lgrey
thickness = 2
</axis>
</axes>


</plot>

</plots>


################################################################
# The remaining content is standard and required.

<image>
# Included from Circos distribution.
<<include etc/image.conf>>                
</image>

# RGB/HSV color definitions, color lists, location of fonts, fill patterns.
# Included from Circos distribution.
<<include etc/colors_fonts_patterns.conf>>

# Debugging, I/O an dother system parameters
# Included from Circos distribution.
<<include etc/housekeeping.conf>>
""".format(args.bed, args.bed2))


# run circos
