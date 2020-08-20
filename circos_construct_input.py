#!/mnt/lustre/software/linuxbrew/colsa/bin/python3

import sys


parse_sam_log = sys.argv[1]
output = "circos_input.txt"

master_phylum_list = ["Ctenophora", "Porifera", "Cnidaria", "Chordata", "Hemichordata", "Echinodermata","Onychophora", "Kinorhyncha", "Priapulida", "Nematoda", "Tardigrada", "Chaetognatha", "Bryozoa", "Placozoa", "Acanthocephala", "Rotifera", "Gnathostomulida", "Gastrotricha", "Platyhelminthes", "Entoprocta", "Nemertea", "Brachiopoda", "Annelida", "Mollusca"]

uniprot_gene_name_dict = {}
input_uniprot_gene_names = "/mnt/lustre/hcgs/joseph7e/mitochondrial_metagenomics/final_references/uniprot_accession_geneNames_and_more.txt"
for line in open(input_uniprot_gene_names):
    if line[0] != '#':
        uniprot, stuff, all = line.rstrip().split(' ')
        genome, gene_name,protein = all.split(':')
        uniprot_gene_name_dict[uniprot] = gene_name

uniprot_locations_dict = {}
input_gene_locations = "/mnt/lustre/hcgs/joseph7e/mitochondrial_metagenomics/final_references/genome_protein_locations.csv"
for line in open(input_gene_locations):
    accession,uniprot,start,stop = line.rstrip().split(',')
    tmp = start
    tmp2 = stop
    if int(start) > int(stop):
        start = tmp2
        stop = tmp
    if accession in uniprot_locations_dict.keys():
        uniprot_locations_dict[accession][uniprot] = [start,stop]
    else:
        uniprot_locations_dict[accession] = {}
        uniprot_locations_dict[accession][uniprot] = [start,stop]
genome_length_dict = {}
input_genome_lengths = "/mnt/lustre/hcgs/joseph7e/mitochondrial_metagenomics/circos/data_files/accession_genome_size.txt"
for line in open(input_genome_lengths):
    acc, phyla,size = line.rstrip().split()
    genome_length_dict[phyla] = [acc,size]




phylum_gene_start_stop = {}
for phyla in genome_length_dict.keys():
    phylum_gene_start_stop[phyla] = {}
    for u,l in uniprot_locations_dict[genome_length_dict[phyla][0]].items():
        try:
            g = uniprot_gene_name_dict[u]
        except:
            #print (phyla, u, ','.join(uniprot_locations_dict[genome_length_dict[phyla][0]]))
            g = u
        s = l[0]
        s2 = l[1]
        phylum_gene_start_stop[phyla][g] = [s,s2]

### construct karyotype file
gene_phylum_dict = {}
kary_handle = open('karyotype_file.txt','w')



count = 1
for phyla in master_phylum_list:
    color = "green"
    kary_handle.writelines("chr - " + phyla + ' ' + phyla + ' 0 ' + str(genome_length_dict[phyla][1]) + " " + color + '\n')
    count += 1

bands_file = open("bands.txt", 'w')
for phyla in master_phylum_list:
    saved_starts_and_stops = []
    for uni,locations in uniprot_locations_dict[genome_length_dict[phyla][0]].items():
        if uni in uniprot_gene_name_dict.keys():
            gene_name = uniprot_gene_name_dict[uni]
        else:
            gene_name = uni
        start = str(locations[0])
        stop = str(locations[1])

        for pairs in saved_starts_and_stops:
            if int(start) >= pairs[0] and int(start) <= pairs[1]:
                start = str(pairs[1]+1)
            if int(stop) >= pairs[0] and int(stop) <= pairs[1]:
                stop = str(pairs[0]-1)
        saved_starts_and_stops.append([int(start),int(stop)])

        bands_file.writelines("band " + phyla + " " + uni + " " + uni + " " + start + " " + stop + " " + gene_name + "\n" )

## add in band information, these are all the geno locations

kary_handle.close()


hit_dict = {}
phylum_dict = {}
quality_phylum_dict = {}
gene_list = ["ATP6","ATP8","ATP9","COX1","COX2","COX3","CYTB","hypothetical_protein","ND1","ND2","ND3","ND4","ND4L","ND5","ND6","none","heg"]

count = 0
flag = False
for line in open(parse_sam_log):
    line = line.rstrip()
    #print (line)
    if flag == False:
        if line.startswith('#acc'):
            pass
        elif line.startswith('####'):
            flag = True
            # print (line)

        else:
            if line.startswith('__________'):
                hit_dict[current_read] = [current_read,ref_acc,read_phylum,ref_phylum,quality,read_taxonomy,ref_taxonomy,read_gene,ref_gene]
                #print (current_read,ref_acc,read_phylum,ref_phylum,quality,read_gene,ref_gene)
                if ref_gene in ["COI", "cox1", "CO1"]:
                    ref_gene = "COX1"
                if ref_gene in ["cox3", "COIII","CO3"]:
                    ref_gene = "COX3"
                if ref_gene in ["CTYB","Cytb", "CYTb","cytb"]:
                    ref_gene = "CYTB"
                if ref_gene in ["COXII","CO2", "cox2", "COII"]:
                    ref_gene = "COX2"
                if ref_gene == "nad1":
                    ref_gene = "ND1"
                if ref_gene == "nad5":
                    ref_gene = "ND5"
                if ref_gene in ["nad4"]:
                    ref_gene = "ND4"
                if ref_gene in ["ATPase6","AT6"]:
                    ref_gene = "ATP6"

                #print (ref_gene)
                index_read_gene = 100
                index_ref_gene = 100
                i = 0
                for g in gene_list:
                    if read_gene == g:
                        index_read_gene = i
                    if ref_gene == g:
                        index_ref_gene = i
                    i += 1
                if read_phylum in phylum_dict.keys():
                    if ref_phylum in phylum_dict[read_phylum].keys():
                        if read_gene in phylum_dict[read_phylum][ref_phylum].keys():
                            phylum_dict[read_phylum][ref_phylum][read_gene][index_ref_gene] += 1
                        else:
                            phylum_dict[read_phylum][ref_phylum][read_gene] = [0]*len(gene_list)
                            phylum_dict[read_phylum][ref_phylum][read_gene][index_ref_gene] += 1
                    else:
                        phylum_dict[read_phylum][ref_phylum] = {}
                        phylum_dict[read_phylum][ref_phylum][read_gene] = [0]*len(gene_list)
                        phylum_dict[read_phylum][ref_phylum][read_gene][index_ref_gene] += 1

                else:
                    phylum_dict[read_phylum] = {}
                    phylum_dict[read_phylum][ref_phylum] = {}
                    phylum_dict[read_phylum][ref_phylum][read_gene] = [0]*len(gene_list)
                    phylum_dict[read_phylum][ref_phylum][read_gene][index_ref_gene] += 1

                ########### ADD QUALITY LISTS
                reference_gene = gene_list[index_ref_gene]
                if read_phylum in quality_phylum_dict.keys():
                    if ref_phylum in quality_phylum_dict[read_phylum].keys():
                        if read_gene in quality_phylum_dict[read_phylum][ref_phylum].keys():
                            if reference_gene in quality_phylum_dict[read_phylum][ref_phylum][read_gene].keys():
                                quality_phylum_dict[read_phylum][ref_phylum][read_gene][reference_gene].append(int(quality))
                            else:
                                quality_phylum_dict[read_phylum][ref_phylum][read_gene][reference_gene] = [int(quality)]
                        else:
                            quality_phylum_dict[read_phylum][ref_phylum][read_gene] = {}
                            quality_phylum_dict[read_phylum][ref_phylum][read_gene][reference_gene] = [int(quality)]
                    else:
                        quality_phylum_dict[read_phylum][ref_phylum] = {}
                        quality_phylum_dict[read_phylum][ref_phylum][read_gene] = {}
                        quality_phylum_dict[read_phylum][ref_phylum][read_gene][reference_gene] = [int(quality)]

                else:
                    quality_phylum_dict[read_phylum] = {}
                    quality_phylum_dict[read_phylum][ref_phylum] = {}
                    quality_phylum_dict[read_phylum][ref_phylum][read_gene] = {}
                    quality_phylum_dict[read_phylum][ref_phylum][read_gene][reference_gene] = [int(quality)]
                ### END ADDING QUALITY LISTS

                #if read_phylum == "Chordata":

                    #print (read_phylum, 10, 20,ref_phylum,12000,12010)
                count = 0
            else:
                if count == 0:
                    current_read,ref_acc,read_phylum,ref_phylum,quality = line.split(',')
                    count = 1
                elif count == 1:
                    read_taxonomy = line.split(',')
                    count = 2
                elif count == 2:
                    ref_taxonomy = line.split(',')
                    count = 3
                else:
                    tmp,read_gene,ref_gene = line.split(' ')

    # else:
    #     print (line)

# create a dictionary with the counts for each match and do an average
#
# gastrotricha, gastrotricha 100 hits, 100 genes match
# gastrotricha, chordata 100 hits, blah blah blah
#
# this way you can incorporate every single run.

links_file = open("links.txt", 'w')
linkends_file = open("linkends.txt", 'w')

for read_p,k in phylum_dict.items():
    ref_size = int(genome_length_dict[read_p][1])
    for ref_p,c in phylum_dict[read_p].items():
        for read_g, w in phylum_dict[read_p][ref_p].items():
            current_index = 0
            for gene_counter in w:
                if gene_counter != 0:


                    thickness = .05*int(gene_counter)
                    if thickness < 1:
                        thickness = 1
                    thickness_add = "thickness="+str(thickness) + "p"
                    read_begin = 0; read_end = 10
                    ref_begin = ref_size-10; ref_end = ref_size
                    ref_gene_name = gene_list[current_index]
                    qualities = quality_phylum_dict[read_p][ref_p][read_g][ref_gene_name]
                    average_quality = sum(qualities)/len(qualities)
                    link_color = "black"
                    if average_quality < 2:
                        link_color = "red"
                    elif average_quality < 10:
                        link_color = "orange"
                    elif average_quality < 25:
                        link_color = "green"
                    elif average_quality < 40:
                        link_color = "blue"
                    else:
                        link_color = "purple"

                    link_color_add = "color="+link_color
                    link_value = "value="+str(average_quality)

                    if gene_counter != len(qualities):
                        print ("the gene quality list and gene counter were nto the same size")
                        sys.exit()
                    if ref_p != read_p:
                        try:
                            read_begin = phylum_gene_start_stop[read_p][read_g][0]
                            read_end = phylum_gene_start_stop[read_p][read_g][1]
                        except:
                            read_begin = 0
                            read_end = 10
                        try:
                            ref_begin = phylum_gene_start_stop[ref_p][ref_gene_name][0]
                            ref_end = phylum_gene_start_stop[ref_p][ref_gene_name][1]
                        except:
                            ref_begin = 0
                            ref_end = 10

                    # if ref_p == read_p:
                    #     ref_begin = ref_size-10; ref_end = ref_size

                    print (read_p,ref_p,read_g,ref_gene_name,gene_counter)
                    links_file.writelines(read_p + " "  + str(read_begin) + " " +  str(read_end) + " " + ref_p + " "  + str(ref_begin) + " " + str(ref_end) +  " " + thickness_add + "," +link_value + "," + link_color_add +  "\n")
                    linkends_file.writelines(ref_p + " "  + str(ref_begin) + " " + str(ref_end) +  " 0 \n")
                current_index += 1



