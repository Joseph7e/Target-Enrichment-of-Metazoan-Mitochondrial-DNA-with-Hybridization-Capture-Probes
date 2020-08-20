#!/mnt/lustre/software/linuxbrew/colsa/bin/python3

import sys,os

results_url = sys.argv[1]
#mitos_prefix = "http://mitos.bioinf.uni-leipzig.de/"
mitos_prefix = "http://mitos2.bioinf.uni-leipzig.de/" # for mitos2
sample = sys.argv[2]

os.system('mkdir' +  ' mitos_' + sample)

def Wget(url,outputname):
    os.system('wget '+ '"' + url + '"' + ' -O ' + 'mitos_'+sample+'/'+outputname)

Wget(results_url,sample + '_mitos_html_file.txt')

for line in open('mitos_'+sample+'/'+sample + '_mitos_html_file.txt'):
    if ">GFF file<" in line:
        url = line.split('"')[1]
        Wget(mitos_prefix + url,sample+'.gff')
    if ">FAS file<" in line:
        url = line.split('"')[1]
        Wget(mitos_prefix + url,sample+'.fas')
    if ">FAA file<" in line:
        url = line.split('"')[1]
        Wget(mitos_prefix + url,sample+'.faa')
    if ">TBL file<" in line:
        url = line.split('"')[1]
        Wget(mitos_prefix + url,sample+'.tbl')
    if ">BED file<" in line:
        url = line.split('"')[1]
        Wget(mitos_prefix + url,sample+'.bed')
    if ">protein plot<" in line:
        url = line.split('"')[1]
        Wget(mitos_prefix + url,sample+'_protein_plot.pdf')
    if ">ncRNA plot<" in line:
        url = line.split('"')[1]
