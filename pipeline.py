#!/usr/bin/env python3

import numpy as np
from estimate_length import *
from collections import defaultdict
import itertools	
import os
import sys
import time
import random
import subprocess

# URL to get annotation GTF from
gtf_url = 'ftp://ftp.ensembl.org/pub/current_gtf/homo_sapiens/Homo_sapiens.GRCh38.84.chr.gtf.gz'

# Simulate a full workflow.
# Folder where bamfile, bioanalyzer profile, genome and gtf are in
folder_in = 'test_data'
gtf = os.path.join(folder_in, 'Homo_sapiens.GRCh38.84_chr9.gtf.gz')
genome = os.path.join(folder_in, 'Homo_sapiens.GRCh38.dna.chromosome.9.fa')

# Create output directory for storing everything
folder_out = os.path.join(folder_in, 'output')
try:
    os.mkdir(folder_out)
except Exception:
    pass

print(folder_out)

### 1. Download annotation
print ("downloading annotation for chr9...", end=" ", flush=True)
if os.path.isfile(gtf):
    print ('skipping [ file already exists ]')
else:
    start_time = time.time()
    subprocess.call('wget -q ' + gtf_url + ' -O - | zcat | grep "^9\t" | gzip --best > ' + gtf, shell=True)
    print ('done [', round(time.time() - start_time, 2), 'seconds ]')

### 2. Extract utr information from gtf file
print ("extracting 3'UTR information ...", end=" ", flush=True)
if os.path.isfile(os.path.join(folder_out, 'utr_annotation.bed')):
    print ('skipping [ file already exists ]')
else:
    start_time = time.time()
    old_stdout = sys.stdout
    sys.stdout = open(os.path.join(folder_out, 'utr_annotation_temp.bed'), 'w')
    extract_three_prime_utr_information(gtf, bed_name_attributes = ["gene_name"])
    sys.stdout = old_stdout
    print ('done [', round(time.time() - start_time, 2), 'seconds ]')

    ### 2.1 Clean utr from haplotypes and junk chromosomes
    with open(os.path.join(folder_out, 'utr_annotation_temp.bed'), 'r') as fin, open(os.path.join(folder_out, 'utr_annotation_unsorted.bed'), 'w') as fout:
        for line in fin:
            if line.startswith('chrGL') or line.startswith('chrKI'):
                continue
            else:
                fout.write(line)
    os.remove(os.path.join(folder_out, 'utr_annotation_temp.bed'))

    ### 2.2 Sort the utr file alphabetically
    subprocess.call('sort -V test_data/output/utr_annotation_unsorted.bed > test_data/output/utr_annotation.bed', shell=True)
    os.remove(os.path.join(folder_out, 'utr_annotation_unsorted.bed'))

### 3. Extract polyA intervals from genome
print ('extracting polyA intervals from genome ...', end=" ", flush=True)
if os.path.isfile(os.path.join(folder_out, 'pAi.bed')):
    print ('skipping [ file already exists ]')
else:
    start_time = time.time()
    extract_pAi_from_genome(genome, window=10, occurences=7, consecutive=6)
    os.remove('pAi_temp.bed')
    os.rename('pAi.bed', os.path.join(folder_out, 'pAi.bed'))
    print ('done [', round(time.time() - start_time, 2), 'seconds ]')

### 4. Add gene information to polyA intervals
print ('adding gene annotation to pAi intervals ...', end=" ", flush=True)
if os.path.isfile(os.path.join(folder_out, 'pAi_gene.bed')):
    print ('skipping [ file already exists ]')
else:
    start_time = time.time()
    annotate_pAi_with_gene(os.path.join(folder_out, 'pAi.bed'), 
                           os.path.join(folder_out, 'utr_annotation.bed'))
    os.rename('pAi_gene.bed', os.path.join(folder_out, 'pAi_gene.bed'))
    print ('done [', round(time.time() - start_time, 2), 'seconds ]')

### 5. Merge polyA intervals with 3'UTRs into a dictionary
print ("merging polyA intervals with 3'UTR ...", end=" ", flush=True)
start_time = time.time()
pAi_full = merge_pAi_and_utr_intervals(os.path.join(folder_out, 'utr_annotation.bed'),
                                       os.path.join(folder_out, 'pAi_gene.bed'))
print ('done [', round(time.time() - start_time, 2), 'seconds ]')

### 6. Read bioanalyzer information
print ('reading bioanalyzer profile ...', end=" ", flush=True)
start_time = time.time()
bio_size = np.array([])
bio_intensity = np.array([])
with open(os.path.join(folder_in, 'ds_012_50fix_bioanalyzer.txt'), 'r') as f:
    for line in f:
        bio_size = np.append(bio_size, int(line.split()[0]))
        bio_intensity = np.append(bio_intensity, float(line.split()[1]))
f_size, f_prob = discretize_bioanalyzer_profile(bio_size, bio_intensity, 10)
print ('done [', round(time.time() - start_time, 2), 'seconds ]')

### 7. Read bamfile
print ('reading bamfile into memory ...', end=" ", flush=True)
start_time = time.time()
bamfile = defaultdict(list)
with gzip.open(os.path.join(folder_in, 'ds_012_50fix_bamfile.txt.gz'), 'rt') as f:
    for columns in (row.strip().split() for row in f):
        gene = columns[12][8:]
        bamfile[gene].append([columns[3], columns[11], columns[18]])
print ('done [', round(time.time() - start_time, 2), 'seconds ]')

### 8. Collapsing PCR duplicates
#print ('collapsing PCR duplicates ...', end=" ", flush=True)
#start_time = time.time()
#for gene in bamfile:
#    temp_list = bamfile[gene]
#    temp_list.sort()
#    bamfile[gene] = list(temp_list for temp_list,_ in itertools.groupby(temp_list))
#print ('done [', round(time.time() - start_time, 2), 'seconds ]')

### 9. Estimate tail lengths per gene.
# focus on particular genes as examples (single 3'UTRs)
print ('setting up a tail range of', end=" ")
tail_range = tail_length_range(10, 550, 30)
for length in tail_range:
    print (length, end=" ")
print ('\n')   

### 10. Read all single UTR genes with no pAi in the UTRs
genes = []
with open(os.path.join(folder_in, 'single_utr_no_pAi_genes.txt'), 'r') as f:
    for line in f:
        genes.append(line.rstrip())

### 11. iterate over all genes and predict tails
with open (os.path.join(folder_out, 'tail_lengths.txt'), 'w') as results, open (os.path.join(folder_out, 'coverage.txt'), 'w') as cov:
    for index in range(len(genes)):
        gene = genes[index]
        print ('estimating polyA tail length for gene', gene, '...', end=" ", flush=True)
        reads = []
        for item in bamfile[gene]:
            if (int(pAi_full[gene][0]['start']) - int(item[0]) <= max(f_size)):
                reads.append(int(item[0]))
        #reads = [ reads[i] for i in sorted(random.sample(range(len(reads)), 100)) ]
        # Put threshold for number of reads required
        if len(reads) < 100:
            print ('not enough reads for analysis [', len(reads), ']')
            continue
        print (len(reads), 'reads will be used for the analysis ...', end=" ", flush=True)
        start_time = time.time()
        probs = estimate_poly_tail_length(reads, tail_range, pAi_full[gene],
                                          0, f_size, f_prob, False)
        print ('done [', round(time.time() - start_time, 2), 'seconds ]')
        results.write(gene + ',' + str(probs) + '\n')
        cov.write(gene + ',' + str(list(int(pAi_full[gene][0]['start']) - np.array(reads))) + '\n')








sys.exit()

# select the gene here
gene = 'CKS2'
reads = []
for item in bamfile[gene]:
    if (int(pAi_full[gene][0]['start']) - int(item[0]) <= max(f_size)):
        reads.append(int(item[0]))

reads = [ reads[i] for i in sorted(random.sample(range(len(reads)), 1000)) ]

print (len(reads), 'reads will be used for the analysis')

print ('estimating unweighted polyA tail length for gene', gene, '...', end=" ",
       flush=True)
probs = estimate_poly_tail_length(reads, tail_range, pAi_full[gene],
                                      0, f_size, f_prob, False)
print ('done [', round(time.time() - start_time, 2), 'seconds ]')
print ('unweighted probabilities for', gene, 'are')
print (probs)

print ('\n')
print ('estimating weigthed polyA tail length for gene', gene, '...', end=" ",
       flush=True)
probs = estimate_poly_tail_length(reads, tail_range, pAi_full[gene],
                                      0, f_size, f_prob, True)
print ('done [', round(time.time() - start_time, 2), 'seconds ]')
print ('weighted probabilities for', gene, 'are')
print (probs)
