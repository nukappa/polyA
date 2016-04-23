#!/usr/bin/env python3

import numpy as np
from estimate_length import *
import os
import sys
import time

# Simulate a full workflow.
# Folder where bamfile, bioanalyzer profile, genome and gtf are in
folder_in = 'test_data'
gtf = os.path.join(folder_in, 'Homo_sapiens.GRCh38.83_chr9.gtf')
genome = os.path.join(folder_in, 'Homo_sapiens.GRCh38.dna.chromosome.9.fa')

# Create output directory for storing everything
folder_out = os.path.join(folder_in, 'output')
try:
    os.mkdir(folder_out)
except Exception:
    pass

### 1. Extract utr information from gtf file
print ("extracting 3'UTR information ...", end=" ", flush=True)
start_time = time.time()
old_stdout = sys.stdout
sys.stdout = open(os.path.join(folder_out, 'utr_annotation.bed'), 'w')
extract_three_prime_utr_information(gtf)
sys.stdout = old_stdout
print ('done [', round(time.time() - start_time, 2), 'seconds ]')

### 2. Extract polyA intervals from genome
print ('extracting polyA intervals from genome ...', end=" ", flush=True)
if os.path.isfile(os.path.join(folder_out, 'pAi.bed')):
    print ('skipping [ file already present ]')
else:
    start_time = time.time()
    extract_pAi_from_genome(genome, window=10, occurences=7, consecutive=6)
    os.remove('pAi_temp.bed')
    os.rename('pAi.bed', os.path.join(folder_out, 'pAi.bed'))
    print ('done [', round(time.time() - start_time, 2), 'seconds ]')

### 3. Add gene information to polyA intervals
print ('adding gene annotation to pAi intervals ...', end=" ", flush=True)
start_time = time.time()
annotate_pAi_with_gene(os.path.join(folder_out, 'pAi.bed'), 
                       os.path.join(folder_out, 'utr_annotation.bed'))
os.rename('pAi_gene.bed', os.path.join(folder_out, 'pAi_gene.bed'))
print ('done [', round(time.time() - start_time, 2), 'seconds ]')

### 4. Merge polyA intervals with 3'UTRs into a dictionary
print ("merging polyA intervals with 3'UTR ...", end=" ", flush=True)
start_time = time.time()
pAi_full = merge_pAi_and_utr_intervals(os.path.join(folder_out, 'utr_annotation.bed'),
                                       os.path.join(folder_out, 'pAi_gene.bed'))
print ('done [', round(time.time() - start_time, 2), 'seconds ]')

### 5. Read bioanalyzer information

### 6. Read bamfile

### 7. Estimate tail lengths per gene!
