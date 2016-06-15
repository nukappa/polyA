#!/usr/bin/env python3


#########
# about #
#########

__version__ = "0.1.2"
__author__ = ["Nikolaos Karaiskos","Marcel Schilling"]
__credits__ = ["Nikolaos Karaiskos","Mireya Plass Pórtulas","Marcel Schilling","Nikolaus Rajewsky"]
__status__ = "beta"
__licence__ = "GPL"
__email__ = "marcel.schilling@mdc-berlin.de"


###########
# imports #
###########

import unittest
import os
from estimate_length import *
from simulate import *
import numpy as np
import sys
import subprocess


##############
# parameters #
##############

# defines the float precision to be tested
PRECISION = 12

pAi = [{'start' : 500, 'end' : 541, 'strand' : '+', 'is_tail' : False},
       {'start' : 600, 'end' : 621, 'strand' : '+', 'is_tail' : False},
       {'start' : 650, 'end' : 690, 'strand' : '+', 'is_tail' : True}]

reads = [550, 567, 568, 578, 579, 581, 600, 611]

Lrange = tail_length_range(10, 200, 25)


folder_in = 'test_data'
folder_out = os.path.join(folder_in, 'output')
gtf_url = 'ftp://ftp.ensembl.org/pub/current_gtf/homo_sapiens/Homo_sapiens.GRCh38.84.chr.gtf.gz'
gtf = os.path.join(folder_in, 'Homo_sapiens.GRCh38.84_chr9.gtf.gz')

f_size_sim = np.array([400])
f_prob_sim = np.array([1])

tail_range_sim = tail_length_range(40, 50, 1)


##################
# initialization #
##################

bio_size = []
bio_intensity = []


##############
# input data #
##############

with open(os.path.join('test_data', 'ds_012_50fix_bioanalyzer.txt'), 'r') as f:
    for line in f:
        bio_size.append(int(line.split()[0]))
        bio_intensity.append(float(line.split()[1]))

f_size, f_prob = discretize_bioanalyzer_profile(np.array(bio_size), np.array(bio_intensity), 5)

# taken from pipeline.py:
try:
    os.mkdir(folder_out)
except Exception:
    pass


# generate annotation BED (if not existing):
if not(os.path.isfile(os.path.join(folder_out, 'utr_annotation.bed'))):

    # generate GTF (if not existing):
    if not(os.path.isfile(gtf)):
        subprocess.call('wget ' + gtf_url + ' -O - | zcat | grep "^9\t" | gzip --best > ' + gtf, shell=True)

    # taken from pipeline.py:
    old_stdout = sys.stdout
    sys.stdout = open(os.path.join(folder_out, 'utr_annotation_temp.bed'), 'w')
    extract_three_prime_utr_information(gtf, bed_name_attributes = ["gene_name"])
    sys.stdout = old_stdout

    ### 1.1 Clean utr from haplotypes and junk chromosomes
    with open(os.path.join(folder_out, 'utr_annotation_temp.bed'), 'r') as fin, open(os.path.join(folder_out, 'utr_annotation_unsorted.bed'), 'w') as fout:
        for line in fin:
            if line.startswith('chrGL') or line.startswith('chrKI'):
                continue
            else:
                fout.write(line)
    os.remove(os.path.join(folder_out, 'utr_annotation_temp.bed'))

    ### 1.2 Sort the utr file alphabetically
    subprocess.call('sort -V test_data/output/utr_annotation_unsorted.bed > test_data/output/utr_annotation.bed', shell=True)
    os.remove(os.path.join(folder_out, 'utr_annotation_unsorted.bed'))

pAi_sim = defaultdict(list)
with open(os.path.join(folder_out, 'utr_annotation.bed'), 'r') as f:
    for line in f.readlines():
        chr, start_position, end_position, gene, strand, score = line.split('\t')
        pAi_sim[gene].append({'start' : end_position, 'end' : 0,
                              'strand' : strand, 'is_tail' : True})

genes = []
with open(os.path.join(folder_in, 'single_utr_no_pAi_genes.txt'), 'r') as f:
    for line in f:
        genes.append(line.rstrip())


#################
# simulate data #
#################


    reads_sim=simulate_reads(genes,pAi_sim)


##########################
# analyze simulated data #
##########################

    est_pAlen = {}
    for gene in dict.keys(reads_sim):
        if len(reads_sim[gene]) < 100:
            continue
        probs = estimate_poly_tail_length(reads_sim[gene], tail_range_sim,
                                          pAi_sim[gene], 0, f_size_sim,
                                          f_prob_sim, False)
        est_pAlen[gene]=int(round(sum(tail_range_sim*probs))) # expected value


#########
# tests #
#########

class TestStringMethods(unittest.TestCase):

    def test_bioanalyzer_probabilities_summing_to_one(self):
        self.assertEqual(round(sum(f_prob), PRECISION), 1)

    def test_prob_d_given_L_summing_to_one(self):
        Lrange = tail_length_range(10, 250, 5)
        prob_sum = 0
        for length in Lrange:
            prob_sum += prob_d_given_L(650-97, pAi, 2, length, f_size, f_prob, Lrange)
        self.assertEqual(round(prob_sum, PRECISION), 1)

    def test_prob_read_given_pAi_summing_to_one(self):
        prob_sum = 0
        for interval in range(len(pAi)):
            prob_sum += prob_d_given_pAi(650-97, pAi, interval, f_size, f_prob)
        self.assertEqual(round(prob_sum, PRECISION), 1)

    def test_prob_pAi_given_read_summing_to_one(self):
        prob_sum = 0
        for interval in range(len(pAi)):
            prob_sum += prob_pAi_given_d(pAi, interval, 650-97, f_size, f_prob)
        self.assertEqual(round(prob_sum, PRECISION), 1)

    def test_estimate_poly_tail_length_probs_summing_to_one(self):
        self.assertEqual(sum(estimate_poly_tail_length(reads, Lrange, pAi, 2, f_size, f_prob, True)), 1) 


    def test_simulated_data_resulting_in_expected_value_pAlen_correct(self):
        self.assertTrue(all(est_pAlen[gene]==42 for gene in est_pAlen))

#######
# run #
#######

suite = unittest.TestLoader().loadTestsFromTestCase(TestStringMethods)
unittest.TextTestRunner(verbosity=2).run(suite)
