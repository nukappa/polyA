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


################################
# seed random number generator #
################################

# Ensure independence of the test results from random choices
# According to http://stackoverflow.com/questions/11526975/set-random-seed-programwide-in-python#comment15239525_11527011
# the random seed has to be set before any other import that imports the
# random module to avoid seeding with system time.
import numpy as np
np.random.seed(42)


###########
# imports #
###########

import unittest
import os
from estimate_length import *
from simulate import *
import sys
import subprocess
from scipy.stats import power_divergence
from scipy.stats import pearsonr


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

f_size_sim = np.array([399, 400])
f_prob_sim = np.array([.1, .9])
reads_per_gene = 450
pAlen_sim = 42

tail_range_sim = tail_length_range(40, 50, 1)

# Power to use in the Cressie-Read power divergence statistic
# This allows to easily switch between different test statistics when
# comparing a sample to a potentially multinomial discrete distribution.
# See http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.stats.power_divergence.html#scipy.stats.power_divergence
# for more details.
# The following quote from the documentation cited above explains common
# choices of lambda:
# > "pearson"             1     Pearson's chi-squared statistic.
# >                             In this case, the function is
# >                             equivalent to `stats.chisquare`.
# > "log-likelihood"      0     Log-likelihood ratio. Also known as
# >                             the G-test [R256]_.
# > "freeman-tukey"      -1/2   Freeman-Tukey statistic.
# > "mod-log-likelihood" -1     Modified log-likelihood ratio.
# > "neyman"             -2     Neyman's statistic.
# > "cressie-read"        2/3   The power recommended in [R258]_.
power_divergence_lambda = 1

# p-value threshold to declare distributions equal based on power
# deviation test
alpha_distcomp = .05

# Pearson's R threshold to declare distributions equal based on correlation
r_threshold = .9

# p-value threshold to declare distributions equal based on correlation
alpha_cor = alpha_distcomp


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

fragment_sizes_sim, pAoffsets_sim, reads_sim = simulate_reads(genes, pAi_sim,
                                                              f_size_sim,
                                                              f_prob_sim,
                                                              reads_per_gene,
                                                              pAlen_sim)


##########################
# analyze simulated data #
##########################

probs_estimated = {}
for gene in dict.keys(reads_sim):
    if len(reads_sim[gene]) < 100:
        continue
    probs_estimated[gene] = estimate_poly_tail_length(reads_sim[gene],
                                                      tail_range_sim,
                                                      pAi_sim[gene], 0,
                                                      f_size_sim, f_prob_sim,
                                                      False)


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

    def test_number_of_simulated_reads_correct(self):
        self.assertTrue(all(len(reads_sim[gene]) == reads_per_gene for gene in probs_estimated))

    def test_singleUTR_no_pAi_genes_have_UTR(self):
        self.assertTrue(all([len([interval for interval in pAi_sim[gene] if interval['is_tail']])!=0 for gene in genes]))

    def test_singleUTR_no_pAi_genes_dont_have_multiple_UTR(self):
        self.assertTrue(all([len([interval for interval in pAi_sim[gene] if interval['is_tail']])<2 for gene in genes]))

    def test_singleUTR_no_pAi_genes_have_singleUTR(self):
        self.assertTrue(all([len([interval for interval in pAi_sim[gene] if interval['is_tail']])==1 for gene in genes]))

    def test_singleUTR_no_pAi_genes_have_no_pAi(self):
        self.assertEqual(max([len([interval for interval in pAi_sim[gene] if not interval['is_tail']]) for gene in genes]),0)

    def test_simulated_fragment_size_distributions_match_that_to_simulate(self):
        for gene in genes:
            f_min = min(np.append(f_size_sim, fragment_sizes_sim[gene]))
            f_max = max(np.append(f_size_sim, fragment_sizes_sim[gene]))
            probs_to_simulate = np.zeros(f_max - f_min + 1)
            n_simulated = np.zeros(f_max - f_min + 1)
            for f in range(f_min, f_max + 1):
                probs_to_simulate[f - f_min] = ([prob for prob,size
                                                 in zip(f_prob_sim, f_size_sim)
                                                 if size == f] + [0])[0]
                n_simulated[f - f_min] = sum(fragment_sizes_sim[gene] == f)

            # Calculate p-value for the two distributions being different
            p_divergence = power_divergence(n_simulated,
                                            reads_per_gene * probs_to_simulate,
                                            lambda_=power_divergence_lambda)[1]

            # Invert p-value to test if distributions are identical
            # Note that this is statistically not correct as not being able
            # to reject the null hypothesis does not generally prove it but
            # this seems like the best approach possible (plus: It is commonly
            # used in normality test).
            self.assertTrue(all(probs_to_simulate == n_simulated
                                                     / reads_per_gene)
                            or ((1 - p_divergence) <= alpha_distcomp))

    def test_simulated_data_resulting_in_expected_pAlen_distribution(self):
        probs_simulated = np.array([int(length == pAlen_sim)
                                    for length in tail_range_sim])
        for gene in genes:
            r, p_val = pearsonr(probs_estimated[gene], probs_simulated)
            self.assertTrue(all(probs_estimated[gene] == probs_simulated) or
                            (r >= r_threshold and p_val <= alpha_cor))


#######
# run #
#######

suite = unittest.TestLoader().loadTestsFromTestCase(TestStringMethods)

# exit 0 only if all tests pass (see http://stackoverflow.com/a/24972157/2451238)
sys.exit(not unittest.TextTestRunner(verbosity=2).run(suite).wasSuccessful())
