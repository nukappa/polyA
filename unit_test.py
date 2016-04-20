import unittest
import os
from estimate_polyA_tail_length import *
from extract_annotation import *
import numpy as np

# defines the float precision to be tested
PRECISION = 12

pAi = [{'start' : 500, 'end' : 541, 'strand' : '+', 'is_tail' : False},                                                                                                                                                                                                               
       {'start' : 600, 'end' : 621, 'strand' : '+', 'is_tail' : False},                                                                                                                                                                                                               
       {'start' : 650, 'end' : 690, 'strand' : '+', 'is_tail' : True}]    

reads = [550, 567, 568, 578, 579, 581, 600, 611]

bio_size = []
bio_intensity = []

with open(os.path.join('test_data', 'ds_012_50fix_bioanalyzer.txt'), 'r') as f:
    for line in f:
        bio_size.append(int(line.split()[0]))
        bio_intensity.append(float(line.split()[1]))

f_size, f_prob = discretize_bioanalyzer_profile(np.array(bio_size), np.array(bio_intensity), 5)

Lrange = tail_length_range(10, 200, 25)

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

suite = unittest.TestLoader().loadTestsFromTestCase(TestStringMethods)
unittest.TextTestRunner(verbosity=2).run(suite)
