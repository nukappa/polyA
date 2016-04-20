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

bio_size = []
bio_intensity = []

with open(os.path.join('test_data', 'ds_012_50fix_bioanalyzer.txt'), 'r') as f:
    for line in f:
        bio_size.append(int(line.split()[0]))
        bio_intensity.append(float(line.split()[1]))

size, probability = discretize_bioanalyzer_profile(bio_size, bio_intensity, 5)
size_np, probability_np = discretize_bioanalyzer_profile_np(np.array(bio_size), np.array(bio_intensity), 5)

class TestStringMethods(unittest.TestCase):

    def test_bioanalyzer_probabilities_summing_to_one(self):
        self.assertEqual(round(sum(probability), PRECISION), 1)

    def test_prob_d_given_L_summing_to_one(self):
        Lrange = tail_length_range(10, 250, 5)
        prob_sum = 0
        for length in Lrange:
            prob_sum += prob_d_given_L(650-97, pAi, 2, length, size, probability, Lrange)
        self.assertEqual(round(prob_sum, PRECISION), 1)

    def test_prob_read_given_pAi_summing_to_one(self):
        prob_sum = 0

        for interval in range(len(pAi)):
            prob_sum += prob_d_given_pAi(650-97, pAi, interval, size, probability)
        self.assertEqual(round(prob_sum, PRECISION), 1)

    def test_bioanalyzer_probabilities_summing_to_one_np(self):
        self.assertEqual(round(sum(probability_np), PRECISION), 1)

    def test_migration_to_np(self):
        self.assertEqual(prob_d_given_pAi(650-97, pAi, 2, size_np, probability_np), 
                         prob_d_given_pAi_np(650-97, pAi, 2, size_np, probability_np))
        self.assertEqual(prob_pAi_given_d(pAi, 2, 650-97, size_np, probability_np),
                         prob_pAi_given_d_np(pAi, 2, 650-97, size_np, probability_np))

suite = unittest.TestLoader().loadTestsFromTestCase(TestStringMethods)
unittest.TextTestRunner(verbosity=2).run(suite)
