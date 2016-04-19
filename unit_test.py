import unittest
import os
from estimate_polyA_tail_length import *
from extract_annotation import *

# defines the float precision to be tested
PRECISION = 12

class TestStringMethods(unittest.TestCase):

    def test_bioanalyzer_probabilities_summing_to_one(self):
        bio_size = []
        bio_intensity = []
        
        with open(os.path.join('test_data', 'ds_012_50fix_bioanalyzer.txt'), 'r') as f:
            for line in f:
                bio_size.append(int(line.split()[0]))
                bio_intensity.append(float(line.split()[1]))

        size, probability = discretize_bioanalyzer_profile(bio_size, bio_intensity, 5)
        self.assertEqual(round(sum(probability), PRECISION), 1)

    def test_prob_d_given_L_summing_to_one(self):
        bio_size = []
        bio_intensity = []

        with open(os.path.join('test_data', 'ds_012_50fix_bioanalyzer.txt'), 'r') as f:
            for line in f:
                bio_size.append(int(line.split()[0]))
                bio_intensity.append(float(line.split()[1]))
        
        size, probability = discretize_bioanalyzer_profile(bio_size, bio_intensity, 5)
        Lrange = tail_length_range(10, 250, 5)
        prob_sum = 0
        for length in Lrange:
            prob_sum += prob_d_given_L(237, length, size, probability, Lrange)
        self.assertEqual(round(prob_sum, PRECISION), 1)

    def test_prob_read_given_pAi_summing_to_one(self):
        bio_size = []
        bio_intensity = []

        with open(os.path.join('test_data', 'ds_012_50fix_bioanalyzer.txt'), 'r') as f:
            for line in f:
                bio_size.append(int(line.split()[0]))
                bio_intensity.append(float(line.split()[1]))

        size, probability = discretize_bioanalyzer_profile(bio_size, bio_intensity, 5)
        prob_sum = 0
        pAi = [{'name' : 0, 'start' : 500, 'end' : 541, 'strand' : '+', 'is_polyA' : False},
               {'name' : 1, 'start' : 600, 'end' : 621, 'strand' : '+', 'is_polyA' : False},
               {'name' : 2, 'start' : 650, 'end' : 690, 'strand' : '+', 'is_polyA' : True}]

        for interval in range(len(pAi)):
            prob_sum += prob_d_given_pAi(650-97, pAi, interval, size, probability)
        self.assertEqual(round(prob_sum, PRECISION), 1)

suite = unittest.TestLoader().loadTestsFromTestCase(TestStringMethods)
unittest.TextTestRunner(verbosity=2).run(suite)
