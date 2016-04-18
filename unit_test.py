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

suite = unittest.TestLoader().loadTestsFromTestCase(TestStringMethods)
unittest.TextTestRunner(verbosity=2).run(suite)
