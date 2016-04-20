import os
import numpy as np
from estimate_polyA_tail_length import *

pAi = [{'start' : 500, 'end' : 541, 'strand' : '+', 'is_tail' : False},                                                                                                                                                                                                               
       {'start' : 600, 'end' : 621, 'strand' : '+', 'is_tail' : False},                                                                                                                                                                                                               
       {'start' : 650, 'end' : 690, 'strand' : '+', 'is_tail' : True}]    

bio_size = []
bio_intensity = []

with open(os.path.join('test_data', 'ds_012_50fix_bioanalyzer.txt'), 'r') as f:
    for line in f:
        bio_size.append(int(line.split()[0]))
        bio_intensity.append(float(line.split()[1]))

size_np, probability_np = discretize_bioanalyzer_profile_np(np.array(bio_size), np.array(bio_intensity), 5)

print (prob_pAi_given_d(pAi, 2, 650-97, size_np, probability_np) == prob_pAi_given_d_np(pAi, 2, 650-97, size_np, probability_np))
