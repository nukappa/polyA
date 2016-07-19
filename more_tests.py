import os
import numpy as np
from estimate_length import *
import time
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

pAi = [{'start' : 650, 'end' : 0, 'strand' : '+', 'is_tail' : True}]

bio_size = []
bio_intensity = []

with open(os.path.join('test_data', 'ds_012_50fix_bioanalyzer.txt'), 'r') as f:
    for line in f:
        bio_size.append(int(line.split()[0]))
        bio_intensity.append(float(line.split()[1]))


f_size, f_prob = discretize_bioanalyzer_profile(np.array(bio_size), np.array(bio_intensity), 10)

Lrange = tail_length_range(10, 500, 20)

read = 450
for length in Lrange:
    pAi[0]['end'] = pAi[0]['start'] + length
    print (prob_d_given_pAi(read, pAi, 0, f_size, f_prob))

sys.exit()

print ('computing length probabilities')
start_time = time.time()
print (estimate_poly_tail_length([550, 567, 568, 578, 579, 581, 600, 611], Lrange, pAi, 0, f_size, f_prob, True))
print (time.time() - start_time, 'seconds elapsed')
