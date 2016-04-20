import os
import numpy as np
from estimate_polyA_tail_length import *
import time

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

Lrange = tail_length_range(10, 200, 50)
Lrange_np = tail_length_range_np(10, 200, 50)

print ('comparing individual reads probabilities')
start_time = time.time()
for Length in Lrange:
    print (prob_d_given_L_weighted(650-97, pAi, 2, Length, size_np, probability_np, Lrange))
print (time.time() - start_time, 'seconds elapsed')
start_time = time.time()
for Length in Lrange_np:
    print (prob_d_given_L_weighted_np(650-97, pAi, 2, Length, size_np, probability_np, Lrange_np))
print (time.time() - start_time, 'seconds elapsed')

print ('\n' + 'comparing final function timing')
start_time = time.time()
print (estimate_poly_tail_length([550], Lrange, pAi, 2, size_np, probability_np, True))
print (time.time() - start_time, 'seconds elapsed')

start_time = time.time()
print (estimate_poly_tail_length_np([550], Lrange, pAi, 2, size_np, probability_np, True))
print (time.time() - start_time, 'seconds elapsed')
