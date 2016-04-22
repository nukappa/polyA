import os
import numpy as np
from estimate_polyA_tail_length import *
from extract_annotation import *
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

f_size, f_prob = discretize_bioanalyzer_profile(np.array(bio_size), np.array(bio_intensity), 5)

Lrange = tail_length_range(10, 500, 500)

print ('computing length probabilities')
start_time = time.time()
print (estimate_poly_tail_length([550, 567, 568, 578, 579, 581, 600, 611], Lrange, pAi, 2, f_size, f_prob, True))
print (time.time() - start_time, 'seconds elapsed')
