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

import gzip
import os
import time


#############
# functions #
#############

### FUNCTION BELOW WILL PROBABLY BE DEPRECATED
def compute_distances_of_reads_to_all_pAi(read_coordinate, pAi):
    """Computes the distances of the beginning of the read until the beginning
       of all given pAis. Returns pAi itself but with the distances inserted."""
    for interval in range(len(pAi)):
        pAi[interval]['distance'] = read_coordinate-pAi[interval]['start']
    return pAi

def discretize_bioanalyzer_profile(size, intensity, bin_size):
    """Discretizes a given bioanalyzer profile intensity=f(size) by putting 
       fragment sizes into bins of given bin_size. The intensities are 
       then transformed into probabilities."""
    size = [round(s/bin_size)*bin_size for s in size]
    probability = []

    for bin in range(min(size), max(size)+1, bin_size):
        indices = [i for i, x in enumerate(size) if x == bin]
        total_bin_intensity = 0

        for j in indices:
            total_bin_intensity += intensity[j]

        probability.append(total_bin_intensity/sum(intensity))

    return list(sorted(set(size))), probability

def step_function(x):
    """The 'Heaviside function'. For x=0 it returns zero, which is more
       appropriate in the current context."""
    return 1 * (x > 0)

def tail_length_range(start, end, step):
    return list(range(start, end, step))

def prob_d_given_pAi(read_coordinate, pAi, interval, f, prob_f):
    """Computes the conditional probability P(d|pAi) that a read to
       originate from the particular pAi, given a bioanalyzer profile."""
    nominator = 0
    for fragment in range(len(f)):
        nominator += (prob_f[fragment] * step_function(f[fragment] - (pAi[interval]['start']-read_coordinate)) * 
                      step_function(pAi[interval]['end'] - read_coordinate - f[fragment]))

    # normalization factor for sum(prob)=1
    norm_factor = 0
    for interval2 in range(len(pAi)):
        temp_sum = 0
        for fragment in range(len(f)):
            temp_sum += (prob_f[fragment] * step_function(f[fragment] - (pAi[interval2]['start']-read_coordinate)) *
                         step_function(pAi[interval2]['end'] - read_coordinate - f[fragment]))
        norm_factor += temp_sum

    return nominator/norm_factor

def prob_pAi_given_d(pAi, interval, read_coordinate, f, prob_f):
    """Computes the conditional probability P(pAi|d) for a pAi to give 
       rise to the read d. Prior probabilities for each pAi are taken to
       be homogeneous, namely 1/N, N=number of pAis."""
    nominator = prob_d_given_pAi(read_coordinate, pAi, interval, f, prob_f)
    denominator = 0
    for interval2 in range(len(pAi)):
        denominator += prob_d_given_pAi(read_coordinate, pAi, interval2, f, prob_f)
    return nominator/denominator

def prob_d_given_L(d, L, f, prob_f, length_range):
    """Computes the conditional probability P(d|L) given a bioanalyzer 
       profile and an interval for L."""
    nominator = 0
    for fragment in range(len(f)):
        nominator += prob_f[fragment] * 1/L * step_function(L-f[fragment]+d)
    
    # compute the norm_factor for sum(prob)=1
    norm_factor = 0
    for length in length_range:
        temp_sum = 0
        for fragment in range(len(f)):
            temp_sum += prob_f[fragment] * 1/length * step_function(length-f[fragment]+d)
        norm_factor += temp_sum
        
    return nominator/norm_factor

def assign_probabilities_to_pAi(pAi):
    """Given a set of distances of a read to all pAi as computed by 
       compute_distances_of_reads_to_all_pAi, this function assigns 
       probabilities for the read originating from each of the pAis."""
    return 0

def estimate_poly_tail_length(distances, tail, weights):
    return 0


########
# main #
########

if __name__ == '__main__':
    start_time = time.time()
    bio_size = []
    bio_intensity = []
    with open(os.path.join('test_data', 'ds_012_50fix_bioanalyzer.txt'), 'r') as f:
        for line in f:
            bio_size.append(int(line.split()[0]))
            bio_intensity.append(float(line.split()[1]))

    size, probability = discretize_bioanalyzer_profile(bio_size, bio_intensity, 5)

    Lrange = tail_length_range(10, 150, 10)
    print ('computing probabilities of read with distance 97 originating from polyA', 
           'tail with length in range(10, 150, 10)')
    for length in Lrange:
        print ('for length', length, 'probability is', prob_d_given_L(97, length, size, probability, Lrange))

    print ('computing probabilities for observing read originating from pAi')

    pAi = [{'name' : 0, 'start' : 500, 'end' : 541, 'strand' : '+', 'is_polyA' : False},
           {'name' : 1, 'start' : 600, 'end' : 621, 'strand' : '+', 'is_polyA' : False},
           {'name' : 2, 'start' : 650, 'end' : 690, 'strand' : '+', 'is_polyA' : True}]
    for interval in range(len(pAi)):
        print ('for pAi', interval, pAi[interval]['is_polyA'], 'being a polyA, probability is', 
                prob_d_given_pAi(650-97, pAi, interval, size, probability))

    print ('computing probability for pAi giving rise to read d')
    for interval in range(len(pAi)):
        print ('for pAi', interval, 'probability is', prob_pAi_given_d(pAi, interval, 650-97, size, probability))

    print (time.time() - start_time, 'seconds elapsed')
