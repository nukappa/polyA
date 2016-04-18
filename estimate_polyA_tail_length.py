#!/usr/bin/env python3


#########
# about #
#########

__version__ = "0.1.1.2"
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
    print(time.time() - start_time, 'seconds elapsed')
