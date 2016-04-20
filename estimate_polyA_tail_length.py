#!/usr/bin/env python3


#########
# about #
#########

__version__ = "0.1.3"
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
import numpy as np
import time


#############
# functions #
#############

def discretize_bioanalyzer_profile(size, intensity, bin_size):
    """Discretizes a given bioanalyzer profile intensity=f(size) by putting 
       fragment sizes into bins of given bin_size. The intensities are 
       then transformed into probabilities."""
    bins = np.arange(min(size), max(size)+1, bin_size)
    size = np.digitize(size, bins) * bin_size + min(size) - bin_size
    probability = np.array([sum(intensity[size == x])/sum(intensity) for x in np.unique(size)])
    return np.unique(size), probability

def step_function(x):
    """The 'Heaviside function'. For x=0 it returns zero, which is more
       appropriate in the current context."""
    return 1 * (x > 0)

def tail_length_range(start, end, step):
    return np.arange(start, end, step)

def prob_d_given_pAi(read_coordinate, pAi, interval, f, prob_f):
    """Computes the conditional probability P(d|pAi) that a read to
       originate from the particular pAi, given a bioanalyzer profile."""
    nominator = sum(prob_f * step_function(f - pAi[interval]['start'] + read_coordinate) * 
                    step_function(pAi[interval]['end'] - read_coordinate - f))
    # normalization factor for sum(prob)=1
    norm_factor = sum([sum(prob_f * step_function(f - pAi[i]['start'] + read_coordinate) * 
                       step_function(pAi[i]['end'] - read_coordinate - f)) for i in range(len(pAi))])
    return nominator/norm_factor

def prob_pAi_given_d(pAi, interval, read_coordinate, f, prob_f):
    """Computes the conditional probability P(pAi|d) for a pAi to give 
       rise to the read d. Prior probabilities for each pAi are taken to
       be homogeneous, namely 1/N, N=number of pAis."""
    nominator = prob_d_given_pAi(read_coordinate, pAi, interval, f, prob_f)
    denominator = sum([prob_d_given_pAi(read_coordinate, pAi, intrv, f, prob_f) for intrv in range(len(pAi))])
    return nominator/denominator

def prob_d_given_L(read_coordinate, pAi, interval, Length, f, prob_f, length_range):
    """Computes the conditional probability P(d|L) given the genomic coordinate
       of the read, a set of pAis, which of the pAis is the polyA tail, a length
       value, a bioanalyzer and a range for L."""
    nominator = sum(prob_f * 1/Length * step_function(pAi[interval]['start'] + Length - read_coordinate - f))
    norm_factor = sum([sum(prob_f * 1/length * 
                       step_function(pAi[interval]['start'] + length - read_coordinate - f)) for length in length_range])
    return nominator/norm_factor

def prob_d_given_L_weighted(read_coordinate, pAi, interval, Length, f, prob_f, length_range):
    """Computes the conditional probability P(d|L) given the genomic coordinate
       of the read, a set of pAis, which of the pAis is the polyA tail, a length
       value, a bioanalyzer and a range for L."""
    pAi[interval]['end'] = pAi[interval]['start'] + Length
    nominator = sum(prob_f * 1/Length * step_function(pAi[interval]['end'] - read_coordinate - f) *
                    prob_d_given_pAi(read_coordinate, pAi, interval, f, prob_f))

    # compute the norm_factor for sum(prob)=1
    norm_factor = 0
    for length in length_range:
        pAi[interval]['end'] = pAi[interval]['start'] + length
        norm_factor += sum(prob_f * 1/length * step_function(pAi[interval]['end'] - read_coordinate - f) * 
                           prob_d_given_pAi(read_coordinate, pAi, interval, f, prob_f))
    return nominator/norm_factor

def estimate_poly_tail_length(reads, tail_range, pAi, interval, f, prob_f, weighted):
    """Takes a set of reads (list of read_coordinates), a range of polyA tail
       lengths, a set of internal priming intervals and a bioanalyzer profile.
       Homogeneous prior probabilities for the p(L) are assumed."""
    L_probs = []
    nominator = np.ones(len(tail_range))
    for read in reads:
        read_probs = []
        for L in tail_range:
            if weighted:
                read_probs.append(prob_d_given_L_weighted(read, pAi, interval, L, f, prob_f, tail_range))
            else:
                read_probs.append(prob_d_given_L(read, pAi, interval, L, f, prob_f, tail_range))
        nominator *= read_probs
    return nominator/sum(nominator)


########
# main #
########

if __name__ == '__main__':
    start_time = time.time()
    print ('\n',  time.time() - start_time, 'seconds elapsed')
