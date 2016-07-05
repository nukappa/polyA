#!/usr/bin/env python3


#########
# about #
#########

__version__ = "0.1.3"
__author__ = ["Marcel Schilling"]
__credits__ = ["Nikolaos Karaiskos","Mireya Plass Pórtulas","Marcel Schilling","Nikolaus Rajewsky"]
__status__ = "beta"
__licence__ = "GPL"
__email__ = "marcel.schilling@mdc-berlin.de"


###########
# imports #
###########

import numpy as np


#############
# functions #
#############

def simulate_reads(genes,pAi,f_size,f_prob,reads_per_gene=100,pAlen=42):
    """Simulates reads based on fixed poly(A)-tail length distribution."""
    fragment_sizes = {}
    pAoffsets = {}
    reads = {}
    f_cum = np.cumsum(f_prob)
    for gene in genes:
        intervals = []
        for interval in pAi[gene]:
            if(interval['is_tail']):
                intervals.append(interval)
        if len(intervals) == 0:
            continue
        interval=intervals[0] # pick first 3' UTR isoform for now
        fragment_sizes[gene] = np.zeros(reads_per_gene, dtype=int)
        n_simulated = 0
        for size, prob in zip(f_size, f_prob):
            n_to_simulate = int(prob * reads_per_gene)
            fragment_sizes[gene][n_simulated:(n_simulated
                                              + n_to_simulate)] = size
            n_simulated += n_to_simulate
        for read in range(n_simulated, reads_per_gene):
            r=np.random.random()
            fragment_sizes[gene][read] = f_size[[i for i in range(len(f_cum))
                                                 if r <= f_cum[i]][0]]
        pAoffsets[gene] = np.random.random_integers(0, pAlen,
                                                    size=reads_per_gene)
        reads[gene] = (int(interval['start']) + pAoffsets[gene]
                       - fragment_sizes[gene] - 1)
    return(fragment_sizes, pAoffsets, reads)
