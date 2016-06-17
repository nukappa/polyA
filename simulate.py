#!/usr/bin/env python3


#########
# about #
#########

__version__ = "0.1.2"
__author__ = ["Marcel Schilling"]
__credits__ = ["Nikolaos Karaiskos","Mireya Plass Pórtulas","Marcel Schilling","Nikolaus Rajewsky"]
__status__ = "beta"
__licence__ = "GPL"
__email__ = "marcel.schilling@mdc-berlin.de"


###########
# imports #
###########

import numpy as np
import random


#############
# functions #
#############

def simulate_reads(genes,pAi,f_size,f_prob,reads_per_gene=100,pAlen=42):
    """Simulates reads based on fixed poly(A)-tail length distribution."""
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
        reads[gene]=[]
        # list collapsing as in http://stackoverflow.com/a/952952/2451238
        for fragment_length in random.sample([item for list in [[size] * int(round(prob*reads_per_gene)) for size,prob in zip(f_size,f_prob)] for item in list], reads_per_gene):
            reads[gene].append(int(interval['start'])+random.randint(0,pAlen)-fragment_length-1)
    return(reads)
