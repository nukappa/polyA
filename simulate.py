#!/usr/bin/env python3


#########
# about #
#########

__version__ = "0.1"
__author__ = ["Marcel Schilling"]
__credits__ = ["Nikolaos Karaiskos","Mireya Plass Pórtulas","Marcel Schilling","Nikolaus Rajewsky"]
__status__ = "beta"
__licence__ = "GPL"
__email__ = "marcel.schilling@mdc-berlin.de"


#############
# functions #
#############

def simulate_reads(genes,pAi,reads_per_gene=100,pAlen=42,fragment_length=400):
    """Simulates reads based on fixed poly(A)-tail length distribution."""
    reads = {}
    for gene in genes:
        intervals = []
        for interval in pAi[gene]:
            if(interval['is_tail']):
                intervals.append(interval)
        if len(intervals) == 0:
            continue
        interval=intervals[0] # pick first 3' UTR isoform for now
        reads[gene]=[]
        for read in range(reads_per_gene):
            reads[gene].append(int(interval['start'])-1+pAlen-fragment_length)
    return(reads)
