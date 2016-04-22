#!/usr/bin/env python3


#########
# about #
#########

__version__ = "0.1.1.1"
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
from collections import defaultdict
import sys

#############
# functions #
#############

# Test if a file is gzip compressed or not
def is_gzip_file(filename):
    try:
        # This will raise OSError for uncompressed files & has no side
        # effects for compressed files:
        gzip.GzipFile(filename).peek(1)
        return True
    except OSError:
        return False

# Open a file (gzip compressed or not)
def open_file(filename):
    if (is_gzip_file(filename)):
        return gzip.open(filename,'rt')
    else:
        return open(filename,'rt')

# Should write it to a file
def identify_polyA_intervals(genome, conditions):
    return 0


# Read annotation from GTF file
# For now, this will output BED to STDOUT, but this might be changed to
# either file output or returing the GTF data as a python object
def extract_three_prime_utr_information(gtf_file,
                                        feature_utr3 = "three_prime_utr"):

    # The following parameters define the parsing of the input GTF file.
    # They were chosen according to the standard described in
    # http://genome.ucsc.edu/FAQ/FAQformat.html#format4.

    # Lines starting with the following character will be skipped:
    comment_char = '#'

    # Lines will be split into fields by the following character:
    field_separator = '\t'

    # The attributes field will be split into type/value pairs by the
    # following string:
    attributes_separator = "; "

    # The attribute with the following index will be used as gene ID:
    attribute_index_gene_id = 0

    # Attribute type/value pairs will be split by the following character:
    attribute_separator = ' '

    # The following index will be used to get the value of an attribute
    # type/value pair:
    attribute_value_index = 1

    # The following character will be removed from the beginning and end of
    # attribute values:
    attribute_value_quote = '"'

    # Read GTF input line by line
    with open_file(gtf_file) as gtf:
        for line in gtf:

            # Skip comment lines
            if (line[0] == comment_char):
                continue

            # Split lines into fields
            seqname, source, feature, start, end, score, strand, frame, attributes = line.rstrip().split(field_separator)

            # Skip lines not defining 3' UTRs
            if (feature != feature_utr3):
                continue

            # Convert from 1-based closed to 0-based open intervals
            start = (int(start) - 1)

            # Get gene ID from attributes list
            # I guess there is a more pythonic way to do that generating a
            # dictionary but this works.
            # Split accross several lines for clarity.
            attributes = attributes.split(attributes_separator)
            gene_id = attributes[attribute_index_gene_id].split(attribute_separator)
            gene_id = gene_id[attribute_value_index].strip(attribute_value_quote)

            # output BED line
            print("{seqname}\t{start}\t{end}\t{gene_id}\t{strand}\t{score}".format(**locals()))
    return 0

def merge_pAi_and_utr_intervals(pAi_bed, utr_bed):
    """Merges pAi intervals with 3'UTRs into a big dictionary, suitable
       for downstream analysis. Requires gene annotated pAi bed file."""
    genes_dict = defaultdict(list)
    pAi_full = defaultdict(list)
    
    with open("test_data/utr_annotation.bed", 'r') as f:
        current_gene = 'none'
        for line in f.readlines():
            chr, start_position, end_position, gene, strand, score = line.split('\t')
            if current_gene != gene:
                min_position = 10**15
                max_position = 1
            min_position = min(min_position, int(start_position))
            max_position = max(max_position, int(end_position))
            genes_dict[gene] = [min_position, max_position]
            current_gene = gene

            pAi_full[gene].append({'start' : start_position, 'end' : end_position,
                                   'strand' : strand, 'is_tail' : True})

    with open("test_data/pAi.bed", 'r') as f:
        for line in f.readlines():
            chr, start_position, end_position, gene, strand, score = line.split('\t')
            if int(start_position) >= genes_dict[gene][0] and int(end_position) <= genes_dict[gene][1]:
                pAi_full[gene].append({'start' : start_position, 'end' : end_position,
                                       'strand' : strand, 'is_tail' : False})
    return 0

def annotate_pAi_with_gene(pAi_bed, utr_bed):
    with open('test_data/utr_annotation.bed', 'r') as utr, open('test_data/pAi.bed', 'r') as pAi, open('pAi_gene.bed', 'w') as pAi_out:
        cur_chr, cur_start, cur_end, cur_gene, cur_strand, cur_score = utr.readline().split('\t')
        for line in utr:
            chr, start, end, gene, strand, score = line.split('\t')
            print (gene, cur_gene)
            sys.exit()
            while gene == cur_gene:
                print (gene, cur_gene)
                chr, start, end, gene, strand, score = line.split('\t')
                cur_start = min(cur_start, start)
                cur_end = max(cur_end, end)
                line = next(utr)
                print ('stuck here')
            for line2 in pAi:
                pAi_chr, pAi_start, pAi_end, pAi_gene, pAi_strand = line2.split('\t')
                if pAi_start >= cur_start and pAi_end <= cur_end:
                    print ('found overlap')
                if pAi_start > cur_end: 
                    print ('breaking for loop')
                    break
            cur_chr, cur_start, cur_end, cur_gene, cur_strand, cur_score = line.split('\t')

########
# main #
########

# Only run the following code if this module is run directly
if __name__ == '__main__':

    # GTF file to read annotation from
    gtf = "test_data/Homo_sapiens.GRCh38.83_chr9.gtf"

    # read GTF file and print BED output to STDOUT
    #extract_three_prime_utr_information(gtf)

    # annotate pAi with gene
    annotate_pAi_with_gene(1, 2)
