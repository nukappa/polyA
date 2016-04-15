#!/usr/bin/env python3


#########
# about #
#########

__version__ = "0.1"
__author__ = "Marcel Schilling"
__credits__ = ["Nikolaos Karaiskos","Mireya Plass Pórtulas","Marcel Schilling","Nikolaus Rajewsky"]
__status__ = "beta"
__licence__ = "GPL"
__email__ = "marcel.schilling@mdc-berlin.de"


###########
# imports #
###########

import gzip


##############
# parameters #
##############

# input parameters #

# GTF file to read annotation from.
# For now, use fixed GTF path for testing.
# Note that thanks to Andrei & Marta (but not to Python) this coded can
# deal with both, gzip compressed & uncompressed input files,
# regardless of the file name.
gtf_file = "test_data/Homo_sapiens.GRCh38.83_chr9.gtf"


# file reading parameters #

# The following string will be used to identify if the input GTF was
# gzip compressed or not:
#gzip_magic_string = '\x1f\x8b'

# GTF parsing parameters #

# The following parameters define the parsing of the input GTF file.
# They were chosen according to the standard described in
# http://genome.ucsc.edu/FAQ/FAQformat.html#format4.

# Lines starting with the following character will be skipped:
comment_char = '#'

# Lines will be split into fields by the following character:
field_separator = '\t'

# Lines defining 3' UTRs will be identified by the following feature
# name:
feature_utr3 = "three_prime_utr"

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


#############
# functions #
#############

# open a file (gzip compressed or not)
def open_file(filename):
    try:
        # This will raise OSError for uncompressed files & has no side
        # effects for compressed files:
        gzip.GzipFile(filename).peek(1)
        return gzip.open(filename,'rt')
    except OSError:
        return open(filename,'rt')


########
# main #
########

# Open GTF file to read & iterate over lines
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
