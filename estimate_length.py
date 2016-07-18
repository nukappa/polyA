#!/usr/bin/env python3


#########
# about #
#########

__version__ = "0.1.5.1"
__author__ = ["Nikolaos Karaiskos","Marcel Schilling"]
__credits__ = ["Nikolaos Karaiskos","Mireya Plass Pórtulas","Marcel Schilling","Nikolaus Rajewsky"]
__status__ = "beta"
__licence__ = "GPL"
__email__ = "marcel.schilling@mdc-berlin.de"


###########
# imports #
###########

import gzip
import numpy as np
from collections import defaultdict
import time
import decimal
from scipy.interpolate import interp1d


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


# Read annotation from GTF file
# For now, this will output BED to STDOUT, but this might be changed to
# either file output or returing the GTF data as a python object
def extract_three_prime_utr_information(gtf_file,
                                        bed_name_attributes = ["gene_id",
                                                               "gene_name"],
                                        bed_name_separator = "|",
                                        feature_utr3 = "three_prime_utr",
                                        feature_gene = "gene",
                                        feature_transcript = "transcript",
                                        feature_exon = "exon"):

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

    # Attribute type/value pairs will be split by the following
    # character:
    attribute_separator = ' "'

    # The following index will be used to get the value of an attribute
    # type/value pair:
    attribute_value_index = 1

    # The following character will be removed from the beginning and end
    # of attribute values:
    attribute_value_quote = '"'

    # The following character will be used to separte fields in the BED
    # output:
    bed_separator = '\t'

    # This set will be used to store all 3' UTR BED entries for the
    # current gene.
    three_prime_utrs = set()

    # Read GTF input line by line
    with open_file(gtf_file) as gtf:
        for line in gtf:

            # Skip comment lines
            if (line[0] == comment_char):
                continue

            # Split lines into fields
            (seqname, source, feature, start, end, score, strand, frame,
                attributes) = line.rstrip().split(field_separator)

            # Output BED line for each (different) 3' UTR isoform of the
            # previous gene & re-initialize 3' UTR set for current gene
            if (feature == feature_gene):
                for three_prime_utr in three_prime_utrs:
                    print(bed_separator.join(str(field) for field in
                                             three_prime_utr))
                three_prime_utrs = set()
                continue

            # Re-initialize extension length for new transcripts
            if (feature == feature_transcript):
                extension_length=0
                continue

            # Convert from 1-based closed to 0-based open intervals
            start = (int(start) - 1)
            end = int(end)

            # Store coordinates of last exon of the current transcript
            if (feature == feature_exon):
                exon=dict(start = start, end = end)
                continue

            # Skip lines not defining 3' UTRs
            if (feature != feature_utr3):
                continue

            # Split attributes into type/value pairs
            attributes = attributes.split(attributes_separator)
            attributes = [attribute.rstrip(attributes_separator) for attribute
                          in attributes]
            attributes = [attribute.split(attribute_separator) for attribute in
                          attributes]
            attributes = dict((key, value[:-1]) for (key, value) in attributes)

            # Construct BED name field from specified GTF attributes
            gene=bed_name_separator.join(attributes[attribute] for attribute in
                                         bed_name_attributes)

            # Count 3' UTR nucleotides in upstream exons
            if (strand == "+" and exon["end"] != end) or \
               (strand == "-" and exon["start"] != start):
                extension_length -= end - start

            # Determine number of extra nucleotides in last exon
            # compared to 3' UTR only
            else:
                if (strand == "+"):
                    extension_length += start - exon["start"]
                else:
                    extension_length += exon["end"] - end

            # Append current 3' UTR to 3' UTRs of current gene (if not
            # seen before) replacing the coordinates by those of the
            # last exon and the score by the number of nucleotides added
            # to the 3' UTR by this extension (negative for spliced 3'
            # UTRs)
            three_prime_utrs.add((seqname, exon["start"], exon["end"], gene,
                                  strand, extension_length))

    # Output BED line for each (different) 3' UTR isoform of the last
    # gene
    for three_prime_utr in three_prime_utrs:
        print(bed_separator.join(str(field) for field in three_prime_utr))


def merge_pAi_and_utr_intervals(utr_bed, pAi_bed):
    """Merges pAi intervals with 3'UTRs into a big dictionary, suitable
       for downstream analysis. Requires gene annotated pAi bed file."""
    pAi_full = defaultdict(list)

    with open(utr_bed, 'r') as f:
        for line in f.readlines():
            chr, start_position, end_position, gene, strand, score = line.split('\t')
            pAi_full[gene].append({'start' : end_position, 'end' : 0,
                                   'strand' : strand, 'is_tail' : True})
    with open(pAi_bed, 'r') as f:
        for line in f.readlines():
            chr, start_position, end_position, gene, strand = line.split('\t')
            pAi_full[gene.strip()].append({'start' : start_position, 'end' : end_position,
                                   'strand' : strand.strip(' \n'), 'is_tail' : False})
    return pAi_full

def extract_pAi_from_genome(genome, window, occurences, consecutive):
    with open(genome, 'r') as f, open('pAi_temp.bed' ,'w') as pAi:
        lines = (line.rstrip('\n') for line in f)
        for line in lines:
            if '>' in line:
                chromosome = str(line.split()[0][1:])
                genomic_coordinate = 0
                prefix = ''
                continue
            line = prefix + line
            c = 0
            while c <= len(line)-window:
                segment = line[c:(c+window)]
                if consecutive*'A' in segment:
                    pAi.write('%s\t%i\t%i\t%s\t%s\n' %(chromosome, genomic_coordinate,
                                                       genomic_coordinate+window, '.', '+'))
                    c += 1
                    genomic_coordinate += 1
                    continue
                if segment.count('A') >= occurences:
                    pAi.write('%s\t%i\t%i\t%s\t%s\n' %(chromosome, genomic_coordinate,
                                                       genomic_coordinate+window, '.', '+'))
                    c += 1
                    genomic_coordinate += 1
                    continue
                if consecutive*'T' in segment:
                    pAi.write('%s\t%i\t%i\t%s\t%s\n' %(chromosome, genomic_coordinate,
                                                       genomic_coordinate+window, '.', '-'))
                    c += 1
                    genomic_coordinate += 1
                    continue
                if segment.count('T') >= occurences:
                    pAi.write('%s\t%i\t%i\t%s\t%s\n' %(chromosome, genomic_coordinate,
                                                       genomic_coordinate+window, '.', '-'))
                    c += 1
                    genomic_coordinate += 1
                    continue
                c += 1
                genomic_coordinate += 1
            prefix = line[c:]

    with open('pAi_temp.bed', 'r') as fi, open('pAi.bed', 'w') as fo:
        previous_line = fi.readline().split('\t')
        for line in fi:
            line = line.split('\t')
            if (int(line[1]) > int(previous_line[1]) and int(line[1]) <= int(previous_line[2]) 
                and line[4] == previous_line[4]):
                previous_line[2] = line[2]
            else:
                fo.write('\t'.join(previous_line))
                previous_line = line

def annotate_pAi_with_gene(pAi_bed, utr_bed):
    """Annotates the pAi bed file with gene names as they appear in the utr
       annotation file. It needs improvement since for every gene it reads 
       the pAi bed file again from the beginning. Stems from the fact that 
       utr bed file is not sorted."""
    with open(utr_bed, 'r') as utr, open(pAi_bed, 'r') as pAi, open('pAi_gene.bed', 'w') as pAi_out:
        cur_chr, cur_start, cur_end, cur_gene, cur_strand, cur_score = utr.readline().split('\t')
        for line in utr:
            chr, start, end, gene, strand, score = line.split('\t')
            if gene == cur_gene:
                cur_start = min(cur_start, start)
                cur_end = max(cur_end, end)
                continue
            for line2 in pAi:
                pAi_chr, pAi_start, pAi_end, pAi_gene, pAi_strand = line2.split('\t')
                if pAi_chr != cur_chr:
                    continue
                if int(pAi_start) >= int(cur_start) and int(pAi_end) <= int(cur_end):
                    if pAi_strand.strip(' \n') == str(cur_strand):
                        pAi_out.write('%s\t%i\t%i\t%s\t%s' %(pAi_chr, int(pAi_start), 
                                                             int(pAi_end), cur_gene, pAi_strand))
                if int(pAi_start) > int(cur_end): 
                    break
            cur_chr, cur_start, cur_end, cur_gene, cur_strand, cur_score = line.split('\t')

### Will be deprecated in the future. Interpolate from scipy performs much better.
def discretize_bioanalyzer_profile_old(size, intensity, bin_size):
    """Discretizes a given bioanalyzer profile intensity=f(size) by putting 
       fragment sizes into bins of given bin_size. The intensities are 
       then transformed into probabilities."""
    bins = np.arange(min(size), max(size)+1, bin_size)
    size = np.digitize(size, bins) * bin_size + min(size) - bin_size
    probability = np.array([sum(intensity[size == x])/sum(intensity) for x in np.unique(size)])
    return np.unique(size), probability

def discretize_bioanalyzer_profile(size, intensity, bin_size):
    f = interp1d(size, intensity)
    new_size = np.linspace(min(size), max(size), 
                           num=round(max(size-min(size))/bin_size))
    new_size = np.round(new_size).astype(int)
    probability = f(new_size)/sum(f(new_size))
    return new_size, probability    

def step_function(x):
    """The 'Heaviside function'. For x=0 it returns zero, which is more
       appropriate in the current context."""
    return 1 * (x > 0)

def tail_length_range(start, end, step):
    return np.arange(start, end, step)

def prob_d_given_pAi(read_coordinate, pAi, interval, f, prob_f):
    """Computes the conditional probability P(d|pAi) that a read to
       originate from the particular pAi, given a bioanalyzer profile."""
    nominator = sum(prob_f * step_function(f - int(pAi[interval]['start']) + read_coordinate) *
                    step_function(int(pAi[interval]['end']) - read_coordinate - f) * 
                    1/(int(pAi[interval]['end']) - int(pAi[interval]['start'])))
    
    # normalization factor for sum(prob)=1
    norm_factor = sum([sum(prob_f * step_function(f - int(pAi[i]['start']) + read_coordinate) *
                       step_function(int(pAi[i]['end']) - read_coordinate - f) * 
                       1/(int(pAi[i]['end']) - int(pAi[i]['start']))) for i in range(len(pAi))])
    if norm_factor == 0:  
         return 0
    else:
        return nominator/norm_factor

def prob_pAi_given_d(pAi, interval, read_coordinate, f, prob_f):
    """Computes the conditional probability P(pAi|d) for a pAi to give 
       rise to the read d. Prior probabilities for each pAi are taken to
       be homogeneous, namely 1/N, N=number of pAis."""
    nominator = prob_d_given_pAi(read_coordinate, pAi, interval, f, prob_f)
    denominator = sum([prob_d_given_pAi(read_coordinate, pAi, intrv, f, prob_f) for 
                       intrv in range(len(pAi))])
    return nominator/denominator

def prob_d_given_L(read_coordinate, pAi, interval, Length, f, prob_f, length_range):
    """Computes the conditional probability P(d|L) given the genomic coordinate
       of the read, a set of pAis, which of the pAis is the polyA tail, a length
       value, a bioanalyzer and a range for L."""
    nominator = sum(prob_f * 1/Length * step_function(int(pAi[interval]['start']) +
                                                      Length - read_coordinate - f))
    norm_factor = sum([sum(prob_f * 1/length *
                       step_function(int(pAi[interval]['start']) + length - 
                                     read_coordinate - f)) for length in length_range])
    return nominator/norm_factor

def prob_d_given_L_weighted(read_coordinate, pAi, interval, Length, f, prob_f,
                            length_range):
    """Computes the conditional probability P(d|L) given the genomic coordinate
       of the read, a set of pAis, which of the pAis is the polyA tail, a length
       value, a bioanalyzer and a range for L."""
    pAi[interval]['end'] = int(pAi[interval]['start']) + Length
    nominator = sum(prob_f * 1/Length * step_function(int(pAi[interval]['end']) - 
                                                      read_coordinate - f) *
                    prob_d_given_pAi(read_coordinate, pAi, interval, f, prob_f))

    # compute the norm_factor for sum(prob)=1
    norm_factor = 0
    for length in length_range:
        pAi[interval]['end'] = int(pAi[interval]['start']) + length
        norm_factor += sum(prob_f * 1/length * step_function(int(pAi[interval]['end']) 
                                                             - read_coordinate - f) *
                           prob_d_given_pAi(read_coordinate, pAi, interval, f, prob_f))
    return nominator/norm_factor


def estimate_poly_tail_length(reads, tail_range, pAi, interval, f, prob_f,
                              weighted):
    """Takes a set of reads (list of read_coordinates), a range of polyA tail
       lengths, a set of internal priming intervals and a bioanalyzer profile.
       Homogeneous prior probabilities for the p(L) are assumed."""
    L_probs = []
    nominator = np.zeros(len(tail_range))
    possible = np.ones(len(tail_range), dtype=bool)
    read_probs = np.zeros(len(tail_range))
    for read in reads:
        for index, L in [(index, L) for index, L, consider
                      in zip(range(len(tail_range)), tail_range, possible)
                      if consider]:
            if weighted:
                read_probs[index] = prob_d_given_L_weighted(read, pAi,
                                                            interval, L, f,
                                                            prob_f,
                                                            tail_range)
            else:
                read_probs[index] = prob_d_given_L(read, pAi, interval, L, f,
                                                   prob_f, tail_range)
            possible[index] = read_probs[index] > 0
        nominator[possible] += np.log(read_probs[possible])
    nominator = [decimal.Decimal(value).exp()
                 if nonzero else decimal.Decimal(0)
                 for value, nonzero in zip(nominator, possible)]
    norm_factor = sum(nominator)
    return [float(value/norm_factor) for value in nominator]


########
# main #
########

# Only run the following code if this module is run directly
if __name__ == '__main__':
    print (0)
