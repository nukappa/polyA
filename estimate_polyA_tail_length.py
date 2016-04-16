import os
import time

def compute_distances_of_reads_to_all_pAi(bamfil):
    return 0

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

def assign_probabilities_to_pAi(distances, pAi):
    return 0

def estimate_poly_tail_length(distances, tail, weights):
    return 0


if __name__ == '__main__':
    start_time = time.time()
    bio_size = []
    bio_intensity = []

    with open(os.path.join('test_data', 'ds_012_50fix_bioanalyzer.txt'), 'r') as f:
        for line in f:
            bio_size.append(int(line.split()[0]))
            bio_intensity.append(float(line.split()[1]))

    print(discretize_bioanalyzer_profile(bio_size, bio_intensity, 5))
    print(time.time() - start_time)
