import time

def extract_pAi_from_genome(genome, window, occurences, consecutive):
    with open('test_data/test.fa', 'r') as f, open('pAi_temp.bed' ,'w') as pAi:
        lines = (line.rstrip('\n') for line in f)
        for line in lines:
            if '>' in line:
                chromosome = 'chr' + line.split()[0][1:]
                genomic_coordinate = 0
                prefix = ''
                continue
            line = prefix + line
            c = 0
            while c <= len(line)-window:
                segment = line[c:(c+window)]
                if consecutive*'A' in segment:
                    pAi.write('%s\t %i\t %i\t %s\n' %(chromosome, genomic_coordinate, genomic_coordinate+window, '+'))
                    c += 1
                    genomic_coordinate += 1
                    continue
                if segment.count('A') >= occurences:
                    pAi.write('%s\t %i\t %i\t %s\n' %(chromosome, genomic_coordinate, genomic_coordinate+window, '+'))
                    c += 1
                    genomic_coordinate += 1
                    continue
                if consecutive*'T' in segment:
                    pAi.write('%s\t %i\t %i\t %s\n' %(chromosome, genomic_coordinate, genomic_coordinate+window, '-'))
                    c += 1
                    genomic_coordinate += 1
                    continue
                if segment.count('T') >= occurences:
                    pAi.write('%s\t %i\t %i\t %s\n' %(chromosome, genomic_coordinate, genomic_coordinate+window, '-'))
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
            if int(line[1]) > int(previous_line[1]) and int(line[1]) <= int(previous_line[2]):
                previous_line[2] = line[2]
            else:
                fo.write('\t'.join(previous_line))
                previous_line = line
            

start_time = time.time()
extract_pAi_from_genome('genome', 9, 7, 6)
print (time.time() - start_time)

