import sys

# corrects reads based on replacement map, avoiding infrequent kmers, and returns indices of corrections as well as new reads
def correct_reads(reads, replacement_map, frequent_kmers, k):
    indices = []
    new_reads = []
    # loops through each read
    for i, read in enumerate(reads):
        start = 0
        while start+k <= len(read):
            modified = False # indicates if we should increment index by 1 or k to check next kmer
            # checks if current window from start to start+k is infrequent kmer
            if read[start:k+start] in replacement_map:
                curr_read = read[start:k+start]
                replaced = False # indicates if proceeding loop finds optimal replacement without creating new infrequent kmer
                for replacement in replacement_map[curr_read]:
                    read = read[:start] + replacement + read[k+start:]
                    # checks if new infrequent kmers formed by checking a shifting window of size k on both sides
                    counter = 1
                    flag = False # indicates if new infrequent kmer was formed
                    while start-counter >= 0 and start+k+counter <= len(reads[i]) and counter < k:
                        if read[start-counter:k+start-counter] not in frequent_kmers or \
                            read[start+counter:k+start+counter] not in frequent_kmers:
                            flag = True
                            break
                        counter += 1
                    if flag:
                        continue # continues loop with next possible replacement
                    else:
                        replaced = True
                        modified = True
                        break # successfully found a replacement or checked all possible replacements
                # if every replacement caused new infrequent kmer, use most frequent one by default
                if not replaced:
                    if replacement_map[curr_read] != []:
                        read = read[:start] + replacement_map[curr_read][0] + read[k+start:]
                        modified = True
                        if str(i) not in indices:
                            indices.append(str(i))
                    else:
                        read = read[:start] + curr_read + read[k+start:]
                else:
                    if str(i) not in indices:
                        indices.append(str(i))
            if modified:
                start += k ## changes output if switch to 1
            else:
                start += 1
        new_reads.append(read)
    return indices, new_reads

# calculates replacements for each infrequent kmer
def get_replacements(infrequent_kmers, frequent_kmers, kmer_count, d):
    replacement_map = {} # {infrequent kmer: [descending list of frequent_kmer replacements with lowest hamming dist]}
    for infreq in infrequent_kmers:
        least_dist = float('inf')
        possibilities = []
        for kmer in frequent_kmers:
            dist = hamming(infreq, kmer)
            if dist <= d:
                if dist < least_dist:
                    least_dist = dist
                    possibilities = [(kmer, kmer_count[kmer])]
                elif dist == least_dist:
                    possibilities.append((kmer, kmer_count[kmer]))
        possibilities = sorted(possibilities, key=lambda x: x[1], reverse=True)
        possibilities = [tuple[0] for tuple in possibilities]
        replacement_map[infreq] = possibilities
    return replacement_map

# computes hamming distance
def hamming(str1, str2):
    counter = 0
    for i in range(len(str1)):
        if str1[i] != str2[i]:
            counter += 1
    return counter

# calculates all kmers from reads, computing counts of each kmer, and adding them to infrequent or frequent set based on t
def get_kmers(reads, k, t):
    infrequent = set()
    frequent = set()
    kmer_count = {} # dict of {kmer:count}
    for read in reads:
        start = 0
        end = k
        while end <= len(read):
            kmer = read[start: end]
            if kmer in kmer_count:
                kmer_count[kmer] += 1
                if kmer_count[kmer] >= t and kmer not in frequent:
                    frequent.add(kmer)
                    infrequent.remove(kmer) 
            else:
                kmer_count[kmer] = 1
                infrequent.add(kmer)
            start += 1
            end += 1
    return infrequent, frequent, kmer_count

# computes contamination sequences of each read
def correction():
    try:
        reads_file = sys.argv[1]
        k = int(sys.argv[2])
        t = int(sys.argv[3]) 
        d = int(sys.argv[4])
        output_file = sys.argv[5]

    except:
        print("ERROR: Should be 4 inputs to correction algorithm.")
        return
    
    with open(reads_file, 'r') as file:
        entry = file.readline().replace('\n', "").upper()
        reads = [] # list of entries in reads file(strings)
        while entry != '':
            reads.append(entry)
            entry = file.readline().replace('\n', "").upper()

    if k <= 0:
        print("ERROR: k value should be positive.")
        return
    
    if t < 0:
        print("ERROR: t value should be non-negative.")
        return

    if d < 0:
        print("ERROR: d value should be non-negative.")
        return       

    infrequent_kmers, frequent_kmers, kmer_count = get_kmers(reads, k, t)
    replacement_map = get_replacements(infrequent_kmers, frequent_kmers, kmer_count, d)
    replaced_indices, new_reads = correct_reads(reads, replacement_map, frequent_kmers, k)

    with open(output_file, 'w') as file:
        for read in new_reads:
            file.write(read+'\n')

if __name__ == "__main__":
    correction()