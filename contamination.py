import sys

# calculates seeds from alphabet which score greater than T against any seed from the query
def get_seeds(read_kmer, vector_kmers):
    seeds = [] # list of vector_kmer tuples
    for kmer in vector_kmers:
        if kmer[0] == read_kmer:
            seeds.append(kmer)
    return seeds

# calculates all kmers from vector
def get_vector_kmers(vector, k):
    kmers = []
    start = 0
    end = k
    while end <= len(vector):
        kmer = (vector[start: end], start) # saves tuple of (kmer, index position #)
        kmers.append(kmer)
        start += 1
        end += 1
    return kmers

# computes contamination sequences of each read, and returns reads without any contamination
def contamination():
    try:
        reads_file = sys.argv[1]
        vector_file = sys.argv[2]
        k = int(sys.argv[3])
        output_file = sys.argv[4]

    except:
        print("ERROR: Should be 3 inputs to de-contamination algorithm.")
        return
    
    with open(reads_file, 'r') as file:
        entry = file.readline().replace('\n', "").upper()
        reads = [] # list of entries in reads file(strings)
        while entry != '':
            reads.append(entry)
            entry = file.readline().replace('\n', "").upper()

    with open(vector_file, 'r') as file:
        vector = file.readline().replace('\n', "").upper()
    
    if k <= 0:
        print("ERROR: k value should be positive.")
        return

    vector_kmers = get_vector_kmers(vector, k)
    contaminated = []
    new_reads = []
    for i in range(len(reads)):
        seeds_left = get_seeds(reads[i][:k], vector_kmers)
        seeds_right = get_seeds(reads[i][-k:], vector_kmers)
        if seeds_left or seeds_right:
            contaminated.append(i)

    for i in range(len(reads)):
        if i not in contaminated:
            new_reads.append(reads[i])

    with open(output_file, 'w') as file:
        for read in new_reads:
            file.write(read+'\n')


if __name__ == "__main__":
    contamination()