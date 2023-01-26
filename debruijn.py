import sys
import networkx as nx

# merges two strings such that the first and second share a substring with a unique prefix and suffix
def merge_two_strings(seq1, seq2):
    n = 1
    max_n = 1
    while n <= len(seq1) and n <= len(seq2):
        if seq1[-n:] == seq2[:n]:
            max_n = n
        n += 1
    if len(seq1) > len(seq2):
        new_val = seq1 + seq2[max_n:]
    else:
        new_val = seq1[:-max_n] + seq2
    return new_val

# computes simplified debruijn graph of input reads of length k
def debruijn():
    try:
        reads_file = sys.argv[1]

    except:
        print("ERROR: Should be 1 input to debruijn algorithm.")
        return
    
    with open(reads_file, 'r') as file:
        entry = file.readline().replace('\n', "").upper()
        reads = [] # list of entries in reads file(strings)
        while entry != '':
            reads.append(entry)
            entry = file.readline().replace('\n', "").upper()

    # create debruijn graph using networkx with nodes of k-1 mers and edges implicit kmers based on reads
    G = nx.MultiDiGraph(directed=True)
    for read in reads:
        first_val = read[:-1]
        second_val = read[1:]
        if first_val not in G:
            G.add_node(first_val)
        if second_val not in G:
            G.add_node(second_val)
        G.add_edge(first_val, second_val)

    # simplify singleton paths
    while True:
        # no nodes found yet, so done until node found
        done = 0
        G0 = G.copy()
        for node in G:
            if len(list(G.successors(node))) == 1:
                next_node = list(G.successors(node))[0]
                if len(list(G.predecessors(next_node))) == 1:
                    # merge
                    pred = G.predecessors(node)
                    succ = G.successors(next_node)
                    curr_mid = node[1:]
                    next_mid = next_node[:-1]
                    if len(curr_mid) > len(next_mid):
                        new_val = node[0] + curr_mid + next_node[-1]
                    else:
                        new_val = node[0] + next_mid + next_node[-1]
                    G0.add_node(new_val)
                    for in_edge in pred:
                        G0.add_edge(in_edge, new_val)
                    for out_edge in succ:
                        G0.add_edge(new_val, out_edge)
                    G0.remove_node(next_node)
                    G0.remove_node(node)
                    # one node merged, need to restart loop
                    done = 1
                    break
            else:
                # no nodes merged, done
                done = 2
        G = G0
        if done == 0 or done == 2:
            break
        else:
            continue
    
    # save dot graph
    nx.nx_pydot.write_dot(G, 'debruijn.dot')
    with open('debruijn.dot', 'w') as file:
        file.write('digraph {\n')
        for node in G:
            edges = list(G.successors(node))
            file.write(node + ';\n')
            for edge in edges:
                label = merge_two_strings(node, edge)
                file.write(node + ' -> ' + edge + '  [label=' + label + '];\n')
        file.write('}')

    # find optimal sequence(s)
    source_nodes = []
    dest_nodes = []
    for node in G:
        if len(list(G.predecessors(node))) == 0:
            source_nodes.append(node)
        if len(list(G.successors(node))) == 0:
            dest_nodes.append(node)
    longest = 0
    seqs = []
    for source in source_nodes:
        for dest in dest_nodes:
            for path in nx.all_simple_paths(G, source=source, target=dest):
                for i in range(1, len(path)):
                    node = path[0]
                    next_node = path[i]
                    new_val = merge_two_strings(node, next_node)
                    path[0] = new_val
                    if len(new_val) > longest:
                        longest = len(new_val)
                        seqs = [new_val]
                    elif len(new_val) == longest:
                        seqs.append(new_val)
    seqs.sort()
    for path in seqs:
        print(path)

if __name__ == "__main__":
    debruijn()