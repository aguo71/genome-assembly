SCRIPT COMMAND CALLS - 1:
contamination.sh can be ran by: sh contamination.sh <reads.txt> <vector.txt> <k>
    this script implements de-contamination algorithm for genome assembly
correction.sh can be ran by: sh correction.sh <reads.txt> <k> <t> <d>
    this script implements correction algorithm for genome assembly
extract.sh can be ran by: sh extract.sh
    this script extracts 50-mers from the hiv-1_genome.txt file and mutates them with a global 1% error rate
plot.sh can be ran by: sh plot.sh
    this script calculates and plots correction algorithm efficiency using k from {6, ..., 25} and t = {4, 6, 8, 10, 12}
plot2.sh can be ran by: sh plot2.sh
    this script calculates and plots correction algorithm efficiency using k from {6, ..., 25} and t = {2, 4, 6}

PROGRAM METHODOLOGY/FILES:
Overall, contamination.py and correction.py implement the two respectively-named project algorithms. Extract.py is used to extract 50-mers
from the hiv-1_genome and mutates them with global 1% error. correction.py also writes the corresponding inputs k and t to plot.txt for Global.py, and 
writes the corrected sequence into corrected.txt also for Global.py. It then writes the true reads and mutated reads into true_reads.txt and mutated_reads.txt
respectively. Global.py calculates ungapped global alignment score between two sequences, and writes this quantity to plot.txt. plot_graph.py reads from plot.txt
and plots a -log quantity, as well as calculates some important quantities of the input data.


ALGORITHM DESIGN EXPLANATION:
When an infrequent kmer needs to be replaced, my program tries to find the most frequent kmer (with lowest distance at most d away) 
to replace the current infrequent kmer, but if the algorithm notices after that replacement there are still infrequent kmers in the k window 
around, I try replacing the original infrequent kmer with the next most frequent kmer instead. I try this until
I find a replacement which removes all nearby infrequent kmers, or I default to using the most frequent kmer (with the minimum distance away). 
If there are no possible replacements for the given infrequent kmer, it is kept in the sequence.
Then, I jump ahead k steps in the read if a replacement was made, or just 1 step ahead otherwise.

The rationale for this is that errors are unlikely, so within a window of k, there should only be a single 
substitution error. If a kmer results in another error being spotted nearby (within k), those errors are too close in 
proximity to be likely given the error rate. This means we don't want to replace multiple times in the same window, so I try to find 
the most optimal replacement which removes all nearby infrequent kmers as well (and prevents new ones from 
forming). If this is not possible, the algorithm defaults to using the most frequent kmer because it is not sure
what the optimal replacement is. It jumps ahead k steps to avoid doing duplicate work after replacing a kmer since 
the process of replacing already considers the nearby window of plus or minus k since computational time/efficiency is also important.

----
SCRIPT COMMAND CALLS - 2:

The debruijn script can be called as such: sh debruijn.sh <reads.txt>, where <reads.txt> is a user input. 

This will automatically generate a debruijn.dot file which is a 
representation of a graph, which can be visualized using the shell command: dot -Tpng debruijn.dot -o debruijn.png.

The assembly script can be called as such: sh assembly.sh <reads.txt> <vector file> <k for contamination> <output file>  <k for correction>
<t> <d> <output file 2>. 

Each of these is a user input, where reads.txt is the input reads file, vector file is file which contains vector sequence, k for contamination is the
parameter k for contamination algorithm, output file is where the output of contamination will be stored, k for correction is the parameter for the correction algorithm, 
t is the parameter for correction, d is the parameter for correction, output file 2 is where the corrected reads will be stored. This will result in the final 
assembled sequence to be printed out.

No auxiliary files created, apart from debruijn.dot which is the output of debruijn program, and debruijn.png which is the visualization of the graph. Output.txt and output2.txt were my intermediate
output files for assembly.

No known bugs.

My correction algorithm was designed with a bias towards performance/speed rather than accuracy, but I think this slight dip in accuracy results in a
magnified effect on the assembled sequence. Thus, if my overall assembly process results in a sub-optimal sequence, it is likely due to the design choices of
correction from my pr2 choices which were made without knowledge of this context that it might be used in.


FUTURE STEPS:
In the future, should try breaking the corrected reads into kmers before running the assembly algorithm to improve performance.