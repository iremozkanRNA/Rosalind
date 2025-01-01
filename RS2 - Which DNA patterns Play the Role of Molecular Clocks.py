# What is the expected number of occurrences of a 9-mer in 500 random DNA strings, each of length 1000? Assume that the sequences are formed by selecting each nucleotide (A, C, G, T) with the same probability (0.25).

import random

def random_dna_sequence(lenght):
    return "".join(random.choices("ACGT", k=lenght))

def random_multiple_dna_sequences(num_sequences, sequence_length):
    return [random_dna_sequence(sequence_length) for _ in range(num_sequences)]

dna_sequences = random_multiple_dna_sequences(500,1000)

prob_nt =0.25
lenght_9mer=9
len_dna_seq=1000
num_dna_seq=500

prob_9mer = prob_nt ** lenght_9mer
#The number of occurences of a 9 elength sub-string in a string of length 1000
possible_positions_in_one_string = len_dna_seq - lenght_9mer +1
expected_occurences_in_one_string = possible_positions_in_one_string * prob_9mer
expected_occurences_in_all_strings = num_dna_seq * expected_occurences_in_one_string
print(expected_occurences_in_all_strings)

# Probability Calculation: The probability of any specific 9-mer is calculated as 0.25**9
# Positions Calculation: The number of possible positions for a 9-mer in a single string is (1000-9+1).
# Expected Occurrences in One String: Multiply the probability of the 9-mer by the number of possible positions.
# Total Expected Occurrences: Multiply the expected occurrences in one string by the total number of strings (500).

# Implanted Motif Problem: Find all (k, d)-motifs in a collection of strings.
#   Input: A collection of strings Dna, and integers k and d.
#   Output: All (k, d)-motifs in Dna.

def HammingDistance(p,q):
    if len(p) != len(q):
        return ValueError("The sequences must of equal lenght")
    else:
        return sum(1 for i in range(len(p)) if p[i] != q[i])
    
def Neighbors(pattern,d):
    """ the set of all k-mers whose Hamming distance from Pattern does not exceed d."""
    if d == 0:
        return {Pattern}
    if len(pattern) == 1:
        return {"A", "C", "G", "T"}
    neighbors = set()
    suffixneighbors = Neighbors(pattern[1:],d)
    for suffix in suffixneighbors:
        if HammingDistance(pattern[1:],suffix) < d:
            for x in "ACTG":
                neighbors.add(x +suffix)
        else:
            neighbors.add(pattern[0]+suffix)
    return neighbors

def motif_enumeration(DNA, k, d):
    """Find all (k,d) motifs in a collection of DNA strings"""
    Patterns = set()
    # iterate over each DNA string in DNA
    for dna_string in DNA:
        #iterate over each k-mer in the dna_string
        for i in range(len(dna_string -k+1)):
            kmer = dna_string[i:i+k]
            kmer_neighbors = Neighbors(kmer,d)
            # check if each kmer_neighbors appear in all DNAstrings with at most d mismatches
            for neighbor in kmer_neighbors:
                if all(any(HammingDistance(neighbor,dna_string[j:j+k]) <=d 
                           for j in range(len(dna_string)-k+1))
                           for dna_string in DNA):
                    Patterns.add(neighbor)
    return Patterns


# Exercise Break:
# Compute the entropy of the NF-κB motif matrix (reproduced below). 

from Bio.Seq import Seq
from Bio import motifs
import math
# Step 1: Define the motifs as DNA sequences
instances = [
    Seq("TCGGGGGTTTTT"),
    Seq("CCGGTGACTTAC"),
    Seq("ACGGGGATTTTC"),
    Seq("TTGGGGACTTTT"),
    Seq("AAGGGGACTTCC"),
    Seq("TTGGGGACTTCC"),
    Seq("TCGGGTATAACC"),
    Seq("TCGGGGATTCAT"),
    Seq("TCGGGGATTCCT"),
    Seq("TAGGGGAACTAC")
]
# Step 2: Create a motif object
motif = motifs.create(instances)
# Step 3: Calculate entropy for each column in the motif matrix
def calculate_entropy(motif):
    total_entropy = 0
    for position in range(len(motif.consensus)):
        counts = [motif.counts[nucleotide][position] for nucleotide in "ACTG"]
        total = sum(counts)
        column_entropy = 0
        for i in counts:
            if i > 0:
                frequency = i/total
                column_entropy -= frequency*math.log2(frequency)
        total_entropy += column_entropy
    return total_entropy

entropy = calculate_entropy(motif)
print(f"total entropy is {entropy}")
