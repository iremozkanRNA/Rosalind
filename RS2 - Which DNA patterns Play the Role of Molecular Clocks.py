# What is the expected number of occurrences of a 9-mer in 500 random DNA strings, each of length 1000? Assume that the sequences are formed by selecting each nucleotide (A, C, G, T) with the same probability (0.25).

import random

def random_dna_sequence(length):
    return "".join(random.choices("ACGT", k=length))

def random_multiple_dna_sequences(num_sequences, sequence_length):
    return [random_dna_sequence(sequence_length) for _ in range(num_sequences)]

dna_sequences = random_multiple_dna_sequences(500,1000)

prob_nt =0.25
lenght_9mer=9
len_dna_seq=1000
num_dna_seq=500
# Probability Calculation: The probability of any specific 9-mer is calculated as 0.25**9
prob_9mer = prob_nt ** lenght_9mer
#The number of occurences of a 9 length sub-string in a string of length 1000
# Positions Calculation: The number of possible positions for a 9-mer in a single string is (1000-9+1).
possible_positions_in_one_string = len_dna_seq - lenght_9mer +1
# Expected Occurrences in One String: Multiply the probability of the 9-mer by the number of possible positions.
expected_occurences_in_one_string = possible_positions_in_one_string * prob_9mer
# Total Expected Occurrences: Multiply the expected occurrences in one string by the total number of strings (500).
expected_occurences_in_all_strings = num_dna_seq * expected_occurences_in_one_string
print(expected_occurences_in_all_strings)


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


# Equivalent Motif Finding Problem: 
# Given a collection of strings, find a collection of k-mers (one from each string) that minimizes the distance between all possible patterns and all possible collections of k-mers.
# Input: A collection of strings Dna and an integer k.
# Output: A k-mer Pattern and a collection of k-mers, one from each string in Dna, minimizing d(Pattern, Motifs) among all possible choices of Pattern and Motifs.

 # Find a k-mer Pattern that minimizes d(Pattern, Dna) over all k-mers Pattern,
 # the same task that the Equivalent Motif Finding Problem is trying to achieve. 
 # Such a k-mer is called a median string for Dna.

# Median String Problem: Find a median string.

# Input: A collection of strings Dna and an integer k.
# Output: A k-mer Pattern that minimizes d(Pattern, Dna) among all possible choices of k-mers.

# Insert your median_string function here, along with any subroutines you need
def HammingDistance(p,q):
    if len(p) != len(q):
        raise ValueError("Strings must be of equal length")
    else:
        return sum(1 for i in range(len(p)) if p[i] != q[i])
    
def DistanceBtwPatternsandStrings(pattern, DNA):
    """Calculate the total distance between a pattern and a collection of DNA strings."""
    total_distance = 0
    for dna_string in DNA:
        min_distance = float('inf')
        for i in range(len(dna_string) - len(pattern) + 1):
            kmer = dna_string[i:i+len(pattern)]
            current_distance = HammingDistance(pattern,kmer)
            if current_distance < min_distance:
                min_distance = current_distance
        total_distance += min_distance
    return total_distance

from itertools import product

def median_string(dna, k):
    # Generate all possible kmers in a string
    patterns = ["".join(x) for x in product("ACTG", repeat=k)]
    distance = float('inf')
    # Iterate over each pattern
    for pattern in patterns:
        if distance > DistanceBtwPatternsandStrings(pattern, dna):
            distance = DistanceBtwPatternsandStrings(pattern, dna)
            Median = pattern
    return Median


# Exercise Break: Compute Pr("​​​​​​​TCGTGGATTTCC" | Profile), where Profile is the matrix shown above.

# Step 1: Define the Profile matrix
profile = {
    "A": [0.2, 0.2, 0.0, 0.0, 0.0, 0.0, 0.9, 0.1, 0.1, 0.1, 0.3, 0.0],
    "C": [0.1, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.4, 0.1, 0.2, 0.4, 0.6],
    "G": [0.0, 0.0, 0.1, 0.1, 0.9, 0.9, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0],
    "T": [0.7, 0.2, 0.0, 0.0, 0.1, 0.1, 0.0, 0.5, 0.8, 0.7, 0.3, 0.4]
}

# Step 2: Define the sequence
sequence = "TCGTGGATTTCC"

# Step 3: Calculate the probability of the sequence given the profile matrix
def probability(sequence):
    prob = 1 
    for i in range(len(sequence)):
        prob *= profile[sequence[i]][i]
    return prob

probability = probability(sequence)
print(probability)

# Profile-most Probable k-mer Problem: Find a Profile-most probable k-mer in a string.
# Input: A string Text, an integer k, and a 4 × k matrix Profile.
# Output: A Profile-most probable k-mer in Text.

#profile=[{'A': 0.2, 'C': 0.4, 'G': 0.3, 'T': 0.1},
#{'A': 0.2, 'C': 0.3, 'G': 0.3, 'T': 0.2},
#{'A': 0.3, 'C': 0.1, 'G': 0.5, 'T': 0.1},
#{'A': 0.2, 'C': 0.5, 'G': 0.2, 'T': 0.1},
#{'A': 0.3, 'C': 0.1, 'G': 0.4, 'T': 0.2}]


def reformat_profile(profile_file):
    """
    A helper function to reformat the space-delimited profile matrices as dictionaries.
    Input: profile_file (str) - Path to the file containing the profile matrix.
    Output: A list of dictionaries where each dictionary represents a column in the profile matrix.
    """
    profile_dicts, letters, profile_lists = {}, 'ACGT', []
    for line_number, line in enumerate(open(profile_file)):
        letter_probs = line.strip().split(' ')
        letter = letters[line_number]
        for column_number, column in enumerate(letter_probs):
            profile_dicts.setdefault(column_number, {})
            profile_dicts[column_number][letter] = float(column)
    for column in sorted(list(profile_dicts.keys())):
        profile_lists.append(profile_dicts[column])
    return profile_lists


def calculate_kmer_probability(kmer, profile):
    """
    Calculate the probability of a k-mer based on a given profile matrix.
    Input:
        kmer (str) - The k-mer whose probability is to be calculated.
        profile (list of dict) - The profile matrix as a list of dictionaries.
    Output:
        prob (float) - The calculated probability of the k-mer.
    """
    prob = 1.0
    for i, nucleotide in enumerate(kmer):
        prob *= profile[i][nucleotide]
    return prob



def profile_most_probable_kmer(text, k, profile):
    
    max_prob = -1
    most_probable_kmer = text[:k]  # Default to the first k-mer in case of ties or no valid probabilities

    # Iterate through all possible k-mers in the string
    for i in range(len(text) - k + 1):
        kmer = text[i:i + k]
        prob = calculate_kmer_probability(kmer, profile)
        if prob > max_prob:
            max_prob = prob
            most_probable_kmer = kmer

    return most_probable_kmer

# GreedyMotifSearch() constructs Profile(Motifs) and 
# selects the Profile-most probable k-mer from Dnai based on this profile matrix.

# Implement GreedyMotifSearch().
# Input: Integers k and t, followed by a space-separated collection of strings Dna.
# Output: A collection of strings BestMotifs resulting from 
# applying GreedyMotifSearch(Dna, k, t). If at any step you find more than one Profile-most 
# probable k-mer in a given string, use the one occurring first.



def GreedyMotifSearch(Dna, k, t):
    """
    Finds the best motifs in a collection of DNA strings using a greedy approach.
    
    Parameters:
        Dna (list): A list of DNA strings.
        k (int): The length of the motif.
        t (int): The number of DNA strings.
        
    Returns:
        list: A list of the best motifs (one from each string in Dna).
    """
    # Initialize BestMotifs with the first k-mers from each string
    BestMotifs = [dna[:k] for dna in Dna]
    
    # Iterate over all possible k-mers in the first string
    for i in range(len(Dna[0]) - k + 1):
        Motif1 = Dna[0][i:i + k]  # Select a k-mer from the first string
        Motifs = [Motif1]  # Initialize Motifs with this first k-mer
        
        # Build motifs iteratively for the rest of the strings
        for j in range(1, t):
            Profile = build_profile(Motifs)
            Motif_j = profile_most_probable_kmer(Dna[j], k, Profile)
            Motifs.append(Motif_j)
        
        # Update BestMotifs if the new set of motifs has a lower score
        if score(Motifs) < score(BestMotifs):
            BestMotifs = Motifs
    
    return BestMotifs


def build_profile(motifs):
    """
    Builds a profile matrix from a list of motifs.
    
    Parameters:
        motifs (list): A list of motifs (strings).
        
    Returns:
        dict: A profile matrix as a dictionary with keys 'A', 'C', 'G', 'T'.
    """
    k = len(motifs[0])
    t = len(motifs)
    
    # Initialize profile matrix with zeros
    profile = {nucleotide: [0] * k for nucleotide in "ACGT"}
    
    # Count occurrences of each nucleotide at each position
    for motif in motifs:
        for i, nucleotide in enumerate(motif):
            profile[nucleotide][i] += 1
    
    # Convert counts to probabilities by dividing by t (number of motifs)
    for nucleotide in "ACGT":
        profile[nucleotide] = [count / t for count in profile[nucleotide]]
    
    return profile


def profile_most_probable_kmer(dna, k, profile):
    """
    Finds the Profile-most probable k-mer in a DNA string.
    
    Parameters:
        dna (str): A DNA string.
        k (int): The length of the k-mer.
        profile (dict): The profile matrix as a dictionary.
        
    Returns:
        str: The most probable k-mer based on the profile matrix.
    """
    max_prob = -1
    most_probable_kmer = dna[:k]  # Default to the first k-mer
    
    # Iterate through all possible k-mers in the DNA string
    for i in range(len(dna) - k + 1):
        kmer = dna[i:i + k]
        prob = 1
        
        # Calculate the probability of this k-mer based on the profile
        for j, nucleotide in enumerate(kmer):
            prob *= profile[nucleotide][j]
        
        # Update most probable k-mer if this one has a higher probability
        if prob > max_prob:
            max_prob = prob
            most_probable_kmer = kmer
    
    return most_probable_kmer


def score(motifs):
    """
    Calculates the score of a set of motifs by summing mismatches to the consensus sequence.
    
    Parameters:
        motifs (list): A list of motifs (strings).
        
    Returns:
        int: The total number of mismatches across all positions.
    """
    consensus = ""
    k = len(motifs[0])
    
    # Build consensus sequence by finding most frequent nucleotide at each position
    for i in range(k):
        counts = {nucleotide: 0 for nucleotide in "ACGT"}
        for motif in motifs:
            counts[motif[i]] += 1
        consensus += max(counts, key=counts.get)
    
    # Calculate total mismatches between consensus and all motifs
    total_score = 0
    for motif in motifs:
        total_score += sum(1 for i in range(k) if motif[i] != consensus[i])
    
    return total_score


# Example Usage
if __name__ == "__main__":
    Dna = [
        "AGCTGACCTG",
        "CCTGAGCTGA",
        "GACTGAGCTA",
        "TGACTGACCT"
    ]
    
    k = 3  # Length of motif
    t = len(Dna)  # Number of DNA strings
    
    result = GreedyMotifSearch(Dna, k, t)
    
    print("Best Motifs:")
    print(result)

# Implement GreedyMotifSearch() with pseudocounts.

# Input: Integers k and t, followed by a space-separated collection of strings Dna.
# Output: A collection of strings BestMotifs resulting from applying GreedyMotifSearch(Dna, k, t)
# with pseudocounts. If at any step you find more than one Profile-most probable k-mer in a given string,
# use the one occurring first.


import sys

# Please do not remove package declarations because these are used by the autograder.
def build_profile(Motifs):
    k = len(Motifs[0])
    t = len(Motifs)
    # Initialize profile with pseudocounts
    Profile = {nucleotide: [1] * k for nucleotide in "ACTG"}
    
    # Count occurrences of each nucleotide at each position
    for motif in Motifs:
        for i, nucleotide in enumerate(motif):
            Profile[nucleotide][i] += 1
            
    # Convert counts to probabilities by dividing by (t + 4) (Laplace correction)
    for nucleotide in "ACTG":
        Profile[nucleotide] = [count / (t + 4) for count in Profile[nucleotide]]
        
    return Profile

def score(Motifs):
    """Calculates the score of a set of motifs by summing mismatches to the consensus sequence."""
    k = len(Motifs[0])
    consensus = ""
    
    # Build consensus sequence by finding most frequent nucleotide at each position
    for i in range(k):
        counts = {nucleotide: 0 for nucleotide in "ACTG"}  # Initialize counts to 0
        for motif in Motifs:
            counts[motif[i]] += 1
        consensus += max(counts, key=counts.get)
    
    # Calculate total mismatches between consensus and all motifs
    total_score = 0
    for motif in Motifs:
        total_score += sum(1 for j in range(k) if motif[j] != consensus[j])
    
    return total_score

def profile_most_probable_kmer(Dna, k, Profile):
    """Finds the Profile-most probable k-mer in a DNA string."""
    Max_prob = -1
    most_probable_kmer = Dna[:k]
    
    # Iterate through all possible k-mers in the DNA string
    for i in range(len(Dna) - k + 1):
        kmer = Dna[i:i+k]
        Probs = 1
        
        # Calculate the probability of this k-mer based on the profile
        for j, nucleotide in enumerate(kmer):
            Probs *= Profile[nucleotide][j]
        
        # Update most probable k-mer if this one has a higher probability
        if Probs > Max_prob:
            Max_prob = Probs
            most_probable_kmer = kmer
            
    return most_probable_kmer

def greedy_motif_search_pseudocounts(Dna, k, t):
    """Finds the best motifs using a greedy search with pseudocounts."""
    BestMotifs = [dna[:k] for dna in Dna]  # Initialize BestMotifs with first k-mer from each DNA string
    
    # Iterate through all possible k-mers in the first DNA string
    for i in range(len(Dna[0]) - k + 1):
        Motif1 = Dna[0][i:i+k]
        Motifs = [Motif1]
        
        # Build motifs iteratively using the profile matrix
        for j in range(1, t):
            Profile = build_profile(Motifs)  # Build profile from current motifs
            Motif_j = profile_most_probable_kmer(Dna[j], k, Profile)  # Find most probable k-mer
            Motifs.append(Motif_j)  # Append new motif to Motifs
            
        # Update BestMotifs if current motifs have a lower score
        if score(Motifs) < score(BestMotifs):
            BestMotifs = Motifs
            
    return BestMotifs