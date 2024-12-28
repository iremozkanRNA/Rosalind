# Pattern Count code
def PatternCount(Text, Pattern):
  count = 0
  for i in range(len(Text) - len(Pattern)):
    if Text[i:i+len(Pattern)] == Pattern:
      count += 1
  return count

PatternCount("ACAACTATGCATACTATCGGGAACTATCCT", "ACTAT")

# The most frequency k-mer in a string
# Frequent Words Problem: Find the most frequent k-mers in a string.

def frequencyTable(Text,k):
    """Builds a frequency table for all k-mers in the given text"""
    freqMap={}
    for i in range(len(Text)-k+1):
        Pattern=Text[i:i+k]
        if Pattern in freqMap:
            freqMap[Pattern] += 1
        else:
            freqMap[Pattern] = 1
    return freqMap        

def frequent_words(Text,k):
    """Finds the most frequent k-mers in a given text"""
    frequentPatterns =[]
    freqMap=frequencyTable(Text,k)
    maxValue=max(freqMap.values())
    for pattern in freqMap:
        if freqMap[pattern]==maxValue:
            frequentPatterns.append(pattern)
    return frequentPatterns

# Find the reverse complement of a DNA string.

# Insert your reverse_complement function here, along with any subroutines you need
def reverse_complement(sequence):
    """Calculate the reverse complement of a DNA pattern."""
    #Generate the dictionary for nt complements
    complement = {"A":"T", "T":"A", "G":"C", "C":"G"}
    #Generate the complement of the DNA string
    complemented = "".join(complement[base] for base in sequence)
    #Rever the complemented string
    reversed = complemented[::-1]
    return reversed

# Find all occurrences of a pattern in a string.
# Input: Two strings, Pattern and Genome.
# Output: A collection of integers specifying all starting positions where Pattern appears as a substring of Genome.

def pattern_matching(pattern, genome):
    """Find all occurences of a pattern in a genome"""
    positions = []
    for i in range(len(genome)-len(pattern)+1):
        if genome[i:i+len(pattern)] == pattern:
            positions.append(i)
    return positions

# Return a space-separated list of starting positions (in increasing order) where CTTGATCAT appears as a substring in the Vibrio cholerae genome.
os.getcwd()
'/Users/merinakzo/Documents'
os.chdir("/Users/merinakzo/Downloads")
vc = open("Vibrio Cholerae.txt")
def pattern_matching(pattern,genome):
    positions = []
    for i in range(len(genome)-len(pattern)+1):
        if genome[i:i+len(pattern)] == pattern:
            positions.append(i)
    return positions

# Clump Finding Problem: Find patterns forming clumps in a string.
# Input: A string Genome, and integers k, L, and t.
# Output: All distinct k-mers forming (L, t)-clumps in Genome.

# L: Slide a window of fixed length L along the genome, looking for a clump of k-mers.
# L=500 reflects the typical length of ori in bacterial genomes.

# t: a k-mer appearing at least t times
# TGCA forms a (L,t)-(25,3)-clump in the following Genome:
# gatcagcataagggtccCTGCAATGCATGACAAGCCTGCAGTtgttttac

def find_clumps(Text, k, L, t):
    """ Find patterns forming clumps in a genome
    Parameters:
    - Text: The genome (string).
    - k: Length of the pattern (k-mer).
    - L: Length of the window.
    - t: Minimum number of times a k-mer must appear in the window.
    
    Returns:
    - A list of all distinct k-mers forming clumps."""
    Patterns = []
    n = len(Text)
    for i in range(n-L + 1):
        Window = Text[i:i+L]
        freqMap =frequencyTable(Window, k)
        for key in freqMap:
            if freqMap[key] >= t:
                Patterns.append(key)            
    return list(set(Patterns)) # set removes duplicates


# Find the minimum skew of a DNA string.
# Input: A DNA string Genome.
# Output: All integer(s) i minimizing Skew(Prefix_i(Text)) over all values of i (from 0 to |Genome|).

# Skew is the difference between the total number of Gs and Cs in the first i nucleotides of a genome.
# Skew(Genome)= #G - #C
def skew(Genome):
    skew = [0]
    for i in range(len(Genome)):
        if Genome[i] =="G":
            skew.append(skew[i]+1)
        elif Genome[i] == "C":
            skew.append(skew[i]-1)
        else:
            skew.append(skew[i])
    return skew

# Peculiar Statistics of the Forward and Reverse Half-Strands
# Input: A DNA string Genome.
# Output: All integer(s) i minimizing Skew(Prefix_i(Text)) over all values of i (from 0 to |Genome|).
def minimum_skew(Genome):
    """Find positions in a genome where the skew diagram attains a minimum."""
    skew=[0]
    for i in range(len(Genome)):
        if Genome[i] =="G":
            skew.append(skew[i]+1)
        elif Genome[i] == "C":
            skew.append(skew[i] -1)
        else:
            skew.append(skew[i])
    minskew = min(skew)
    return [i for i in range(len(skew)) if skew[i]== minskew]

# Find the hamming distance between two strings.
# CGAAT and CGGAC have two mismatches. 
# The number of mismatches between strings p and q is called the Hamming distance between these strings and is denoted HammingDistance(p, q).
# Compute the Hamming distance between two strings.
# Input: Two strings of equal length.
# Output: The Hamming distance between these strings.

def hamming_distance(p,q):
    """Calculate the Hamming distance btw two strings"""
    return sum(1 for i in range(len(p)) if p[i] != q[i])

# Find all approximate occurrences of a pattern in a string.
# Input: Strings Pattern and Text along with an integer d.
# Output: All starting positions where Pattern appears as a substring of Text with at most d mismatches.

def approximate_pattern_matching(pattern, Text, d):
    """Find all starting positions where Pattern appears as a substring of Text with at most d mismatches."""
    positions=[]
    for i in range(len(Text)-len(pattern)+1):
        substring=Text[i:i+len(pattern)]
        if hamming_distance(pattern,substring) <= d:
            positions.append(i)
    return positions

def approximate_pattern_count(Text, pattern,d):
    """Find the number of times that a pattern appears in a string with at most d mismatches."""
    count=0
    for i in range(len(Text) - len(pattern)+1):
        substring=Text[i:i+len(pattern)]
        if hamming_distance(pattern,substring) <= d:
            count += 1
    return count

# Find the most frequent k-mers with mismatches in a string.
# Input: A string Text as well as integers k and d.
# Output: All most frequent k-mers with up to d mismatches in Text.

