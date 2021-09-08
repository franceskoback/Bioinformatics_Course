import random

"""
The functions in Motif.py will return 0 for an entire motif probability even if only
one of the positions has a 0 probability of existing in the consensus string.
This doesn't make sense because a motif that differs from the consensus string
at every position will also get a total probability of 0.
In order to improve this unfair scoring, bioinformaticians often substitute zeroes
with small numbers called pseudocounts.
"""
# Input:  String Text, an integer k, and profile matrix Profile
# Output: String of most probable pattern
def ProfileMostProbablePattern(Text, k, Profile):
    n = len(Text)
    maximum = -1
    probable_pattern = ''
    for i, letter in enumerate(Text):
        for i in range(n-k+1):
            pattern = Text[i:i+k]
            probability = Pr(pattern,Profile)
            if (probability > maximum):
                maximum = probability
                probable_pattern = pattern
    if maximum == -1:
        return Text[0:0+k]
    else:
        return probable_pattern

# Input:  A set of kmers Motifs
# Output: CountWithPseudocounts(Motifs)
def CountWithPseudocounts(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    count = {}
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
            count[symbol].append(1)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1
    return count

"""
ProfileWithPseudocounts(Motifs) that takes a list of strings Motifs as input and
returns the profile matrix of Motifs with pseudocounts as a dictionary of lists
"""

# Input:  A set of kmers Motifs
# Output: ProfileWithPseudocounts(Motifs)
def ProfileWithPseudocounts(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    profile = {}
    count = CountWithPseudocounts(Motifs)
    for key, motif_lists in sorted(count.items()):
        profile[key] = motif_lists
        for motif_list, number in enumerate(motif_lists):
            motif_lists[motif_list] = number/(float(t+4))
    return profile

# motif1 = "AACGTA"
# motif2 = "CCCGTT"
# motif3 = "CACCTT"
# motif4 = "GGATTA"
# motif5 = "TTCCGG"
# motifs = [motif1, motif2, motif3, motif4, motif5]
#
# print(ProfileWithPseudocounts(motifs))

"""
 Write a function GreedyMotifSearchWithPseudocounts(Dna, k, t) that takes a list
 of strings Dna followed by integers k and t and returns the result of running
 GreedyMotifSearch, where each profile matrix is generated with pseudocounts
 """
# Input:  A list of kmers Dna, and integers k and t (where t is the number of kmers in Dna)
# Output: GreedyMotifSearch(Dna, k, t)
def GreedyMotifSearchWithPseudocounts(Dna, k, t):
    BestMotifs = []
    # search through DNA string
    for i in range(0, t):
        # starts by setting BestMotifs equal to the first k-mer from each string in Dna
        BestMotifs.append(Dna[i][0:k])
    n = len(Dna[0])
    # ranges over all possible k-mers in Dna[0], trying each one as Motifs[0]
    for i in range(n-k+1):
        Motifs = []
        Motifs.append(Dna[0][i:i+k])
        for j in range(1, t):
            # builds a profile matrix Profile for this lone k-mer, and sets Motifs[1] equal to the Profile-most probable k-mer in Dna[1]
            P = ProfileWithPseudocounts(Motifs[0:j])
            # sets Motifs[i] equal to the Profile-most probable k-mer from Dna[i] based on this profile matrix
            Motifs.append(ProfileMostProbablePattern(Dna[j], k, P))
        # GreedyMotifSearch checks whether Motifs outscores the current best scoring collection of motifs, BestMotifs
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs

# Input:  A set of kmers Motifs
# Output: A consensus string of Motifs.
def Consensus(Motifs):
    k = len(Motifs[0])
    profile = ProfileWithPseudocounts(Motifs)
    consensus = ""
    for j in range(k):
        maximum = 0
        frequentSymbol = ""
        for symbol in "ACGT":
            if profile[symbol][j] > maximum:
                maximum = profile[symbol][j]
                frequentSymbol = symbol
        consensus += frequentSymbol
    return consensus

# Input:  A set of k-mers Motifs
# Output: The score of these k-mers.
def Score(Motifs):
    count = 0
    consensus = Consensus(Motifs)
    for motif in Motifs:
        for index, letter in enumerate(motif):
            if letter != consensus[index]:
                count += 1
    return count

# Input:  String Text and profile matrix Profile
# Output: Probability value
def Pr(Text, Profile):
    p = 1
    # loop through each index(char) in text
    for index,char in enumerate(Text):
        for key, profile_lists in sorted(Profile.items()):
            if char == key:
                p *= profile_lists[index]
    return p

# k = 3
# t = 5
# Dna = ["GGCGTTCAGGCA", "AAGAATCAGTCA", "CAAGGAGTTCGC", "CACGTCAATCAC", "CAATAATATTCG"]
# print(GreedyMotifSearchWithPseudocounts(Dna, k, t))

# Input:  A profile matrix Profile and a list of strings Dna
# Output: Profile-most probable k-mer from each row of Dna
def Motifs(Profile, Dna, k):
    probable_kmer = []
    for each_string in Dna:
        probable_kmer.append(ProfileMostProbablePattern(each_string, k, Profile))
    return probable_kmer

# Profile = {'A': [0.8, 0.0, 0.0, 0.2],
#            'C': [0.0, 0.6, 0.2, 0.0],
#            'G': [0.2, 0.2, 0.8, 0.0],
#            'T': [0.0, 0.2, 0.0, 0.8]}
#
# Dnas = ["TTACCTTAAC", "GATGTCTGTC", "ACGGCGTTAG", "CCCTAACGAG", "CGTCAGAGGT"]
#
# print(Motifs(Profile, Dnas))

# Input:  A list of strings Dna, and integers k and t
# Output: RandomMotifs(Dna, k, t)
# HINT:   You might not actually need to use t since t = len(Dna), but you may find it convenient
def RandomMotifs(Dna, k, t):
    n = len(Dna[0])
    random_motifs = []
    for each_string in Dna:
        random_start = random.randint(0,n-k)
        random_motifs.append(each_string[random_start:random_start+k])
    return random_motifs
#
# Dnas = ["TTACCTTAAC", "GATGTCTGTC", "ACGGCGTTAG", "CCCTAACGAG", "CGTCAGAGGT"]
# k = 3
# t = len(Dnas)
# print(RandomMotifs(Dnas, k, t))

# Input:  Positive integers k and t, followed by a list of strings Dna
# Output: return a list of random kmer motifs
def RandomizedMotifSearch(Dna, k, t):
    random_motifs = RandomMotifs(Dna, k, t)
    best_motifs = random_motifs
    while True:
        profile = ProfileWithPseudocounts(random_motifs)
        check_random = Motifs(profile, Dna, k)
        if Score(check_random) > Score(best_motifs):
            best_motifs = check_random
        else:
            return best_motifs

"""
Input:  Positive integers k and t, followed by a list of strings Dna
Output: A list containing BestMotifs resulting from running RandomizedMotifSearch(Dna, k, t) 1000 times.
"""
def RepeatedRandomizedMotifSearch(Dna, k, t):
    BestScore = float('inf')      # start the "best score" as infinity
    BestMotifs = []               # output variable
    for i in range(1000):         # run RandomizedMotifSearch 1000 times
        Motifs = RandomizedMotifSearch(Dna, k, t)
        CurrScore = Score(Motifs)
        if CurrScore < BestScore: # if this set of Motifs is better than BestMotifs, swap
            BestScore = CurrScore
            BestMotifs = Motifs
    for i in BestMotifs:
        print(i)
    return BestMotifs
# k = 15
# t = 20
# Dna = ["AGTTCACGGAACACCTATTCTGGATCGAGGGAGCTAATGTATGGAGGGTACGCAAGGGATACATAATATGGACTCAAATTATCCTATGAGACTCCAAGGGCACCGAAGAAGCTTCTGATACTCACAGCTAAACAACGGAGCTACAGGATAATGTAGGAGATGGCTCGGTATTCAAGACTCAATAGACAGTTCACGGAACACC",
#         "TATTCTGGATCGAGGGAGCTAATGTATGGAGGGTACGCAAGGGATACATAATATGGACTCAAATTATCCTATGAGACTCCAAGGGCACCGAAGAAGCTTCTGATACTCACAGCTAAACAACGGAGCTACAGGATAATGTAGGAGATGGCTCGGTATTCAAGACTCAATAGACAGTTCACGGAACACATGCTCCGCCGACGCC",
#         "GACTGCGATGACACACATTCATGATCTGCCGACGCCCAGTGCGCCGCCGGCTGCTATATCGCTACAGTAGTATATGCCTCAGCTACCGAGAAAGCGCGCGTTGTTTCGACAGAGCTCGGGAGTCCCTACGACCACGGCAACTCTCAGGTACCTACGTTGAACGCTTCGATGGGCTACAGTGACTAGCGATGCTTCAGAAGGG",
#         "CGGGAACTGTCTCAGATTCTGCATTTACGTCCCACCTTATAATGTGATGTACTATTACTGGTATCATATCCCGCAAATAGAGATTCGTTTGTTCGAACAATGTGCGCTAGAGATGCAGAGTGGACGGTCGTCCATCATCGGAAGGTTAGGTTCTGATCCGCCACTGTGGGTAGGCCTTCACTAGTATGAGGTCCCGACGCGC",
#         "CACAACTAGGAACCGGATTCGTACATGGACCGTACCGCGTTCCTCAAAAGGATACCTGATCAACTGATCTTCGCCAAGACTGTTGTTGGCCGTGACACTTAACTGGAAATGTTCGTTGCTCATGCGAGACGTCGCCGACGCATCATCACACGTTGACAAGAAGCTGGCAGAGCTGTGGCTAGGCACAATGGCAGATATGGGC",
#         "CGTAAATGCGTGAGGGGCCTGAAGGTAATCTTCTACATATGGCACCTTCTGATATGCGTCTTCGCCCAACAGGGCCTAAGGCATGCTATGATATGCCGACGCCCGACCTTTCCCCGTGTGACCAATCAAAATATTAACAAATGTGGCTCTACGGATATGAATGAGAGTAGCCAGTGAGTGACTAATATGATTGCATGCGGTT",
#         "GTCAAGATCAGTTACATTGACCTAGTCACGATACGGCCGTGATGTCCAGTTGAGTTCGGGGGGGATTTTAACCAGGGTCTTGTATCCAGGTGGTTTATTTATCTAGGACACTCCGACCCCTGAAGCAAGTGCTGGCGGGAGAATGCGTAGTGATTCATGAGTCGCCGCGACTCTTAGAAAAATAAAAACGGTGTAAGCAAAC",
#         "AAGGAGTCTAAGTGACAGAGTATCCAGTAAAGCCGTACACTTTCTAAGTCGCGAGGTTACAGCTACAGCCATGGGTTGAACTGAATGAAGATGGCAAACGGCGTTTGATCACATTTTCACGGATTACCTACCACACAGATAGGTGATGCTCCGCCGACGCGTATCGACCACGCAGCTGACAAAGTCCGCCGCTATGTTCGAG",
#         "GGTCTATAGTTACCCTTCCAATCGATTTCTCACGTAGTATGAGTGAACGACGCATCTTTGCTTAGATGAAATGCATATCTGATTCATATCATACATAGACACGGACTTGGAGTATGACAACACGAATAGATTCCCTCCTTATGTGATTAAATCCGGTTTGGCTGTGCCAATAGCGCGAACAGTCTTCTACTGTTACAACCAC",
#         "AACACCGCGTCTCTTGTAAATGAATGAACTAATAGAATAATCTCGACGGGTTTAGCTAGTCGCCGACGCCGGCAGGTGATTATCTATCCCCTGGCAGCCCCGTAGGACGATAGTATGAGTTCTGTGCACTAGTCCAACAAACTCCGTGCCATAACAGGATGGTTCCGGCTTCGACCCCCAACATATTACCTTGGGGCCTCAC",
#         "ACGCACGAAGGCGATTATGTCAGATATACTGTCAGACGGACATTCCGGGATTAATAACATATGGCGCCAATTACTCCGGATAATCTTCGCGTCCACAAGTAAACATCGACTGACTGTTGCATGAGTCTTTGACGCACACGTCATGCATATCCAGACTAAGCCTCCTCGAAACTTGAACGTGGGATTGAGTGGACGGATGTAT",
#         "TCAGGGTCGCAGGGTACGTCAGTTTGAAACGTGTGAGATCCTTCGGACATAGGGAGTCGCCGACGAAAGAAGAGCCAAATTCAGCGTCTAAGATGTACATCGCGTAAGCTAGTAAGATCAAAATTTGCCTTAGTGCTCGGTCGGATCGCGCTCTCTGCGGTACATGGGTTTCCTAGGACGGATTCTGCACTGCAATAAGTCA",
#         "TAGCCGAGGGCATCAGGAGATGAGTCGCCAGAGCCTAGCGAACATCCAGCCCATTAACTAGGAAACCTTTACACCGTTTGGGGGCCGATCTGCTGACACGAGCAAGAGCGCCCATTGTCTCTTCTCCAGGTGTTGTCAGGCGACCTACGGGCGAACTGGACCAAAACTACCGAACTCTTCGCTCAACATGATCGTTGACATT",
#         "AGTCTCCTAGACGTGAATTGTCGCACGCACATTATTGAGTTATGAGGCCAGGATGACTGCAAACTCCGATCAGACAACTGGCTATCGTAGTATAATCCATCGCCGACGCAAAACCAATAGCGGAGAGTACTTTAACCCGACTTTGTTAGGCGTAGCTACAAGGGCAGCTGAGACGTAGAGATGAGCGTGCATGACGCGAGAT",
#         "ATTCAAGTCGTGGTAACAGTTTGTAGTGACGCCGTCAATGAGTTAGCGACGCCCTTAAAGAAGAAAGCGATTCATAGCGCAAGGCCATTAGGGCGCCCTATCATCAATCTCTTATACTCTTGTTTTCCTGCGAAGGCGAACCCCGTGAGTACGGCACCGAGGTTTAGCTAGTATGGTCAATGATTCGCGAGAACATCGGTAG",
#         "TTGAGTCGCCGACTGGCTAACAAACAGTAACGTAGCTCATGTGTTTGTCCAGAAAATAGTCACCGTACGTCGCTCATGCAACAATGTTGTGGTGATGGCAGTGGTCGACAATGTGCAGCGGCACCGCATCTTTGACTGGCGGTGATTACGCCCATTCTAGGGAGGGGTCGTACAAGCGGCTATGATCATCATGCCCGCTGCG",
#         "ATCAATGTGAACCCAAACTACGGACCTGTACCTACTCGAGGAAATAGCACACAGGACCATCACATGAGCACCCGACGCGCCTCACCAACAAGCACTAGTAATGGTAAAAGTTGTACTGGACTTATATCGAGCGTCCTAAACACAAGCGCAAGGCGTGTTTCCCACGGGAGGGCCATAGCCATCAACTTTTCGCCTTCTCGGA",
#         "ATAGGTCCCCCCGGGAATGAGTCGCATTCGCTGCCAAACCTAAATGTCATGACAAATGTGGGGACGGTCATTATGTACAACAACGGATGATCAATGTCTATATACTCAGCTGCTCCGAACGAAGCGCCGTAAAACGTTGGATATGTAGCCTTCTTCGTTTAGGTATTCAGCCCATGATCAGGCTAAAGACGAAAATATGGTT",
#         "GAGGAATCTGCTTCGGTCATAGGAAGGATATCGTGAGTATAGTTGTAAGCAGTGATTCACACGTTACCGATGAGTCGATTACGCGGGAAGCTAGCCGAACAAATAACCGATGGCGGCACCACAGAGTTAACCCAGCAGTCCACCGTACGCATTCTGAACAATGCAGAAGCCGATCACGCTATCAAGCCACGTATCATGTTGG",
#         "ATCTTCTTTCAACTGGGAACTCATTTGTTATCTCTGCTGACTGTGTCACTGCGTCATATCTGATGCAGGAGCTATGGTACGGGTTGTACACATGAGTCGCCGATTATGGTCATTCGGAGATGGTTTCCTCTGGGGACAACACTGGTCTCCCGGGATAGCTGAATGATGGTTTCAGCGGACAATAAGGTTGGATCATGGCGGC"]
# RepeatedRandomizedMotifSearch(Dna, k, t)
"""
The function should divide each value in Probabilities by the sum of all values
in  Probabilities, then return the resulting dictionary
"""

# Input: A dictionary Probabilities, where keys are k-mers and values are the
# probabilities of these k-mers (which do not necessarily sum up to 1)
# Output: A normalized dictionary where the probability of each k-mer was
# divided by the sum of all k-mers' probabilities
def Normalize(Probabilities):
    probability_sum = sum(Probabilities.values())
    for key, value in sorted(Probabilities.items()):
        Probabilities[key] = value / probability_sum
    return Probabilities


# Probabilities = {'A': 0.15, 'B': 0.6, 'C': 0.225, 'D': 0.225, 'E': 0.3}
# print(Normalize(Probabilities))

"""
This function takes a dictionary Probabilities whose keys are k-mers and whose
values are the probabilities of these k-mers. The function should return a
randomly chosen k-mer key with respect to the values in Probabilities
"""

# Input:  A dictionary Probabilities whose keys are k-mers and whose values are the probabilities of these kmers
# Output: A randomly chosen k-mer with respect to the values in Probabilities
def WeightedDie(Probabilities):
    random_float = random.uniform(0,1)
    input_keys = []
    input_values = []
    for key, value in sorted(Probabilities.items()):
        input_keys.append(key)
        input_values.append(value)

    current = 0
    range_of_values = []
    for i, value in enumerate(input_values):
        current += value
        range_of_values.append(current)

    for i, value in enumerate(range_of_values):
        if value >= random_float:
            return(input_keys[i])

# Probabilities = {'AA': 0.2, 'AT': 0.4, 'CC': 0.1, 'GG': 0.1, 'TT': 0.2}
# print(WeightedDie(Probabilities))

"""
Now that we can simulate a weighted die roll over a collection of probabilities
of strings, we need to make this function into a subroutine of a larger function
that randomly chooses a k-mer from a string Text based on a profile matrix profile
"""

# Input:  A string Text, a profile matrix Profile, and an integer k
# Output: ProfileGeneratedString(Text, profile, k)
def ProfileGeneratedString(Text, profile, k):
    n = len(Text)
    probabilities = {}
    for i in range(0,n-k+1):
        probabilities[Text[i:i+k]] = Pr(Text[i:i+k], profile)
    probabilities = Normalize(probabilities)
    return WeightedDie(probabilities)

"""
RandomizedMotifSearch may change all t strings in Motifs in a single iteration.
This strategy may prove reckless, since some correct motifs (captured in Motifs)
may potentially be discarded at the next iteration.
GibbsSampler is a more cautious iterative algorithm that discards a single k-mer
from the current set of motifs at each iteration and decides to either keep it
or replace it with a new one.
"""

def GibbsSampler(Dna, k, t, N):
    # randomly select k-mers Motifs = (Motif1, …, Motift) in each string from Dna
    random_motifs = RandomMotifs(Dna, k, t)
    # ﻿BestMotifs ← Motifs
    best_motifs = random_motifs
    # for j ← 1 to N
    for j in range(N):
    #     i ← randomly generated integer between 1 and t
        i = random.randint(0,t)
    #     Profile ← profile matrix formed from all strings in random_motifs except for random_motifs
        profile = ProfileWithPseudocounts(random_motifs)
    #     Motif ← Profile-randomly generated k-mer in the i-th string
        check_random = Motifs(profile, Dna, k)
        if Score(check_random) > Score(best_motifs):
            best_motifs = check_random
        else:
            return best_motifs

k = 8
t = 5
N = 100

Dna = ["CGCCCCTCTCGGGGGTGTTCAGTAACCGGCCA",
    "GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG",
    "TAGTACCGAGACCGAAAGAAGTATACAGGCGT",
    "TAGATCAAGTTTCAGGTGCACGTCGGTGAACC",
    "AATCCACCAGCTCCACGTGCAATGTTGGCCTA"]

output=GibbsSampler(Dna, k, t, N)

with open('outputfile2.txt', 'w') as f:
        for item in output:
            f.write("%s\n" % item)