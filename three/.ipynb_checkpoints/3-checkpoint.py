import numpy as np
from numpy import inf

def main():
    print("hello world")
    Motifs1 = [
"TCGGGGGTTTTT",
"CCGGTGACTTAC",
"ACGGGGATTTTC",
"TTGGGGACTTTT",
"AAGGGGACTTCC",
"TTGGGGACTTCC",
"TCGGGGATTCAT",
"TCGGGGATTCCT",
"TAGGGGAACTAC",
"TCGGGTATAACC"
]
    k=5
    d=1
    inputlist=['ACAACATTAAGGTTATATGCCATAC',
'CATCCAAAGCACACCAAGGATCTAA',
'CTCATCATACTCCGCATCACGTACA',
'CATTCTCCACAGATTCGCTTACCTA',
'AATAGAAGGTCGCTGCATTCCAATC',
'CCTTCGCCTCCATTCTCTTCGGCAG']

    #patterns = MotifEnumeration(inputlist, k, d)
    #len(patterns)
    #print(patterns)

    #print(Entropy(Motifs1))
    
    #Dna=['CCGGGGAGGCGTGATCTAATGGTAGTGAGGGGACCAGTGGTA','GTTCTGCGCGAGGTGATGGAAAAATACCATGTGAAGTGTTGC','GTGACGAATCTGTACCATTATCGGGCCAGGATTAATGGTTCT','CACTAAACGATAGTGATGACGTCGTGTCTGTTCCCCATGTGG','ATTGATTTCGTAGGCTGTACTTCTGCGCATGTGACGGGGGAG','TTTAGGGCGAGAATATTCAGGGGGACCAACTCAATAGTGAGG','TATGCGATGCGTATTAGACCCAGCACATCGGTGACGGAACAC','TGCACGGGTGTGGTAGAGGTGACGCTTGGATATTGCAACTCG','GTGACGGTTAGCAGCGAGGCTCGTAGTTTATCAGTACAAAAT','ACGAGAGAACACGTGACGCGCCTACATTACGACATCTTCCGG']
    #k=6;
    text="ACTCGCGATGTCGACGGCTGTGTAGTATGTAGTAAAGTTTAAACACTCCCCTTAACTATTGTAACATTGACCAGCCGCGCCTGTGGGGGTCTTCACTAGAACCCCTCGTTCAAGAAGAGACTTAGAGATTTCTGTAATCATTGGCGCTGTAAAGCGTGCTGTCGCCTCGACGTAAGGATTAAGTAAACGAGTGTTCCCCGCTACTACGATCAAAGCGAGTACGGAGCAAAATGACTGCGCGGGCTTGGATCACTCAAAAGGGCGCCACAGGGATCGTTGGGTTGCAGTGGCTAGTTTAGTCGCAGTAGCGACCCGATAGTGTGAGCATTTACGGACTCTCCTTCCTAAGTTAATCCGCATGATTCCATGGTTGAATCATTGCTGATCGGTGTCATGGATCCGTACGTTAGATGAGTTTGGGAAAGCCACCTCTATTAGACGGGCAAAGACTAAATTGCTTTGAACGCATAGGTGCGATAATGGCGCTTTGGTCCAGCCTCTAAACGTGGATACGGTTGCGGGGCTGCCCGAAGAGTAGAGGAGGATGAATGGAGTGGTGCATGGCCGTGATTTTTTATGCAGTCGAGTTTTATTAATTACGGAAAGGCCTCTCGCCAAAGAATATATGCCCCGTTGGCAGGGGTTAGGAACTCCCCACATCAGAGTTGAGCTGTCGCTATTCGTTTAGTCTCACCTACTCTCCGTCTCTCTATCTCGAAGGTTAGCAGTACTTATCTGAGGGGTGATTTAATTGACCTTGGGGCAGATCTACGATCTCAAGCAGAGCAGTGGAGTTTTATTAAGCAGTGCCTCGCATGATATTAGAGAAGTAATTCCGTCACGTGTTTCTCATTAGCCCCCCCTTTTCAATGGCTAGCGACCGGAATTCTACCAGCCGGAATGACTTATGTGATAACGGTCCCCCTCAAAGGAGCGTGGCGCGTAGGATTGTTTTGTTAAGGCATGTACTAGACGGTTGGCTCCAAGCCAGGGGGCCTGA"
    #k=13;
    #profile= [[0.25, 0.276, 0.355, 0.263, 0.25, 0.211, 0.25, 0.224, 0.316, 0.276, 0.303, 0.224, 0.237],[0.25, 0.276, 0.25, 0.289, 0.25, 0.289, 0.263, 0.303, 0.237, 0.263, 0.276, 0.224, 0.276],[0.263, 0.224, 0.184, 0.211, 0.25, 0.184, 0.276, 0.184, 0.263, 0.25, 0.211, 0.342, 0.263],[0.237, 0.224, 0.211, 0.237, 0.25, 0.316, 0.211, 0.289, 0.184, 0.211, 0.211, 0.211, 0.224]] 
    #print(ProfileMostProbableKmer(text, k, profile))
    dna=["GGCGTTCAGGCA","AAGAATCAGTCA","CAAGGAGTTCGC","CACGTCAATCAC","CAATAATATTCG"]
    #print(len(dna))
    k=3
    t=5
    #print(GreedyMotifSearch(dna, k, t))
    mat = ProfileMat(['GGC'])
    #print(FindPortableKmer('GGCGTTCAGGCA', 3,mat))
    print(MedianString(['CTCGATGAGTAGGAAAGTAGTTTCACTGGGCGAACCACCCCGGCGCTAATCCTAGTGCCC','GCAATCCTACCCGAGGCCACATATCAGTAGGAACTAGAACCACCACGGGTGGCTAGTTTC','GGTGTTGAACCACGGGGTTAGTTTCATCTATTGTAGGAATCGGCTTCAAATCCTACACAG'],7))
          #filename="/Users/franceskoback/Documents/bioinformatics_course/three/textfile2.txt"
    #Pattern='CCTAC'
    #with open(filename) as f:
       # Dna = open(filename).read().split()
        #print(Dna)
    
    #print(Dna)
    #print(DistanceBetweenPatternAndStrings(Pattern, Dna))
        #res = GreedyMotifSearch(dnaSet, 12, 25)
       # for i in res:
        #    print(i)
    #print([HammingDistance(motif, homogeneous*len(motif)) for homogeneous in 'ACGT'])
    #print(score(['GGC','AAG']))
    #BestMotifs = [seq[0:k] for seq in dna]
    #score(BestMotifs)
    
    
   


    

def MotifEnumeration(Dna, k, d):
        PatternsSet =set()
        tf=True
        for listelem in Dna:
            for i in range(0,len(listelem)-k+1):
                Pattern=listelem[i:i+k]
                Patternstars=Neighbors(Pattern,d)
                for patternstar in Patternstars:
                    tf= True
                    for listelm2 in Dna:
                        count=ApproximatePatternCount(listelm2,patternstar,d)
                        if count==0: 
                            tf=False
                    
                    if tf==True and patternstar not in PatternsSet:
                        PatternsSet.add(patternstar)
        if len(PatternsSet)==0:
            return("nothing")
        return PatternsSet

def HammingDistance(p, q):
    count=0
    if len(p) != len(q):
        return
    for i in range(0,len(p)):
        if p[i] != q[i]:
            count+=1
    return count
    
def Neighbors(Pattern,d):
    nucleotide_options=["A","C","G","T"]
    if d == 0:
        return Pattern
    if len(Pattern) == 1:
        return nucleotide_options
    neighborhood = []
    RecursiveOutput = Neighbors(Pattern[1:],d)
    for portion in RecursiveOutput:
        if HammingDistance(Pattern[1:],portion) < d:
            for n in nucleotide_options:
                neighborhood.append(n+portion)
        else:
            neighborhood.append(Pattern[:1]+portion)
    return neighborhood

def ApproximatePatternMatching(Text, Pattern, d):
    positions = [] # initializing list of positions
    for i in range(0,len(Text)-len(Pattern)+1):
        atthispos=HammingDistance(Text[i:i+len(Pattern)],Pattern)
        if atthispos<=d:
            positions.append(i)
    return positions

def ApproximatePatternCount(Text, Pattern, d):
    output=ApproximatePatternMatching(Text,Pattern,d)
    return len(output)


k = 3
d = 1
inputlist = ['ATTTGGC', 'TGCCTTA', 'CGGTATC', 'GAAAATT']
#ATA ATT GTT TTT




def Entropy(Motifs):
    Profile = {}
    # check that all motifs are the same length
    ListLength = len(Motifs)
    L1 = len(Motifs[1]) # length of the first motif
    print('There are {} motifs of length {}'.format(ListLength, L1))
    for i in range(len(Motifs)):
        if len(Motifs[i]) != L1:
            ShortMotif = Motifs[i]
            ShortMotifLen = len(ShortMotif)
            print('Oops, Motif {} is {} nucleotides instead of {}!'.format(ShortMotif, ShortMotifLen, L1))
            break
        
    # fill all positions with frequency of 0
    for nucleotide in 'ACGT':
        values = [0] * L1
        Profile[nucleotide] = values
        
    # iterate through each position in the motif matrix, counting nucleotide frequencies
    TotalEntropy = 0
    for key, values in Profile.items():
        for Motif in Motifs:
            for i in range(len(Motif)):
                if Motif[i] == key:
                    Profile[key][i] += 1
        
        # convert nucleotide frequencies to probabilities
        for i in range(len(values)):
            Profile[key][i] = Profile[key][i] / float(ListLength)
        
        # calculate total entropy (Sum of (Prob_value * log2 Prob_n))
        import math
        for value in values:
            if value > 0:
                TotalEntropy += abs(value * math.log(value, 2))
            else: continue
            
    return(TotalEntropy)

def MedianString(Dna, k):
    testlist= [];
    for string in Dna:
        for i in range(len(string)-k+1):
            if string[i:i+k] not in testlist:
                testlist.append(string[i:i+k])
    dictionary={}            
    for pattern in testlist:
        totaldist=1000
        for dna in Dna:
            rowdistance=1000
            for i in range(len(dna)-k+1):
                if HammingDistance(pattern,dna[i:i+k])<rowdistance:
                    rowdistance=HammingDistance(pattern,dna[i:i+k])
            totaldist+=rowdistance
        dictionary[pattern]=totaldist
    temp = min(dictionary.values()) 
    res = [key for key in dictionary if dictionary[key] == temp] 
                
    return res

def ProfileMostProbableKmer(text, k, profile):
    dictionary={}
    for i in range(len(text)-k+1):
        pattern= text[i:i+k]
        count=1
        for j in range(len(pattern)):
            if pattern[j]=="A":
                count=count*profile[0][j]
            if pattern[j]=="C":
                count=count*profile[1][j]
            if pattern[j]=="G":
                count=count*profile[2][j]
            if pattern[j]=="T":
                count=count*profile[3][j]
        dictionary[pattern]=count
    temp = max(dictionary.values()) 
    res = [key for key in dictionary if dictionary[key] == temp] 
    #print(dictionary)
    return res[0]

def GreedyMotifSearch(dna, k, t):
    BestMotifs = [seq[0:k] for seq in dna]
    seq = dna[0]
    for i in range(len(seq) - k):
        motifs = [seq[i:i+k]]
        for j in range(1, len(dna)):
            mat = ProfileMatWithPesudoCount(motifs)
            tempMotif = ProfileMostProbableKmer(dna[j], k, mat)
            motifs.append(tempMotif)
        if score(motifs,k) < score(BestMotifs,k):
            BestMotifs = motifs
    return BestMotifs

def ProfileMat(dna):
    seqNum = float(len(dna))
    nucleotide = ['A', 'C', 'G', 'T']
    matrix = []
    for i in range(len(dna[0])):
        base_i = [seq[i] for seq in dna]
        colProfile = [float(base_i.count(n))/seqNum for n in nucleotide]
        matrix.append(colProfile)
    return [list(i) for i in zip(*matrix)]


def ProfileMatWithPesudoCount(dna):
    seqNum = float(len(dna))
    nucleotide = ['A', 'C', 'G', 'T']
    matrix = []
    for i in range(len(dna[0])):
        base_i = [seq[i] for seq in dna]
        colProfile = [float(base_i.count(n) + 1)/seqNum + 4 for n in nucleotide]
        matrix.append(colProfile)
    return [list(i) for i in zip(*matrix)]

def ProfileDict(matrix, j):
    return {'A': matrix[0][j], 'C': matrix[1][j], 'G': matrix[2][j], 'T': matrix[3][j]}

def score(motifs,k):
    score=0
    for i in range(k):
        motif = ''.join([motifs[j][i] for j in range(len(motifs))])
        score += min([HammingDistance(motif, homogeneous*len(motif)) for homogeneous in 'ACGT'])
        #print([HammingDistance(motif, homogeneous*len(motif)) for homogeneous in 'ACGT'])
        #print(score)
    return score
    
def DistanceBetweenPatternAndStrings(Pattern, Dna):
    k = len(Pattern)
    Distance = 0
    for string in  Dna: 
        hammingDistance = float("inf")
        for i in range(0, len(string)-k+1):
            if hammingDistance > HammingDistance(Pattern, string[i:i+k]):
                hammingDistance = HammingDistance(Pattern, string[i:i+k])
        Distance = Distance + hammingDistance
    return Distance
    
if __name__=="__main__":
    main()

    
    