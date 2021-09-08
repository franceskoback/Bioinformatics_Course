import random
from itertools import *

def main():
    print("hello world")
    #Dna1=["CGCCCCTCTCGGGGGTGTTCAGTAACCGGCCA","GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG","TAGTACCGAGACCGAAAGAAGTATACAGGCGT","TAGATCAAGTTTCAGGTGCACGTCGGTGAACC","AATCCACCAGCTCCACGTGCAATGTTGGCCTA"]
    Dna=["ATGAGGTC","GCCCTAGA","AAATAGAT","TTGTGCTA"]


    
    
    #k=15
    #t=20
    #N=2000
    k=3
    t=4
    N=1000
    my_list=RandomizedMotifSearch(Dna, k, t)
    #with open('outputfile2.txt', 'w') as f:
        #for item in my_list:
           #f.write("%s\n" % item)
    print(my_list)
    filename="/Users/franceskoback/Documents/bioinformatics_course/four/untitled.txt"
    #with open(filename) as f:
        #Dna=open(filename).read().split()
        #print(Dna)
    #my_list= RandomizedMotifSearch(Dna1, k, t)
    #with open('outputfile.txt', 'w') as f:
        #for item in my_list:
            #f.write("%s\n" % item)
    #my_list=GibbsSampler(Dna, k, t, N) 
    #with open('outputfile2.txt', 'w') as f:
        #for item in my_list:
            #f.write("%s\n" % item)
    p1 = float ( (600.-15) / (600.-15+1)   )
    p2 = 1 - p1
    counter = 0
    for seq in combinations(range(10),2):
        counter +=1
    p3= pow(p2,2) * pow(p1,8) * counter
    print(p3)
    
    
def GibbsSampler(Dna, k, t, N):
    #CurrentMotifs=RandomMotifs(Dna,k) #randomly select k-mers Motifs = (Motif1, …, Motift) in each string from Dna
    CurrentMotifs=RandomizedMotifSearch(Dna,k,t)
    BestMotifs=CurrentMotifs #BestMotifs ← Motifs
    for j in range(1,N): #for j ← 1 to N
        i=random.randint(0,t) #i ← Random(t)
        Profile= ProfileMatWithPseudoCount(Dna)
            #Profile ← profile matrix constructed from all strings in Motifs except for Motifi
        Motif=Motifs(Profile,Dna,k)#Motifi ← Profile-randomly generated k-mer in the i-th sequence
        if score(Motif,k)<score(BestMotifs,k):#if Score(Motifs) < Score(BestMotifs)
            BestMotifs=Motif #BestMotifs ← Motifs
    return BestMotifs
    
def RandomizedMotifSearch(Dna, k, t):
    count=0
    #CurrentMotifs=RandomMotifs(Dna,k)
    CurrentMotifs=["GTC","CCC","ATA","GCT"]
    BestMotifs=CurrentMotifs
    while count<1:
        CurrentMotifs=RandomMotifs(Dna,k)
        count=count+1# randomly select k-mers Motifs = (Motif1, …, Motift) in each string from Dna
        while True:
            count=count+1
            Profile= ProfileMatWithPseudoCount(CurrentMotifs) #Profile ← Profile(Motifs)
            CurrentMotifs=Motifs(Profile,Dna,k)
            if score(CurrentMotifs,k)<score(BestMotifs,k):
                BestMotifs=CurrentMotifs
            else:
                return BestMotifs
    return BestMotifs
        
def HammingDistance(p, q):
    count=0
    if len(p) != len(q):
        return
    for i in range(0,len(p)):
        if p[i] != q[i]:
            count+=1
    return count

def Motifs(profile,dna,k):
    motifs=[]
    for i in range (len(dna)):
        newmotif=ProfileMostProbableKmer(dna[i],k,profile)
        motifs.append(newmotif)
    return motifs

def RandomMotifs(dna,k):
    length=len(dna)
    motifs=[]
    for i in range(length):
        randomnum=random.randint(0, len(dna[0])-k)
        motifs.append(dna[i][randomnum:randomnum+k])
    return motifs
        
def ProfileMatWithPseudoCount(dna):
    seqNum = float(len(dna))
    nucleotide = ['A', 'C', 'G', 'T']
    matrix = []
    for i in range(len(dna[0])):
        base_i = [seq[i] for seq in dna]
        colProfile = [float(base_i.count(n) + 1)/seqNum + 4 for n in nucleotide]
        matrix.append(colProfile)
    return [list(i) for i in zip(*matrix)]

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

def score(motifs,k):
    score=0
    for i in range(k):
        motif = ''.join([motifs[j][i] for j in range(len(motifs))])
        score += min([HammingDistance(motif, homogeneous*len(motif)) for homogeneous in 'ACGT'])
        #print([HammingDistance(motif, homogeneous*len(motif)) for homogeneous in 'ACGT'])
        #print(score)
    return score

if __name__=="__main__":
    main()