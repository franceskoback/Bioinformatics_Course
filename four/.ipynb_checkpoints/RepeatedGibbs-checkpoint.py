# first, import the random package
import random
import copy


# Write a GibbsSampler() function here along with any subroutines that you need.
# GibbsSampler() should return a list of strings.
def GibbsSampler(Dna, k, t, N):
    # your code here
    motifs= []
    best_motifs = []
    for i in range(len(Dna)):
        rd = random.randint(0,len(Dna[0]) - k)
        motifs.append(Dna[i][rd:rd + k])
        best_motifs.append(Dna[i][rd:rd + k])
    for i in range(N):
        rd = random.randint(0, t - 1)


        motifs.pop(rd)
        pro_mat = ProfileMatrix_pseudo(motifs, 1)
        probs = []
        cumulative_prob = 0
        for j in range(len(Dna[0])-k+1):
            pattern = Dna[rd][j:j+k]
            prob = 1
            for z in range(len(pattern)):
                prob *= pro_mat[pattern[z]][z]
            cumulative_prob += prob
            probs.append([prob, pattern])
        for j in range(len(probs)):
            probs[j][0] = probs[j][0] / cumulative_prob
        probs = sorted(probs, key = lambda x:x[0])
        for j in range(1,len(probs)):
            probs[j][0] = probs[j-1][0] + probs[j][0]
        rand = random.uniform(0,1)
        #print("len probs: ", len(probs))
        #print(probs)
        #print("random numb: " +str(rand))
        ##print("ith : " + str(rd))
        #print(motifs)
        for j in range(len(probs)):
            if rand < probs[j][0]:
                motifs.insert(rd, probs[j][1])
                #print(probs[j][1])
                #print(motifs)
                #print("add")
                #motifs = motifs[:rd]+[probs[j][1]] + motifs[rd:]
                #print("motif length ", len(motifs))
                #print(motifs)
                break
        #print(motifs)
        if Score(motifs) < Score(best_motifs):
            #print("changed best motif")
            best_motifs = copy.deepcopy(motifs)
    
    #print(Score(['AGACG','AGACG','ACACG','ACACG','AGACG','AGATG','AGAGG','ATACG','ACACT','GGAGG']))
    return best_motifs
                      
    

# Fill in your RepeatedGibbsSampler() function here. You should simply call GibbsSampler() n times and return the best Motifs. This function should return a list of strings.
def RepeatedGibbsSampler(Dna, k, t, N, n):
    # your code here
    best_motifs= []
    for i in range(len(Dna)):
        rd = random.randint(0,len(Dna[0]) - k)
        best_motifs.append(Dna[i][rd:rd + k])
    for i in range(n):
        good_motifs = GibbsSampler(Dna,k,t,N)
        print("curr score:", Score(good_motifs))
        print("best score:", Score(best_motifs))
        if Score(good_motifs) < Score(best_motifs):
            best_motifs = copy.deepcopy(good_motifs)
        
    return best_motifs


def ProfileMatrix(dna):
    d = {'A':0, 'C':1, 'G':2, 'T':3}
    #res = [[0 for y in range(len(dna[0]))] for x in range(4)]
    pro_mat = {'A':[0 for x in range(len(dna[0]))],'C':[0 for x in range(len(dna[0]))],'G':[0 for x in range(len(dna[0]))],'T':[0 for x in range(len(dna[0]))]}
    
    
    for i in range(len(dna)):
        for j in range(len(dna[0])):
            #res[d[dna[i][j]]][j] += (1/len(dna))
            pro_mat[dna[i][j]][j] += (1/len(dna))
    return pro_mat     

def ProfileMatrix_pseudo(dna, pseudocount):
    pro_mat = {'A':[pseudocount for x in range(len(dna[0]))],'C':[pseudocount for x in range(len(dna[0]))],'G':[pseudocount for x in range(len(dna[0]))],'T':[pseudocount for x in range(len(dna[0]))]}
    
    
    for i in range(len(dna)):
        for j in range(len(dna[0])):
            pro_mat[dna[i][j]][j] += 1
    for nucleotide in pro_mat:
        for j in range(len(pro_mat[nucleotide])):
            pro_mat[nucleotide][j] = pro_mat[nucleotide][j] / (len(dna) + pseudocount * len(dna))
    return pro_mat


def ProfileMostProbableKmer(text, k, profile):
    
    maximum = -1
    res = ""
    for i in range(len(text) - k + 1):
        pattern = text[i:i + k]
        prob = 1
        for j in range(k):
            prob *= profile[pattern[j]][j]
        if prob > maximum:
            maximum = prob
            res = pattern
    return res
    
def Score(motifs):
    count = 0
    consensus = Consensus(motifs)
    for i in range(len(motifs)):
        count += HammingDistance(consensus, motifs[i])
    return count
    
def Consensus(motifs):
    res = ""
    
    for col in range(len(motifs[0])):
        d = {'A':0, 'C':0, 'G':0, 'T':0}
        for row in range(len(motifs)):
            d[motifs[row][col]] += 1
        most = max(d, key = lambda i:d[i])
        res += most
    #print("consensus: " + str(res))
    return res

def HammingDistance(p, q):

    # your code here

    ans = 0
    for i in range(len(p)):
        if p[i] != q[i]:
            ans += 1
    return ans

def read_file(input_file):
    '''
    >>> k,t,N,Dna = read_file('untitled.txt')
    '''
    f = open(input_file)
    data = [item.strip() for item in f.readlines()]
    k,t,N = map(int,data[0].split(' '))
    f.close()
    return (k,t,N,data[1:])
    
#print(GibbsSampler(['TTGACTAAGAAGACGCCTATGGGGTGGCTCATAGCTAGACAACTTTGAC','TAAGAAGACGCCTATGGGGTGGCTCATAGCTAGACAACGGATTTTTGAC','TCAGCCAGACCAGGCGGAAGAGTGGGTCTGACATTCGCTTGCACACGCA','AGCCCTGCATGTCGCGTTACCCGGGACACGCTCCAACGACTGCCAACCG','AGAAACTTTAGCGCTCGCAACTACCACCGGGGCCAGACGGGTTGCGGCT','GGGACTGTAAGATGTCAGCATGCCGTAATATGTTGATTTTTCGATGGCA','CTAGAGGTCAGGCAGTCTACACCGGAACAAGTCCCTGTGCACTGCGAGA','GTGGATACGGGATCGGCTCGTCAGGATTTGAATCTAATGTGATCCCACT','ATCTCTTAAACGGCACGTTTACCGAACACTGGCTTTTGAGCGTGGGAAT','ATCTTGGAGGGCTCAACGCCGCCGATGTGCGGGCCTGGAACTCATTCCC'], 5,10,1500))

Dna1=["CGCCCCTCTCGGGGGTGTTCAGTAACCGGCCA","GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG",
"TAGTACCGAGACCGAAAGAAGTATACAGGCGT","TAGATCAAGTTTCAGGTGCACGTCGGTGAACC",
"AATCCACCAGCTCCACGTGCAATGTTGGCCTA"]
#print(RepeatedGibbsSampler(['TTGACTAAGAAGACGCCTATGGGGTGGCTCATAGCTAGACAACTTTGAC','TAAGAAGACGCCTATGGGGTGGCTCATAGCTAGACAACGGATTTTTGAC','TCAGCCAGACCAGGCGGAAGAGTGGGTCTGACATTCGCTTGCACACGCA','AGCCCTGCATGTCGCGTTACCCGGGACACGCTCCAACGACTGCCAACCG','AGAAACTTTAGCGCTCGCAACTACCACCGGGGCCAGACGGGTTGCGGCT','GGGACTGTAAGATGTCAGCATGCCGTAATATGTTGATTTTTCGATGGCA','CTAGAGGTCAGGCAGTCTACACCGGAACAAGTCCCTGTGCACTGCGAGA','GTGGATACGGGATCGGCTCGTCAGGATTTGAATCTAATGTGATCCCACT','ATCTCTTAAACGGCACGTTTACCGAACACTGGCTTTTGAGCGTGGGAAT','ATCTTGGAGGGCTCAACGCCGCCGATGTGCGGGCCTGGAACTCATTCCC'], 5,10,1500, 20))
#print(RepeatedGibbsSampler(Dna1,8,5,1500,200))

filename="untitled.txt"

k,t,N,Dna = read_file(filename)
n=20
my_list=RepeatedGibbsSampler(Dna,k,t,N,n)
with open('outputfile2.txt', 'w') as f:
    for item in my_list:
        f.write("%s\n" % item)
print(my_list)
    