def PatternCount(Text, Pattern):
    count=0
    pattern_length=len(Pattern)
    for i in range(0,len(Text)-len(Pattern)+1):
        if Text[i:i+pattern_length] == Pattern:
            count=count+1
    return count

print(PatternCount("ACTGTACGATGATGTGTGTCAAAG", "TGT"))

def FrequentWords(Text, k):
        FrequentPatterns=[]
        Count=[]
        for i in range(0,(len(Text)-k)):
            Pattern= Text[i:i+k]
            Count.append(PatternCount(Text, Pattern))
        maxCount=max(Count)
        for i in range(0,(len(Text)-k)):
            if Count[i] == maxCount and not Text[i:i+k] in FrequentPatterns:
                FrequentPatterns.append(Text[i:i+k])
        return FrequentPatterns

print(FrequentWords("CGGAGGACTCTAGGTAACGCTTATCAGGTCCATAGGACATTCA", 3))
    
def ReverseComplement(Pattern):
    length_Pattern=len(Pattern)
    Complement=""
    for i in range(length_Pattern-1,-1,-1):
        if Pattern[i]=="A":
            Complement=Complement+"T"
        elif Pattern[i]=="T":
            Complement=Complement+"A"
        elif Pattern[i]=="C":
            Complement=Complement+"G"
        elif Pattern[i]=="G":
            Complement=Complement+"C"
    return Complement
    
print(ReverseComplement("CCAGATC"))

#ATAT
#GATATATGCATATACTT

# fill in your PatternMatching() function along with any subroutines that you need.
def PatternMatching(Pattern, Genome):
    index_list=[]
    pattern_length=len(Pattern)
    for i in range(0,len(Genome)-len(Pattern)+1):
        if Genome[i:i+pattern_length] == Pattern:
            index_list.append(i)
    return index_list

mydata=open("Vibrio_cholerae_genome.txt")
genome=mydata.read()
data=PatternMatching("CGC","ATGACTTCGCTGTTACGCGC")
s=' '.join([str(pos) for pos in data])
print(s)


from collections import defaultdict

def ClumpFinding(genome, k, L, t):
    dictionary = defaultdict(list)
    output = set()

    for i in range(0,len(genome) - k + 1):
        portion = genome[i:i + k]
        
        while dictionary[portion] and k-dictionary[portion][0]+i> L:
            dictionary[portion].pop(0)

        dictionary[portion].append(i)
        if len(dictionary[portion]) == t:
            output.add(portion)
    return output


#data=ClumpFinding("CGGACTCGACAGATGTGAAGAACGACAATGTGAAGACTCGACACGACAGAGTGAAGAGAAGAGGAAACATTGTAA",5,50,4)



mydata=open("E_coli.txt")
genome=mydata.read()
data=ClumpFinding(genome,9,500,3)
print(len(data))
#s=' '.join([str(pos) for pos in data])
#print(s)

def PatternToNumber(pattern):
    int_string=""
    for i in pattern:
        if i=="A":
            int_string=int_string+"0"
        elif i=="C":
            int_string=int_string+"1"
        elif i=="G":
            int_string=int_string+"2"
        elif i=="T":
            int_string=int_string+"3"
    return(int(int_string,4))
            
#print(PatternToNumber("ATGCAA"))
import numpy as np
def ComputingFrequencies(Text, k):
    FrequencyArray=[]
    for i in range(0,4**k):
        FrequencyArray.append(0)
    for i in range(0,len(Text)-k+1):
        Pattern=Text[i:i+k]
        j=PatternToNumber(Pattern)
        FrequencyArray[j]=FrequencyArray[j] + 1
    returnstring=""
    for i in FrequencyArray:
        returnstring=returnstring+str(i)
    return returnstring

names=ComputingFrequencies("ACGCGGCTCTGAAA",2)

#print(names)

def PatternToNumber(Pattern):
        if not Pattern:
            return 0
        symbol=Pattern[-1]
        Prefix=Pattern[0:len(Pattern)-1]
        return 4*PatternToNumber(Prefix) + SymbolToNumber(symbol)
    
def SymbolToNumber(symbol):
    if symbol=="A":
        Number=0
    elif symbol=="C":
        Number=1
    elif symbol=="G":
        Number=2
    elif symbol=="T":
        Number=3
    return Number 

print(PatternToNumber("ACCGGAGCAACCCTCT"))

def NumberToPattern(index, k):
    if k ==1:
        return NumberToSymbol(index)
    prefixIndex= int(index/4)
    r= index % 4
    symbol=NumberToSymbol(r)
    PrefixPattern= NumberToPattern(prefixIndex,k-1)
    return PrefixPattern + symbol

def NumberToSymbol(Number):
    symbol=""
    if Number==0:
        symbol="A"
    elif Number==1:
        symbol="C"
    elif Number==2:
        symbol="G"
    elif Number==3:
        symbol="T"
    return symbol

print(NumberToPattern(7794,8))

