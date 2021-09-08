def skew(dna):
    output=[]
    output.append(0)
    for i in dna:
        if i=="C":
            output.append(-1+output[-1])
        elif i=="G":
            output.append(1+output[-1])
        elif i=="A" or i=="T":
            output.append(output[-1])
    return output

#output=skew("GAGCCACCGCGATA")
#s=' '.join([str(pos) for pos in output])
#print(s)



def MinimumSkew(Genome):
    skew_data=skew(Genome)
    result_string=""
    for i in range(1,len(skew_data)):
        if skew_data[i]==min(skew_data):
            result_string=result_string+" "+str(i)
    return result_string[1:]
    
#skewtwo MUCH Faster than the combination of minimum skew and skew
def skewtwo(dna):
    output=[]
    output.append(0)
    minsofar=0
    skewindex=[]
    count=0
    for i in dna:
        count +=1 
        if i=="C":
            output.append(-1+output[-1])
        elif i=="G":
            output.append(1+output[-1])
        elif i=="A" or i=="T":
            output.append(output[-1])
        if output[-1]<minsofar:
            minsofar=output[-1]
            skewindex=[count]
        elif output[-1]==minsofar:
            skewindex.append(count)
    return skewindex

    
    
    
output=skewtwo("GCATACACTTCCCAGTAGGTACTG")
#mydata=open("dataset_7_6.txt")
#genome=mydata.read()
#output=skewtwo(genome)
print(output)
#s=' '.join([str(pos) for pos in output])
#print(s)

def HammingDistance(p, q):
    count=0
    if len(p) != len(q):
        return
    for i in range(0,len(p)):
        if p[i] != q[i]:
            count+=1
    return count

    
#output=HammingDistance("CTTGAAGTGGACCTCTAGTTCCTCTACAAAGAACAGGTTGACCTGTCGCGAAG","ATGCCTTACCTAGATGCAATGACGGACGTATTCCTTTTGCCTCAACGGCTCCT")
print(output)

# Input:  Strings Pattern and Text along with an integer d
# Output: A list containing all starting positions where Pattern appears
# as a substring of Text with at most d mismatches
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

# Insert your Hamming distance function on the following line.

output=ApproximatePatternCount("CATGCCATTCGCATTGTCCCAGTGA","CCC",2)

print(output)
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
            
def FrequentWordsWithMismatches(Text, k, d):
    FrequencyArray=[]
    CodeList={}
    for i in range(0,4**k):
        CodeList[i]=0
    for i in range(0,len(Text)-k+1):
        NeighborsList=Neighbors(Text[i:i+k],d)
        for pattern in NeighborsList:
            number=PatternToNumber(pattern)
            CodeList[number]+=1
    maxcount=max(CodeList.values())
    
    for number in range(0,4**k):
        if CodeList[number]==maxcount:
            Pattern=NumberToPattern(number,k)
            if Pattern not in FrequencyArray: 
                FrequencyArray.append(Pattern)
    
    return FrequencyArray


def FrequentWordsWithMismatchesAndReverseComplements(Text, k, d):
    FrequencyArray=[]
    CodeList=[0] * (4**k)
    for i in range(0,len(CodeList)):
        pattern=NumberToPattern(i,k)
        CodeList[i]=ApproximatePatternCount(Text,pattern,d)+ApproximatePatternCount(Text,ReverseComplement(pattern),d)
        maxcount=max(CodeList)
        
    for number in range(0,len(CodeList)):
        if CodeList[number]==maxcount:
            Pattern=NumberToPattern(number,k)
            if Pattern not in FrequencyArray: 
                FrequencyArray.append(Pattern)
    return FrequencyArray
        
                
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



    
output=Neighbors("ACGT",3)

#fh = open("hello.txt","w")
#fh.write('\n'.join(Neighbors("GGATTACA", 2)))
#fh.close()

print(len(output))

#output=FrequentWordsWithMismatches("ACGTTGCATGTCGCATGATGCATGAGAGCT",4,1)
#print(output)

#output = FrequentWordsWithMismatchesAndReverseComplements("TATCTTACCTTACCTTTATGACGGATATTATCTTACACACCGCGACGATATCGGAGAGATATTATACACCTTCTTGAGACTTCTTCGACCTTCGCTTTATCTTGACGTATGACTTCTTACACGAACCGCGACACCGACACCGTATACCGACACACGACGACTATTATACACACACCTTGACTTACGACTTCTTACACGAACCGTATCGACGATATTATCGTATTATTATCTTCTTGATATGAGA",7,3)
#print(output)

#print(ReverseComplementSuperPro("ACGC"))