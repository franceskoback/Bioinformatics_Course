# put your python code here
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 19:34:27 2019

@author: davicin
"""
import random
import copy


def profileMostProbableKmer(text, k, profile): 
    '''
    Find a Profile-most probable k-mer in a string.
    
    :param str text: A pattern of nucleotides
    :param int k: Length of nucleotides in the sequence to find
    :param dict profile: Profile matrix, represented as a dictionary of lists
    :return: A profile-most probable k-mer in text
    :rtype: str  
    '''
    maxPr = float('-inf') # Fix to pick first ocurrence as most probable in w.c.e.
    prMost = ''
    for i in range(len(text) - k + 1):
        pattern = text[i:i+k]
        pr = 1
        for j, nucleotide in enumerate(pattern):
            pr *= profile[nucleotide][j]
            
        if maxPr < pr:
            maxPr = pr
            prMost = pattern
    return prMost

def scoreMotifs(motifs):
    '''
    Find the total score of a certain number of motifs (not entrophy, just differences)
    
    :param list[str] motifs: List of motifs in dna
    :return: Score of differences between motifs and concensus string (difference)
    :rtype: int 
    '''
    
    # New implementation. To avoid redundant logic
    nStr = len(motifs) # total of motifs
    lenMotif = len(motifs[0]) # number of nucleotides in the motifs
    
    score = 0
    for i in range (lenMotif):
        nucleotidesStep = [strand[i] for strand in motifs]
        countMaxNucl = nucleotidesStep.count(mostFrequent(nucleotidesStep))
        score += (nStr - countMaxNucl)
    
    return score
        
def mostFrequent(eList): 
    '''
    Find the most frequent element in a list
    
    :param list[] eList: List of elements of any type
    :return: Most frequent element in list
    :rtype: int,str, idk lol 
    '''
    return max(set(eList), key = eList.count) 


def greedyMotifSearchWithPseudocounts(dna, k, t, pseudocount):
    '''
    A collection bestMotifs of k-mers, one from each string in dna, 
    minimizing Score(bestMotifs) from all motif collections of length t,
    using pseudocounts for the creation of profiles from Motifs.
    
    :param list[str] dna: List of upstream regions of genes
    :param int k: Length of nucleotides in a sequence
    :param int t: Number of strings in dna (can be computed, lol)
    :param int pseudocount: Integer to add for the computation of profiles
    :return: A collection of strings BestMotifs (each string for each dna(i))
    :rtype: list[str]  
    '''
    
    bestMotifs = [text[:k] for text in dna]
    baseStrand = dna[0]
    otherStrands = dna[1:t]
    
    for i in range(len(baseStrand) - k + 1):
        kmer = baseStrand[i:i+k]
        motifs = list()
        motifs.append(kmer) #motif(1)
        for strand in otherStrands:
            profileStep = createProfileFromMotifsWithPseudocounts(motifs,pseudocount)
            motifStep = profileMostProbableKmer(strand, k, profileStep)
            motifs.append(motifStep)
        if scoreMotifs(motifs) < scoreMotifs(bestMotifs):
            bestMotifs = motifs
    return bestMotifs
    
def createProfileFromMotifsWithPseudocounts(motifs,pseudocount):
    '''
    Create dict-type var representing the profile of a certain number of motifs,
    using pseudocounts
    
    :param list[str] motifs: List of motifs in dna
    :param int pseudocount: Integer to add for the computation of profiles
    :return: Profile matrix, represented as a dictionary of lists
    :rtype: dict 
    '''
    
    nStr = len(motifs) + (pseudocount*4) # sum of motifs and pseudocount added(per each nucleotide)
    lenMotif = len(motifs[0]) # number of nucleotides in the motifs
    profileDict = {'A':list(),'C':list(),'G':list(),'T':list()}
    for i in range(lenMotif):
        nucleotidesStep = [strand[i] for strand in motifs]
        auxDict = {'A': 0,'C': 0,'G': 0,'T': 0}
        for nucl in nucleotidesStep:
            auxDict[nucl] += 1
        
        auxDict = {k: (v + pseudocount)/nStr for k, v in auxDict.items()}  
        for k, v in profileDict.items():
            profileDict[k].append(auxDict[k])
    return profileDict




def randomizedMotifSearch(dna, k, t, pseudocount): #t var useless
    '''
    Generate a set of motifs with a random init, in an iteration of profile-motifs 
    until the score cannot be improved.
    :param list[str] dna: List of upstream regions of genes
    :param int k: Length of nucleotides in a sequence
    :param int t: Number of strings in dna
    :param int pseudocount: Integer to add for the computation of profiles
    :return: A collection of strings BestMotifs (each string for each dna(i))
    :rtype: list[str]   
    '''
    
    randomLimit = len(dna[0]) - k
    motifs = list()
    for text in dna:
        randomInit = random.randint(0,randomLimit)
        motif = text[randomInit : randomInit + k]
        motifs.append(motif)
    bestMotifs = copy.deepcopy(motifs)
    while True:
        profile = createProfileFromMotifsWithPseudocounts(motifs, pseudocount) 
        motifs = list()
        for text in dna:
            motif = profileMostProbableKmer(text, k, profile)
            motifs.append(motif)      
        if scoreMotifs(motifs) < scoreMotifs(bestMotifs):
            bestMotifs = motifs
        else:
            return bestMotifs
        

def runThousandRandomMotifSearch(dna, k, t, pse):
    """
    Run the randomizedMotifSearch a thousand times and keep best motifs
    :param list[str] dna: List of upstream regions of genes
    :param int k: Length of nucleotides in a sequence
    :param int t: Number of strings in dna
    :param int pseudocount: Integer to add for the computation of profiles
    :return: A collection of strings BestMotifs (each string for each dna(i))
    :rtype: list[str]
    """
    motifs = randomizedMotifSearch(dna, k, t, pse)
    bestMotifs = motifs
    
    for i in range(1000):
        motifs = randomizedMotifSearch(dna, k, t, pse)
        if scoreMotifs(motifs) < scoreMotifs(bestMotifs):
                bestMotifs = motifs
                
    return bestMotifs

# Thx Chitrasen Mohanty
if __name__ == "__main__":
    k,t = [int(a) for a in input().strip().split(" ")]
    dna = []
    for _ in range(t):
        dna.append(input())
        
    ans = runThousandRandomMotifSearch(dna,k,t,1) #Hardcoded pseudocount
    for a in ans:
        print(a)
