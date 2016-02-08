# -*- coding: utf-8 -*-
"""
GENEFINDER
@author: Lydia Zuehsow
"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))


def get_complement(nucleotide):
    """ 
    # Returns the complementary nucleotide
    # nucleotide: a nucleotide (A, C, G, or T) represented as a string
    # returns: the complementary nucleotide"""

    if (nucleotide is 'A'):
        return 'T'
    elif (nucleotide is 'C'):
        return 'G'
    elif (nucleotide is 'T'):
        return 'A'
    elif (nucleotide is 'G'):
        return 'C'
    else: return



def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence
        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    reversecomplement = ''
    reversedna = dna[::-1]
    for i in range(0, len(dna)):
        reversecomplement = reversecomplement + get_complement(reversedna[i])
    return reversecomplement

def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.
        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    """

    index = 0
    EndCodon = ['TAG', 'TGA', 'TAA']

    while index < len(dna):
        codon = dna[index:index+3]
        if codon in EndCodon:
            return dna[0:index]
        else:
            index = index + 3
    return dna    

def find_all_ORFs_oneframe(dna, index):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.
        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    """

    StartCodon = 'ATG'
    ORF_list = []

    while index < len(dna):
        codon = dna[index:index+3]
        if codon == StartCodon:
            ORF_list.append(rest_of_ORF(dna[index:len(dna)]))
            index = index + len(dna)
        index = index + 3
    return ORF_list

def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.
        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """

    complete_ORF_list = []

    for n in range(0, 3):
        if find_all_ORFs_oneframe(dna, n) != []:
            complete_ORF_list = complete_ORF_list + find_all_ORFs_oneframe(dna, n)
    return complete_ORF_list

def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.
        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """

    Really_Complete_ORF_List = find_all_ORFs(dna) + find_all_ORFs(get_reverse_complement(dna))
    return Really_Complete_ORF_List

def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """

    Really_Complete_ORF_List = find_all_ORFs_both_strands(dna)
    if not Really_Complete_ORF_List:
        Max_ORF = 'Z'
    else:
        Max_ORF = max(Really_Complete_ORF_List, key = len)
        Max_ORF = str(Max_ORF)
    return Max_ORF
longest_ORF("ATGCGAATGTAGCATCAAA")

def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence
        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF 
        Max_ORF = [3,2] & [1]
        Good test: dna = ['ATGCGAATGTAGCATCAAA','TAAAAAGCGGTGAATCATC']
        Good Test 2: dna = ['AAACGGCGTTTTAAACAGA', 'AAACAGGAGTTACGTATAC']
        """
    Max_ORF_List = []

    for i in range(0, num_trials):
        dna = shuffle_string(dna)
        #print "dna: {}".format(dna)
        Max_ORF_List.append(longest_ORF(dna))
        #print Max_ORF_List
        #print '     '

    Longest_Maximum_ORF = max(Max_ORF_List, key = len)
    if Longest_Maximum_ORF == 'Z':
        #print 'No ORFs were found.'
        return
    else:
        return Longest_Maximum_ORF
    
#longest_ORF_noncoding("ATGCGAATGTAGCATCAAA", 5)

def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).
        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment
        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    index = 0
    amino_acid = ''

    while index < len(dna):
        if index+3 <= len(dna):
            codon = dna[index:index+3]
            amino_acid = amino_acid + aa_table[codon]
            index = index + 3
        else:
            break
    return amino_acid

def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna
        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """

    amino_acid_list = []

    threshold = len(longest_ORF_noncoding(dna,3))
    print threshold

    Really_Complete_ORF_List = find_all_ORFs_both_strands(dna)
    #print Really_Complete_ORF_List

    for i in range(0,len(Really_Complete_ORF_List)):
        temp_ORF = Really_Complete_ORF_List[i]
        #print len(temp_ORF)
        if len(temp_ORF) > threshold:
            amino_acid = coding_strand_to_AA(temp_ORF)
            amino_acid_list.append(amino_acid)
    return amino_acid_list

if __name__ == "__main__":

    from load import load_contigs
    contigs = load_contigs()
    name = contigs[7][0]
    dna = contigs[7][1]
    print name

    cleaned_dna = ''

    for n in dna:
        if n in 'ATCG':
            cleaned_dna += n
        else:
            cleaned_dna += 'G'

    print gene_finder(cleaned_dna)


    # import doctest
    # print doctest.run_docstring_examples(get_complement('A'), globals())
    # print doctest.run_docstring_examples(get_complement('A'), globals())