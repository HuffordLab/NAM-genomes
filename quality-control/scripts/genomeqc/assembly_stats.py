#!/usr/bin/python

"""
This script calculates all the basic length metrics for the given input genome assembly.

Usage: python assembly_stats.py <genome assembly file (fasta format)> <output file name>  

"""

from Bio import SeqIO
import sys
import statistics
import numpy as np



inputfile = sys.argv[1]
outputfile = sys.argv[2]

sys.stdout=open(outputfile,"w")

records = list(SeqIO.parse(inputfile, "fasta"))

number_of_scaffolds = len(records)

print("number of scaffolds - ",number_of_scaffolds)



len_seq = [len(rec) for rec in records]

total_size_scaffolds = sum(len_seq)

print("total size of sequences - ", total_size_scaffolds)

number_seq_greater_25k = sorted(i for i in len_seq if i>25000)

print("number of sequences greater than 25k - ", sum(number_seq_greater_25k))

#print((sum(number_seq_greater_25k)/(estimated_genome_size*1000000))*100)

print("largest sequence - ", max(len_seq))

print("smallest sequence - ", min(len_seq))

number_seq_greater_1k = len(sorted(i for i in len_seq if i>1000))

print("number of sequences greater than 1k - ", number_seq_greater_1k)

print("percentage of sequences greater than 1k - ", (number_seq_greater_1k/number_of_scaffolds)*100)

number_seq_greater_10k = len(sorted(i for i in len_seq if i>10000))

print("number of sequences greater than 10k - ", number_seq_greater_10k)

print("percentage of sequences greater than 10k - ", (number_seq_greater_10k/number_of_scaffolds)*100)

number_seq_greater_100k = len(sorted(i for i in len_seq if i>100000))

print("number of sequences greater than 100k - ", number_seq_greater_100k)

print("percentage of sequences greater than 100k - ", (number_seq_greater_100k/number_of_scaffolds)*100)

number_seq_greater_1M = len(sorted(i for i in len_seq if i>1000000))

print("number of sequences greater than 1M - ", number_seq_greater_1M)

print("percentage of sequences greater than 1M - ", (number_seq_greater_1M/number_of_scaffolds)*100)

number_seq_greater_10M = len(sorted(i for i in len_seq if i>10000000))

print("number of sequences greater than 10M - ", number_seq_greater_10M)

print("percentage of sequences greater than 10M - ", (number_seq_greater_10M/number_of_scaffolds)*100)


#calculates N50 and L50 values
sorted_len = sorted(len_seq, reverse=True)

sum_sorted_length = sum(sorted_len)

testSum = 0

N50 = 0

N50con = 0

L50 = 0

i = 0

for con in sorted_len:

    testSum += con

    N50con += 1

    i += 1

    if sum_sorted_length/2.0 <= testSum:

       N50 = con

       L50 = i

       break



print ("N50 value - ", N50)

print ("L50 value - ", L50)


#calculates A,C,G,T,N percentages
counterA1 = 0
counterA2 = 0

for record in records:

    counterA1 += record.seq.count('A') 
    counterA2 += record.seq.count('a')


print ("%A - ", ((counterA1+counterA2)/total_size_scaffolds)*100)


counterC1 = 0
counterC2 = 0

for record in records:
    
    
    counterC1 += record.seq.count('C')
    counterC2 += record.seq.count('c')


print ("%C - ",((counterC1+counterC2)/total_size_scaffolds)*100)

counterG1 = 0
counterG2 = 0

for record in records:

    counterG1 += record.seq.count('G')
    counterG1 += record.seq.count('g')

print ("%G - ",((counterG1+counterG2)/total_size_scaffolds)*100)



counterT1 = 0
counterT2 = 0

for record in records:

    counterT1 += record.seq.count('T')
    counterT2 += record.seq.count('t')

print ("%T - ",((counterT1+counterT2)/total_size_scaffolds)*100)



counterN = 0

for record in records:

    counterN += record.seq.count('N')


print ("N count - ", counterN)
print ("%N count - ", (counterN/total_size_scaffolds)*100)

sys.stdout.close()
