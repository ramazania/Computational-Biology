"""
Ali Ramazani
Winter 2024
CS 362 - Computational Biology
Final Project: Hirschberg's Algorithm

This file contains the implementation of Hirschberg's algorithm for global sequence alignment.
It takes two sequences and returns the optimal alignment of the two sequences in linear space.
"""

from Bio import SeqIO
import numpy as np
import sys

def read_fasta(file_path):
    sequences = []
    with open(file_path, "r") as file:
        for record in SeqIO.parse(file, "fasta"):
            sequences.append(str(record.seq))
    return sequences

def read_scoring_file(file_path):
    with open(file_path, "r") as file:
        # skip the first line
        file.readline()
        # read all lines and store them in a list
        lines = file.readlines()
        # extract match, mismatch, and gap from the list
        for line in lines:
            line = line.strip().split(',')
            match, mismatch, gap = int(line[0]), int(line[1]), int(line[2]) 
        return match, mismatch, gap

#This function returns the last column of the Needleman-Wunsch score matrix in linear space
def linear_space_alignment_score(s1, s2, match=2, mismatch=-1, gap=-2):
    n, m = len(s1), len(s2)
    prev_score = np.zeros(m+1)
    curr_score = np.zeros(m+1)
    for j in range(m+1):
        prev_score[j] = j * gap
    for i in range(1, n+1):
        curr_score[0] = i * gap
        for j in range(1, m+1):
            character_match = s1[i-1] == s2[j-1]
            if character_match:
                score = match
            else:
                score = mismatch
            curr_score[j] = max(prev_score[j-1] + score, prev_score[j] + gap, curr_score[j-1] + gap)
        prev_score = curr_score.copy()
    return curr_score


def hirschberg_alignment(seq1, seq2, match=2, mismatch=-1, gap=-2):
    n = len(seq1)
    m = len(seq2)
    middle = n // 2
    optimal_nodes = []
    seq1_left = seq1[:middle]
    left_score = linear_space_alignment_score(seq1_left, seq2)
    print(left_score)
    seq1_right = seq1[middle:]
    right_score = linear_space_alignment_score(seq1_right[::-1], seq2[::-1])
    print(right_score[::-1])
    sum_of_score = left_score + right_score[::-1]
    print(sum_of_score)
    max_score = max(sum_of_score)
    optimal_node = np.argmax(sum_of_score)
    optimal_nodes.append((middle, max_score, optimal_node))
    
    return optimal_nodes
    


seq1 = "AGTACGCA"
seq2 = "TATGC"

#print(linear_space_alignment_score(seq2, seq1))
print(hirschberg_alignment(seq1, seq2))


