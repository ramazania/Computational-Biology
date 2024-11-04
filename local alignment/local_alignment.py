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

def local_alignment(seq1, seq2, match, mismatch, gap):
    n = len(seq1)
    m = len(seq2)

    # Initialize the scoring matrix
    score_matrix = np.zeros((n + 1, m + 1))

    # Fill the scoring matrix
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            # calculate the match score
            current_char_seq1 = seq1[i - 1].upper()
            current_char_seq2 = seq2[j - 1].upper()
            # check for match/mismatch
            characters_match = current_char_seq1 == current_char_seq2
            if characters_match:
                match_score = score_matrix[i - 1][j - 1] + match
            else:
                match_score = score_matrix[i - 1][j - 1] + mismatch
            # calculate the delete and insert scores
            delete_score = score_matrix[i - 1][j] + gap
            insert_score = score_matrix[i][j - 1] + gap
            # Update the current cell in the scoring matrix with the maximum of match, delete, and insert scores
            score_matrix[i][j] = max(0, match_score, delete_score, insert_score)

    # Find the maximum score and its position
    max_score = 0
    max_position = (0, 0)
    for i in range(len(score_matrix)):
        for j in range(len(score_matrix[0])):
            if score_matrix[i][j] > max_score:
                max_score = score_matrix[i][j]
                max_position = (i, j)
   
    # Traceback to find the alignment
    alignment_seq1 = ""
    alignment_seq2 = ""
    i, j = max_position
    while i > 0 and j > 0 and score_matrix[i][j] != 0:
        # Gap in sequence 2
        if score_matrix[i][j] == score_matrix[i - 1][j] + gap:
            alignment_seq1 = seq1[i - 1] + alignment_seq1
            alignment_seq2 = "-" + alignment_seq2
            i -= 1
        # Gap in sequence 1
        elif score_matrix[i][j] == score_matrix[i][j - 1] + gap:
            alignment_seq1 = "-" + alignment_seq1
            alignment_seq2 = seq2[j - 1] + alignment_seq2
            j -= 1
        # Match or mismatch
        else:
            alignment_seq1 = seq1[i - 1] + alignment_seq1
            alignment_seq2 = seq2[j - 1] + alignment_seq2
            i -= 1
            j -= 1
            
    alignment_score = max_score
    return alignment_score, alignment_seq1, alignment_seq2

def main():
    # Get command line arguments
    fasta_file1 = sys.argv[1]
    fasta_file2 = sys.argv[2]
    scoring_file = sys.argv[3]

    # Read sequences from FASTA files
    sequence1, sequence2 = read_fasta(fasta_file1), read_fasta(fasta_file2)

    # Read scoring parameters from scoring file
    match, mismatch, gap = read_scoring_file(scoring_file)

    # Perform local alignment
    alignment_score, alignment_seq1, alignment_seq2 = local_alignment(sequence1[0], sequence2[0], match, mismatch, gap)

    # Print the result
    print("Alignment Score:", alignment_score)
    print("Optimal Local Alignment:")
    print(alignment_seq1)
    print(alignment_seq2)

if __name__ == "__main__":
    main()
