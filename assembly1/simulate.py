import random
import sys

def simulate_dna_sequencing(fasta_file, coverage, read_length, error_rate):
    # Read the DNA sequence from the FASTA file
    with open(fasta_file, 'r') as file:
        # skip the first line
        file.readline()
        dna_sequence = file.read().replace('\n', '')
        
    # Calculate the number of reads (N) using Lander-Waterman statistics
    num_reads = int((coverage * len(dna_sequence)) / read_length)
    
    # Generate simulated reads
    reads = []
    for i in range(num_reads):
        # Randomly choose the starting index in the genome
        start_index = random.randint(0, len(dna_sequence) - read_length)
        read = dna_sequence[start_index:start_index + read_length]
        # Introduce errors based on the specified error rate
        for i in range(read_length):
            if random.uniform(0, 1) < error_rate:
                # Flip a coin to choose a different nucleotide
                nucleotides = ['A', 'C', 'G', 'T']
                nucleotides.remove(read[i])
                read = read[:i] + random.choice(nucleotides) + read[i+1:]
        reads.append(read)
    return reads

def write_reads_to_file(reads, output_file):
    # Write the simulated reads to a file
    with open(output_file, 'w') as file:
        for read in reads:
            file.write(read + '\n')

if __name__ == "__main__":
    # Read command-line arguments
    fasta_file = sys.argv[1]
    coverage = int(sys.argv[2])
    read_length = int(sys.argv[3])
    error_rate = float(sys.argv[4])

    simulated_reads = simulate_dna_sequencing(fasta_file, coverage, read_length, error_rate)
    write_reads_to_file(simulated_reads, 'reads.txt')
    #write_reads_to_file(simulated_reads, 'sample_c12_r_50_e0.00.txt')
    #write_reads_to_file(simulated_reads, 'sample_c12_r_50_e0.01.txt')
