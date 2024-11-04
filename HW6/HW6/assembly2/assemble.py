from collections import defaultdict
import sys

def construct_de_bruijn_graph(reads, k):
    # Create a defaultdict to store the De Bruijn graph. Dictionary with values as list:
    graph = defaultdict(list)
    
    # Iterate through each read in the list of reads
    for read in reads:
        # Generate k-mers and add edges to the graph
        for i in range(len(read) - k + 1):
            kmer = read[i:i+k-1]
            next_kmer = read[i+1:i+k]
            # Add the next k-mer as an outgoing edge from the current k-mer
            graph[kmer].append(next_kmer)
    
    return graph

def assemble_contigs(graph):
    # Initialize a list to store assembled contigs
    contigs = []
    # Initialize a set to keep track of visited nodes during DFS traversal
    visited = set()

    def dfs(node, path):
        # Extend the current path with the last character of the node
        path += node[-1]
        visited.add(node)

        # If it's a leaf node, add the path as a contig
        if len(graph[node]) == 0:
            contigs.append(path)
            return

        # Recursively explore neighbors
        for neighbor in graph[node]:
            if neighbor not in visited:
                dfs(neighbor, path)

    # Perform DFS traversal for each unvisited node in the graph
    for node in graph:
        if node not in visited:
            dfs(node, "")

    return contigs

def calculate_N50_L50(contigs):
    # Sort contigs by length in descending order
    contigs.sort(key=len, reverse=True)
    
    # Total length of all contigs
    total_length = sum(len(contig) for contig in contigs)
    # Half of the total length
    half_length = total_length / 2
    
    running_sum = 0
    n50 = 0
    l50 = 0

    # Iterate through sorted contigs to calculate N50 and L50
    for contig in contigs:
        running_sum += len(contig)

        if running_sum >= half_length and n50 == 0:
            # Length of the contig at N50
            n50 = len(contig)
            # Number of contigs at L50
            l50 = contigs.index(contig) + 1

    return n50, l50

def write_contigs_to_file(contigs, output_file):
    # Write assembled contigs to a file
    with open(output_file, 'w') as file:
        for contig in contigs:
            file.write(contig + '\n')

if __name__ == "__main__":
    # Reads file provided as a command-line argument
    reads_file = sys.argv[1]
    # k-mer size provided as a command-line argument
    k = int(sys.argv[2])

    # Read DNA sequencing reads from the file
    with open(reads_file, 'r') as file:
        reads = [line.strip() for line in file]

    # Construct the De Bruijn graph
    de_bruijn_graph = construct_de_bruijn_graph(reads, k)

    # Assemble contigs
    contigs = assemble_contigs(de_bruijn_graph)

    # Write contigs to contigs.txt
    write_contigs_to_file(contigs, 'contigs.txt')

    # Calculate and print N50 and L50 scores
    n50, l50 = calculate_N50_L50(contigs)
    print(f"N50 for assembly: {n50}")
    print(f"L50 for assembly: {l50}")
