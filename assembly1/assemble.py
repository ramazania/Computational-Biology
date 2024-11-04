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

def count_vertices_and_edges(graph):
    # Count the number of vertices and unique edges in the graph
    num_vertices = len(graph)
    #set(edges) to remove duplicate entries
    num_edges = sum(len(set(edges)) for edges in graph.values())
    return num_vertices, num_edges

def generate_dot_file(graph, dot_filename):
    # Generate a DOT file for visualization of the De Bruijn graph
    with open(dot_filename, 'w') as dot_file:
        dot_file.write('digraph mygraph{\n')
        for vertex, edges in graph.items():
            for edge in edges:
                dot_file.write(f'"{vertex}"->"{edge}"\n')
        dot_file.write('}\n')

if __name__ == "__main__":
    # Read command-line arguments
    reads_file = sys.argv[1]
    k = int(sys.argv[2])

    # Read the list of reads from the input file
    with open(reads_file, 'r') as file:
        reads = []
        for line in file:
            reads.append(line.strip())

    # Construct the De Bruijn graph
    de_bruijn_graph = construct_de_bruijn_graph(reads, k)

    # Count vertices and edges in the De Bruijn graph
    num_vertices, num_edges = count_vertices_and_edges(de_bruijn_graph)

    # Print the results to the screen
    print(f"Number of vertices in graph: {num_vertices}")
    print(f"Number of edges in graph: {num_edges}")

    # Generate DOT file for visualization if the graph has 30 or fewer vertices
    if num_vertices <= 30:
        dot_filename = 'graph.dot'
        generate_dot_file(de_bruijn_graph, dot_filename)
        print(f"DOT file '{dot_filename}' generated for visualization.")
