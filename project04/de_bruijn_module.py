"""De Bruijn Graph Genome Assembly Module.

This module provides classes and functions for constructing De Bruijn graphs
from sequencing reads and performing genome assembly using Eulerian path
traversal.
"""

from collections import defaultdict
import random
import numpy as np


class DeBruijnGraph:
    """Main class for De Bruijn graphs and genome assembly.

    This class builds De Bruijn graphs from sequencing reads and performs
    genome assembly by finding Eulerian paths through connected components
    of the graph.

    Attributes:
        graph (defaultdict): Adjacency list representation of graph edges.
            Keys are (k-1)-mers, values are lists of adjacent (k-1)-mers.
        k (int): The k-mer size used for graph construction.

    Example:
        >>> reads = ["ATGGCGTACG", "GCGTACGTTA", "ACGTTACCAT"]
        >>> dbg = DeBruijnGraph(reads, k=6)
        >>> contigs = dbg.assemble_contigs(seed=42)
        >>> len(contigs) > 0
        True
    """

    def __init__(self, reads, k, seed=42):
        """Initialize De Bruijn graph from sequencing reads.

        Args:
            reads (list): List of DNA sequence strings.
            k (int): K-mer size for graph construction.

        Example:
            >>> reads = ["ATGGCG", "GCGTGC", "TGCAAC"]
            >>> dbg = DeBruijnGraph(reads, k=4)
            >>> len(dbg.graph) > 0
            True
        """
        self.graph = defaultdict(list)
        self.k = k
        self.build_graph_from_reads(reads, k)
        self.seed = seed

    def add_edge(self, left, right):
        """Add a directed edge to the graph.

        Args:
            left (str): Source (k-1)-mer node.
            right (str): Destination (k-1)-mer node.

        Example:
            >>> dbg = DeBruijnGraph([], k=4)
            >>> dbg.add_edge("ATG", "TGG")
            >>> "TGG" in dbg.graph["ATG"]
            True
        """
        pass

    def remove_edge(self, left, right):
        """Remove a directed edge from the graph.

        Args:
            left (str): Source (k-1)-mer node.
            right (str): Destination (k-1)-mer node.

        Example:
            >>> dbg = DeBruijnGraph([], k=4)
            >>> dbg.add_edge("ATG", "TGG")
            >>> dbg.remove_edge("ATG", "TGG")
            >>> len(dbg.graph["ATG"])
            0
        """
        return self.graph[left].remove(right)

    def build_graph_from_reads(self, reads, k):
        """Build De Bruijn graph from multiple sequencing reads.

        Extracts all k-mers from all reads and adds edges between
        consecutive (k-1)-mers within each k-mer.

        Args:
            reads (list): List of DNA sequence strings.
            k (int): K-mer length for graph construction.

        Example:
            >>> reads = ["ATGGC", "TGGCA"]
            >>> dbg = DeBruijnGraph([], k=4)
            >>> dbg.build_graph_from_reads(reads, 4)
            >>> "ATG" in dbg.graph
            True
        """
        for read in reads:
            for i in range(len(read) - k + 2):
                
                if i > 0:
                    prev_kmer = current_kmer
                current_kmer = read[i:i+k-1]
                
                if i > 0:
                    self.graph[current_kmer].append(prev_kmer)
        print(self.graph)

    def eulerian_walk(self, node, graph, seed=None):
        """Perform recursive Eulerian walk on a graph component.

        This is a recursive function that follows all edges from a node
        to traverse the graph, building a path in reverse order.

        Args:
            node (str): Current node to traverse from.
            graph (defaultdict): Graph or subgraph to traverse.
            seed (int, optional): Seed for random edge selection.

        Returns:
            list: List of (k-1)-mers traversed (in reverse order).

        Example:
            >>> reads = ["ATGGCG"]
            >>> dbg = DeBruijnGraph(reads, k=4)
            >>> graph_copy = defaultdict(list, dbg.graph)
            >>> tour = dbg.eulerian_walk("ATG", graph_copy, seed=42)
            >>> len(tour) > 0
            True
        """

        # Recursion Method
        # While list is not empty
        if len(graph[node]) == 0:
            return [node] # <- we should talk about this eventually

        # Walk
        #  Randomly choose a value from the list
        next_node = np.random.choice(graph[node])

        # Remove the value from that list
        self.remove_edge(node, next_node)

        # Get that key, and continue to randomly choose a value until we reach an empty list

        x = self.eulerian_walk(next_node, graph, seed=seed)
        x.append(node)
        return x

        # Alternative Approach (Loop method)
        #y = [node]
        #while len(graph[node]) > 0:
        #    next_node = np.random.choice(graph[node])
        #    self.remove_edge(node, next_node)
        #    y.append(next_node)
        #    node = next_node

    def assemble_contigs(self, seed=None):
        """Assemble all contigs from the De Bruijn graph.

        Finds all connected components and generates an Eulerian path
        for each component, producing multiple assembled contigs.

        Args:
            seed (int, optional): Random seed for reproducible assembly.

        Returns:
            list: List of assembled contig sequences (DNA strings).

        Example:
            >>> reads = ["ATGGCGTACG", "GCGTACGTTA", "ACGTTACCAT"]
            >>> dbg = DeBruijnGraph(reads, k=6)
            >>> contigs = dbg.assemble_contigs(seed=42)
            >>> all(isinstance(c, str) for c in contigs)
            True
        """
        contigs = []
        # If all the edges has to be walked once that means our end game, then the values in our default dict needs to be empty
        # while edges are present
        while any(edge for edge in self.graph.values()):
            # Get a starting node
            starting_node = self.get_starting_node()
            # Traverse the path - > contig
            walk = self.eulerian_walk(starting_node, self.graph, self.seed)
            # Convert contig into a sequence
            contig = self.tour_to_sequence(walk)
            # Append it to a list
            contigs.append(contig)

        return contigs


    def tour_to_sequence(self, tour):
        """Convert a tour of (k-1)-mers into a DNA sequence.

        Args:
            tour (list): List of (k-1)-mer strings in order.

        Returns:
            str: Assembled DNA sequence.

        Example:
            >>> dbg = DeBruijnGraph([], k=4)
            >>> tour = ['ATG', 'TGG', 'GGC', 'GCG']
            >>> dbg.tour_to_sequence(tour)
            'ATGGCG'
        """
        seq = [tour[0]]
        for kmer in tour[1:]:
            seq.append(kmer[-1])

        return ''.join(seq)


    def get_assembly_stats(self, contigs):
        """Calculate assembly statistics for assembled contigs.

        Args:
            contigs (list): List of contig sequences.

        Returns:
            dict: Dictionary containing assembly statistics:
                - num_contigs: Total number of contigs
                - total_length: Total assembled sequence length
                - longest_contig: Length of longest contig
                - shortest_contig: Length of shortest contig
                - mean_length: Mean contig length
                - n50: N50 statistic

        Example:
            >>> contigs = ["ATGGCG", "TTTAAA", "CCCCCCCCCC"]
            >>> dbg = DeBruijnGraph([], k=4)
            >>> stats = dbg.get_assembly_stats(contigs)
            >>> stats['num_contigs']
            3
        """
        stats = dict()
        # Number contigs == length of contigs list
        stats['num_contigs'] = len(contigs)

        # Total Assembled seq length == total number of bp generated across contigs
        bp_count = 0
        for contig in contigs:
            bp_count += len(contig)
        stats['total_length'] = bp_count

        # Length of longest contig
        stats['longest_contig'] = max(len(contig) for contig in contigs)

        # Length of shortest contig
        stats['shortest_contig'] = min(len(contig) for contig in contigs)

        # Mean contig length
        stats['mean_length'] = sum(len(contig) for contig in contigs) / len(contigs)

        # If you sort all contigs by size, the N50 is the length of the contig where the sum of lengths starting from the longest equals 50% of the total assembly size.
        sorted_contigs = sorted(contigs, key=len, reverse=True)
        bp_half = round(bp_count / 2, 0)
        sum_length = 0
        for contig in sorted_contigs:
            sum_length += len(contig)
            if sum_length >= bp_half:
                break
        stats['n50'] = sum_length

        return stats


    def write_fasta(self, contigs, filename, organism="Mus musculus"):
        """Write assembled contigs to a FASTA file.

        Args:
            contigs (list): List of contig sequences.
            filename (str): Output FASTA filename.

        Example:
            >>> contigs = ["ATGGCG", "TTTAAA"]
            >>> dbg = DeBruijnGraph([], k=4)
            >>> dbg.write_fasta(contigs, "output.fasta")
        """
        with open(filename, 'w') as outfile:
            for counter, contig in enumerate(contigs):
                header = f'>Contig {counter + 1} [organism={organism}] [moltype=DNA]'
                outfile.write(f'{header}\n{contig}\n')


    def get_starting_node(self):
        # Find a starting node (a node where nothing is being directed to it)
        ## Ensure that a key does not exist as a value (computationally expensive)

        # Could we start randomly and not worry about
        # Cast the list as a set
        # Changing initialization of defaultdict
        # Iterate through the keys in default dict
        # get the values and typecast it to a set
        # {key: set(self.graph[key]) for key in self.graph}

        # Create a dictionary of sets
        graph_set = {}
        for key, value in self.graph.items():
            graph_set[key] = set(value)

        # Initialized a list to hold starting node candidates
        starting_node_candidates = []

        # The starting vertex needs out-degree = in-degree + 1
        for node in self.graph.keys():
            in_degree = 0
            # Out degree is the list of nodes
            out_degree = len(self.graph[node])

            for other_node in graph_set.keys():
                if other_node == node:
                    continue
                else:
                    if node in graph_set[other_node]:
                        in_degree += self.graph[other_node].count(node)
            if out_degree == in_degree + 1:
                starting_node_candidates.append(node)

        return np.random.choice(starting_node_candidates)




