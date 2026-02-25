import time
from datetime import datetime
from de_bruijn_module import DeBruijnGraph
import gzip

def read_fastq(filename, n_rows=10000):
    """Read sequences from a FASTQ file.

    Args:
        filename (str): Path to FASTQ file.

    Returns:
        list: List of DNA sequence strings (quality scores ignored).

    Example:
        >>> reads = read_fastq("test_reads.fastq")
        >>> len(reads)
        10
    """
    sequences = []
    f = gzip.open(filename, 'rt') #use gzip here

    with f:
        line_count = 0
        for line in f:
            line_count += 1
            if line_count % 4 == 2:  # Sequence line in FASTQ format
                sequences.append(line.strip())
            if len(sequences) == n_rows:
                break
    return sequences


def write_statistics_file(
    stats_file,
    input_file,
    num_reads,
    avg_read_length,
    k_mer_size,
    random_seed,
    num_nodes,
    num_edges,
    stats,
    contig_lengths,
    timing,
    coverage_estimate,
    assembly_fraction
):
    """Write comprehensive assembly statistics to text file.
    
    Args:
        stats_file (str): Output filename.
        input_file (str): Input FASTQ filename.
        num_reads (int): Number of reads processed.
        avg_read_length (float): Average read length.
        k_mer_size (int): K-mer size used.
        random_seed (int): Random seed used.
        num_nodes (int): Number of graph nodes.
        num_edges (int): Number of graph edges.
        stats (dict): Assembly statistics.
        contig_lengths (list): Sorted list of contig lengths.
        timing (dict): Timing information.
        coverage_estimate (float): Estimated sequencing coverage.
        assembly_fraction (float): Assembly size as % of genome.
    """
    with open(stats_file, 'w') as f:
        f.write("="*80 + "\n")
        f.write("MOUSE GENOME ASSEMBLY STATISTICS\n")
        f.write("="*80 + "\n\n")
        
        f.write(f"Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Input File: {input_file}\n")
        f.write(f"K-mer Size: {k_mer_size}\n")
        f.write(f"Random Seed: {random_seed}\n\n")
        
        f.write("-"*80 + "\n")
        f.write("INPUT DATA\n")
        f.write("-"*80 + "\n")
        f.write(f"Number of reads:         {num_reads:,}\n")
        f.write(f"Average read length:     {avg_read_length:.1f} bp\n")
        f.write(f"Total sequencing data:   {num_reads * avg_read_length:,.0f} bp\n")
        f.write(f"Estimated coverage:      {coverage_estimate:.1f}x\n")
        f.write(f"Read time:               {timing['read_time']:.2f} seconds\n\n")
        
        f.write("-"*80 + "\n")
        f.write("DE BRUIJN GRAPH CONSTRUCTION\n")
        f.write("-"*80 + "\n")
        f.write(f"Graph nodes:             {num_nodes:,}\n")
        f.write(f"Graph edges:             {num_edges:,}\n")
        f.write(f"Average out-degree:      {num_edges/num_nodes:.2f}\n")
        f.write(f"Construction time:       {timing['graph_time']:.2f} seconds\n\n")
        
        f.write("-"*80 + "\n")
        f.write("ASSEMBLY RESULTS\n")
        f.write("-"*80 + "\n")
        f.write(f"Number of contigs:       {stats['num_contigs']:,}\n")
        f.write(f"Total assembly length:   {stats['total_length']:,} bp\n")
        f.write(f"Assembly vs. genome:     {assembly_fraction:.2f}%\n")
        f.write(f"Longest contig:          {stats['longest_contig']:,} bp\n")
        f.write(f"Shortest contig:         {stats['shortest_contig']:,} bp\n")
        f.write(f"Mean contig length:      {stats['mean_length']:,.1f} bp\n")
        f.write(f"N50:                     {stats['n50']:,} bp\n")
        f.write(f"Assembly time:           {timing['assembly_time']:.2f} seconds\n\n")
        
        f.write("-"*80 + "\n")
        f.write("TOP 20 LONGEST CONTIGS\n")
        f.write("-"*80 + "\n")
        for i, length in enumerate(contig_lengths[:20], 1):
            f.write(f"{i:3d}. {length:10,} bp\n")
        f.write("\n")
        
        f.write("-"*80 + "\n")
        f.write("CONTIG LENGTH DISTRIBUTION\n")
        f.write("-"*80 + "\n")
        bins = [
            (">100kb", sum(1 for x in contig_lengths if x > 100000)),
            (">50kb", sum(1 for x in contig_lengths if x > 50000)),
            (">10kb", sum(1 for x in contig_lengths if x > 10000)),
            (">5kb", sum(1 for x in contig_lengths if x > 5000)),
            (">1kb", sum(1 for x in contig_lengths if x > 1000)),
            (">500bp", sum(1 for x in contig_lengths if x > 500)),
        ]
        for bin_name, count in bins:
            f.write(f"Contigs {bin_name:8s}:     {count:,}\n")
        f.write("\n")
        
        f.write("-"*80 + "\n")
        f.write("TIMING SUMMARY\n")
        f.write("-"*80 + "\n")
        total_time = timing['total_time']
        f.write(f"Read time:               {timing['read_time']:8.2f} seconds "
                f"({timing['read_time']/total_time*100:5.1f}%)\n")
        f.write(f"Graph construction:      {timing['graph_time']:8.2f} seconds "
                f"({timing['graph_time']/total_time*100:5.1f}%)\n")
        f.write(f"Assembly:                {timing['assembly_time']:8.2f} seconds "
                f"({timing['assembly_time']/total_time*100:5.1f}%)\n")
        f.write(f"Total time:              {total_time:8.2f} seconds "
                f"({total_time/60:.2f} minutes)\n\n")
        
        f.write("="*80 + "\n")
        f.write("END OF REPORT\n")
        f.write("="*80 + "\n")



def assemble_mouse_genome(
    input_file="data/mouse_SE_150bp.fq.gz",
    output_fasta="mouse_assembly.fasta", 
    stats_file="mouse_assembly_stats.txt",
    k_mer_size = 51,
    random_seed = None
):
    """Main driver function for mouse genome assembly.
    
    This function orchestrates the complete assembly pipeline:
    1. Load FASTQ reads
    2. Build De Bruijn graph
    3. Assemble contigs
    4. Calculate statistics
    5. Write output files
    
    Returns:
        dict: Dictionary containing assembly results including:
            - dbg: DeBruijnGraph object
            - contigs: List of assembled sequences
            - stats: Assembly statistics dictionary
            - timing: Performance timing information
    """
    
    timing = {}
    
    start_time = time.time()
    reads = read_fastq(input_file)
    timing['read_time'] = time.time() - start_time
    
    num_reads = len(reads)
    total_bases = sum(len(read) for read in reads)
    avg_read_length = total_bases / num_reads if num_reads > 0 else 0
    
    start_time = time.time()
    dbg = DeBruijnGraph(reads, k=k_mer_size)
    timing['graph_time'] = time.time() - start_time
    
    # Calculate graph statistics
    num_nodes = len(dbg.graph)
    num_edges = sum(len(neighbors) for neighbors in dbg.graph.values())
    avg_degree = num_edges / num_nodes if num_nodes > 0 else 0
    
    start_time = time.time()
    contigs = dbg.assemble_contigs(seed=random_seed)
    timing['assembly_time'] = time.time() - start_time
    
    stats = dbg.get_assembly_stats(contigs)
    
    # Display distribution of contig lengths
    contig_lengths = sorted([len(c) for c in contigs], reverse=True)
    
    # Calculate coverage estimate
    genome_size_estimate = 2700000000  # Mouse genome ~2.7 Gbp
    coverage_estimate = (num_reads * avg_read_length) / genome_size_estimate
    assembly_fraction = (stats['total_length'] / genome_size_estimate) * 100
    
    # Write assembled contigs to FASTA
    dbg.write_fasta(contigs, output_fasta)
    
    # Write detailed statistics
    write_statistics_file(
        stats_file,
        input_file,
        num_reads,
        avg_read_length,
        k_mer_size,
        random_seed,
        num_nodes,
        num_edges,
        stats,
        contig_lengths,
        timing,
        coverage_estimate,
        assembly_fraction
    )

    timing['total_time'] = (timing['read_time'] + timing['graph_time'] + 
                           timing['assembly_time'])
    
    return {
        'dbg': dbg,
        'contigs': contigs,
        'stats': stats,
        'timing': timing,
        'graph_stats': {
            'num_nodes': num_nodes,
            'num_edges': num_edges,
            'avg_degree': avg_degree
        },
        'coverage': coverage_estimate
    }

if __name__ == "__main__":
    assemble_mouse_genome()
