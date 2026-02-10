# Introduction
This program utilizes a Markov Model algorithm to attempt to identify a possible binding site in a *.bam file of ChIPseq sequences. In the version we have created here, we are analyzing a subset of 3 million ChIPseq reads from the archive entry SRR9090854, but the program can be changed to analyze other *.bam (if *.bam.bai file is also provided) if desired.
The Gibbs Motif Finder is designed to find statistically significant nucleotide sequences of length-k in a genome that has been read-separated. Our program implements this algorithm using two strategies, a "pythonic" one that should be understandable to the class and a "fast" version that is suited to the scope of the problem (3 million reads). We find that a ratio of above 1:1 iterations per read is necessary for finding distinct motifs, hence defaulting to the "fast" option.

This program uses several files:
- project03.ipynb, a notebook demonstrating functionality at a high level
- project03.py, the main script where the Gibbs function lives
- utils.py, a class of static methods that use numba for the fast method
- sequence_database.py, a python object for storing the flat array of the genome for fast iteration.

The algorithm behaves as expected, though the convergence condition needs refinement as it does not use entropy and instead reilies on a fixed threshold for the difference in matrices.


# Pseudocode

```python
function GibbsMotifFinder(DNA, k-length)
    random pick of k-length sequences from each line of DNA as Motifs
    build frequency matrix for random k-mers
    
    for j ← 1 to 10000 or Motifs stops changing
        i ← Random(N) where N is number of DNA entries
        PFM ← PFM constructed from all Motifs except for Motif_i
        PWM ← PWM constructed from all Motifs except for Motif_i
        Seq_i ← retrieve sequence at position i
        score all motifs in seq_i
        pick best score with argmax(probs(score))
        update PFM with new nucleotide sequences
    return PFM
```

# Successes
Our group met many times to discuss possible implementations and decide on a general algorithm for the program to follow. We made sure all group members understood the intended flow of the program. We also did peer programming which helped everyone participate in the inception of the key parts of the data processing. Eric provided a lot of code that read the *.bam file quickly and used encoded sequences to speed up processing of the sequences. This allowed us to get data faster than requiring many minutes of initialization and saved iteration time. Because of this, we spent less time waiting for the program to run while troubleshooting.

# Struggles
Because our program was so fast, our initial approach was to read in the whole file and process it each time, instead of subsetting it. We realized that we were not getting a lot of information this way due to the sheer number of our initial motif scores, and that we should instead run our program many times on many randomized subsets.

# Personal Reflections
## Eric Arnold
Group leader's reflection on the project
- Implementation of Gibbs Sampling methods.
- more experience with numpy
- 

# Struggles
- What is the rule of thumb for the number of iterations for a dataset of this size? 10,000 seems way too few.
- 

# Personal Reflections
## Group Leader
Eric: I decided to develop the "fast" pipeline as well in addition to the vanilla pythonic one we all worked on. I did this because I knew that we would not be able to get through iterations on 3 million rows of data in a reasonable amount of time (1 iteration of our pythonic loop took 36s). 6-million iterations with the fast method took roughly 7 hours. We tried to keep the versions somewhat consistient as far as the results are concerned and incorporated tests to ensure that the output was equal.

There are many elements of this project that remain to be explored and could be revised. Chiefly:
- The initialization state and whether to include reverse complements in these motifs; we could try this
- The use of background frequencies for score calcluation
- The specifics of how score distributions for the fwd and reverse motifs are combined and used in the selection function
- The conversion of the score logits into probabilities (we used softmax with a temperature of .1).
- Subsampling and different subsampling strategies. I explored an iterative method of freezing/unfrezzing rows and gradually adding to the PFM calculation that might bias the algorithm towards a certain region of the search space, but might enable faster convergence.
- Termination condition. We used a threshold for the difference in the PWMs between iterations, but the threshold for this is non-normalized and would change with the length of k-mer used. I noticed at the last minute that motif_ops.ic_content was the required function, so this could be implemented in the revision.
- Aggregation of the parallel output -- to choose the best matrix or to sum them

Overall, this was an exciting project and opened up my thinking about unupervised tokenization learning strategies. The application of the MCMC method makes me wonder about the fundamental differences between genomic grammar and natural language.

## Meghana Ravi
Other members' reflections on the project

## Victoria Van Berlo
This was the most challenging project yet, but definitely conceptually intriguing because it was in the realm of bioinformatics and the possible result of our program was unknown at conception. The number of functions and the extra code supplied by Eric made it a blessing and a curse, where the code ran quickly, but it was difficult to keep track of all the functions. Our group met many times which helped everyone stay on track and troubleshoot problems and redirect program flow as needed.

# Generative AI Appendix
A lot of prompts for debugging and inquiring about numpy arrays and functions. Here are couple of highlights:


