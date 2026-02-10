# Introduction
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

## Other member
Other members' reflections on the project

# Generative AI Appendix
A lot of prompts for debugging and inquiring about numpy arrays and functions. Here are couple of highlights:


