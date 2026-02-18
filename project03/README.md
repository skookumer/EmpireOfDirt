# Introduction
This program utilizes a Markov Model algorithm to attempt to identify a possible binding site in a *.bam file of ChIPseq sequences. In the version we have created here, we are analyzing a subset of 3 million ChIPseq reads from the archive entry SRR9090854, but the program can be changed to analyze other *.bam (if *.bam.bai file is also provided) if desired.

The Gibbs Motif Finder is designed to find statistically significant nucleotide sequences of length-k in a genome that has been read-separated. Our program implements this algorithm using two strategies, a "pythonic" one that should be understandable to the class and a "fast" version that is suited to the scope of the problem (3 million reads). We find that a ratio of above 1:1 iterations per read is necessary for finding distinct motifs, hence defaulting to the "fast" option.

Version 2 contains the following files:
- project03_v2.ipynb, a notebook which contains the main program and cells to visualize results
- seq_tools.py, new functions made to support project03_v2.ipynb
- motif_ops.py, one change was made to this file to support project03_v2.ipynb

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
During the second week we remade the program to work in more simple code that was easy to follow and run. This allowed all group members to easily understand the algorithm and participate in bug fixing and suggesting improvements.

# Struggles
Because our program was so fast, our initial approach was to read in the whole file and process it each time, instead of subsetting it. We realized that we were not getting a lot of information this way due to the sheer number of our initial motif scores, and that we should instead run our program many times on many randomized subsets. 
During the second week we really struggled with getting MACS3 working which limited the amount of time we had to spend on the second version of our program.
One difficulty we ran into was trying to improve our convergence citerion by requiring the IC difference to remain below a certain threshold for five consecutive iterations. We expected this to capture stabilization better. However, this approach did not work as expected and produced fluctuating IC values as iterations progressed.  

# Personal Reflections
## Group Leader
Eric: This was a fun project that tested our abilities. We met numerous times over the last two weeks to discuss the implementation and to refine our approach. The collaborative coding sessions were important to keeping everyone on the same page. Numerous questions remain to be explored about the convergence of the algorithm and the interpretation of the results. One question I have is how the reads should contribute to the input of the file. It seems like we should treat each read like a uniform kernel in a density estimation, but I'm still not clear on how these should be treated after the peaks are handled by MACS3. I expect MACS3 to be doing this under the hood, but the function we used jsut returns flat sequences, which might trnansform the results in an undesired manner. So there is still work to be done as far as understanding exactly how the algorithm should be implemented.

## Meghana Ravi
This project was definitely the one I struggled most with so far. I struggled with understanding the details of the concept in the beginning and how to apply them when writing code. My team was really supportive, Eric and Victoria helped me understand details I was confused about and I learned a lot of new things about coding throughout the process thanks to them. As we discussed and implemented the concepts I definitely started understanding better and I found it to be a very challenging but rewarding process. The extra week to work on this project gave me more time to understand the algorithm and how it works, but I also struggled a little with the dataset and trying to filter it. I think ultimately though, these two weeks have helped me really understand the concept and working of the algorithm.

## Victoria Van Berlo
This was the most challenging project yet, but definitely conceptually intriguing because it was in the realm of bioinformatics and the possible result of our program was unknown at conception. The number of functions and the extra code supplied by Eric made it a blessing and a curse, where the code ran quickly, but it was difficult to keep track of all the functions. Our group met many times which helped everyone stay on track and troubleshoot problems and redirect program flow as needed. 
After having the second week to solve the input problem, I spent the extra time really trying to understand the algorithm and we came up with a new working program from start to finish where I really feel like I understand now. 

# Generative AI Appendix
A lot of prompts for debugging and inquiring about numpy arrays and functions. Here are couple of highlights:

Gemini-3

Prompt:
```python

    def initialize_motif(seq, k):                       #random pick of k-length sequences for each line of DNA as Motifs
        idx = random.randint(0, len(seq) - k)
        seq = utils.seq_to_array(seq)
        motif = seq[seq[idx:idx+k]]                                 #get the motif
        background = np.concatenate([seq[:idx], seq[idx + k:]])     #get the background
        output = {"k-mer": motif, "idx": idx, "bg": background}
        return output

ok, here's the code that gets everything
```
pointed out the indexing issue going on with motif

prompt: I just want to replace that one list comp with an njit function, that's all

prompt:
```python
    @njit
    def convert_sequences(all_bytes):
        mdata = np.empty(len(all_bytes))
        for i in range(len(all_bytes)):
            mdata[i] = base_to_index[all_bytes[i]]
```

prompt: am i basically creating an indptr
prompt: 
```python
        start = 0
        for i in range(n_rows):
            end = indptr[i]
            seq = seqs[start:end]
```
ok, this is what I ahve
prompt: no, indptr[i] is the length of the first sequence, so we should end at the first sequence

Sonnet-4.5

prompt: Does gibbs sampling have an upper limit for sequences
prompt: I'm finding that my algorithm isn't converging and is essentially staying at baseline
prompt: I'm using softmax to convert scores to probabilities and am and also scoring reverse sequences



