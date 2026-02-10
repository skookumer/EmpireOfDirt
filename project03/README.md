# Introduction
This program utilizes a Markov Model algorithm to attempt to identify a possible binding site in a *.bam file of ChIPseq sequences. In the version we have created here, we are analyzing a subset of 3 million ChIPseq reads from the archive entry SRR9090854, but the program can be changed to analyze other *.bam (if *.bam.bai file is also provided) if desired.

# Pseudocode

```python
function GibbsMotifFinder(DNA, k-length)
    random pick of k-length sequences from each line of DNA as Motifs
    build frequency matrix for random k-mers
    build non-motif position array
    
    for j ← 1 to 10000 or Motifs stops changing
        i ← Random(N) where N is number of DNA entries
        PWM ← PWM constructed from all Motifs except for Motifi
        Motifi ← select position m from PWM-scored k-mers in DNAi in probabilistic fashion from score distribution
    return PFM
```

# Successes
Our group met many times to discuss possible implementations and decide on a general algorithm for the program to follow. We made sure all group members understood the intended flow of the program. We also did peer programming which helped everyone participate in the inception of the key parts of the data processing. Eric provided a lot of code that read the *.bam file quickly and used encoded sequences to speed up processing of the sequences. This allowed us to get data faster than requiring many minutes of initialization and saved iteration time. Because of this, we spent less time waiting for the program to run while troubleshooting.

# Struggles
Because our program was so fast, our initial approach was to read in the whole file and process it each time, instead of subsetting it. We realized that we were not getting a lot of information this way due to the sheer number of our initial motif scores, and that we should instead run our program many times on many randomized subsets.

# Personal Reflections
## Eric Arnold
Group leader's reflection on the project

## Meghana Ravi
Other members' reflections on the project

## Victoria Van Berlo
This was the most challenging project yet, but definitely conceptually intriguing because it was in the realm of bioinformatics and the possible result of our program was unknown at conception. The number of functions and the extra code supplied by Eric made it a blessing and a curse, where the code ran quickly, but it was difficult to keep track of all the functions. Our group met many times which helped everyone stay on track and troubleshoot problems and redirect program flow as needed.

# Generative AI Appendix
As per the syllabus
