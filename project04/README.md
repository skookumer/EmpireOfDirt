### Description

Please see project04_mouse.ipynb for the notebook stuff to see that the program is working. Leave the .fq.gz file as-is in the data directory for it to work. Functions are separated out into a utils script and a debruijn module script for ease of reading and editing. 

### Pseudocode

```python
10,000 Fragments
remove the redundant ones

CHoose first fragment
Init k-mer size
iterate through the fragments
    choose fragment
    split into kmers
    if key in graph:
        add edge to previous node
    else:
        add edge to new node

conduct eularian walk
    choose starting node
    walk dict keys until some condition (no outgoing edges)
    return the contig

assemble contigs
    combine the walks

get stats

```

### Successes
1. We were successfully able to implement the algorithm for the given toy data.

### Struggles
1. We were struggling to identifying our starting point for graph traversal.
2. Transitioning the logic from graphs to sequences(contigs).
3. It was also hard for us to come up with a more appropriate stopping condition (syntatically) for looping until all the edges are consumed. 

### Reflections

## Group lead (Eric Arnold):
The algorithm works for the toy problem. I feel like I would need to examine this with a bigger sequence and/or a natural language sequence to fully grasp what's going on. Recursion was a fun topic to talk about especially in retrospect of all the CS stuff I've done. Thinking in terms of recursion is making a lot of sense and we wrote the function in both formats. I was a little confused about the goal of the project because it seems that we don't really have a coherent genomic sequence at this point -- we just have a bunch of continguous fragments. So it seems like the extra step required to get a reference genome is coming in module 5. I was a bit confused by this because it seems like we didn't go all the way in this one.

## Group members:
Aaronie Jersha Jenyfred: Chantera and Eric were amazing to work with. We met several times during the week and did group coding sessions which were really useful in terms of brainstorming different ideas. With some prior knowledge about the biological aspect of the De Bruijn graph, the concepts were easier to understand but creating an algorithmic pipeline was hard. I think after a point I started focusing more on the syntax over the algorithm especially while working on the function def_assemble_contigs, which made things a bit more complicated than it needed to be. But we successfully implemented the De Bruijn algorithm for the toy data. Also, understanding the usage of default dictionaries was interesting.

### Generative AI appendix:

