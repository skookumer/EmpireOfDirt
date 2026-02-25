### Description


### Pseudocode

```python
10,000 Fragments
remove the redundant ones

CHoose first fragment
Init k-mer size
iterate through the fragments
    choose fragment
    split into kmers (make the debrujin graph)
    compare the two graphs
    If no overlapping:
        some condition?

```

### Successes
1. We were successfully able to implement the algorithm for the given toy data.

### Struggles
1. We were struggling to identifying our starting point for graph traversal.
2. Transitioning the logic from graphs to sequences(contigs).
3. It was also hard for us to come up with a more appropriate stopping condition (syntatically) for looping until all the edges are consumed. 

### Reflections

## Group lead (Eric Arnold):

## Group members:
Aaronie Jersha Jenyfred: Chantera and Eric were amazing to work with. We met several times during the week and did group coding sessions which were really useful in terms of brainstorming different ideas. With some prior knowledge about the biological aspect of the De Bruijn graph, the concepts were easier to understand but creating an algorithmic pipeline was hard. I think after a point I started focusing on the syntax over the algorithm especially for the function def_assemble_contigs, which made things a bit more complicated than it needed to be. But we successfully implemented the De Bruijn algorithm for the toy data. Also, understanding the usage of default dictionaries was interesting.

### Generative AI appendix:

