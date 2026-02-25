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
4. Handling memory complexity when dealiing with large datasets

### Reflections

## Group lead (Eric Arnold):

## Group members:
Aaronie Jersha Jenyfred: Chantera and Eric were amazing to work with. We met several times during the week and did group coding sessions which were really useful in terms of brainstorming different ideas. With some prior knowledge about the biological aspect of the De Bruijn graph, the concepts were easier to understand but creating an algorithmic pipeline was hard. I think after a point I started focusing more on the syntax over the algorithm especially while working on the function def_assemble_contigs, which made things a bit more complicated than it needed to be. But we successfully implemented the De Bruijn algorithm for the toy data. Also, understanding the usage of default dictionaries was interesting.

Chantera: I enjoyed working with Eric and Jersha on this project. Both members challenged me to think more conceptually about the algorithm and its biological application. Discussing the concepts with my team helped deepen my understanding of the algorithm, although implementation proved to be more difficult. Specifically, I struggled with implementing recursion and using `defaultdicts`. It was also difficult to differentiate between what should be stored in the class's state vs when it should just be a local copy in the function. Since the algorithm removes edges, I had to be mindful about mutating self.graph vs the local copy of the graph. Finally, working with real read fragments highlighted that they do not statisfy balance conditions presentied by Eulerian path as k-mer multiplicity and uneven coverage violates that balance. Our implementation extracts path using the semi-balance conditions which we deemed more biologically appropriate.

### Generative AI appendix:
Generative AI was used to review reflection for grammatical clarity and minor syntax corrections.
