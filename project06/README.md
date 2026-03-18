# Description

Welcome to our phylogenetic tree builder. Please see the notebook project06.ipynb for our demonstration. The other scripts are helper functions. We decided to test multiple distance metrics against one another and are showing a couple of those in the code. A critical question during this project was how to convert the smith-waterman score into a metric. Ultimately, this looked like the complement of the smith waterman score normalized by the minimum of the maximum possible scores for both sequences (i.e. the sequences are aligned against themselves).

$1 - SW(A, B) / min(SW(A), SW(B))$

We also discussed the differences between the recursive approach to the problem shown in class and the version shown in the notebook. It became clear that the notebook option was the way to go because our scores probably wouldn't satisfy the triangle inequality.

# Pseudocode

```python

Open file
    read sequences
    return sequences as strings

Neighbor-Joining Algorithm

let D be a matrix D x n
let D_names be an array with names

intialize the tree (empty dict/node obj)
initialize n leaf nodes (one per sequence)
initalize D (smith waterman distances normalized)

While active_nodes > 2 do
    row_sum, col_sum = D rowsum; colsum
    update Q <- Q_ij = (n-2) * D_ij - row_sum - col_sum         #update Q

    target = argmin(Q)                                          #find the min of Q

    Calculate branch lengths:
        $d_i = \frac{d(i,j) + (r_i - r_j)/(n-2)}{2}$
        $d_j = d(i,j) - d_i$
        where $r_i = \sum_{k=1}^n d(i,k)$

    find branch_length(k|ij)                                    #create internal node k #maybe a little confusing here

    for x in remaining_nodes do                                 #for loop may not be necessary here
        D(k, x) <- (D(i, x) + D(j, x) - D(i, j)) / 2            #update the distance matrix
        D(x, k) <- (D(i, x) + D(j, x) - D(i, j)) / 2            #symmetrical?
        #make the row/column here, then block update

    pop i, j from D; D_names

    add k to tree
        e.g. nested dict {D_names[i]: {k: branch_length}}
        e.g. obj i.add_child(k, branch_legnth)

    n -= 1

add branch ij with length d(i,j)
return tree

Convert to newick string e.g. ((B:0.2,(C:0.3,D:0.4)E:0.5)F:0.1)A
    walk tree recursively
    append node names to output list + distance

Plot using Phylo

```

# Successes
A successful implementation of the algorithm; conversion of newick strings; overall introduction to phylogenetic trees.

# Struggles
Not necessarily a struggle, but the matrix manipulations required some thought. Recursively building the newick string required two passes, and the final version ended up being much simpler than our initial draft. 

# Reflections

Eric: I was initially very confused about the difference between the algorithm discussed in class and the one in the notebook. It was difficult to know what the fundamental differences between the two were and why the iterative approach were without fully grasping the math or seeing the pseudocode for the recursive approach. Our coding sessions went over a lot across two days, and probably seemed fast as Sneha mentions. I think the issue is that familiarity with matrix manipulation was key for this project. A point of consternation was the removal of nodes i and j from the adjacency matrix and the addition of k while maintaining alignment and proper dimentions. I already had experience with tree objects, walking them, and recursive functions, so these were not a fundamental challenge at least for me. I think the project was a good way of pushing the class further into a number of directions simultaneously.

Sneha: This project was a challenge for me because the group worked at a much faster pace than I was used to in previous group settings. I found it difficult at times to keep up with the speed of development and discussions, which made me feel somewhat overwhelmed, especially early on. However, this also pushed me to become more proactive in reviewing the code independently and making sure I truly understood each component. One of the most challenging parts for me was understanding the recursion used to walk the tree and construct the Newick string. Initially, I struggled to grasp how the function was traversing the tree and building the output from the bottom up. The idea that each node relies on the fully computed results of its children before contributing its own portion of the string was not immediately intuitive. It was interesting to play around with the plots and look ar different scoring methods to see how they were affected.

# AI appendix
We used AI to generate some plots during the discovery phase, and for some syntax lookup.
