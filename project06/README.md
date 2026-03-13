# Description

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

3. Calculate branch lengths:
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
Plot using SeqIO? (some biopython package)

```

# Successes

# Struggles

# AI appendix
