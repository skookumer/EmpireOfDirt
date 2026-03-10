### Description

Please see demonstration.ipynb for a demonstration of our code. We did three separate implementations after discussing the concept of the algorithm; these are separate scripts in the folder first_code_runthrough_sm.py, project5.py, and fast_align.py. See these scripts specifically if interested in commentary and details about the implementation.

### Pseudocode

```python
function initialize_matrix
    set first row and col to zeros
    score using max function for adjacent cells

function traceback
    start_idx = argmax(score_matrix)
    for cell in adjacent_cells:
        if score matches predicted
          save movement step
          move on to next cell
        elif score == 0:
          break loop
```

### Successes
Three unique implementations of the algorithm.

### Struggles
Squashing bugs; comparing output


### Reflections

Eric Arnold: with the recursive approach, it was good practice to debug the whole thing to see where it was going wrong. For instance, the numerical precision for my score matrix was woefully inadequate (int-8) and only reared its head when trying to align longer sequences and comparing to the other implementations of the function. It was also good practice transferring the array methods into numpy. While pure arrays would be optimal, the function ultimately makes use of a Jit dynamic array (JitList) to store results. The recursive approach only makes sense if you want to get every possible aligned sequence, but it breaks catastrophically when input sequences become complex an can have multiple gaps in different places to get the optimal score, such as with homologous or repeated sections, as discovered in testing.

### GenAI appendix
Used to generate some of the examples, troubleshoot syntax with numpy
