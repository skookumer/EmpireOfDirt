# Project 11

# Description


# Pseudocode

```python

class profile-HMM
    N is the number of states (1 = Match, 2 = Insert, 3 = Delete)
    L is the max sequence length
    
    let vocab be the 20 amino acids and "-" and "Z" and "Not" (nothing)   #assuming Z is the non-existing amino acid
    let transitions be a stack of all transition probabilities          #dim: N x N x L
    let emissions matrix be a stack of all emission probabilities       #dim: N x vocab x L

    #Emission matrix example:
    M  0.1  0.23...
    I  0.5, 0.25...
    D  0    0...       1.0
        P    R    Y    Not


    function _forward_table
        let seq be the input sequence
        





    
```

# toy problem:
## execute baum_welch on the sequence

seq = VG--H

### _forward_table

transitions[0] = #our initial state matrix from the stack
I: 0.1
M: 0.8
D: 0.0
emissions[0] =
     G      V    NOT
I:  0.80  0.20  0.00  
M:  0.5   0.5   0.00 
D:  0.00  0.00  0.00
      [   ,   ,   ,   ,  ]
fwd = [   ,   ,   ,   ,  ]
      [   ,   ,   ,   ,  ]
  1. For our first emission [V]G--H, multiply the initial state probs by the emission probs. transitions[0] x emissions[0][V]
        [0.1 * 0.2 = .02,   ]
  fwd = [0.8 * 0.5 = .4,   ,]
        [0.0 * 0.0 = 0  ,   ]
  2.
      transitions[1] =
          I   M   D
      I: 0.2 0.7 0.1
      M: 0.3 0.5 0.2
      D: 0.0 0.3 0.7

      emissions[1] =
           G      V    NOT
      I:  0.30  0.50  0.20  
      M:  0.60  0.30  0.10
      D:  0.30  0.40  0.30
  
     For our second emission V[G]--H multiply the previous scores by transition column and the current value of the emission matrix:
    fwd[0][1] = sum(prev_score_col * current_emission_col) * current emission_i_j
          [.02, (.02 * .7 + .4 * .5 + 0 * .3) * 0.3 = .0642]
    fwd = [.08,.4 *   ,]
          [0.0,   ,   ,]
    
### _backward_table
  1.
      transitions[1] =
          I   M   D
      I: 0.2 0.7 0.1
      M: 0.3 0.5 0.2
      D: 0.0 0.3 0.7

      emissions[1] =
           G      V    NOT
      I:  0.30  0.50  0.20  
      M:  0.60  0.30  0.10
      D:  0.30  0.40  0.30
          [   ,   ,   ,   ,  ]
    bwd = [   ,   ,   ,   ,  ]
          [   ,   ,   ,   ,  ]
-

# Reflections

## Eric Arnold

## Thu Thu Han

## Stefanie Moreno

# AI Appendix
