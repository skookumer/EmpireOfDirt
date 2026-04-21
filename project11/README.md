# Project 11

# Description

Here is our pseudocode and HMM module for project 11, which includes updated versions of backwards, forwards, baum-welch, viterbi, and an implementation of MSA. Please see demo.ipynb for a demonstration that the code works (albeit perhaps incorrectly) with the test sequences. The HMM module is capable of being run in both modes, but some functionality is exclusive to one, i.e. printing the model only works in normal HMM mode. The fix for this is to implement inheritance, but that's for a later time. All new code is either at the top or bottom of the object. Initialization methods are at the top and updated profile methods are at the bottom.

What distinguishes our profile HMM from the regular HMM is A) the use of matrices for every index and B) masking the matrices to reflect the inability to transition from D -> I or I -> D and the emissions in M and I states for gap. We chose to incorporate gap "-" but to silence it by zeroing out the probabilities for those positions. We think this works correctly but relying on the zeroed probabilities in the computation of gamma is a potential problem since there is no hard-guarantee that problematic transitions do not happen, other than the masking. In other words, everything is handled probabilistically.

# Pseudocode

## MSA initialization

```python
Function: MSA alignment
Description: Construct Match, Insert, and Delete states from a multiple sequence alignment.

1.	Determine the number of alignment columns.

2.	For each column:
        Count the number of matches and gap.
    	If the proportion of gaps is less than 50%, mark this column as a Match column.
    	Otherwise, mark this column as an Insertion column.

3.	Create the following states in order:
        Begin state
        For each column k:
    	    Match state Mk
    	    Insert state Ik
    	    Delete state Dk
    	End state

4.	Initialize emission distributions:
        For each Match state Mk:
    	    Count residues in column k (excluding gaps).
    	    Convert counts to probabilities.
    	For each Insert state Ik:
        	Initialize with background frequencies or uniform distribution.
    	For each Delete state Dk:
        	No emissions.

5.	Initialize transition structure:
    	From Begin: transitions to M1 and D1.
    	For each k from 1 to K−1:
        	Mk → Mk+1, Ik, Dk+1
        	Ik → Ik, Mk+1
        	Dk → Dk+1, Mk+1
    	From final states Mk, Ik, Dk: transitions to End.

6.	Initialize transition probabilities:
    	For Match and Insert states:
            estimate from alignment paths.
    	For Delete states:
            set all outgoing transitions to 1 (deterministic).

7.	Return the full Profile HMM structure.

```

## Baum-Welch + fwd + bwd
```python

class profile-HMM
    N is the number of states (1 = Match, 2 = Insert, 3 = Delete)
    L is the max sequence length

    let vocab be the 20 amino acids and "-" and "Z"                     #assuming Z is the non-existing amino acid
    let transitions be a stack of all transition probabilities          #dim: N x N x L
    let emissions matrix be a stack of all emission probabilities       #dim: N x vocab x L

    #Emission matrix example:
    M  0.1  0.23...
    I  0.5, 0.25...
    D  0    0...       1.0
        P    R    Y    

    let avg_len be the average length of the input sequences # normailzing the column len of the matrix : avg len = sum of seq leq / total no of seqs


    function _forward_table
        let seq be the input sequence # loop through one seq by one seq
        let forward matrix be the unfilled matrix # initialize the matrix with 0s 
        initialize probability from the initial prob array # get first character from obs seq 

        for index in sequence
            get emission column from emission matrix at index # get emmission prob for first char index 

            for state in states:
                previous scores is forward matrix [index - 1] 
                transitions is the matrix of transition probabilities at index
                emission is the value in emission column at state

                let scores be prev scores x transitions
                forward matrix [state, index] = product(scores) * emission # calculation of forward table

    function _backward_table
        let seq be the reversed input sequence sliced by avg_len # reversing sliced obs
        let backward matrix be the unfilled matrix # initilaize the backward matrix with 0s
        initialize first column of backward to be prob = 1 

        for index in sequence
            get emission column from emission matrix at avg_len - index

            for state in states:
                previous scores is backward matrix [index - 1]
                transitions is the matrix of transition probabilities at avg_len - index
                scores is previous scores x transitions x emission column
                backward [state, index] = product(scores)

    function baum_welch_profile

        set iterations
        set threshold for convergence
        initialize accumulation matrices, num = avg_len 
        

        while not converged and less than set iterations do

            initial accumulation
            emission accumulation   (3 x 22 x avg_len) # accumalation matrix for emissions at each obs in avg_len
            transition accumulation (3 x 3 x avg_len)  # accumalation matrix for transitions at each obs in avg_len

            for seq in input sequences do
                let gamma be the posterior probability from forward/backward
                let p_seq be the sum of the last column of forward
                gamma = fwd + bwd - p_seq
                add gamma column [0] to initial accumulation

                for index in seq range of (avg_len) do
                    # adjust the emissin matrix associtate with current index and current seq
                    add gamma column [index] to emission accumulation [index] 

                    if index < avg_len - 1 do # -1 to avoid index error
                        add to trans accumulation:
                            forward [index] *
                            transitions [index] *
                            emissions [index] *
                            backward [index + 1] /
                            p_seq #to normalize


        normalize accumulations by the sum of accumulated probabilities
        check convergence  
        reassign matrices for next iteration


```python
    
```

# Reflections

## Eric Arnold

This project was a lot to cram into 6 days relative to the previous projects. I can see why Marcus wanted us to only submit pseudocode. Nevertheless, we felt that to truly know how the pseudocode should go, we should actually implement the algorithm and see where it breaks. As expected, the overall idea that a profile HMM is just an HMM with a set of matrices for each index was true, however indexing specific matrices in the model was troublesome in the conversion process. A bug we encountered with numpy arrays was the flattening of the masks for 2-d matrices, which we had to mitigate by processing row by row. Another bug was that in the trans_adj calculation, we actually need to point the emission matrix to the next index to correspond with the next observation, which was an unexpected change. Finally, indexing in the backwards calculation needs to be reversed, and we must slice off the last part of the input sequence so that we start at the end corresponding with our truncated forward pass.

An unintuitve aspect of MSA initialization is determining what we caleld the "consensus states" based upon the observations. It seems that by defaulting to match for most indices, where the majority of emissions are not gap, the model is liable to commit type II errors where we should actually reject a match state. We questioned why MSA and Baum-Welch co-exist and have come to the conclusion that this is probably why. MSA provides a good-enough set of initial probabilities but they're probably biased towards match states, which fwd-bwd would register as less probable. It would be interesting to test this hypothesis by examining the proportion of match states assigned over time given MSA initialization for longer sequences.

## Thu Thu Han

This week's project was definitely a learning curve for me and honestly even till now I am not satisfied with my understanding of profile HMMs. Putting Baum Welch, Viterbi, forward and backward together to get a profile HMM within one week was quite a lot of work even though I understand the concepts individually it was challenging to work with them in the process of building profile HMM which made me question my understanding in the first place. I think while our decision to work with matrices two weeks ago did turn out to be in our favor this week, keeping in track of indexes were quite challenging. I think I also filled in my understanding gaps in Balum welch from last week while brain storming on this project. While MSA was a simpler algorithm compared to Balum welch, one of the most challenging aspects was also trying to put the two algorithms together for profile HMM. I am also glad this week we only have to pseudocode and not actually implement them because I believe that would take longer than just one week to complete for me

## Stefanie Moreno

This assignment was challenging in that there was a lot of conflicting information across resources and very different interpretations among team members. Despite meeting many times over the week, I didn’t feel that we were able to fully overcome the difficulties caused by miscommunication, differences in coding proficiency, and the general complexity of Profile HMMs. Even so, working through these obstacles forced me to clarify the underlying concepts much more deeply than I would have otherwise. I had to reconcile multiple definitions of match, insert, and delete states, understand how silent states behave in the algorithms, and translate the theoretical model into a consistent pseudocode framework. Although the group struggled to align on a unified approach, the process helped me strengthen my own understanding of how Profile HMMs are constructed and why their architecture is so sensitive to small conceptual misunderstandings. In the end, the assignment pushed me to be more precise, more patient, and more confident in my ability to navigate ambiguity in computational biology.



# AI Appendix
"Daddy Claude," as he is known, is a friend, but he wasn't much help with the convergence issues, which needed to be resolved manually by printing. He was of help in understanding why the 2d masks were falling apart though and helping fix that bug.
