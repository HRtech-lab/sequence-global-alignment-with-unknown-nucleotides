# A Python Program / Project

**This is a Python project from my early days as a Computer Science student**

_This programm was created for the sixth semester class Bioinformatics
and is one of the final projects for the class_

> #### Description of project
>
>>A Python script that for a parameter k, calculates the universal alignment of 2 sequences, with limitation that the alignment contains at most k unknown nucleotides. Nucleotide sequences are sometimes written in a 5-character alphabet, A, T, G, C, and N where N stands for an undefined nucleotide. A sequence with an N nucleotide is referred to as a degenerate sequence. Generally, a sequence with k unknown nucleotides has 4^k different interpretations.

> #### Implementation of project
> We first load 2 sequences from fasta files or create a 4 nucleotide array to randomly create 2 sequences. We create a matrix D to store the best dynamic programming scores for each pair of sequences, the matrix has dimensions of 1st sequence + 1 * 2nd sequence + 1. We create the sim_mat scaling matrix where we save values 2 for pair, 1 for purine and purine or pyrimidine with pyrimidine and -1 for purine with pyrimidine. These values will help us to replace and match. We also define a gap_score penalty of -2.
The first element of matrix D is zero or blank while the first row and column is filled in with gap_score adding it once more. To calculate matrix D we have 3 rules. If the highest score comes from the element D[i-1, j-1] and the scoring clause for each nucleotide means we have a pair, if it comes from the element D[i-1, j] and the scoring penalty then we have a gap in the first sequence while if it comes from element D[I, j-1] and the scoring penalty then we have a gap in the second sequence.
The next step is the traceback to create the alignment. We start from the last element D[m + 1, n + 1] to D[1,1]. To find the path that led to the optimal score, we first start with the matrix element D[m, n] and repeat the calculations backwards to find out where we came from. In other words, from box D[i, j] we recalculate the results at positions D[i-1, j-1], D[i-1, j] and D[i, j-1]. By doing this we can rebuild the path we came from. Our priority is to take the diagonal route first (if any). Adding a gap to the second priority and adding a gap to the third priority. Again, this choice of priority is completely arbitrary.
By reversing the nucleotides in the alignment we found we complete the universal alignment.


> #### About this project
>
> - The comments to make the code understandable, are within the .py archive
> - This project was written in IDLE, Python’s Integrated Development and Learning Environment.
> - Biological data used from https://www.ncbi.nlm.nih.gov/gene/ (genes 281894, 396218)
> - This program runs for Python version 2.7
> - This repository was created to show the variety of the work I did and experience I gained as a student
>
