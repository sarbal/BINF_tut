## Test yourself SOLUTIONS 
1. Create a dictionary called dna_map that maps DNA bases to their complementary bases. I.e., A -> T, C -> G, etc..

2. Use str.maketrans() to convert this to a "translation table" named dna_table.

3. By passing dna_table to the string method .translate(), then using indexing with a step of -1 ([::-1]), compute the reverse complement of the MCS sequence given below.

```
mcs_seq = 'GAGACCCAAGCTGGCTAGCGTTTAAACTTAAGCTTGGTACCGAGCTCGGATCCACTA' \
          'GTCCAGTGTGGTGGAATTCTGCAGATATCCAGCACAGTGGCGGCCGCTCGAGTCTAG' \
          'AGGGCCCGTTTAAACCCGCTGATCAGCCT'
```

4. Repeat question 3 but using the SeqIO module functions.

  
