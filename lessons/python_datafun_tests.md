## Test yourself SOLUTIONS 
1. Create a dictionary called dna_map that maps DNA bases to their complementary bases. I.e., A -> T, C -> G, etc..

```
dna_map = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}

# For example:
dna_map['A']
```

2. Use str.maketrans() to convert this to a "translation table" named dna_table.

```
dna_table = str.maketrans(dna_map)
print(dna_table)
```


3. By passing dna_table to the string method .translate(), then using indexing with a step of -1 ([::-1]), compute the reverse complement of the MCS sequence given below.

```
mcs_seq = 'GAGACCCAAGCTGGCTAGCGTTTAAACTTAAGCTTGGTACCGAGCTCGGATCCACTA' \
          'GTCCAGTGTGGTGGAATTCTGCAGATATCCAGCACAGTGGCGGCCGCTCGAGTCTAG' \
          'AGGGCCCGTTTAAACCCGCTGATCAGCCT'

### Reverse complement 
mcs_seq_rc = mcs_seq.translate(dna_table)[::-1]
### Complement/translate 
mcs_seq_c = mcs_seq.translate(dna_table)

# Print forward and RC sequences
print('FW:', mcs_seq)
print('C: ', mcs_seq_c)
print('RC:', mcs_seq_rc)
print('RC:', mcs_seq_c[::-1])
```


4. Repeat question 3 but using the SeqIO module functions.

https://biopython.org/docs/1.75/api/Bio.Seq.html

```
from Bio.Seq import Seq
mcs_seq = 'GAGACCCAAGCTGGCTAGCGTTTAAACTTAAGCTTGGTACCGAGCTCGGATCCACTA' \
          'GTCCAGTGTGGTGGAATTCTGCAGATATCCAGCACAGTGGCGGCCGCTCGAGTCTAG' \
          'AGGGCCCGTTTAAACCCGCTGATCAGCCT'
mcs_seq_record  = Seq(mcs_seq)

mcs_seq_record_c = mcs_seq_record.complement()
mcs_seq_record_rc = mcs_seq_record.reverse_complement()

print('FW:', mcs_seq_record)
print('C: ', mcs_seq_record_c)
print('RC:', mcs_seq_record_rc)
```

  
