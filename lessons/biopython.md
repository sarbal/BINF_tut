# Week 5: Biopython and assignment help 
Adapted from https://biopython.org/docs/1.75/api/Bio.Entrez.html
 
## Objectives 
This week's tutorial is on biopython. You will learn: 
- How to access and manipulate data from different databases 

## Setting up
Start a new notebook. Save the file as "yourname_week5.ipynb". 
As before, copy the code into your notebook as chunks. 

## Biopython
Biopython is a set of freely available tools for biological computation written in Python by an international team of developers.
It is a distributed collaborative effort to develop Python libraries and applications which address the needs of current and future work in bioinformatics. 
Quick install: 

```
pip install biopython
```
Then import into your notebook (or console/terminal or script):

```
import Bio
```
There is a lot of functionality in [biopython](https://biopython.org/docs/ ) (too much to cover here) but revolve around sequences and sequence analysis. 
These include BLAST searches, downloading sequences from NCBI, Phylogenetics, Cluster analysis, Graphics, etc. 
The most useful function is accessing NCBI through their e-utils API. 
Some extra tutorials here: https://biopython-tutorial.readthedocs.io/en/latest/notebooks/00%20-%20Tutorial%20-%20Index.html#  *NOTE: these will be important for your assigmment...*   

### Sequences 
```
from Bio.Seq import Seq

my_seq = Seq("CATGTAGACTAG")

# print out some details about it
print("seq %s is %i bases long" % (my_seq, len(my_seq)))
```

### Reading and writing Sequence Files
Use the SeqIO module for reading or writing sequences as SeqRecord objects.
```
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

record = SeqRecord(
    Seq("MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF"),
    id="YP_025292.1",
    name="HokC",
    description="toxic membrane protein, small",
)
print(record)
# As Genbank entry
Bio.SeqIO.write(record, "HokC.gbk", "gb")
# As FASTA file 
Bio.SeqIO.write(record, "HokC.fasta", "fasta")
```
Other formats here: 
- https://biopython.org/wiki/SeqIO

  
