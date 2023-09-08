# Week 4: Snakes on a plane 
## Objectives 
This week's tutorial is on flow control and data wrangling in python. You will learn: 
- Flow control statements
- Functions and vectorised computations
- Biopython

## Setting up
Start a new notebook. Save the file as "yourname_week4.ipynb". 
As before, copy the code into your notebook as chunks. 

## Download 
- [csv file](data.csv)
- [json file](data.json)
    
## if statements
Allow blocks of code to be executed only under specific conditions.
```
x = 2
y = 4

if x==y:
    print('In block 1')
    print('They are equal!')

elif x>y:
    print('In block 2')
    print('x is more than y')

else:
    print('In block 3')
    print('y is more than x')
```
Note the indentation within each code block. Code within the same block must have the same indentation level, since this is how Python code blocks. Although the amount of indentation doesn't actually matter, you should adhere to the PEP 8 standard that code blocks be indented with 4 spaces, not with tabs.

## for loops
The most common type of loop in data analysis is the for loop. A for loop executes a code block once for each value in a collection of values that you specify.

```
# Define list of numbers
num_list = list(range(10))
print('num_list = ', num_list)

# Print these numbers one by one
for n in num_list:
    print('->', n)
```

For loops can loop over any "iterable" object, such as a string.
```
dna = 'AACTTGGTTAAATTATGGG'
for c in dna:
    print('->', c)
```

If we want to fill a list with content, we can write a for loop:
```
# dna_list = []
for c in dna:
    dna_list.append(c)
print('dna_list =', dna_list)
```

Alternatively, we can use a "list comprehension" to do this in one line
```
simple_list = [c for c in name]
print('simple_list = ', simple_list)
```

If we want to use the index of each element in a list, we can create a counter that is incremented in each run through the code block
```
i=0
for c in name:
    print('name[%d] = %s'%(i,c))
    i += 1
```

Alternatively, we can use enumerate(), which takes any iterable as input and outputs a an iterable of index-value pairs. This will assign values to TWO variables in each pass of the loop.
This strategy is said to be more "Pythonic" 
```
for i, c in enumerate(name):
    print("name[%d] = %s"%(i,c))
```

Range is the Pythonic way of iterating over consecutive integers.
```
for x in range(10):
    print(x)
```

## while loops
The while loop keeps going as long as the argument it is passed evaluates to "True".  
```
i = 1
n = 10
while i <= n;
    print(i)
    i = i + 1
```
## Vectorized computations
Vectorized computations give you faster operations by using fast, low-level code to operate on bulk data
```
N = 100
ns = np.arange(N)
ns
```

## Functions
```
def n_factorial(n): 
    # Check that n is of the right type and form
    assert isinstance(n,int),'Input is not an integer'
    assert n >= 0, 'Input needs to be positive'
    assert n <= 1000, 'Intput is too large!'

    # Initialize return variable
    val = 1

    # Loop over n
    for i in range(1,n+1):
        val *= i

    return val
```

## Input 
### Command line input 
Reading in input from the command-line is in the argv from the sys module: 
```
import sys
print 'Number of arguments:', len(sys.argv), 'arguments.'
print 'Argument List:', str(sys.argv)
```
Note, the first argument is the name of the script. 

### Reading/writing files 
Reading in input from a file with the open() function. It takes in two parameters: filename, and method.
There are four different modes for opening a file:
- "r" - Read - Default value. Opens a file for reading, error if the file does not exist
- "a" - Append - Opens a file for appending, creates the file if it does not exist
- "w" - Write - Opens a file for writing, creates the file if it does not exist
- "x" - Create - Creates the specified file, returns an error if the file exists
  
```
f = open("data.txt", "r")
print(f.readline()) # per line
print(f.read()) # whole document 
f.close()
```
In addition you can specify if the file should be handled as binary or text mode
- "t" - Text - Default value. Text mode
- "b" - Binary - Binary mode (e.g. images)

### Reading in files with pandas
```
import pandas as pd
df = pd.read_csv('data.csv')
print(df.to_string())

df = pd.read_json('data.json')
print(df.to_string())
```

As a data frame, you can manipulate/look at your data: 

```
print(df.head(10))
print(df.tail()) 
```

Dropping empty cells
```
new_df = df.dropna()
df.dropna(inplace = True) # replaces within the same dataframe
```

Dropping rows 
```
df.dropna(subset=['Date'], inplace = True)
```

Replacing cells 
```
df.fillna(130, inplace = True) # replace NAs with fill value
df["Calories"].fillna(130, inplace = True) # replace values in the column "Calories" with fill value 
x = df["Calories"].mean()
df["Calories"].fillna(x, inplace = True) # replace with mean
```

Changing format 
```
df['Date'] = pd.to_datetime(df['Date'], format="mixed", yearfirst=True)   
print(df.to_string())
```

Set values that are "wrong" (mistakes in data entry/recording) to a max value.

```
for x in df.index:
  if df.loc[x, "Duration"] > 120:
    df.loc[x, "Duration"] = 120
```

Or you can drop the rows that meet that criteria: 
```
for x in df.index:
  if df.loc[x, "Duration"] > 120:
    df.drop(x, inplace = True)
```

Another issue could be duplicated rows. You can use this check to find them: 
```
print(df.duplicated())
df.drop_duplicates(inplace = True)
```
 
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


### Using E-utils
```
from Bio import Entrez
Entrez.email = "my.email@unsw.edu.au"
Entrez.tool = "my_script.py"
```
Searching:
```  
handle = Entrez.esearch(db="XXX", term=query)
record = Entrez.read(handle)
```

Fetching:
```
handle = Entrez.efetch(db="XXXX", id=input, retmode="xml")
record = Entrez.read(handle)
```

Cross IDs search:
```
handle = Entrez.elink(dbfrom="assembly", db="nucleotide", from_uid=accession)
links = Entrez.read(handle)
```

More here: 
https://biopython.org/docs/1.75/api/Bio.Entrez.html

## Test yourself! 
1. Create a dictionary called dna_map that maps DNA bases to their complementary bases. I.e., A -> T, C -> G, etc..
2. Use str.maketrans() to convert this to a "translation table" named dna_table.
3. By passing dna_table to the string method .translate(), then using indexing with a step of -1 ([::-1]), compute the reverse complement of the MCS sequence given below.
```
mcs_seq = 'GAGACCCAAGCTGGCTAGCGTTTAAACTTAAGCTTGGTACCGAGCTCGGATCCACTA' \
          'GTCCAGTGTGGTGGAATTCTGCAGATATCCAGCACAGTGGCGGCCGCTCGAGTCTAG' \
          'AGGGCCCGTTTAAACCCGCTGATCAGCCT'
```
4. Repeat question 3 but using the SeqIO module functions.
  
 
Solutions: Next week!

Back to the [homepage](../README.md)
