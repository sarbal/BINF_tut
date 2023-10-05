# Week 4: Snakes on a plane 
Adapted from https://www.w3schools.com/python/pandas/pandas_cleaning.asp and 
https://biopython.org/docs/1.75/api/Bio.Entrez.html
 
## Objectives 
This week's tutorial is on flow control and data wrangling in python. You will learn: 
- Flow control statements
- Functions and vectorised computations
- More on pandas for data wrangling
- Biopython

## Setting up
Start a new notebook. Save the file as "yourname_week4.ipynb". 
As before, copy the code into your notebook as chunks. 

## Download 
- [text/tab file](../data/DatasaurusDozen.txt)
- [csv file](../data/data.csv)
- [json file](../data/data.json)
    
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
dna_list = []
for c in dna:
    dna_list.append(c)
print('dna_list =', dna_list)
```

Alternatively, we can use a "list comprehension" to do this in one line
```
simple_list = [c for c in dna]
print('simple_list = ', simple_list)
```

If we want to use the index of each element in a list, we can create a counter that is incremented in each run through the code block
```
i=0
for c in dna:
    print('dna[%d] = %s'%(i,c))
    i += 1
```

Alternatively, we can use enumerate(), which takes any iterable as input and outputs a an iterable of index-value pairs. This will assign values to TWO variables in each pass of the loop.
This strategy is said to be more "Pythonic" 
```
for i, c in enumerate(dna):
    print("dna[%d] = %s"%(i,c))
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
while (i < n):
    print(i)
    i += 1 
```
## Vectorized computations
Vectorized computations give you faster operations by using fast, low-level code to operate on bulk data
```
import numpy as np
N = 100
ns = np.arange(N)
ns
```

## Functions
Functions, as in other languages, are a block of code/statements that return a value and/or perform a specific task. The main idea is to join repeatedly done tasks together and "call" the function, rather than write the same code over multiple times.
For example, if we want to calculate the factorial of a number, we can write a factorial function that takes in an integer as input (n) and returns the factorial of n (labelled val). So now everytime we wish to compute this, we simply call our function. 
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

print(n_factorial(5))
print(n_factorial(100))
print(n_factorial(0))
print(n_factorial(1001))
```

## Input 
### Command line input 
Reading in input from the command-line is in the argv variable from the sys module: 
```
import sys
print('Number of arguments:', len(sys.argv), 'arguments.')
print('Argument List:', str(sys.argv))
```
Note, this only works if you run a python script (there is no terminal) and the first argument (ie index 0) is the name of the script. 
What happens when you run this code in an interactive session or your notebook? Test it! 


### Reading/writing files 
Reading in input from a file with the open() function. It takes in two parameters: filename, and method.
There are four different modes for opening a file:
- "r" - Read - Default value. Opens a file for reading, error if the file does not exist
- "a" - Append - Opens a file for appending, creates the file if it does not exist
- "w" - Write - Opens a file for writing, creates the file if it does not exist
- "x" - Create - Creates the specified file, returns an error if the file exists
The name of the object that refers to the open file is the "stream handle" from the input/ouput (IO). In the code below, we've named our handle "f". There are certain functions that we can use with this, including reading and writing to the file. And also, we need to make sure to close the IO streamonce we've finished our tasks. Note, once we open a stream and read it in, we cannot "go back", the file is read in until it is "exhausted" (run out of input, hence a "stream") and we need to close and open it again.   
```
f = open("DatasaurusDozen.txt", "r")
print(f.readline()) # per line
print(f.read()) # whole document 
f.close()
```
(we will have a play with this data file later I promise) 

In addition you can specify if the file should be handled as binary or text mode
- "t" - Text - Default value. Text mode
- "b" - Binary - Binary mode (e.g. images)

### Reading in files with pandas
```
import pandas as pd
df = pd.read_csv('data.csv')
print(df.to_string())

df_js = pd.read_json('data.json')
print(df_js.to_string())
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
The general use of e-utils is one of these functions: 
- esearch: searching a database with a query
- efetch: fetching a record from a database with an key/ID
- elink: cross-searching database with an ID

Each function is broadly run like so, where the search returns a "handle", and you read the records/data fromt that handle in. Note, as before, you can only read in the stream once, and will have to repeat the function call if you do not store the record.   
```  
handle = Entrez.esearch(db="XXX", term=query)
record = Entrez.read(handle)
```
  
Searching:
Take a look at the website for different databases and queries you can use. 
E.g., Searching for the human taxonomic ID
```  
handle = Entrez.esearch(db="taxonomy", term="Human")
record = Entrez.read(handle)
len(record)
print(record['IdList'])
taxid = record['IdList'][0]
```

Fetching:
Now we have the ID/key for the term we want. We can then retreive that entry from the taxonomy database using that ID.

```
handle = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
record = Entrez.read(handle)
```

The record that we read in from the search handle is a list. We access each result using an index.
e.g., 

```
first_result = record[0]
```

And for multiple results, we can iterate through:
```
for result in record:
    print(result)
```

The results from this search are individual dictionaries, so we can view the data through keys() and values(). 
```
print(first_result.keys())
print(first_result.values())
```

And then retrieve individual results based on keys. For Taxonomic entries, some of these inlcude: 
```
print(first_result['ScientificName'])
print(first_result['Lineage'])
```
 
Cross IDs search:
Now, sometimes we have the key from once database and want the ID from a second. We can do this with the elink function.  
```
handle = Entrez.elink(dbfrom="taxonomy", db="assembly", from_uid=taxid)
links = Entrez.read(handle)
```

How many assembly entries do we get for the human taxon? 
As before, we first access the first entry, and then the set of links, and final the ids. It is a little convoluted because of the list -> dict -> list -> dict structure of the results.  
```
len(links[0]['LinkSetDb'][0]['Link'] ) 
```

Let's take the first result, and access the record in the assembly database. 
```
assembly_id = links[0]['LinkSetDb'][0]['Link'][0]['Id']
handle = Entrez.efetch(db="assembly", id=assembly_id, retmode="xml")
record = Entrez.read(handle)
print(record)
```

As you can see, not much in that! Let's use the taxid ID to search for other database entries, and once again select the first ID. 
```
handle = Entrez.elink(dbfrom="assembly", db="nuccore", from_uid=assembly_id)
links = Entrez.read(handle)
nuc_id = links[0]['LinkSetDb'][0]['Link'][0]['Id']
handle = Entrez.efetch(db="nuccore", id=nuc_id, retmode="xml")
record = Entrez.read(handle)
print(record)
```
What does this give us? A random region of the genome?!


You can refine your queries by selecting certain fields aswell. For example, this search looks for organisms with the human taxid, and filters on "RefSeq" entries and with the word "CYBB" in the record title.  
```
handle = Entrez.esearch(db="nuccore", term="txid"+str(taxid)+"[Organism] AND refseq[filter] AND CYBB[title]")
record = Entrez.read(handle)
nuc_id = record['IdList'][0]
handle = Entrez.efetch(db="nuccore", id=nuc_id, retmode="xml")
record = Entrez.read(handle)
print(record)
```
What does this give us? 


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
5. Write and run a query using the Entrez e-utils to search and fetch the nucletodie sequence for the mouse mitochondrial genome.
  
 
Solutions: Next week!

Back to the [homepage](../README.md)
