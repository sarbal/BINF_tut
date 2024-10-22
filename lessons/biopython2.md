# Week 7: Biopython and assignment help 
Adapted from https://biopython.org/docs/1.75/api/Bio.Entrez.html
 
## Objectives 
This week's tutorial is on biopython. You will learn: 
- How to access and manipulate data from different databases 

## Setting up
Start a new notebook. Save the file as "yourname_week7.ipynb". 
As before, copy the code into your notebook as chunks. 

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
 

You can refine your queries by selecting certain fields aswell. 
For example, this search looks for organisms with the human taxid, and filters on "RefSeq" entries and with the word "CYBB" in the record title.  
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
