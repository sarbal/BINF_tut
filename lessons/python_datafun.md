# Week 4: Snakes on a plane 
Adapted from https://www.w3schools.com/python/pandas/pandas_cleaning.asp and 
https://biopython.org/docs/1.75/api/Bio.Entrez.html
 
## Objectives 
This week's tutorial is on flow control and data wrangling in python. You will learn: 
- Flow control statements
- Functions and vectorised computations
- More on pandas for data wrangling

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


# Q1: Test yourself!
### In your jupyter notebook, create new chunk for this question. In this chunk, write a for loop that prints out the first ten powers of 2 starting from 2^0. 

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
# Q2: Test yourself!
### In your jupyter notebook, create new chunk for this question. In this chunk, write a function that prints out the first N powers of 2. It should take in the number of values "N" as input, and check that N is not larger than 100.   

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


# Q3: Test yourself!
### In your jupyter notebook, create new chunk for this question. In this chunk, 

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


# Q4: Test yourself!
### In your jupyter notebook, create new chunk for this question. In this chunk,  


[Solutions next week]

Back to the [homepage](../README.md)
