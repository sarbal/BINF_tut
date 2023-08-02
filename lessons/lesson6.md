# Week 6: Python basics
Taken/adapted from: https://ucdavis-bioinformatics-training.github.io/2022-Feb-Introduction-To-Python-For-Bioinformatics/python/python2
First things first! Download these files into your working directory: 
- [lesson5](../data/lesson)
- [helper.py](../data/helper.py)
###

Hello, World! is traditionally the first program to write in any new programming language.
Open a new file n your editor. 
Copy and paste (or type) this into your new file:
```
print("Hello, World!")
```
Save the file and name it “helloworld.py”. Python code files typically use the “.py” extension. 
Then click the triangle “play” button to run your code in the terminal.

Python can be run interactively, however, we will be saving the code to files to be run in VSCode.  

## Defintions
Variables
Data types 
Functions
Casting 

### Variables 
- Named storage that can be used to store anything from numbers to characters to much more complex data.
- Variable names can contain only alphanumeric characters (A-Z, a-z, 0-9) and underscore (_) and no other special characters and cannot start with a number. 
- Type is determined by the value assigned
- Assign values with the equals sign =  

#### Basic Data Types: 
Integers, Floating-point numbers, booleans, strings.
- Integers: whole numbers, negative or positive
```
n = 42 
```
- Floats: basically Real numbers
```
n_pi = 3.14
print(gene_exp)
```
- Strings: Any text (including all symbols and spaces), designated by quotes.
```
gene = "TAF1"
print(gene)
```



### Built-in string methods: 
Strings have a long list of built-in methods to return modified versions of the string. 
The methods do not modify the string, they return a new string.
tmpstr = "hello my name is Nik."
allcaps = tmpstr.upper()
print(allcaps)
newstr = tmpstr.replace("Nik","Joe")
print(newstr)

- you can also use the + symbol to concatenate strings
tmpstr2 = "How are you doing?"
print(tmpstr + " " + tmpstr2)

- you can find the position of one string in another. Returns -1 if not found
print(tmpstr2.find("you"))
print(tmpstr2.find("california"))
Here is a list of the various methods you can use on strings.

PRACTICE: Take a look at the string methods and play around with them. See if you can figure out how to use the “split”, “isdigit”, and “index” methods.

- Booleans
“In computer science, the Boolean data type is a data type that has one of two possible values (usually denoted true and false) which is intended to represent the two truth values of logic and Boolean algebra. It is named after George Boole, who first defined an algebraic system of logic in the mid 19th century.” -wikipedia

control = False
treatment = True
print(control)
print(treatment)
Booleans come into play when making comparisons. The simplest comparisons are comparing integers:
```
print(1<1)
print(1<2)
print(2>1)
print(1<=1)
print(2>=1)
```
- the double equals is an equality comparison, a single equals is for assignment.
```
print(1==1)
print(0==1)
```
- You can also do string comparisons
```
gene_id = "BRCA2"
hw = "HELLO!!! WORLD!!!"
print(gene_id == "BRCA2")
print(hw == "hello")
```
### Built-in datatypes and the type function
You can use the “type” function to query the type of a variable:
```
print("The data type of the variable 'control' is:")
print(type(control))
print("The data type of the variable 'hw' is:")
print(type(hw))
print("The data type of the variable 'gene_count' is:")
print(type(gene_count))
```
Here is a list of the built-in datatypes and their corresponding category:
- Text Type: str
- Numeric Types: int, float, complex
- Sequence Types: list, tuple, range
- Mapping Type: dict
- Set Types: set, frozenset
- Boolean Type: bool
- Binary Types: bytes, bytearray, memoryview

### Casting
In python, you can convert certain values for certain datatypes to other datatypes. Three of the commonly used functions for this are ‘bool’, ‘int’, and ‘str’.
- convert a string to a boolean
```
tmpstr = "Hello"
tmpbool = bool(tmpstr)
print(tmpbool)
```
convert a decimal to a string
the print function can only print strings
anything that isn't a string needs to be converted to a string using the "str" method, to be able to print it.
```
x = 45.5
print("The number is: " + str(x))
```

Convert string to an integer
```
mystr = "50"
myint = int(mystr) + 1
print(str(myint))
```

### Arithmetic: Addition/subtraction, multiplication/division, assignment arithmetic (assignment operators).
```
if flow
print("\nArithmetic\n")
a = 21
b = 3
```
```
# Addition
print(a + b)

# "c += 5" is shorthand for "c = c + 5"
c = 83
c += 5
print(c)

# Subtraction
print(b - a)

# Division
print(a/b)

# Exponents
print(4**b)
#or
expb = pow(4,b)
print(expb)

# Remainder 
print(4 % 3)

# Absolute value
av = abs(22-32)
print(av)

# Round, Floor, Ceiling
print(round(3.2))
print(int(3.2))
import math
print(math.ceil(3.2))
print(math.floor(3.7))
print(int(3.7))
```
The math package has many common math functions you can use.

PRACTICE: Add and multiply various datatypes together using conversions and math methods. Play around with other math methods. Try the “atan2” method and try the “log” method with a different base than the default.

## Basic Data Structures: Lists and Dictionaries.
Lists
```
if flow
if flow
if flow
print("\nLists\n")

gene_list = ["DDX11L1","WASH7P","MIR6859-1","MIR1302-2HG","MIR1302-2","FAM138A"]

# get the first element in the list, 0-indexed
gene1 = gene_list[0]
print(gene1)

# getting the last element in a list
last_gene = gene_list[-1]
print(last_gene)
# OR
last_gene = gene_list[5]
print(last_gene)

# getting a range of the list
print(gene_list[-3:])
print(gene_list[1:3])
print(gene_list[:3])

# The same range concept works for strings
mystring = "The Quick Brown Fox"
print(mystring[4:9])

# get the length of a list
print("The length of gene_list is " + str(len(gene_list))

# lists can have elements of any type
gene_exp = [43.2, 45, 60.1, 12, 0.5, 23]
expval = gene_exp[2]
print(expval)

# creating a new variable equal to a list does NOT create a copy
# both variables point to the same list
gene_list2 = gene_list
gene_list2[2] = "DMR3"
print(gene_list)

# use the copy method to make a actual copy of a list
gene_list2 = gene_list.copy()
gene_list2[2] = "DMR5"
print(gene_list)
print(gene_list2)

# use the "in" keyword to check for membership in a list
print("BRCA2" in gene_list)
Built-in list methods
Lists have built-in methods which WILL modify the list in place.

gene_list.append("BRCA2")
print(gene_list)

gene_list.remove("WASH7P")
print(gene_list)

gene_list.reverse()
print(gene_list)
Here is a list of functions and methods for lists.
```
PRACTICE: Play with some list methods. Try out the “pop”, “insert”, and “sort” methods. Try to insert a new gene name in the middle of gene_list.

### Dictionaries
A set of key/value pairs where the keys are unique. Used primarily when you want to use a string as an index.
Dictionary values are pointed to by the keys. Values can be anything from int, float, and bool to lists, tuples, and other dictionaries.
dictionary
```
print("\nDictionaries\n")

gene_exp_dict = {"DDX11L1":43.2,"WASH7P":45,"MIR6859-1":60.1,"MIR1302-2HG":12,"MIR1302-2":0.5,"FAM138A":23}

# You can access a value in the dictionary by referencing the key
print(gene_exp_dict["WASH7P"])

# You can overwrite a value
gene_exp_dict["WASH7P"] = 39
print(gene_exp_dict["WASH7P"])

# You can add a new value
gene_exp_dict["BRCA2"] = 100
print(gene_exp_dict)

# dict built-in methods

# returns all the keys
print(gene_exp_dict.keys())

# returns all the values
print(gene_exp_dict.values())

# returns all the key,value pairs
print(gene_exp_dict.items())

# you can also check if a key exists in a dictionary
print("BRCA2" in gene_exp_dict)

# to make a copy of a dictionary, you need to use the "copy" method
gene_exp_dict_copy = gene_exp_dict.copy()
```

 Create a dictionary for gene expression of multiple genes for multiple samples. Populate the dictionary with at least 3 genes for 5 samples. Write the code to get the gene expression for the 2nd gene, 4th sample.
