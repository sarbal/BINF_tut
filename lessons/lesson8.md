# Week 9: Snakes on a plane 
## Objectives 
This week's tutorial is on flow control and data wrangling in python. You will learn: 
- Flow control statements
- Functions and vectorised computations

## Setting up
Start a new notebook. Save the file as "yourname_week9.ipynb". 
As before, copy the code into your notebook as chunks. 

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

