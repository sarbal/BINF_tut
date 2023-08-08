# Week 8: Data Viz
First things first! Download these files into your working directory: 
- [lesson](../data/lesson)
- [helper.py](../data/helper.py)
###

### Download data 
Into a terminal, type in: 
```
pip -q install palmerpenguins
```

Then, run an interactive python session:
```
import numpy as np
import pandas as pd
# Visulization package 
import seaborn as sns
```

- Pandas is a library for working with tabular data. Based on the R data.frame library.
```
from palmerpenguins import load_penguins
df = load_penguins()
print(type(df))
df
```  
