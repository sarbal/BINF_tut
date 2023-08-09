# Week 8: Data Viz
### Download data 
First things first! Download these palmer penguins data into your working directory. 
Into a terminal, type in: 
```
pip -q install palmerpenguins
```

### Installing and importing libraries 
```
pip install numpy
pip install pandas
pip install seaborn
```
Then, run an interactive python session in a terminal. In your python session:
```
import numpy as np
import pandas as pd
import seaborn as sns
```
- Pandas is a library for working with tabular data. Based on the R data.frame library.
- Seaborn is a visulization package. 

### Importing data 
```
from palmerpenguins import load_penguins
df = load_penguins()
print(type(df))
df
```  

