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
import matplotlib.pyplot as plt
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
### Summary statistics and properties 
```
df.describe()
df.dtypes
df.columns
```
### Indexing/accessing data in a data frame 
```
df.values
df.loc[i] 
df.iloc[i,j]
df[['bill_length_mm','island']]
df.query("year > 2007")

```

### Plotting 
```
sns.histplot(df['bill_length_mm'],kde=True,bins=20)
plt.show() 
```

```
sns.jointplot(data=df, x="bill_length_mm", y="bill_depth_mm")
plt.show() 
```

```
sns.pairplot(df)
plt.show() 
```

```
g = sns.boxplot(x = 'island',
            y ='body_mass_g',
            hue = 'species',
            data = penguins,
            palette=['#FF8C00','#159090','#A034F0'],
            linewidth=0.3)
g.set_xlabel('Island')
g.set_ylabel('Body Mass')
plt.show() 
```


```
g = sns.lmplot(x="flipper_length_mm",
               y="body_mass_g",
               hue="species",
               height=7,
               data=penguins,
               palette=['#FF8C00','#159090','#A034F0'])
g.set_xlabels('Flipper Length')
g.set_ylabels('Body Mass')
plt.show() 
```

```
import matplotlib.pyplot as pltt
fig ,ax = pltt.subplots(figsize=(15,12), ncols=2,nrows=2)
sns.swarmplot(data=df,x='species',y='body_mass_g',ax=ax[0,0],hue='species')
sns.violinplot(data=df,x='species',y='body_mass_g',ax=ax[0,1])
sns.boxplot(data=df,x='species',y='body_mass_g',ax=ax[1,0])
sns.barplot(data=df,x='species',y='body_mass_g',ax=ax[1,1])
pltt.show()
```


See the documentation at https://pandas.pydata.org/docs/user_guide/io.html


### Show the help message for this function.
```
df.to_csv?
```     

### Save to a file, load from file 
```
df.to_csv("mytaxis.csv")
df = pd.read_csv("mytaxis.csv")
df.head()
```    

### 
```
df["color"].value_counts()
df.groupby("color")["tip"].mean()
df.groupby("color").agg({"tip": "mean"})    
df.groupby("color")["tip"].describe()
```


Test yourself!
1. Imagine a gene expression dictionary like the one from the Dictionary section, except arbitrary keys and expression values. Write an if statement (or statements) that do the following:
2. if both SYF2 and FBX04 exist in the dictionary and both are upregulated (i.e. values > 0), then print “GO:1”
3. Otherwise, if ATF2 does not exist in the dictionary or PLK1 exists and is downregulated (i.e. value < 0), then print “GO:2”
4. Otherwise, if HUS1B exists and HUS1B is downregulated or SYF2 does not exist, then print “GO:3”
5. Otherwise, print “GO:4”


Solutions: Next week!

Back to the [homepage](../README.md)


