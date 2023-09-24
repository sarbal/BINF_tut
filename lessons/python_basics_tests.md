## Test yourself SOLUTIONS
1. Create a dictionary for gene expression of multiple genes for multiple samples. Call it gene_exp_dict. 

```
gene_exp_dict = {"RAC1P":[0,8,9,0,1],"HM13":[0,5,6,7,9],"RPSAP9":[3,0,1,9,0],"MIR32":[5,0,3,5,6],"GTF3C2":[0,6,2,5,67],"RBPMS":[0,5,43,2,1],"MTHFSD":[0,9,8,0,4],"SCD5":[0,4,4,5,76], "CALM1P2":[0,2,2,45,6]}
```

2. Populate the dictionary with at least 10 genes for 5 samples.

```
gene_exp_dict["WASH7P"] = [1,5,6,1,2]
print(gene_exp_dict["WASH7P"]) 
gene_exp_dict["HM13"]  = [0,0,0,0,1]
gene_exp_dict["MIR32"]  = [10,21,34,0,41]

print(gene_exp_dict.keys())
print(gene_exp_dict.values())
```

3. Write the code to get the gene expression for the 4th gene, 2nd sample.

```
genes = list(gene_exp_dict.keys())
gene_target = genes[3]
sample_target = 1 
print(gene_exp_dict[gene_target][sample_target])
```

Now these are obviously not the only (nor best!) solutions. 
Feel free to share yours as comments or issues or reviews on the github. 

