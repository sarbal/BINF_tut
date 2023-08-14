
### How to install a package from source 
If the installation fails, you can try an older/different version of the package by downloading it and then installing it from that source:
```
install.packages(filename, repos = NULL, type="source")
```
Where filename would represent the full path and file name of the package you've downloaded: 
```
path_to_file = "X:/sm_2.2-5.5.zip"
```

Or you can do it directly from the URL (make sure the download method is specified and check it works, depending on your available libraries). 
```
install.packages(url, repos = NULL, method="libcurl") # could also try method="wget"
```
And the URLs could be something like:
```
# Windows
url = "https://cran.r-project.org/bin/windows/contrib/3.6/sm_2.2-5.5.zip" 
# MAC
url = "https://cran.r-project.org/bin/macosx/el-capitan/contrib/3.5/sm_2.2-5.4.tgz"
```

### How to install a package through conda
One useful setup is to use "conda", a package manager and environment management system. Anaconda is the Python and R distribution. Go to the intro page to install it. Then, when you need to install a package or update a tool, you run:  
```
conda install r-XXXX
```
