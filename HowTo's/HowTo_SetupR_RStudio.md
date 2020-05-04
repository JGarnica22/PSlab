# How to Set up R and RStudio :registered:

**R** is a powerful open source program for statistics and graphics. It can run on pretty much any computer and has a very active and friendly support community online.<br/>

Instead of working with R directly, we suggest to use **RStudio**. This is an integrated development environment (IDE) for R. It’s basically a nice front-end for R, giving you a console, a scripting window, a graphics window, and an R workspace, making it easier to see what you are doing and the figures and data you generate.  
<br/>


## 1. Install R and RStudio

Installing R on Windows or Mac is very straightforward. The easiest way is to install it through CRAN.<br/>

To download and install R and RStudio follow these steps (bear in mind that <ins>RStudio requires that R is installed beforehand</ins>).  


### Install R <br/>

1. Download R from http://cran.us.r-project.org/ (click on “Download R for Windows” or "Download R for Mac OS X).
2. In the subdirectory "base" chose "Install R for the first time".
3. You will obtain an executable file.
4. Run the executable file. Leave all default settings in the installation options, so click the button 'Next' until the process is complete.  

### Install RStudio <br/>

1. Download RStudio from http://rstudio.org/download/desktop and install it ("Download RStudio Desktop"). 
2. You will obtain an executable file.
3. Run the executable file with the default settings by clicking 'Next', and wait until the installation finishes.    
<br/>

**FROM THIS POINT ON, EVERYTHING WILL BE EXPLAINED ASSUMING YOU ARE WORKING IN RSTUDIO**


## 2. Install R packages <br/>

R packages are collections of functions and data sets developed by the community. A package will include code, documentation for the package and the functions inside, some tests to check everything works as it should, and data sets. You could also find out what the package does, who the author is, what version the documentation belongs to, the date, the type of license its use, and the package dependencies.<br/>

#### _What Are Repositories?_ <br/>

A repository is a place where packages are located so you can install them from there. Three of the most popular repositories for R packages are:

**1. CRAN:** the official repository, it is a network of ftp and web servers maintained by the R community around the world. The R foundation coordinates it, and for a package to be published here, it needs to pass several tests that ensure the package is following CRAN policies.<br/>
**2. Bioconductor:** this is a topic-specific repository, intended for open source software for bioinformatics. As CRAN, it has its own submission and review processes.<br/>
**3. Github:** this is not R specific, Github is probably the most popular repository for open source projects. Its popularity comes from the unlimited space for open source, the integration with git, a version control software, and its ease to share and collaborate with others. But be aware that <ins>there is no review process associated with it</ins>!  
<br/>

### Install from CRAN:<br/>

1. Run the following function with the name of the package to be installed:
````
install.packages("package")
````
To install more than a package at a time, just write them as a character vector in the first argument of the function:
````
install.packages(c("package_1", "package_2"))
````
<br/>

### Install from Bioconductor: <br/>

1. Executing the following script:
````
source("https://bioconductor.org/biocLite.R")
````
2. Run the following function with the name of the package to be installed:
````
biocLite("package")
````
To install more than a package at a time, just write them as a character vector in the first argument of the function:
````
biocLite(c(package_1", "package_2))
````

**Load packages in RStudio**<br/>

After a package is installed, you are not able to work with its functionalities unless you call it before using a command to use the package:
````
library("package")
````
Another good reminder is that using the command without arguments, it will provide you the list of packages installed in different libraries on your computer :
````
library()
````

**Check and remove packages in RStudio**<br/>

To check what packages are installed on your computer, you can use:
````
installed.packages()
````

Uninstalling a package is straightforward with the following function:
````
remove.packages("package")
````

You can check what packages need an update with a call to the function:
````
old.packages()
````

You can also update all packages by using:
````
update.packages()
````
