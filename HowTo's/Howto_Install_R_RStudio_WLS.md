# How to install R and RStudio in WLS Ubuntu :computer:

As any other tools you can install `r-base` and `rstudio` in your Ubuntu terminal in WLS. This will be useful in case you want to have all in the same system avoiding moving files from Linux to Windows, or in order to use some `R` features that not available for windows, such as Bigwig files import in `trackViewer` package.

R is a programming language commonly used for statistical computing and graphical representation of data. RStudio is a set of integrated tools designed to help R developers be more productive.


## Install R :fishing_pole_and_fish:
From here you can go ahead and install R, but if you use the default Ubuntu repository you'll get an old version of R (R 3.2.3, from 2015). You probably want the latest version of R, so add CRAN as a new package repository for Ubuntu. If you have an Ubuntu version 18.04 or older you can download R 3.6, but if you have a newer Ubuntu version you can download R 4.0. You can explore [CRAN website](https://cran.r-project.org/bin/linux/ubuntu/) to choose the best version for you.

You'll need to run these commands as root, so as you type `sudo` before any command it will be run as root or administrator, and thus it will ask for your account's password you created when you initiated Ubuntu. For instance you could use these repositories:

To install **R 3.6 packages**:
````
sudo echo "deb http://cloud.r-project.org/bin/linux/ubuntu xenial/" 
````

To install **R 4.0 packages** (only allowed for ubuntu 19.10 or superior):
````
sudo echo https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/
````
Then update your applications source list
````
sudo tee -a /etc/apt/sources.list
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9
sudo apt-get update
````

Now you're all set to install the latest version of R, which can be done with:
````
sudo apt-get install r-base
````  
<br/>

\* In case you encounter problems with packages or older versions, use a more interactive installation debian tool: `aptitude`. You will be able to spot and solve possible problems, as well as remove unnecessary packages:
````
sudo apt-get install -y aptitude
sudo aptitude install r-base
````

**And that's it! (Once all the dependencies are installed, which can take a while the first time.) Now you're all ready to run R from the Linux command line!** :+1:  

:bulb: Keep in mind that like with other Ubuntu applications, you will need to have all your R files in the Linux subsystem folders.  
<br/>

If you want to start R in your terminal just run:
````
R
````  
<br/>

## Install RStudio :surfer:
First, visit the  [RStudio downloads page](https://rstudio.com/products/rstudio/download/#download) to grab the latest release of `RStudio` for your Debian based Linux distribution (.deb), then you can manually download it (and move it to your Linux folders) or use `wget` to get the version you want, for instance:

Install `wget` package in case you do not already have it:
````
sudo apt -y install wget
````

Then, download `Rstudio` and install it:
````
wget https://download1.rstudio.org/desktop/bionic/amd64/rstudio-1.2.5042-amd64.deb
sudo apt install ./rstudio-1.2.5042-amd64.deb
````

If you encounter any dependency problems, use:
````
sudo apt -f install "your_rstudio_version"
````

Finally, to start enjoying your Ubuntu-RStudio you will need to visualize your data and code in your default browser, to do so run:

````
sudo rstudio-server start
````
Then open your explorer (in Windows system) and visit this page: **http://localhost:8787/**

It will ask you for a user and password and you will be ready to go!

You will need to do these final two steps every time you want to run `RStudio` from Linux WLS, but you will not need to be running R in your Terminal at the moment.
