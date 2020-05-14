# How to install R and rstudio in WLS ubuntu :computer:

As any other tools you can install `r-base` and `rstudio` in your ubuntu terminal in WLS. This will be useful in case you want to have all in the same system avoiding moving files from linux to windows, or to use some `R` features not available for windows such as Bigwig files import in `trackViewer` package.

R is a programming language commonly used for statistical computing and graphical representation of data. RStudio is a set of integrated tools designed to help R developers be more productive.


## Install R :fishing_pole_and_fish:
From here you can go ahead and install R, but if you use the default Ubuntu repository you'll get an old version of R (R 3.2.3, from 2015). You probably want the latest version of R, so add CRAN as a new package repository for Ubuntu. You'll need to run these commands as root, so enter the password you created previously here if requested:

To install **R 3.6 packages**, use:
````
sudo echo "deb http://cloud.r-project.org/bin/linux/ubuntu xenial/" 
````

To install **R 4.0 packages** (only allowed for ubuntu 19.10 or superior), use:
````
sudo echo https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/
````
Then  update your applications source list
````
sudo tee -a /etc/apt/sources.list
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9
sudo apt-get update
````

Now you're all set to install the latest version of R, which can be done with:
````
sudo apt-get install r-base
````
In case you encounter problems with packages or older versions, use a more interactive installation debian tool: `aptitude` and be able to spot and solve possible problems, as well as removing unnecessary packages:
````
sudo apt-get install -y aptitude
sudo aptitude install r-base
````

And that's it! (Once all the dependencies install, which can take a while the first time.) Now you're all ready to run R from the Linux command line! 
Keep in mind that like with other ubuntu applications you will need to have all your R files in the linux subsystem folders.

If you want to start R in your terminal just run:
````
R
````




## Install RStudio :surfer:
First, visit the  [RStudio downloads page](https://rstudio.com/products/rstudio/download/#download) to grab the latest release of `RStudio` for your Debian based Linux distribution, then you can manually download it (and move it to your linux folders) or use wget to get the version you want:

````
--- Ubuntu 18.04 / Linux Mint 19 / Debian 10 ---
sudo apt -y install wget
wget https://download1.rstudio.org/desktop/bionic/amd64/rstudio-1.2.5042-amd64.deb
sudo apt install ./rstudio-1.2.5042-amd64.deb

--- Ubuntu 16.04 / Linux Mint 18 ---
sudo apt -y install wget
wget https://download1.rstudio.org/desktop/xenial/amd64/rstudio-1.2.5042-amd64.deb
sudo apt install ./rstudio-1.2.5042-amd64.deb

--- Debian 9+ ---
wget https://download1.rstudio.org/desktop/debian9/x86_64/rstudio-1.2.5042-amd64.deb
sudo apt install ./rstudio-1.2.5042-amd64.deb
````

If you encounter any dependency problems, use:
````
sudo apt -f install "your_rstudio_version"
````

Finally, to start enjoying your ubuntu-RStudio you will need to visualize your data and code in your default browser, to do so run:

````
sudo rstudio-server start
````
Then open you explorer (in windows system) and visit this page: **http://localhost:8787/**

It will ask you for a user and password and you will be ready to go!

You will need to do these final two steps every time you want to run `RStudio` from linux WLS, but you will not need to be running R in your terminal at the moment.
