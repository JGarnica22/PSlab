# How to Set up your Terminal in Mac

If you working in Mac, you will already have installed the _Terminal.app_. Terminal provides a command line interface to control the UNIX-based operating system that controls Mac operating system (macOS).<br/>

Currently, after macOS Catalina, Mac uses _zsh_ as the default login and interactive shell ("language").  
<br/>

### 1. Learn how to use the Terminal
In order to use the Terminal, you need to first learn the basic commands. You can find a summary of **regular expressions** (Appendix 2) and **basic shell commands** (Appendix 3) in the following Cheat Sheet:  
http://practicalcomputing.org/files/PCfB_Appendices.pdf  
<br/>

### 2. (Optional) Install _Oh My Zsh_
It is optional, but in order to easily manipulate your Terminal configuration you can install _Oh My Zsh_ (https://ohmyz.sh).

To modify your zsh profile (_.zshrc_) you should open it from your default directory (/Users/$USER)
````
open .zshrc
````
This will open an external text editor (TextEdit) will allow you to modify your Terminal aesthetics. I suggest you change the zsh theme to "_bureau_" (ZSH_THEME="bureau") in order to see in each of your command lines:<br/>
`$USER@computer   directory                      [time]`  
<br/>

### 3. Install Miniconda
Analysis of NGS data requires many different tools. Many of them can be installed and run through **Conda** (https://docs.conda.io/en/latest/). For Mac, _miniconda_ is preferably used. 

To install _miniconda_ go to https://conda.io/projects/conda/en/latest/user-guide/install/macos.html and follow instructions.  

**IMPORTANT! In the process of installation, you will be asked where to install Miniconda. By default, it will be installed into your default location (/Users/$USER). If you want to change that, specify another location**.  

For example: /Users/$USER/Applications/miniconda3  
(\*note that the path you specify must exist: all folders might be already present, but the last miniconda3 folder that will be created by the installtion. I recommend you create an _Applications_ folder in you user where you install all your command-line tools) 
<br/>

After installing _miniconda_ it will initialize by default every time you open Terminal. You can know Conda is active if you see _(base)_ at the beginning of your command line. To avoid activation by default run the following line in Terminal:

````
conda config --set auto_activate_base false
````

**REMEMBER! After doing this, every time you open your Terminal and you need to use Conda, you will need to initialize**.  
To initialize Conda, use:
````
conda activate
````

To deactivate Conda, run:
````
conda deactivate
````
<br/>

### 4. Install tools from Conda
With miniconda you will be able to install many different tools required for NGS data analysis (e.g. `deeptools`, `samtools`, `bedtools`).

In order to install, search the tool in https://anaconda.org and you will find which channel ("source") it can be installed from (e.g. bioconda) and the code to do so. If you want to install, for example, `samtools`, the following command should be used:
````
conda install -c bioconda samtools
````  

After you install the tool, you can run it directly when conda is activated:
````
samtools <command> [options]
````  
<br/>

In order to check which tools you have installed within conda, you can use:
````
conda list
````
<br/>

### 5. Install other tools outside from Conda
Tools that cannot be installed from Conda will have their on installation protocol. Following the installation manual is recommended.

#### EXAMPLE: How to install STAR (v. 2.7.3a)

1. Change directory to your _Applications_ folder
````
cd /Users/$USER/Applications
````

2. Download package from source (GitHub) and unzip:
````
wget https://github.com/alexdobin/STAR/archive/2.7.3a.tar.gz
tar -xzf 2.7.3a.tar.gz
rm *.tar.gz
cd STAR-2.7.3a
````

3. Compile tool:
````
cd source
make STARforMacStatic CXX=/path/to/gcc #e.g. /usr/local/Cellar/gcc/9.2.0_2/bin/g++-9  
````

4. In order to be able to run the tool from any directory without having to locate the executable file, you need to add the directory to your $PATH:  

**IMPORTANT! This needs to be done every time you install a tool**  

##### #How to know what is in your $PATH:
````
echo $PATH
````

##### #You must add the /path/to/yourtool into your $PATH. You can do that in your _.zshrc_ profile.
````
open .zshrc`
````
#Add /path/to/yourtool into your $PATH and save. For example:
PATH=$PATH:/Users/patri/Applications/STAR-2.7.3a/bin/MacOSX_x86_64
#Close and open Terminal to update changes

5. Run STAR:  
````
STAR  [options]... --genomeDir /path/to/genome/index/   --readFilesIn R1.fq R2.fq
````



