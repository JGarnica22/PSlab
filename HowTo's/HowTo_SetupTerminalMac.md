# How to Set up your Terminal in Mac

If you working in Mac, you will already have installed the _Terminal.app_. Terminal provides a command line interface to control the UNIX-based operating system that controls Mac operating system (macOS).<br/>

Currently, after macOS Catalina, Mac uses _zsh_ as the default login and interactive shell ("language").<br/>

#
In order to use the Terminal, you need to first learn the basic commands. You can find a summary of **regular expressions** (Appendix 2) and **basic shell commands** (Appendix 3) in the following Cheat Sheet:  
http://practicalcomputing.org/files/PCfB_Appendices.pdf  
<br/>

#
It is optional, but in order to easily manipulate your Terminal configuration you can install **Oh My Zsh** (https://ohmyz.sh).

To modify your zsh profile (_.zshrc_) you should open it from your default directory (/Users/$USER)
````
open .zshrc
````
This will open an external text editor (TextEdit) will allow you to modify your Terminal aesthetics. I suggest you change the zsh theme to "_bureau_" (ZSH_THEME="bureau") in order to see in each of your command lines:<br/>
`$USER@computer   directory                      [time]`  
<br/>

Analysis of NGS data requires many different tools. Many of them can be installed and run through **Conda** (https://docs.conda.io/en/latest/). For Mac, _miniconda_ is preferably used. 

To install _miniconda_ go to https://conda.io/projects/conda/en/latest/user-guide/install/macos.html and follow instructions.  

**IMPORTANT! In the process of installation, you will be asked where to install Miniconda. By default, it will be installed into your default location (/Users/$USER). If you want to change that, specify another location**.
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

