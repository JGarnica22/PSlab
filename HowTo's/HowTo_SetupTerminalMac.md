# How to Set up your Terminal in Mac

If you are running your analyses in Mac, you will already have installed the _Terminal.app_. Terminal provides a command line interface to control the UNIX-based operating system that controls Mac operating system (macOS).<br/>

Currently, after macOS Catalina, Mac uses _zsh_ as the default login and interactive shell ("language").<br/>

In order to use the Terminal, you need to first learn the basic commands. You can find a summary of **regular expressions** (Appendix 2) and **basic shell commands** (Appendix 3) in the following Cheat Sheet:  
http://practicalcomputing.org/files/PCfB_Appendices.pdf

It is optional, but in order to easily manipulate your Terminal configuration you can install **Oh My Zsh** (https://ohmyz.sh).

To modify your zsh profile (_.zshrc_) you should call from your default directory (/Users/$USER)
````
nano .zshrc
````
This will open the nano editor within Terminal and will allow you to modify your Terminal aesthetics. I suggest you change the zsh theme to "bureau" (ZSH_THEME="bureau") in order to see in each of your command lines:<br/>
$USER@computer directory                [time]  

<Ctrl + X> to exit nano and accept Save.
<br/>

Analysis of NGS data requires many tools, which usually can be installed through **Conda** (https://docs.conda.io/en/latest/). For Mac, _miniconda_ is preferably used. To install _miniconda_ go to https://conda.io/projects/conda/en/latest/user-guide/install/macos.html  
<br/>

After installing _miniconda_ it will initialize by default every time you open Terminal. You can know Conda is active if you see _(base)_ at the beginning of your command line.

To avoid activation by default run the following line in Terminal:

````
conda config --set auto_activate_base false
````

**REMEMBER! After doing this, every time you open your Terminal and you need to use Conda, you will need to initialize**.  
To initialize Conda, use:
````
conda activate
````


