# How to Set up your Terminal in Windows Subsystem for Linux (WLS) :penguin:

If you need to work with Terminal language and tools and your computer is running on windows system you need to enable 
Linux subsystem (WLS) in order to get access to tools only available for Linux and macOS. Moreover, scripts include in this reservoir (***PS***) will only be written in UNIX language (Linux and macOS).<br/>

On top of that, The Windows Subsystem for Linux lets developers run a GNU/Linux environment -- including most command-line tools, utilities, and applications -- directly on Windows, unmodified, without the overhead of a virtual machine.

With WLS you can:

1. Choose your favorite GNU/Linux distributions from the Microsoft Store.
2. Run common command-line free software such as `grep`, `sed`, `awk`, or other ELF-64 binaries.
3. Run Bash shell scripts and GNU/Linux command-line applications including:
    - Tools: vim, emacs, tmux
    - Languages: Javascript/node.js, Ruby, Python, C/C++, C# & F#, Rust, Go, etc.
    - Services: sshd, MySQL, Apache, lighttpd
4. Install additional software using own GNU/Linux distribution package manager.
5. Invoke Windows applications using a Unix-like command-line shell.
6. Invoke GNU/Linux applications on Windows.

In this How-to we will do:<br/>
- Set up **WSL** (with optional **ZSH** and **Oh-my-zsh**) for Windows 10
- Install and set up **Anaconda**
- Install some tools from Anaconda
- Install tools from outside Anaconda  


## 1. Setting up WSL :wrench:
First off we have to do some preliminary setup to get WSL working:<br/>
### 1.1. Turning on Developer mode
First head to the developer settings in your Windows 10 settings menu:<br/>
Settings ➡ Update & Security ➡ For developers ➡ Developer mode

Or simply do a search for **“dev”** and click on **“Developer settings”**<br/>
![](https://miro.medium.com/max/1400/1*TQxU2JgHp2eEqw8Qu7Qeow.png)

Here you might need to sign into an admin account or get temporary access to turn on the Developer mode:<br/>

![](https://miro.medium.com/max/700/1*Lbk8X5xuctAOzDcwoo5bPA.png)

Once enabled we need to turn on the Windows Subsystem for Linux feature.<br/>

### 1.2. Turning on Windows Subsystem for Linux
This menu can be accessed as such by going to:<br/>

Control Panel ➡ Programs ➡ Turn Windows features on and off

Or by searching *“windows features”* and selecting *“Turn Windows features on or off”*.

![](https://miro.medium.com/max/700/1*KgnlItWjj4d525gmay_g_A.png)

Next make sure the feature “Windows Subsystem for Linux” is ticked:<br/>

![](https://miro.medium.com/max/413/1*f7vMExOir3iPyfbNGcJ8Tw.png)

### 1.3. Restart your computer

### 1.4. Install Ubuntu :minidisc:
Head to the Microsoft Store and search for **“Ubuntu”**, select the App you prefer and install it. This app will enable to run Linux commands in  your Windows system.

![](https://miro.medium.com/max/1070/1*KspvSBty03M8zl6nl9bisQ.png)

### UNIX language :u6307:
Now you have a fully functional Linux Subsystem with full admin rights, in case you are not experienced, to use the Terminal, you need to first learn the basic commands. You can find a summary of **regular expressions** (Appendix 2) and **basic shell commands** (Appendix 3) in the following Cheat Sheet:  
http://practicalcomputing.org/files/PCfB_Appendices.pdf  
<br/>
**IMPORTANT! As you will be using a subsystem inside windows `Ubuntu` will create a hidden folder only found by search on windows explorer in case you want to see how look like via explorer. It is highly recommended to explore Linux directories and files on `Ubuntu` app. Very importantly, you must not edit directories or files from the Linux subsystem using Windows explorer. Conversely, you can edit windows files from `Ubuntu` app. In case, you need to, for instance, move file from windows system to Linux subsystem you can explore Window directories in Terminal via `/mnt/UNIT/USERS/`.**

## 2. (Optional) Install _Oh My Zsh_ :necktie:
It is optional, but in order to easily manipulate your Terminal configuration you can install `Oh-My-Zsh` following these instructions:<br/>

### 2.1. Install zsh

Open the Ubuntu app installed from the App Store. We will now install zsh:

````
sudo apt-get install zsh
````

After installing it, type `zsh`, zsh will ask you to choose some configuration. We will do this later on while installing `oh-my-zsh`, so choose option `0` to create the config file and prevent this message to show again.

### 2.2. Installing oh-my-zsh
Before all we need to have git installed:
````
sudo apt-get install git
````

Then, use `curl` to install `oh-my-zsh`:

````
sh -c "$(curl -fsSL https://raw.githubusercontent.com/robbyrussell/oh-my-zsh/master/tools/install.sh)"
````

This will clone the repo and replace the existing `~/.zshrc` with a template from `oh-my-zsh`.

![](https://blog.joaograssi.com/content/images/2018/04/1---Install-oh-my-zsh.PNG)

### 2.3. Configuring zsh/oh-my-zsh
First, we need to make sure `zsh` is executed by default for **Bash** on Ubuntu. This is not mandatory, but if not done you need to type zsh every time. For this, edit the `.bashrc` file with nano: `nano ~/.bashrc` and paste this right after the first comments:

````
if test -t 1; then
exec zsh
fi
````

Save it `Ctrl + shift X` and **restart** your Ubuntu shell. You should be on `zsh` by default now.

### 2.4. Changing the Theme of oh-my-zsh
`oh-my-zsh` has several nice Themes. It's worth checking them out, you have some here: https://github.com/ohmyzsh/ohmyzsh/wiki/Themes 

In order to change your theme you need to modify your zsh profile (`.zshrc`). To do so use again nano: `nano ~/.zshrc`

Edit the `~/.zshrc` again with nano: `nano ~/.zshrc`, and change `ZSH_THEME`

````
# Find and change this
ZSH_THEME="robbyrussell"

# To this, for instance
ZSH_THEME="agnoster"
````

Save it and restart your Ubuntu shell again to see the changes.



## 3. Install and set up Anaconda :snake:

Analysis of NGS data requires many different tools. [Anaconda](https://docs.conda.io/en/latest/) is one of the most widely used solutions for package management for R and Python, it comes with 1.500+ popular packages out of the box which is more than enough for 99% of all Data Science related tasks.

### 3.1. Download Anaconda
Start up your WSL and **download** Anaconda:

````
wget https://repo.anaconda.com/archive/Anaconda3-2019.03-Linux-x86_64.sh
````

### 3.2. Installing Anaconda
Next, execute the downloaded file to **install** Anaconda:
````
Anaconda3-2019.03-Linux-x86_64.sh
````

*Note: Use `bash Anaconda3–2019.03-Linux-x86_64.sh` if using ZSH and you’re running into problems.*

**IMPORTANT! In the process of installation, you will be asked where to install Anaconda. By default, it will be installed into your default location (/home/$USER). If you want to change that, specify another location**.

After following the onscreen instructions to install Anaconda, simply remove the installation file:

````
rm Anaconda3-2019.03-Linux-x86_64.sh
````

### 3.3. Updating Anaconda
Now you should be able to start up your Anaconda **environment**:

````
source ~anaconda3/bin/activate
````

*“~anaconda3/bin/activate”* is default place that Anaconda will install itself but if you chose elsewhere simply point to that directory.

Once activated, initiate a full update:
````
conda update --all
````

### 3.4. (Optional) Turn off Anaconda as your defect environtment
After installing `Anaconda` it will initialize by default every time you open Terminal. You can know Conda is active if you see _(base)_ at the beginning of your command line. To avoid activation by default run the following line in Terminal:

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

### 4. Install tools from Conda :cd:
With Anaconda you will be able to install many different tools required for NGS data analysis (e.g. `deeptools`, `samtools`, `bedtools`).

In order to install, search the tool in https://anaconda.org and you will find which channel ("source") it can be installed from (e.g. bioconda) and the code to do so (*Remember that may the code be different for Linux and macOS installation, you must use **Linux***). For WLS you need to previously create an **environtment** for the tools and then include them in the `PATH` if you do not want needing to activate them every time.

If you want to install, for example, `samtools`, the following command should be used:

Create environtment

This need to be done for tools giving this problem when installing:
````
Solving environment: failed with initial frozen solve. Retrying with flexible solve
````
To solve that create and environtment for the tools before installation:

````
conda create -n samtools
````
Then activate this recently created environtment
````
conda activate samtools
````
Now you can install your tool with
````
conda install -c anaconda samtools
````  
Now you can run the tool when conda is activated after activating its environtment:

````
activate samtools
````
Confirm that tools is properly installed by running help command
````
samtools --help
````
In case you do not want to activate tool environtment each time include this tool in your `PATH` by adding this command at the final of your .zshrc file
````
nano ~.zshrc
export PATH=/home/USER/anaconda3/envs/samtools/bin:$PATH 
````


## 5. Install other tools outside from Anaconda :floppy_disk:
Tools that cannot be installed from Conda will have their on **installation protocol**. Following the installation manual is recommended.
<br/>

### EXAMPLE: How to install `STAR` (v. 2.7.3a)

1. Change directory to the folder you are using to store applications
````
cd /Users/$USER/*
````
<br/>

2. **Download** package from source (*GitHub*) and **unzip**:
````
wget https://github.com/alexdobin/STAR/archive/2.7.3a.tar.gz
tar -xzf 2.7.3a.tar.gz
rm  2.7.3a.tar.gz

cd STAR-2.7.3a
````
<br/>

3. **Compile** tool under Linux:
````
cd STAR/source
make STAR  
````

**NOTE: if application do not find make command, you may still need to install all the basic packages to build code and zlib library, then repeat the previous command.**
````
sudo apt-get install build-essential
sudo apt install zlib1g-dev
````

4. In order to be able to run the tool from any directory without having to locate the executable file as done previously with conda, you need to add the directory to your $PATH on your .zshrc file. **IMPORTANT! This would need to be done every time you install a tool**  
````
export PATH=/home/USER/STAR-2.7.3a/bin/Linux_x86_64:$PATH 
````

5. Run `STAR`:  
````
STAR  [options]... --genomeDir /path/to/genome/index/   --readFilesIn R1.fq R2.fq
````
