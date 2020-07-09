# How to use Cluster :video_game:

This protocol describes how to use and create jobs on **Nord III cluster** supercomputer processors. These are based on Intel SandyBridge processors, iDataPlex Compute Racks, a Linux Operating System and an Infiniband interconnection. Nord3 is a part of the old MareNostrum3. Nord3 machine is installed at BSC (Barcelona Supercomputing Center) and it presents the following configuration:

* 84 IBM dx360 M4 compute nodes with 128 GB of RAM per node
* Infiniband FDR10 interconection network
* SLES 11 SP3
* GPFS Common Storage (home, projects, scratch)

For more information and details check [**Nord III User's Guide**](https://www.bsc.es/user-support/nord3.php). You can also contact with support@bsc.es for further information.

## 1. Connect to Nord III
Open your terminal and use Secure Shell (ssh) tools to login. Other incoming connections such as  telnet, ftp, rlogin, rcp, or rsh are not accepted.
Three login blades are available to access to Nord machine {nord1,nord2,nord3}.bsc.es. Our username is **cek26664** and the users to log in are:


* cek26664@nord1.bsc.es
* cek26664@nord2.bsc.es
* cek26664@nord3.bsc.es


Our password is: ********

Use the following command and you will be asked to type the password:

````
ssh cek26664@nord1.bsc.es
````

Once connected to the machine, you will be presented with a UNIX shell prompt and you will normally be in your home ($HOME) directory. The login nodes serve as front ends and are used typically for editing, compiling, preparation and submition of batch executions. It is not permitted the execution of cpu-bound programs on these nodes,  if some execution needs more cputime than the permitted, this needs to be done through the batch queue system (LSF).

## 2. Password Management
In order to change the password, you have to login to a different machine (dt01.bsc.es). This connection must be established from your local machine.

    % ssh -l username dt01.bsc.es

    username@dtransfer1:~> passwd
    Changing password for username.
    Old Password: 
    New Password: 
    Reenter New Password: 
    Password changed.
Mind that that the password change takes about 10 minutes to be effective.

## 3. File transfer
There are two ways to copy files from/to the Cluster:

* Direct scp or sftp to the login nodes
* Using a Data transfer Machine which shares all the GPFS filesystem for transferring large files

We strongly recommend using **[Cyberduck](https://cyberduck.io/download/)** as a data transfer machine.

### 3.1. How to use Cyberduck
After having installed Cyberduck start a new conexion, then indicate connexion as SFTP. Indicate one of the blades as server  **{nord1,nord2,nord3}.bsc.es** and finally insert our user name cek26664 and password.

Once linked you will be able to easily upload, download and modify files and directories in the cluster machine.

## 4. Modules environtment
The cluster provides a dynamic modification of user's environment via **modulefiles**. Each modulefile contain the information needed to configure the shell for an application or compilation. Typically modulefiles instruct the module command to alter or set shell environment variables such as PATH, MANPATH, etc. The modules needed for a tool or certain job must be loaded beforehand in order to proceed with the job.
### 4.1. Modules tool usage
In order to check the module available at the moment use:
````
module av
module avail
````
Modules can be invoked in two ways: by name alone or by name and version. Invoking them by name implies loading the default module version. This is usually the most recent version that has been tested to be stable (recommended) or the only version available. To invoke modules use:
````
module load <modulename(version)>
````
Other important commands to work with modules include:
````
module list 
````
Show all the loaded modules
````
module purge
````
Removes all the loaded modules
````
module unload <modulename> 
````
Removes all environment changes made by module load command
````
module switch <oldmodule> <newmodule> 
````
Unloads the first module (oldmodule) and loads the second module (newmodule)


## 5. GPFS Filesystem storing space
The IBM General Parallel File System (GPFS) is a high-performance shared-disk file system providing fast, reliable data access from all nodes of the cluster to a global filesystem. GPFS allows parallel applications simultaneous access to a set of files (even a single file) from any node.

These are the GPFS filesystems available in the machine from all nodes:

/apps: Over this filesystem will reside the applications and libraries that have already been installed on the machine. Take a look at the directories to know the applications available for general use.

* **/gpfs/home:** This filesystem has the home directories of all the users, and when you log in you start in your home directory by default. Every user will have their own home directory to store own developed sources and their personal data. A default quota will be enforced on all users to limit the amount of data stored there. Also, it is highly discouraged to run jobs from this filesystem. Please **run your jobs on your group’s /gpfs/projects or /gpfs/scratch** instead.

* **/gpfs/projects:** In addition to the home directory, there is a directory in /gpfs/projects for each group of users. For instance, the group bsc01 will have a /gpfs/projects/bsc01 directory ready to use. This space is intended to store data that needs to be shared between the users of the same group or project. A quota per group will be enforced depending on the space assigned by Access Committee.  
* **/gpfs/scratch:** Each user will have a directory over /gpfs/scratch. Its intended use is to store temporary files of your jobs during their execution. A quota per group will be enforced depending on the space assigned.

In order to assess how much space we are provided with in our quota use:
````
bsc_quota
````

**IMPORTANT:** An incremental backup will be performed daily **only for /gpfs/home.**

## 6. Submitt a job to the cluster
### 6.1. Write your bash script
You can write your code in a .txt file and then upload it to the cluster to execute it. Before coding, make sure that the file is saved as .sh (bash), is executable and in UNIX format.

Command to convert file to executable:
````
chmod u+x <scriptname.sh>
````
Command to convert file to UNIX:
````
sed -i -e 's/\r$//' <scriptname.sh>
````

After this you can start writing your code. This must always start indicating that is a bash script with **#!/bin/bash**

Then you need to indicate all the parameters and directives of this job in the cluster using **#BSUB**. Here you can see the parameters which must be always included, for more check out [IBM support webpage.](https://www.ibm.com/support/knowledgecenter/SSETD4_9.1.2/lsf_command_ref/bsub.1.html) or read the bsub command’s manual from any of Nord’s terminals by using:
````
man bsub
````

-------------------------------------------------------------
**Most commond directives to be included at the beginning of the script:**

````
#BSUB -J job_name
````
Specify the name (description) of the job.

    #BSUB -q debug
Specify the queue for the job to be submitted. The debug queue is only intended for small tests, so there is a limit of 1 job per user, using up to 64 cpus (4 nodes), and one hour of wall clock limit. The queue might be reassigned by LSF internal policy, as with the sequential queue.
NOTE: to know the queues availables for sending jobs use: `bsc_queues`

    #BSUB -W HH:MM
Specify how much time the job will be allowed to run. This is a mandatory field. NOTE: take into account that you can not specify the amount of seconds in LSF. You must set it to a value greater than the real execution time for your application and smaller than the time limits granted to the user. Notice that your job will be killed after the elapsed period.

    #BSUB -cwd pathname
The working directory of your job (i.e. where the job will run).If not specified, it is the current working directory at the time the job was submitted.

    #BSUB -e/-eo file
The name of the file to collect the stderr output of the job. You can use %J for job_id. -e option will APPEND the file, -eo will REPLACE the file.

    #BSUB -o/-oo file
The name of the file to collect the standard output (stdout) of the job. -o option will APPEND the file, -oo will REPLACE the file.

    #BSUB -n number
The number of tasks for the job. In MPI executions corresponds to the number of MPI processes and for sequential executions the number of cores.

    #BSUB -x
Use the nodes exclusively. This is the default behaviour except for sequential executions.

    #BSUB -M number


Then you can start writing your code, remember that if you want to add comments to your code use `# `  with a blank space afterwards.

### 6.2. Submit your job
Once your script is completed uploaded to your working directory and then just use this command:
````
bsub < <scriptname.sh>
````
You can also check out the status of your jobs and jobs from your group by using respectively:
````
bjobs
bsc_jobs
````
Finally, if you want remove a job from the queue or cancel de process use:
````
bkill <job_id>
````