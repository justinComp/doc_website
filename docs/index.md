# Simulation Tutorial 

    For Lampe/Letteri/Taylor Computational Group. Created by Justin Kim (zqp6mj)


## WELCOME!

Welcome to the Computational group! You join this project mainly because you like biomaterials, coding or just curious about MD.

Before you start, I want you to download: 

For people on Mac/Windows, everything should be the same!

Download [Visual Molecular Dynamics](https://www.google.com/url?q=https://www.ks.uiuc.edu/Research/vmd/vmd-1.9.3/&sa=D&source=docs&ust=1740806505460574&usg=AOvVaw3kcYOKcLUVKk_HkS8Yabzt). 

    VMD is mainly used to create our peptide systems using TCL Scripts

Download [Visual Studio Code](https://code.visualstudio.com/download). 

    VSCode is mainly used as a text editor for PDB, PSF, TXT files and IDE for many of the programming languages including Python, TCL, and jupyter notebook
    If you need help with VSCode and how to work around it just let me (Justin) know 

Download [Globus](https://www.globus.org/)

    This will allow you to transfer files from your local computer and Rivanna

Download [Avogadro](https://two.avogadro.cc/index.html)
    
    Like Chem Draw, this program makes our molecules into PDB files 

Go to [CHARMM.GUI](https://charmm-gui.org/) and ask for an account.
    
    Make our PDB files CHARMM Force Field Friendly ands make our psf files.


Also we will be using [Rivanna](https://www.rc.virginia.edu/userinfo/hpc/) for our simulations. They have really tutorials and have friendly help ticket helpers. We will specifically be using [OpenOnDemand](https://www.rc.virginia.edu/userinfo/hpc/login/), which allows us to access the virtual desktop and virtual terminals for our simulations. 

picture of open on demand

Looking at the top of the OpenOnDemand website, there are 4 main tabs you will use: Jobs, Clusters, Interactive Apps, and Interactive Sessions.

* Jobs are used to see what slurm scripts you submitted for a process (usually an analysis script or a simulation)
* Cluster is a local terminal on the head node. This is where you should be submitting jobs.
* Interactive Apps will show you more applications you can open on the super computer. The only important ones are jupyter notebook and desktop.
* My interactive sessions are used to see what interactive applications you have open.

The terminal uses unix commands that I highly recommend you learn. Here are some good sites to learn [1](https://files.fosswire.com/2007/08/fwunixref.pdf) [2](https://learning.rc.virginia.edu/notes/unix-tutorial/).

As you learn more about Rivanna, just remember that you have 3 main file directories: Home, Project and Scratch.

**YOU WANT TO DO SIMULATIONS ON SCRATCH AND SCRATCH ONLY!!!** But keep in mind that files are removed after 90 days, so start putting the files you want to keep in the projects file.

Projects is a shared file between all of us to put our final data and scripts to share.

Home is usually for quick scripts and small files.

Programming Languages
Some of the programming languages you should know to an extent are Python and TCL.

Here are some [TCL tutorials](https://www.tcl-lang.org/man/tcl8.5/tutorial/tcltutorial.html) that you can use to practice.

Here are some [Python tutorials](https://www.w3schools.com/python/) that you can use to practice.

Although it helps to have a programming background, if it is hard to understand that is ok! Learning is the process.

If you don’t understand the code, feel free to just search it up or ask one of us!


See Introduction [here](more_information.md)

