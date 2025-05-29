# Simulation Tutorial

*For the Lampe/Letteri/Taylor Computational Group*
*Created by Justin Kim*

---

## Welcome

Welcome to the Computational Group. Whether you're here because you're interested in biomaterials, coding, or simply curious about molecular dynamics (MD), we’re excited to have you join the team.

Before diving into simulations, please follow the instructions below to install the necessary tools and familiarize yourself with our computing environment.

---

## Required Downloads

**[Visual Molecular Dynamics (VMD)](https://www.ks.uiuc.edu/Research/vmd/vmd-1.9.3/)**
VMD is primarily used to build and visualize peptide systems using TCL scripts.

**[Visual Studio Code (VSCode)](https://code.visualstudio.com/download)**
VSCode is our go-to text editor for `.pdb`, `.psf`, `.txt` files, and scripting in Python, TCL, or Jupyter Notebooks. If you need help with it, feel free to ask Justin.

**[Globus](https://www.globus.org/)**
Globus allows for efficient file transfers between your local computer and Rivanna.

**[Avogadro](https://two.avogadro.cc/index.html)**
This tool is useful for constructing molecular structures and exporting them as `.pdb` files.

**[CHARMM-GUI](https://charmm-gui.org/)**
Create an account on this website. CHARMM-GUI is used to convert `.pdb` files into CHARMM-compatible files and generate `.psf` files for simulations.

---

## Using Rivanna (UVA High-Performance Computing)

We run our simulations on **[Rivanna](https://www.rc.virginia.edu/userinfo/hpc/)**, UVA’s High-Performance Computing cluster. You’ll access it via **[Open OnDemand](https://www.rc.virginia.edu/userinfo/hpc/login/)**, a web-based interface.

### Key Tabs in Open OnDemand

* **Jobs** – View and monitor submitted SLURM job scripts.
* **Clusters** – Terminal access to the head node, where you submit jobs.
* **Interactive Apps** – Open virtual tools such as Jupyter Notebooks and Virtual Desktop.
* **My Interactive Sessions** – View and manage your active sessions.

It is important to be familiar with basic Unix terminal commands. These guides are helpful:

* [Fosswire Unix Cheat Sheet (PDF)](https://files.fosswire.com/2007/08/fwunixref.pdf)
* [UVA RC Unix Tutorial](https://learning.rc.virginia.edu/notes/unix-tutorial/)

---

## File System Overview

Rivanna provides three main directories. Understanding what each is used for will help you manage your data effectively.

* **`/scratch`** – Run all simulations here. Files in this directory are automatically deleted after 90 days, so be sure to back up anything important.
* **`/project`** – Shared among the group. Use this to store final results, key data, and reusable scripts.
* **`/home`** – Your personal directory. Suitable for small scripts or configuration files. Do not run simulations here.

**Important:**
Simulations must be run in `/scratch` to ensure performance and avoid overloading the system.
Files you want to keep should be moved to `/project`.

---

## Programming Languages

You will encounter the following languages during your work in the group:

### TCL

Used primarily for scripting in VMD.
Tutorial: [TCL Basics](https://www.tcl-lang.org/man/tcl8.5/tutorial/tcltutorial.html)

### Python

Used for data analysis, plotting, and automation.
Tutorial: [Python Basics](https://www.w3schools.com/python/)

It’s okay if you’re not familiar with either language. The most important thing is a willingness to learn. Don’t hesitate to ask questions or look things up.

---

## Next Steps

Continue to the [Introduction](more_information.md) for project-specific protocols and the next part of your onboarding.

---

