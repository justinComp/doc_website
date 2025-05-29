# VMD/NAMD System Preparation Tutorial

The necessary VMD files will be downloadable [here](https://drive.google.com/drive/folders/1tE-We6xeah_VM3ksenBdn4Opb3q_biOV?usp=sharing).

Download the most **recent version** of the script called `makefiles.tcl` (check the date to confirm it’s the latest version). You will also need the `.rtf` and `.str` files in your working directory.

---

### Required Files in Your Working Directory

Make sure the following files are in the **same folder** (your working directory):

* `makefiles.tcl` (download from Drive)
* Your `.rtf` file
* Your `.str` file
* Your `.pdb` file
* Your `.psf` file

---

### Editing the Script

Open `makefiles.tcl` in a text editor (e.g., VS Code).

Use the makegrid file is used to build a peptide grid in VMD using our L and D peptide systems.

1. **`speciesL`**

   - **Type:** String  
   - **What it is:** The name of your L-peptide structure (e.g., `Lkyfil`).

2. **`speciesD`**

   - **Type:** String  
   - **What it is:** The name of your D-peptide structure (e.g., `Dkyfil`).

3. **`outputname`**

   - **Type:** String  
   - **What it is:** The base name for the final output files (e.g., `ion1to1kyfil`).

4. **`ratioL`**

   - **Type:** Fraction  
   - **What it is:** The fraction of the total peptide count that will be L-peptides (e.g., `32/64` for 1:1).

5. **`ratioD`**

   - **Type:** Fraction  
   - **What it is:** The fraction of D-peptides (e.g., `32/64` for 1:1).

6. **`totalpep`**

   - **Type:** Integer  
   - **What it is:** The total number of peptides in your simulation (e.g., `64`).

7. **`gridX`**

   - **Type:** Integer  
   - **What it is:** Number of peptides along the x-axis in the 3D grid.

8. **`gridY`**

   - **Type:** Integer  
   - **What it is:** Number of peptides along the y-axis.

9. **`gridZ`**

   - **Type:** Integer  
   - **What it is:** Number of peptides along the z-axis.  
   - **Tip:** `gridX * gridY * gridZ` should equal `totalpep`.

10. **`unit_buffer`**

    - **Type:** Float (in Å)  
    - **What it is:** The spacing between peptides in the grid.  
    - **Recommended:** 32 Å depending on molarity.

11. **`solvatePad`**

    - **Type:** Float (in Å)  
    - **What it is:** The padding around the entire peptide box during solvation.  
    - **Recommended default:** 5 Å.

If something does not work, use **Ctrl + F** and search for the variables above and hard code it into the script (not sure why some grids fail and others work).

---


### Running the Script in VMD

1. **Open VMD**

2. Navigate to the Tk Console:
   `Extensions > Tk Console`

3. Use the following commands to move into your working directory:

```tcl
cd "folder_name"
```

Helpful navigation tips:

* `cd ../` — move back one folder level
* `pwd` — print your current working directory

4. Load the script:

```tcl
source makefiles.tcl
```

5. Run the script:

```tcl
makegrid LnewIons Lkyfil Lkyfil 64/64 0/64 64 4 4 4 32 5
```

Script may take a while as it loads the peptides. If it seems to fail, **restart VMD** and reload the script after hard coding values.
---


### What This Script Does

This script will generate several intermediate `.pdb` and `.psf` files. The final outputs you care about are:

* `ionized.pdb`
* `ionized.psf`

Rename them to keep track of your peptide system. A good naming convention is:

```
ion<ratio><peptidename>.pdb
```

For example:
`ion1to1kyfil.pdb` → Ionized file of KYFIL peptides with a 1:1 ratio of L and D stereochemistry.

---

### Additional Output

You will also see a text file:

```
listofLD$molname.txt
```

This file contains:

* Center of mass
* Water sphere radius
* Cell basis vectors
* Center of the system
* List of peptide names in the order they were created

This is important parameters that you need to copy and paste into the namd parameter files.

---

The necessary NAMD files will be downloadable [here](https://drive.google.com/drive/folders/145Km5R9Yqdb9bxbG-fP4cMt73OEK4t5i?usp=sharing).

### Required NAMD Files in Your Working Directory

Make sure the following files are in the **same folder** (your working directory):

* Your `.namd` file
* Your `.prm` file
* Your `.str` file
* Your `.pdb` file
* Your `.psf` file

---

### Editing .namd File

Find these variables and change it to your specific files.

```namd
structure          ionized.psf
coordinates        ionized.pdb

...

outputName         OutputName ;

...

# Periodic Boundary conditions, BE SURE TO CHANGE THIS FOR EACH PDB WILL BE DIFFERENT EVERY TIME YOU MAKE A NEW SYSTEM
cellBasisVector1 128.29100036621094 0 0
cellBasisVector2 0 128.47200202941895 0
cellBasisVector3 0 0 128.02500247955322
cellOrigin 47.1793327331543 48.19822311401367 48.779972076416016

...

run 				12500000    ;# 50 ns total trajectory therefore 2E+8 fs, but with 2 fs per step, 2E+8/2 = 50000000/2    steps

```

### Editing .slurm File

Then submit this **SLURM job script** to run the full simulation:

```bash
#!/bin/bash
#SBATCH --nodes=7
#SBATCH --ntasks-per-node=36
#SBATCH --mail-user=zqp6mj@virginia.edu
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --time=1-00:00:00
#SBATCH --partition=parallel
#SBATCH -A  [ALLOCATIONS]
#SBATCH -o [OUTFILENAMEOFNAMDFILE].out

module load goolf namd
mpiexec namd2 [NAMEOFNAMDFILE].namd
```