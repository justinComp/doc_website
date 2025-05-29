# VMD/NAMD System Preparation Tutorial

The necessary files will be downloadable [here](https://drive.google.com/drive/folders/1tE-We6xeah_VM3ksenBdn4Opb3q_biOV?usp=sharing).

Download the most **recent version** of the script called `makefiles.tcl` (check the date to confirm it’s the latest version). You will also need the `.rtf` and `.str` files in your working directory.

---

### Required Files in Your Working Directory

Make sure the following files are in the **same folder** (your working directory):

* `makefiles.tcl` (download from Box)
* Your `.rtf` file
* Your `.str` file
* Your `.pdb` file
* Your `.psf` file

---

### Editing the Script

Open `makefiles.tcl` in a text editor (e.g., VS Code).

Use the makegrid file is used to build a peptide grid in VMD using our L and D peptide systems.

1. **`speciesL`**

   * *Type:* String
   * *What it is:* The name of your L-peptide structure (e.g., `Lkyfil`).

2. **`speciesD`**

   * *Type:* String
   * *What it is:* The name of your D-peptide structure (e.g., `Dkyfil`).

3. **`outputname`**

   * *Type:* String
   * *What it is:* The base name for the final output files (e.g., `ion1to1kyfil`).

4. **`ratioL`**

   * *Type:* Fraction
   * *What it is:* The fraction of the total peptide count that will be L-peptides (e.g., `32/64` for 1:1).

5. **`ratioD`**

   * *Type:* Fraction
   * *What it is:* The fraction of D-peptides (e.g., `32/64` for 1:1).

6. **`totalpep`**

   * *Type:* Integer
   * *What it is:* The total number of peptides in your simulation (e.g., `64`).

7. **`gridX`**

   * *Type:* Integer
   * *What it is:* Number of peptides along the x-axis in the 3D grid.

8. **`gridY`**

   * *Type:* Integer
   * *What it is:* Number of peptides along the y-axis.

9. **`gridZ`**

   * *Type:* Integer
   * *What it is:* Number of peptides along the z-axis.
   * *Tip:* `gridX * gridY * gridZ` should equal `totalpep`.

10. **`unit_buffer`**

    * *Type:* Float (in Å)
    * *What it is:* The spacing between peptides in the grid.
    * *Recommended:* 32 Å depending on molarity.

11. **`solvatePad`**

    * *Type:* Float (in Å)
    * *What it is:* The padding around the entire peptide box during solvation.
    * *Recommended default:* `5` Å.

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

4. Run the script:

```tcl
source makefiles.tcl
```
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


