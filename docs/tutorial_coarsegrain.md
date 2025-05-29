Here is a cleaned and professional version of your Coarse Grain with GROMACS tutorial, formatted clearly for reproducibility:

---

# Coarse-Grained Simulation with GROMACS: KYFIL Peptide

## Getting Started

Download the necessary files from this [Google Drive folder](https://drive.google.com/drive/folders/13CHhxHLBMIX4VO-UOdFUi3loL8YrpLb8?usp=sharing).

### Load Required Modules on Rivanna

```bash
module load gcc/11.4.0 openmpi/4.1.4 gromacs/2023.2
module spider gromacs/2023.2
pip install --user vermouth
export PATH=$PATH:/home/YOUR_COMPUTING_ID/.local/bin
```

---

## Step 1: Create Coarse-Grained Peptide

Convert the atomistic peptide PDB file to a coarse-grained model using `martinize2`.

```bash
martinize2 -f KYFIL.pdb -x KYFIL_CG.pdb -o single_KYFIL.top -ff martini3001 -ss C -cter NH2-ter
```

**Flags:**

* `-ff`: Force field (`martini3001`)
* `-ss`: Secondary structure (use `C` for coil)
* `-cter`: Terminus (NH2 amidation at C-terminus)

---

## Step 2: Convert PDB to GRO

```bash
gmx_mpi editconf -f KYFIL_CG.pdb -o KYFIL_CG.gro -d 0.5
```

Adds a 0.5 nm buffer around the peptide.

---

## Step 3: Create Peptide Grid

```bash
gmx_mpi genconf -f KYFIL_CG.gro -nbox 4 4 4 -rot yes -dist 1 2 1 -o KYFIL_64_box_cg.gro
```

**Notes:**

* 64 peptides total (4 × 4 × 4)
* Randomly rotates each peptide
* Spacing defined by `-dist` in nm (x y z)

---

## Step 4: Solvate the Box

```bash
gmx_mpi solvate -cp KYFIL_64_box_cg.gro -cs water.gro -o Grid_KYFIL_CG_water.gro -box 13 13 13
```

Adds coarse-grained water (1 bead = 4 water molecules).

Example output:

* \~18,674 water molecules
* Volume ≈ 2197 nm³
* Density ≈ 2614 g/L

---

## Step 5: Prepare Topology File

1. Copy `single_KYFIL.top` to `system.top`
2. Open `system.top` and edit it to include:

```text
#include "martini_v3.0.0.itp"
#include "martini_v3.0.0_solvents_v1.itp"
#include "martini_v3.0.0_ions_v1.itp"
#include "molecule_0.itp"

[ system ]
Title of the system

[ molecules ]
molecule_0    64
W             [# water molecules]
```

Replace `[ # water molecules ]` with the number from:

```bash
grep -c W Grid_KYFIL_CG_water.gro
```

---

## Step 6: Add Ions

Generate the input file for ion addition:

```bash
gmx_mpi grompp -f ions.mdp -c Grid_KYFIL_CG_water.gro -p system.top -o ions.tpr
```

Add NaCl to 0.15 M concentration:

```bash
gmx_mpi genion -s ions.tpr -pname NA -nname CL -neutral -conc 0.15 -o Grid_KYFIL_CG_ion.gro
```

Choose group `13` (water) when prompted.

Get updated atom counts:

```bash
grep -c W Grid_KYFIL_CG_ion.gro
grep -c CL Grid_KYFIL_CG_ion.gro
grep -c NA Grid_KYFIL_CG_ion.gro
```

Update `system.top` again:

```text
[ molecules ]
molecule_0    64
W             18150
CL            326
NA            198
```

---

## Step 7: Energy Minimization

Prepare:

```bash
gmx_mpi grompp -f minimization.mdp -c Grid_KYFIL_CG_ion.gro -p system.top -o minimization.tpr -maxwarn 2
```

Submit using a SLURM script:

```bash
#!/bin/bash
#SBATCH --nodes=7
#SBATCH --ntasks-per-node=36
#SBATCH --mail-user=YOUR_COMPUTING_ID@virginia.edu
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --time=3-00:00:00
#SBATCH --partition=parallel
#SBATCH -A YOUR_ALLOCATION
#SBATCH -o minimize.out

module purge
module load gcc/11.4.0 openmpi/4.1.4 gromacs/2023.2
gmx mdrun -v -deffnm minimization
```

---

## Step 8: Run Dynamics

Preprocess for dynamics:

```bash
gmx_mpi grompp -f dynamic.mdp -c minimization.gro -p system.top -o dynamic.tpr
```

Submit dynamics run with:

```bash
#!/bin/bash
#SBATCH --nodes=7
#SBATCH --ntasks-per-node=36
#SBATCH --mail-user=YOUR_COMPUTING_ID@virginia.edu
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --time=3-00:00:00
#SBATCH --partition=parallel
#SBATCH -A YOUR_ALLOCATION
#SBATCH -o dynamic.out

module purge
module load gcc/11.4.0 openmpi/4.1.4 gromacs/2023.2
gmx mdrun -v -deffnm dynamic
```

---

Your simulation is now set up and running.

Let me know if you'd like help automating this with a script or integrating multiple peptide types.
