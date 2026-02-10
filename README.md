# Multi-Conformation Continuum Electrostatics

<p align="center">
  <img src="docs/images/mcce_logo1.png" alt="MCCE Logo" style="max-width: 100%; height: auto;">
</p>

# Welcome to __MCCE4-Alpha__!  
Please see our CHANGELOG at the bottom for the latest updates!

## __Quick Introduction__

Given the structure of a macromolucule (in a PDB file), __MCCE4__ can predict the following:

- __pKₐ values__
- __Protonation states__
- __Electrostatic properties__ of biomolecules

In this program, protein side chain motions are simulated explicitly while the dielectric effect of solvent and bulk protein material is modeled by continuum electrostatics.

## __Installation__
#### If you have sudo access or would like a system-wide installation of the needed softwares:
The file `MCCE_bin/sudo_install.txt` has the necessary information for you or your sys admin to install the packages. To display the file, run this command:
```
 cat ./MCCE_bin/sudo_install.txt
```

## __Quick Installation__
#### "Quick Install" script `MCCE_bin/quick_install.sh`:
__Note: The quick install script will not modify an existing conda environment named 'mc4'.__ 
If you want to re-create it, run this command before running the script:
```
 conda env remove -n mc4
```

  1. Clone this repo, then cd into it with this command:
  ```
   git clone https://github.com/GunnerLab/MCCE4-Alpha.git; cd MCCE4-Alpha;
  ```
  
  2. Run the `quick_install.sh` script to download MCCE PBE solver (NGPB) image file and create a conda environment for MCCE4 (this may need several passes if you need to install dependencies such as miniconda and apptainer):
  ```
   bash ./MCCE_bin/quick_install.sh
  ``` 

### What this script does:
  - Checks for required `conda`; Stops if not found so you can install it (commands provided).
  - Create a conda environment for MCCE4 named 'mc4' (using 'mc4.yml').
  - Checks for required `apptainer`; If a system Apptainer installation is not found & an 'unprivilege' version cannot be installed, Apptainer is conda-installed in 'mc4'.
  - Downloads the generic image for NGPB in MCCE4-Alpha/bin.
  - Adds export commands to the PATH variable in ~/.bashrc for:
    * 'MCCE4-Alpha/bin' and 'MCCE4-Alpha/MCCE_bin'
    * the unprivilege version of Apptainer if installed by the script

## __MCCE4-Alpha CLI__ (Optional)
To build and configure the `mc4` CLI tool, follow these steps:

1.  **Install Miniconda** (Skip if Conda is already installed):
    Navigate to your home directory and install Miniconda.
```bash
cd ~/ ; wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && bash ~/Miniconda3-latest-Linux-x86_64.sh -b
```
2. Source path to Miniconda:
```bash
 source ~/miniconda3/bin/activate
 ```
3. Navigate to the root directroy of `MCCE4-Alpha` and run the `setup.sh` script. This will install unprivileged Apptainer and build the `mcce4-alpha.sif` file located in `bin/`.
```bash
./setup.sh
```
>[!Important] You will be prompted to update your .bashrc file. Accept so the mc4 executable can be availble to use for new terminal session.

4.  Source the `.bashrc` file.
```bash
source ~/.bashrc
```
5.  Activate the `mc4` environment that the setup script created.
```bash
conda activate mc4
```
**Usage:**
```bash
mc4 <command>
```
**Example**
```bash
mc4 which python
```

**Expected Output**
```bash
🚀 Running in Apptainer Production Mode...
INFO:    fuse2fs not found, will not be able to mount EXT3 filesystems
/opt/conda/envs/mc4/bin/python
```
> [!NOTE] 
> **Troubleshooting: `apptainer: command not found`**
> The `mc4` CLI wrapper script functions by executing the command: `apptainer exec <path/to/mcce4-alpha.sif> <command> ...`
> * **The Cause:** If Apptainer was installed via the **Conda fallback method**, the `apptainer` binary is hidden inside the `mc4` environment. The wrapper script fails because it cannot find `apptainer` in your system PATH.
> * **The Fix:** You must run `conda activate mc4` to expose the binary before using the CLI.

### Development Mode
To run the CLI in development mode, use the `-d` flag. This mode mounts your local `MCCE4-Alpha` source code into the container, allowing you to test changes immediately without rebuilding the image. Since the repository is bound to the container, ***developers can create new files or modify existing ones***, and these changes will be reflected inside the container.

**Usage:**
```bash
mc4 -d <command>
```

**Example**
```bash
mc4 -d which python
```
**Expected Output:**
```bash
🔧 Running in Apptainer Development Mode...
INFO:    fuse2fs not found, will not be able to mount EXT3 filesystems
/opt/conda/envs/mc4/bin/python
```

## Environment update (01-08-2026):
If your 'mc4' environment predates 01-08-2026, run these commands to update it:
  ```bash
  CLONE=$(dirname $(dirname "$(python3 -c "import os, sys; print(os.path.realpath(sys.argv[1]))" "$(which ms_protonation)")"));
  conda env update -n mc4 -f "$CLONE/mc4.yml
  ```

__🚀 Run Your First Job:__ [Quick Start](https://gunnerlab.github.io/mcce4_tutorial/docs/guide/quick_start/)  

__📖 MCCE4-Alpha Tutorial:__ [Full Documentation](https://gunnerlab.github.io/mcce4_tutorial/)

Comprehensive documentation covering:
- Installation
- Guide: Detailed explanations of all settings
- Example Projects 

## MCCE4-Tools 🔧  
Please also check out the companion repository __MCCE4-Tools__. 

🧰 __Explore Now:__ [MCCE4-Tools GitHub](https://github.com/GunnerLab/MCCE4-Tools)
---

# CHANGELOG:
<!--- NOTE TO EDITOR: Use tis line to indicate that the uses rmust/should update their clone"
  - __Apply changes: cd to your clone, then run `git pull`__
-->
_This section will reflect important changes and will provide you with information on how to apply them; For example, if new python packages are added to the environment file (mc4.yml), then the entry pertaining to that change will list the command(s) to update your environment._ 

* 2026-01-08:
  - Updated python dependencies in mc4.yml
  - __Apply changes: run these commands:__
  ```
  CLONE=$(dirname $(dirname "$(readlink -f "$(which mcce)")")); echo "$CLONE"
  conda env update -n mc4 -f "$CLONE/mc4.yml
  ```

* 2025-11-25:
  - step1.py: Added error trapping on atom.loadline call
  - mfe.py: Updated & moved to MCCE_bin
  - __Apply changes: cd to your clone, then run `git pull`__

* 2025-11-11:
  - Fixed deleterious typo in bin/pdbs_interfaces.py
  - __Apply changes: cd to your clone, then run `git pull`__

* 2025-10-30:
  - Updated README: Added CHANGELOG, link to sudo_install.txt
  - Added topologies for SO4 and PO4 in param/.
  - Updated bin/step3.py with longer timeout value
  - Updated MCCE_bin/quick_install.sh
  - __Apply changes: cd to your clone, then run `git pull`__

---

## Help us improve MCCE4
This is a testing version of MCCE4 development. 
Please let us know about questions, comments or report any issues you encounter [here](https://github.com/GunnerLab/MCCE4-Alpha/issues).
Thank You and we hope you enjoy using MCCE4!  

## MCCE Wiki
[Learn about MCCE, installation, available tools, and research done with MCCE.](https://mccewiki.levich.net) (under construction)

---

Copyright (C) 2024 GunnerLab
This software is distributed under the terms the terms of the MIT licence
