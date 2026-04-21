# Multi-Conformation Continuum Electrostatics

<p align="center">
  <img src="docs/images/mcce_logo1.png" alt="MCCE Logo" style="max-width: 100%; height: auto;">
</p>

# Welcome to MCCE4-Alpha!  
Please see our CHANGELOG at the bottom for the latest updates!

# [__📖 MCCE4-Alpha Tutorial__](https://gunnerlab.github.io/mcce4_tutorial/)

Comprehensive documentation covering:
- Installation
- Guide: Detailed explanations of all settings
- Example Projects

## __Installation__
#### If you have sudo access or would like a system-wide installation of the needed softwares:
The file `MCCE_bin/sudo_install.txt` has the necessary information for you or your sys admin to install the packages. To display the file, run this command:
```
 cat ./MCCE_bin/sudo_install.txt
```

## __Quick Installation__
#### "Quick Install" script `MCCE_bin/quick_install.sh`:
__Note: The quick install script will modify an existing conda environment named 'mc4'.__ 
If you need to re-create it, to troubleshoot an installation issue, for example, run this command before running the script:
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

## OPTIONAL: MCCE4-Alpha CLI (`mc4`) 

The `mc4` command-line tool runs MCCE4-Alpha inside an [Apptainer](https://apptainer.org/) container, giving you a reproducible environment with all solvers (DELPHI, APBS, NGPB) pre-configured. No `sudo` required.

#### Linux

```bash
# 1. Clone the repository
git clone https://github.com/GunnerLab/MCCE4-Alpha.git
cd MCCE4-Alpha

# 2. Run the setup script
bash setup.sh

# 3. Source your shell config (or open a new terminal)
source ~/.bashrc

# 4. You're ready
mc4 step1.py prot.pdb
```

#### macOS

```bash
# 1. Clone the repository
git clone https://github.com/GunnerLab/MCCE4-Alpha.git
cd MCCE4-Alpha

# 2. Run the macOS setup script
bash setup_mac.sh

# 3. Source your shell config (or open a new terminal)
source ~/.zshrc    # or ~/.bashrc

# 4. You're ready
mc4 step1.py prot.pdb
```

---

### What the Setup Scripts Do

#### `setup.sh` (Linux)

1. **Installs Apptainer** — checks for a system installation first; if unavailable, installs it via Conda into user-space (no `sudo` needed).
2. **Downloads the NGPB solver** — fetches the pre-built `NextGenPB.sif` (~1.6 GB) from GitHub Releases.
3. **Builds the container** — creates `mcce4-alpha.sif` using the definition file at `bin/mcce4-alpha.def`. This image contains the `mc4` Conda environment, compiled MCCE4 binaries, and all three PB solvers.
4. **Updates your PATH** — adds `bin/` to your `~/.bashrc` so `mc4` is available in new terminals.

| Flag | Effect |
|------|--------|
| *(none)* | Smart re-run — skips steps that are already complete |
| `--rebuild` | Force a full container rebuild |
| `--build-ngpb` | Build NGPB from source (`bin/recipe_MCCE.def`) instead of downloading the pre-built image |

#### `setup_mac.sh` (macOS)

macOS cannot run Apptainer natively, so the script sets up a lightweight Linux VM via [Lima](https://lima-vm.io/):

1. **Installs Lima & QEMU** via Homebrew.
2. **Creates an Ubuntu 24.04 VM** (`mcce4`) with Apptainer pre-installed and your home directory mounted read-write.
3. **Downloads the NGPB solver** inside the VM.
4. **Builds the container** inside the VM.
5. **Updates your PATH** — adds `bin/` to your `~/.zshrc` (or `~/.bashrc`).

Your local code edits on macOS are instantly visible inside the VM — no rebuild needed.

| Flag | Effect |
|------|--------|
| *(none)* | Smart re-run — skips completed steps |
| `--rebuild` | Delete the VM and start fresh |
| `--build-ngpb` | Build NGPB from source instead of downloading |

---

### Usage

```bash
mc4 <command> [args...]
```

**Examples:**

```bash
mc4 step1.py <pdbfile>         # Run step 1
mc4 step2.py                   # Run step 2
mc4 step3.py                   # Run step 3 (PB solver)
mc4 step4.py                   # Run step 4
mc4 getpdb 4pti                # Fetch a PDB file
mc4 which python               # Check which Python the container uses
mc4 --shell                    # Open an interactive shell inside the container
```

**Expected output** (production mode):

```
/opt/conda/envs/mc4/bin/python
```

---

### Development Mode

Use the `-d` flag to bind-mount your local MCCE4-Alpha source code into the container. Changes to files on your host are reflected immediately — no rebuild required.

```bash
mc4 -d <command> [args...]
```

**Example:**

```bash
mc4 -d step1.py prot.pdb
```

```
🔧 Running in Apptainer Development Mode...
INFO:    fuse2fs not found, will not be able to mount EXT3 filesystems
Preprocessing input pdb file, identifying ligands ...
```

> **Note:** Dev mode is available on Linux. On macOS, code changes are already live via Lima's mount, so `-d` is not needed.

---

### How `mc4` Auto-Detects Its Mode

| Platform | Condition | Mode |
|----------|-----------|------|
| Linux | `mcce4-alpha.sif` exists + `apptainer` on PATH | **Container** — `apptainer exec` |
| macOS | Lima VM `mcce4` exists | **Lima** — routes through `limactl shell` |
| Either | No container or VM found | **Native fallback** — runs commands directly via the `mc4` Conda env |

> In native fallback mode, DELPHI and APBS work but NGPB requires the container.

---

### Troubleshooting

#### `apptainer: command not found`


**Fix:** Activate the environment before using the CLI:

```bash
conda activate mc4
mc4 step1.py prot.pdb
```

#### Container build killed (OOM on HPC login nodes)

If `setup.sh` dies at the `mksquashfs` step, the OOM killer likely terminated it. This is common on shared HPC login nodes with memory limits.

**Options:**

```bash
# Option A: Build on a compute node with more RAM
srun --mem=8G --time=30:00 bash setup.sh --rebuild

# Option B: Point temp files to a larger filesystem
export APPTAINER_TMPDIR=/scratch/$USER/tmp
mkdir -p $APPTAINER_TMPDIR
bash setup.sh --rebuild
```

#### NGPB download failed

If the ~1.6 GB NextGenPB download times out, you can retry or build from source:

```bash
# Retry download
bash setup.sh --rebuild

# Or build from source (~20-30 min, requires bin/recipe_MCCE.def)
bash setup.sh --build-ngpb
```

DELPHI and APBS will still work even if NGPB is unavailable.

#### Lima VM issues (macOS)

```bash
# Check VM status
limactl list

# Restart the VM
limactl stop mcce4 && limactl start mcce4

# Full reset
bash setup_mac.sh --rebuild
```

## Environment update (01-08-2026):
If your 'mc4' environment predates 01-08-2026, run these commands to update it:
  ```bash
  CLONE=$(dirname $(dirname "$(python3 -c "import os, sys; print(os.path.realpath(sys.argv[1]))" "$(which ms_protonation)")"));
  conda env update -n mc4 -f "$CLONE/mc4.yml
  ```

## CHANGELOG:
<!--- NOTE TO EDITOR: Use tis line to indicate that the user must/should update their clone"
  - __Apply changes: cd to your clone, then run `git pull`__
-->
_This section will reflect important changes and will provide you with information on how to apply them; For example, if new python packages are added to the environment file (mc4.yml), then the entry pertaining to that change will list the command(s) to update your environment._ 

* 2026-04-15:
  - Feature Merge: MCCE4 Topology Agent
  - __Apply changes: cd to your clone, then run `git pull`__
  
* 2026-04-09:
  - Feature Merge: MCCE4 GUI
  - Feature Merge: Protein Dipole Vectors and Visualization
  - __Apply changes: cd to your clone, then run `git pull`__
 
* 2026-02-20:
  - Feature Merge: Integrated Apptainer/Singularity containerization for the MCCE4-Alpha CLI tool.
  - Automated environment setup and image building via setup.sh to ensure cross-platform portability.
  - __Apply changes: cd to your clone, then run `git pull`__
  - __Apply changes: run `./setup.sh` to build the new container image__

* 2026-01-30:
  - Updated submit_mcce4.sh and driver_mcce4.sh to pass environment
  - __Apply changes: cd to your clone, then run `git pull`__

* 2026-01-26:
  - Added `numba` in env file:
  - __Apply changes: cd to your clone, then run `git pull`__
  - __Apply changes: run these commands:__
  ```
  CLONE=$(dirname $(dirname "$(python3 -c "import os, sys; print(os.path.realpath(sys.argv[1]))" "$(which mcce)")"));
  conda env update -n mc4 -f $CLONE/mc4.yml
  ```

* 2026-01-20:
  - Comprehensive update of the tutorial site
  - Minimized README file


* 2026-01-08:
  - Updated python dependencies in mc4.yml
  - __Apply changes: run these commands:__
  ``` 
  CLONE=$(dirname $(dirname "$(python3 -c "import os, sys; print(os.path.realpath(sys.argv[1]))" "$(which mcce)")"));
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
