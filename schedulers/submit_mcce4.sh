#!/bin/bash

# Parameter/Options for SLURM (Simple Linux Utility for Resource Management)
#SBATCH --job-name=submit_mcce4
#SBATCH --output=submit_mcce4.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=12G                 # Adjust memory if needed
#SBATCH --time=24:00:00
#SBATCH --export=ALL

#=============================================================================
#-----------------------------------------------------------------------------
# >>> Automated Parameters (best not change)
APPTNR=$(command -v apptainer) || { echo "apptainer not found"; exit 1; }
MCCE=$(command -v mcce)        || { echo "mcce not found"; exit 1; }
PYEX=$(python3 -c "import sys; print(sys.executable)")
PYENV="${PYEX%python3}"
MCBIN=$(cd "$(dirname "$MCCE")" && pwd)
MCCE_HOME=$(cd "$(dirname "$MCCE")/.." && pwd)
APPTAINER_BIN=$(dirname "$APPTNR")
[[ -x "$MCCE_HOME/bin/mcce" ]] || { echo "[ERROR] mcce not found in $MCCE_HOME/bin"; exit 1; }
# <<<

# Set INPUT & MCCE4 Parameters
input_pdb="prot.pdb"               # PATH to input PDB if you have soft-linked your PDB as "prot.pdb", e.g.:  ln -s 4lzt.pdb prot.pdb  
USER_PARAM="./user_param"          # PATH to "user_param" directory containing additonal topology files (local files). This directory must be called "user_param" (default: MCCE_HOME/param)
EXTRA="./user_param/extra.tpl"     # PATH to an different "extra.tpl" file (local file). (default: MCCE_HOME/extra.tpl)
TMP="/tmp"                         # PATH to temporary directory for storing PBE calculation files during step3
CPUS=1                             # Number of CPU cores to use for parallelizable MCCE calculations
EPS=4                              # Protein dielectric constant

# Step control flags
step1="t"               # STEP1: pre-run, pdb-> mcce pdb  (DO_PREMCCE)
step2="t"               # STEP2: make rotamers            (DO_ROTAMERS)
step3="t"               # STEP3: Energy calculations      (DO_ENERGY)
step4="t"               # STEP4: Monte Carlo Sampling     (DO_MONTE)
step_clean="t"          # Clean PBE data from TMP         (BACKUP CLEAN) : Set to f if step3 --debug option is used

# Optional step controls
center="t"              # Center protein structure before MCCE run      : Set to f to skip centering and use input PDB as-is
stepM="f"               # Generate Partial Membranes                    : If true, user MUST satisfy condidtions of stepM.sh, which can be be obtained on MCCE4/inhouse/stepM.sh
stepA="f"               # Run a custom script between step1 and step2   : If true, user MUST satisfy condidtions of their custom script
stepB="f"               # Run a custom script between step2 and step3   : If true, user MUST satisfy condidtions of their custom script
stepC="f"               # Run a custom script between step3 and step4   : If true, user MUST satisfy condidtions of their custom script

# MCCE Simulation
STEP1="$PYEX $MCBIN/step1.py -d $EPS --dry"
STEP2="$PYEX $MCBIN/step2.py -d $EPS -l 1"
STEP3="$PYEX $MCBIN/step3.py -d $EPS -s ngpb -p $CPUS -t $TMP"
STEP4="$PYEX $MCBIN/step4.py --xts -i 7 -n 1"

# Optional MCCE script locations
STEPM="/path/to/stepM.sh"         # Optional StepM: Bash script
STEPA="/path/to/stepA_script.py"  # Optional StepA: Python script to run between step1 and step2.
STEPB="/path/to/stepB_script.py"  # Optional StepB: Python script to run between step2 and step3.
STEPC="/path/to/stepC_script.py"  # Optional StepC: Python script to run between step3 and step4.

# NOTE: User is responsible to precheck if custom scripts work properly and efficiently
# NO USER INPUT NECCESARY BELOW THIS LINE
#------------------------------------------------------------------------------
#==============================================================================

# Initialize Apptainer to ensure job uses user-installed Apptainer and avoid systemd/cgroups (DBus) issues on compute nodes
export PATH="$APPTAINER_BIN:$PATH"
export APPTAINER_CONFIG_FILE="$HOME/.apptainer/apptainer.conf"
mkdir -p "$HOME/.apptainer"
cat > "$APPTAINER_CONFIG_FILE" <<'EOF'
systemd cgroups = no
EOF

# Remove any existing instance of mc_bin from PATH and prepend mc_bin to PATH
PATH=$(echo "$PATH" | tr ':' '\n' | grep -vx "$MCCE_HOME/MCCE_bin" | paste -sd ':' -)
export PATH="$MCCE_HOME/MCCE_bin:$PATH"

# Remove any existing instance of mc_bin from PATH and prepend mc_bin to PATH
PATH=$(echo "$PATH" | tr ':' '\n' | grep -vx "$MCBIN" | paste -sd ':' -)
export PATH="$MCBIN:$PATH"

# Remove any existing instance of PYENV from PATH and prepend PYENV to PATH
PATH=$(echo "$PATH" | tr ':' '\n' | grep -vx "$PYENV" | paste -sd ':' -)
export PATH="$PYENV:$PATH"

echo "============================================================"
echo "MCCE4 SUBMIT SHELL JOB ENVIRONMENT (startup diagnostics)"
echo "------------------------------------------------------------"
echo "Date:             $(date)"
echo "Host:             $(hostname)"
echo "User:             $(whoami)"
echo
echo -e "Apptainer:        $(which apptainer)"
echo -e "Config File:      $APPTAINER_CONFIG_FILE"
echo -e "MCCE_HOME:        $MCCE_HOME"
echo -e "MCBIN:            $MCBIN"
echo -e "Driver:           $MCBIN/driver_mcce4.sh"
echo -e "PYEX, PYENV:      $PYEX"
echo -e "PATH:             $PATH"
echo "============================================================"
echo

# Export environment for downstream script
export PYEX
export input_pdb MCCE_HOME MCBIN USER_PARAM EXTRA TMP
export step1 step2 step3 step4 step_clean
export center stepM stepA stepB stepC
export STEP1 STEP2 STEP3 STEP4
export STEPM STEPA STEPB STEPC

# Inititiate MCCE_HOME PATH and call driver_mcce4.sh
"$MCBIN"/driver_mcce4.sh

# ==============================================================================
# Script Name   : submit_mcce4.sh
# Purpose       : Automate and control the full MCCE4 simulation pipeline including optional custom preprocessing steps.
#
# Description   :
#   This script manages the sequential execution of MCCE4 simulation steps (1 to 4), with optional hooks (stepM, stepA, stepB, stepC)
#   that allow the user to insert custom membrane generation and intermediate processing scripts.
#   It records the timing and success/failure of each step in a detailed log file (`mcce_timing.log`).
#   The script supports flexible control through flags to enable/disable specific MCCE steps or custom steps.
#
# Main Features :
#   - Step 1: Convert standard PDB to MCCE-compatible format
#   - Step 2: Generate rotamers
#   - Step 3: Perform energy calculations
#   - Step 4: Run Monte Carlo sampling
#   - Step clean: Clean temporary pbe_data
#   - Optional StepM: Add membrane-specific conformers (e.g., using IPECE)
#   - Optional StepA/B/C: Insert custom preprocessing scripts between core steps
#   - Intelligent skip logic and output checking to prevent redundant work
#   - Runtime logging with timestamps for each phase
#
# Input Requirements :
#   - input_pdb       : A protein PDB file named `prot.pdb`
#   - MCCE_HOME       : Path to MCCE4 installation directory
#   - USER_PARAM      : Directory with user-defined MCCE parameters (optional)
#   - EXTRA           : Custom extra.tpl file (optional)
#   - TMP             : Path to store temporary pbe_data (default: /tmp)
#   - CPUS            : Multiprocessing for step3
#   - Optional scripts for stepM/A/B/C must exist and be executable if enabled
#
# Output Files :
#   - step1_out.pdb, step2_out.pdb, head3.lst, pK.out
#   - Timing report: mcce_timing.log
#   - Logs for each step: step1.log, step2.log, etc.
#
# Usage :
#   Set control flags (`step1`, `step2`, etc.) to "t" or "f" to enable/disable each step.
#   Set paths to optional scripts as needed.
#   Submit this script to a SLURM cluster or run locally if sbatch is not used.
#
# Author        : Gehan A. Ranepura
# Date Created  : July 2025
# ==============================================================================
