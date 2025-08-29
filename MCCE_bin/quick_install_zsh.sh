#!/bin/zsh
echo "MCCE4 "Quick Install" Script for macOS"
echo ""

set -e

# Use the -h switch to get help on usage.
help_msg() {
  echo ""
  echo "Usage:"
  echo ""
  echo "sh ./MCCE_bin/$(basename "$0") -h"
  echo "sh ./MCCE_bin/$(basename "$0") [env_name]"
  echo ""
  echo " -h                 : Display this help message & exit."
  echo " env_name (optional): Custom name for the conda environment to create (default is 'mc4')."
  echo "                    : Required if you already have a conda environment named 'mc4'."
  echo ""
  exit 0
}

# Process options
while getopts "h" opt; do
  case $opt in
    h) help_msg ;;
    # Exit on invalid options:
    \?) printf "Invalid option: -%s\n" "$OPTARG" >&2
      exit 1
  esac
done
# Shift arguments to access possible env_name after options processing
shift $(( OPTIND-1 ))
# Get arguments from command line
env_name=""
if [ -n "$1" ];
then
  env_name="$1"
  echo "User env: $env_name"
fi

echo "Checking for a download tool (curl or wget)..."
DOWNLOAD_CMD=""
if command -v wget > /dev/null 2>&1;
then
    DOWNLOAD_CMD="wget"
elif command -v curl > /dev/null 2>&1;
then
    DOWNLOAD_CMD="curl"
else
    echo "Error: Neither 'wget' nor 'curl' could be found."
    echo "Please install one of them to download the NGPB container image."
    exit 1
fi
echo "Using '$DOWNLOAD_CMD' for downloads."
echo ""

# Get the script directory path
MCCE_bin=$(cd -- $(dirname -- $0) ; pwd -P);
REPO_PATH=$(dirname "$MCCE_bin")
REPO_bin="$REPO_PATH/bin"
echo "MCCE4/bin folder: $REPO_bin"
echo "MCCE_bin folder: $MCCE_bin"
echo ""
echo "Getting NextGenPB generic container image with $DOWNLOAD_CMD..."
# Need to check if file exists first, otherwise wget will save it as: 'path/NextGenPB.sif.1'
#
# Name of the compiled image when Makefile is used:
sif_mc4="$REPO_bin/NextGenPB_MCCE4.sif"
# Name of the downloaded generic image & source url:
sif_file="$REPO_bin/NextGenPB.sif"
sif_url="https://github.com/concept-lab/NextGenPB/releases/download/NextGenPB_v1.0.0/NextGenPB.sif"

if [ ! -f "$sif_file" ];
then
  echo "Downloading NGPB generic image from $sif_url..."
  if [ "$DOWNLOAD_CMD" == "curl" ];
  then
      curl -L -o "$sif_file" "$sif_url" || { echo "Failed to download NGPB image with curl"; exit 1; }
  else # wget
      wget -O "$sif_file" "$sif_url" || { echo "Failed to download NGPB image with wget"; exit 1; }
  fi
  echo "Soft-linking the generic image as 'NextGenPB_MCCE4.sif'"
  ln -sf "$sif_file" "$REPO_bin/NextGenPB_MCCE4.sif"

else
  echo "NGPB image file already exists: $sif_file"
  echo "Delete it before re-running the script if you need to replace it."
fi

echo ""
echo "ATTENTION: MCCE4-Alpha comes with precompiled `mcce` and `delphi` executable files"
echo "to simplify the installation. These files and the NGPB image MAY NOT WORK on your"
echo "system, but you should try them... as they might!"
echo ""
echo "Please refer to the Installation guide to obtain compiled versions for your system."
echo ""
echo "Conda environment creation..."
if ! command -v conda > /dev/null 2>&1;
then
    DO_CONDA=0
    echo "'conda' could not be found. To use MCCE4, you will need to create (or use) a dedicated"
    echo "environment with the required dependencies in 'requirements.txt'."
    echo "Skipping conda env creation."
else
    DO_CONDA=1
fi

if [ $DO_CONDA == 1 ];
then
  if [ -z "$env_name" ];
  then
    echo "Creating a dedicated conda env with default name: mc4"
    conda env create -f "$REPO_PATH/mc4.yml"
    # assign default name for activation:
    env_name="mc4"
  else
    echo "Creating a dedicated conda env with user provided name: $env_name"
    conda env create -f "$REPO_PATH/mc4.yml" -n "$env_name"
  fi
fi

echo ""
# export lines for the bin folders:
export1="export PATH=\"$REPO_bin:\$PATH\""
export2="export PATH=\"$MCCE_bin:\$PATH\""

# On macOS, login shells source .bash_profile.
rc_file="$HOME/.bash_profile"

# Check if rc file exists
no_rc_msg=""
if [ ! -f "$rc_file" ];
then
  echo "RC file ($rc_file) not found."
  # create a message for latter use:
  no_rc_msg="Create the RC file with this command: ' cd ~; touch $rc_file '"
  echo ""
else
  if grep -qF "$export1" "$rc_file";
  then
    echo "'MCCE4-Alpha/bin' already in PATH."
  fi
  if grep -qF "$export2" "$rc_file"; then
    echo "'MCCE4-Alpha/MCCE_bin' already in PATH."
  fi
fi

echo ""
echo "Almost done!"
echo ""
if [ -n "$no_rc_msg" ]; then
    echo "$no_rc_msg"
    echo ""
fi
echo "Add these export path lines to your $rc_file (if not already there), then save the file:"
echo ""
echo " $export1"
echo " $export2"
echo ""
echo "To effect the changes to your PATH variable, please start a new"
echo "terminal session or run:"
echo ""
echo " source $rc_file"
echo ""

if [ $DO_CONDA == 0 ];
then
    echo "To complete this quick installation, activate a python 3.10 environment"
    echo "that complies with requirement.txt."
else
    echo "To complete this quick installation, activate your environment:"
    echo ""
    echo "  conda activate $env_name"
fi
echo ""
echo "You can then run these commands:"
echo ""
echo "  which getpdb  # 1. check tool is found via PATH"
echo "  getpdb -h     # 2. check several python requirements"
echo ""
echo ""
echo "Thank you for using MCCE4!"
echo "Thank you for opening an issue if you encountered problems with the script:"
echo "https://github.com/GunnerLab/MCCE4-Alpha/issues"
echo ""
