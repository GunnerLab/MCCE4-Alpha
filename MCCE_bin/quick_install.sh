#!/bin/bash

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
CLONE_PATH="$(dirname "$SCRIPT_DIR")"
REPO_bin="$CLONE_PATH/bin"
CLONE="$(basename "$CLONE_PATH")"

if [[ "$CLONE" != "MCCE4" && "$CLONE" != "MCCE4-Alpha" ]]; then
  echo ""
  echo "'$CLONE' is not recognized as a valid MCCE4 clone for this quick install script."
  exit 1
fi

echo "$CLONE 'Quick Install' Script"

if [[ "$(uname -s)" != "Linux" ]]; then
  echo ""
  echo "$CLONE only runs on Linux platforms."
  exit 1
fi

echo ""
echo "If you encounter any problems with this script, please, open an issue here:"
echo "https://github.com/GunnerLab/$CLONE/issues/new?title=quick_install.sh&template=bug_report.md"
echo ""

if [[ "$CLONE" =~ Alpha$ ]]; then
  ENV_NAME="mc4"
  YML="mc4.yml"
else
  ENV_NAME="mc4dev"
  YML="mc4dev.yml"
fi

# Note: legacy code block below; there used to be an option for an alternate env name.
# Use the -h switch to get help on usage.
help_msg() {
  echo ""
  echo "Usage:"
  echo "bash ./MCCE_bin/$(basename "$0")"
  echo ""
  echo "Help option:"
  echo " -h  : Display this help message & exit."
  echo "     : Note: If you already have a conda environment named '$ENV_NAME', it won't be modified;"
  echo "     : yet, it will need to comply with the requirements in $CLONE_PATH/$YML."
  echo "     : To have it recreated, first remove it by running this command (conda required):"
  echo "     :  conda env remove -n $ENV_NAME"
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

# Check for curl download tool needed for apptainer installation:
if ! command -v curl >/dev/null 2>&1; then
  echo ""
  echo "STOP: Download tool 'curl' could not be found."
  echo "TODO: Please, install it with the following commands, then re-run this quick_install script:"
  echo ""
  echo "  cd ~;"
  echo "  sudo apt update; "
  echo "  sudo apt install curl; "
  echo "  cd $CLONE_PATH; "
  echo ""
  exit 1
fi

echo ""
echo "Checking for conda..."
DO_CONDA=0
if command -v conda >/dev/null 2>&1; then
  echo "  Conda found."
  DO_CONDA=1
else
  echo "STOP: conda is required, but not found on your system."
  echo "TODO: Please, install miniconda with the following commands (with lines 2-4 obtained "
  echo "      from this site: https://www.anaconda.com/docs/getting-started/miniconda/install), "
  echo "      then re-run this quick install script:"
  echo ""
  # Note: in the third line below, -b :: batch mode: yes to all questions
  echo "  cd ~; "
  echo "  curl -OL https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh; "
  echo "  bash ~/Miniconda3-latest-Linux-x86_64.sh -b; "
  echo "  ~/miniconda3/bin/conda init bash; "
  echo "  source ~/.bashrc; "
  echo "  cd $CLONE_PATH "
  echo ""
  exit 1
fi

rc_file="$HOME/.bashrc"
if [[ ! -f "$rc_file" ]]; then
  touch "$rc_file"
fi

# export lines for the bin folders:
modified=0
export1="export PATH=\"$REPO_bin:\$PATH\""
export2="export PATH=\"$SCRIPT_DIR:\$PATH\""

echo ""
echo "Checking PATH variables..."
# Check if rc file needs updating; export lines with this repo path
# cannot be commented out since they're needed; if so, there added back.
if grep -q "^${export1}" "$rc_file"; then
  echo "  '$CLONE/bin' already in PATH."
else
  echo "$export1" >> "$rc_file"
  modified=1
fi
if grep -q "^{$export2}" "$rc_file"; then
  echo "  '$CLONE/MCCE_bin' already in PATH."
else
  echo "$export2" >> "$rc_file"
  modified=1
fi
if [[ "$modified" -eq 1 ]]; then
    source "$rc_file"
    cd "$CLONE_PATH"
fi

echo ""
echo "Checking for apptainer (container app for NGPB)..."
APPTAINER_FOUND=0
APPTAINER_CONDA=0
export0=""

if command -v apptainer >/dev/null 2>&1; then
  APPTAINER_FOUND=1
  aptfp=$(which apptainer)
  echo "  Container app 'apptainer' found: $aptfp"
else
  echo "  Container app 'apptainer' is not found on your system."
  if ! command -v rpm2cpio >/dev/null 2>&1; then
    APPTAINER_CONDA=1
    echo "FYI: The 'unprivileged' installation of Apptainer requires rpm2cpio and cpio programs, which you do not have."
    echo "     Apptainer will be installed with conda in the '$ENV_NAME' environment."
    echo ""
  else
    echo "Installing unprivileged version of apptainer..."
    curl -s https://raw.githubusercontent.com/apptainer/apptainer/main/tools/install-unprivileged.sh | bash -s - "$HOME/apptainer"
    APPTAINER_FOUND=1
    export0="export PATH=\"$HOME/apptainer/bin:\$PATH\""
    if ! grep -qF "$export0" "$rc_file"; then
      echo "$export0" >> "$rc_file"
      source "$rc_file"
      cd "$CLONE_PATH"
    fi
  fi
fi

# Check if user namespaces are enabled (required for unprivileged Apptainer)
if command -v unshare >/dev/null 2>&1 && ! unshare -U true >/dev/null 2>&1; then
  echo ""
  echo "WARNING: User namespaces are not enabled! Apptainer will likely fail."
  echo "       If running inside Docker, ensure you use '--security-opt seccomp=unconfined'."
  echo "       Otherwise, contact your system administrator to enable user namespaces."
  echo ""
fi

if [[ "$DO_CONDA" -eq 1 ]]; then
  echo ""
  echo "Conda environment creation..."
  ENVS="$(conda env list)"
  # check if env with default name already exists:
  if echo "$ENVS" | grep -q "^${ENV_NAME}\s\+"; then
    echo "  Updating existing '$ENV_NAME' env..."
    conda env update -n $ENV_NAME -f "$CLONE_PATH/$YML"
  else
    echo "  Creating a dedicated conda env with default name: $ENV_NAME"
    conda env create -f "$CLONE_PATH/$YML"
  fi
  if [[ "$APPTAINER_CONDA" -eq 1 ]]; then
    conda install -n $ENV_NAME -c conda-forge apptainer -y
  fi
fi

RELEASE="NextGenPB_v1.0.0"
# Name of the downloaded generic image & source url:
sif_file="$REPO_bin/NextGenPB.sif"
# Name and path of the compiled image when Makefile is used or as link target to generic sif:
mc4_sif_name="NextGenPB_MCCE4.sif"
sif_mc4="$REPO_bin/$mc4_sif_name"
if [[ "$APPTAINER_FOUND" -eq 1 ]]; then
  echo ""
  echo "Checking for NGPB image file of release $RELEASE..."
  if [[ -f "$sif_mc4" ]]; then
    echo "  NextGenPB_MCCE4.sif is already installed."
    # check if "$sif_mc4" is soft-linked & to "$sif_file":
    if [[ -L "$sif_mc4" ]]; then
      resolved_target="$(readlink -f "$sif_mc4")"
      if [[ "$resolved_target" = "$sif_file" ]]; then
        echo "  $sif_mc4 is soft-linked to the generic image $sif_file. If you want to replace"
        echo "  the generic image, delete the file, then re-run this quick install script."
        chmod +x "$sif_file"
      fi
    fi
  else
    if [[ -f "$sif_file" ]]; then
      echo "  Found NGPB generic image file, $sif_file. If you want to replace it,"
      echo "  delete the file, then re-run this quick install script."
      echo ""
    else
      echo "Getting NextGenPB generic container image from $sif_url (~1.6GB)..."
      sif_url="https://github.com/concept-lab/NextGenPB/releases/download/$RELEASE/NextGenPB.sif"
      curl -L -o "$sif_file" "$sif_url" || { echo "Failed to download NGPB image with curl"; exit 1; }
    fi
    echo "  Soft-linking the generic image as $mc4_sif_name."
    chmod +x "$sif_file"
    ln -sf "$sif_file" "$sif_mc4"
  fi
fi

if [[ "$DO_CONDA" -eq 1 ]]; then
  echo ""
  echo "Now you can run these commands:"
  echo ""
  echo "  conda activate $ENV_NAME"
  echo "  which getpdb      # :: check tool is found via PATH"
  echo "  getpdb -h         # :: check several python requirements"
  echo ""
fi
