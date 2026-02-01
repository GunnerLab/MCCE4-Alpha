#!/bin/bash

# Exit immediately if user hits Ctrl+C
trap "echo '❌ Script interrupted by user'; exit 1" SIGINT

PROJECT_ROOT_DIR=$(dirname "$(readlink -f "$0")")
DOCKER_DIR=$(readlink -f "$PROJECT_ROOT_DIR/Docker")
BIN_DIR=$(readlink -f "$PROJECT_ROOT_DIR/bin")

CONDA_YML="$PROJECT_ROOT_DIR/mc4.yml"
CONDA_ENV_NAME="mc4"

SHELL_CONFIG="$HOME/.bashrc"
if [[ "$SHELL" == *"zsh"* ]]; then
    SHELL_CONFIG="$HOME/.zshrc"
fi

update_shell_config() {
    local dir_path="$1"
    local tool_name="$2"
    local export_line="export PATH=\"$dir_path:\$PATH\""

    if grep -qF "$dir_path" "$SHELL_CONFIG"; then
        echo "✅ The path $dir_path is already configured in $SHELL_CONFIG."
    else
        echo "⚠️  The path $dir_path is not in your PATH."
        echo "   Would you like to append it to $SHELL_CONFIG now? (y/n)"
        read -r -p "   > " response
        
        if [[ "$response" =~ ^([yY][eE][sS]|[yY])$ ]]; then
            echo "" >> "$SHELL_CONFIG"
            echo "# $tool_name" >> "$SHELL_CONFIG"
            echo "$export_line" >> "$SHELL_CONFIG"
            echo "✅ Successfully added to $SHELL_CONFIG"
            export PATH="$dir_path:$PATH"
        else
            echo "🆗 Skipped. You can manually add the line below to your $SHELL_CONFIG:"
            echo "   $export_line"
        fi
    fi
}

setup_conda_env() {
    local force_apptainer="$1"

    echo "🐍 Setting up Conda Environment..."

    # 1. Strict Check: Is 'conda' available right now?
    if ! command -v conda >/dev/null 2>&1; then
        echo "❌ Error: 'conda' command not found in your PATH."
        echo "   Possible reasons:"
        echo "   1. Conda is not installed."
        echo "   2. Conda is installed but not activated."
        echo ""
        echo "   👉 Fix: Run 'source <path/to/conda>/bin/activate' or 'conda init', then retry."
        return 1
    fi

    local envs_list
    envs_list=$(conda env list)
    
    if echo "$envs_list" | grep -q "^${CONDA_ENV_NAME}\s\+"; then
        echo "   Updating existing environment '$CONDA_ENV_NAME'..."
        conda env update -n "$CONDA_ENV_NAME" -f "$CONDA_YML"
    else
        echo "   Creating new environment '$CONDA_ENV_NAME'..."
        conda env create -f "$CONDA_YML"
    fi

    # Fallback: Install Apptainer inside Conda if local install failed
    if [[ "$force_apptainer" -eq 1 ]]; then
        echo "⬇️  Installing Apptainer via Conda (Fallback)..."
        conda install -n "$CONDA_ENV_NAME" -c conda-forge apptainer -y
        echo "✅ Apptainer installed in '$CONDA_ENV_NAME'."
    fi
}

install_apptainer() {
    echo "🔍 Checking for Apptainer..."
    
    if command -v apptainer >/dev/null 2>&1; then
        echo "✅ Apptainer found at: $(which apptainer)"
        return 0
    fi

    echo "⚠️  Apptainer not found."
    
    # --- Attempt 1: Auto-Install Local (Preferred) ---
    if command -v rpm2cpio >/dev/null 2>&1; then
        echo "⬇️  Installing 'unprivileged' Apptainer locally..."
        curl -s https://raw.githubusercontent.com/apptainer/apptainer/main/tools/install-unprivileged.sh | bash -s - "$HOME/apptainer"
        
        # Ask to add this new binary to path
        update_shell_config "$HOME/apptainer/bin" "Apptainer (Unprivileged)"
        return 0
    else
        echo "❌ Missing dependency 'rpm2cpio'. Cannot install locally."
    fi

    # --- Attempt 2: Fallback to Conda ---
    echo "⚠️  Switching to Conda installation method..."
    # If this fails (returns 1), the whole function returns 1
    setup_conda_env 1
}

echo "🔧 Starting MCCE4-Alpha Setup..."

if ! install_apptainer; then
    echo ""
    echo "⛔ SETUP FAILED: Could not install Apptainer via Local OR Conda methods."
    echo "   Please fix the errors above (e.g., install Conda) and run setup again."
    exit 1
fi

update_shell_config "$BIN_DIR" "MCCE4-Alpha CLI"

echo "🚀 Setup Complete!"