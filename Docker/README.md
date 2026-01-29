# MCCE4-Alpha Docker Environment

This directory contains the Docker configuration for the MCCE4-Alpha project. It provides a consistent build and runtime environment, specifically configured to support **Apptainer** (formerly Singularity) running inside the container.

## Prerequisites

- [Docker Engine](https://docs.docker.com/engine/install/)
- [Docker Compose](https://docs.docker.com/compose/install/)
- [Apptainer](https://apptainer.org/docs/admin/main/installation.html#install-unprivileged-from-pre-built-binaries) *(optional: only if converting to .sif locally)*

## Getting Started

### 1. Build the Docker Image
From the root directory of the project initialize the base environment:
```bash
docker build -t mcce4-alpha -f Docker/Dockerfile .
```
>[!IMPORTANT]
There are two dockerfiles. One is for production ready code, `Dockerfile`, while the other is for development purposes, `Dockerfile.dev`. <br> We need to run this command from the root directory so that the build context can locate the necessary files to copy into the conatiner.

### 2. Convert to Apptainer (.sif)
Flatten the Docker layers into a portable, read-only Apptainer image:
```bash
apptainer -v build mcce4-alpha.sif docker-daemon://mcce4-alpha:latest
```

### 3. Verification
Verify that the *Garage Floor* environment logic (PATH hijacking) is active. The command below should point to the Conda environment, not the system Python.
```bash
# Execute the check
apptainer exec mcce4-alpha.sif which python
# Expected Output:
/home/mc4/.conda/envs/mc4/bin/python
```

### 4. Usage
You should now be able to run executables that are baked into the image by running the following:
```bash
apptainer exec mcce4-alpha.sif <exectuable_name>
```

## Security & Infrastructure
This environment is configured to support **Apptainer-in-Docker**. Standard Docker containers restrict capabilities required by Apptainer. To enable these features securely, the `docker-compose.yml` applies the following settings:

*   **`security_opt: - seccomp:unconfined`**: Allows the `unshare` system call, which Apptainer uses to create User Namespaces (required for unprivileged execution).
*   **`devices: - /dev/fuse`**: Exposes the FUSE device, allowing Apptainer to mount `.sif` images without root privileges.
>[!IMPORTANT]
If running manually via docker run, you must include these flags: `--security-opt seccomp=unconfined --device /dev/fuse`

## Development Workflow
This workflow is ***specifically*** for developers who want to test changes within an isolated container environment. It allows you to experiment with new scripts, tools, or configurations without cluttering your local machine or risking dependency conflicts.

### 1. Initialize Environment
First, generate a `.env` file tailored to yoru specific User and Group IDs.
```bash
cd MCCE4-Alpha/Docker/scripts
./setup_env.sh
```
**Why we do this**: Linux identifies file owners by numeric IDs (UID/GID), not names. Your laptop user is likely `1000`. This script detects your specific IDs and saves them to `Docker/.env`. When Docker builds the image, it uses these numbers to "morph" the internal `mc4` user to match your identity exactly.

### 2. Build and Start the Container
Navigate back to the `Docker` directory and launch the environment.
```bash
cd ../
docker compose up --build -d
```
**What this does**:
- `--build`: Rebuilds the image to ensure your detached UID/GID and any Dockerfile chages are applied
- `d`: Runs the container in the background so you can continue using your terminal

### 3. Enter the Development Environment
Open an interactive bash session inside the isolated container to start testing.
```bash
docker exec -it MCCE4-Alpha /bin/bash
```
### Volume Mounts
>[!NOTE]
The project root is mounted to `/home/mc4/MCCE4-Alpha`. Because of the `setup_env.sh` step, any file you create or edit inside the container will be owned by you on your host machine, and vice-versa. You can write code in your favorite IDE (VS Code, etc.) on your laptop, and immediately run it inside the container.
