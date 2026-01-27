# MCCE4-Alpha Docker Environment

This directory contains the Docker configuration for the MCCE4-Alpha project. It provides a consistent build and runtime environment, specifically configured to support **Apptainer** (formerly Singularity) running inside the container.

## Prerequisites

- [Docker Engine](https://docs.docker.com/engine/install/)
- [Docker Compose](https://docs.docker.com/compose/install/)
- Apptainer (optional: only if converting to .sif locally)

## Getting Started

### 1. Build the Docker Image
Initialize the base environment from the project root:
```bash
docker build -t mcce4-alpha -f Docker/Dockerfile .
```

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

**[!IMPORTANT]** If running manually via docker run, you must include these flags: `--security-opt seccomp=unconfined --device /dev/fuse`

## Development Workflow

### Volume Mounts
The project root (`../`) can be mounted to `/home/mc4/MCCE4-Alpha` inside the container. Any changes made to the source code on your host machine will be immediately visible inside the container.

### Troubleshooting
If you encounter errors like `Failed to create user namespace` when running Apptainer:
1. Ensure you are starting the container with `docker compose up` (which applies the security config).
2. If running manually via `docker run`, you must include `--security-opt seccomp=unconfined --device /dev/fuse`.