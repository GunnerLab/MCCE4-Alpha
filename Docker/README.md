# MCCE4-Alpha Docker Environment

This directory contains the Docker configuration for the MCCE4-Alpha project. It provides a consistent build and runtime environment, specifically configured to support **Apptainer** (formerly Singularity) running inside the container.

## Prerequisites

- [Docker Engine](https://docs.docker.com/engine/install/)
- [Docker Compose](https://docs.docker.com/compose/install/)

## Usage

### 1. Build the Image
Navigate to this directory and build the Docker image:
```bash
cd Docker
docker compose build
```

### 2. Start the Container
Run the container in detached mode:
```bash
docker compose up -d
```

### 3. Access the Shell
Open an interactive shell inside the running container:
```bash
docker exec -it MCCE4-Alpha bash
```

### 4. Stop the Container
When finished, stop and remove the container:
```bash
docker compose down
```

## Configuration & Security

This environment is configured to support **Apptainer-in-Docker**. Standard Docker containers restrict capabilities required by Apptainer. To enable these features securely, the `docker-compose.yml` applies the following settings:

*   **`security_opt: - seccomp:unconfined`**: Allows the `unshare` system call, which Apptainer uses to create User Namespaces (required for unprivileged execution).
*   **`devices: - /dev/fuse`**: Exposes the FUSE device, allowing Apptainer to mount `.sif` images without root privileges.

## Volume Mounts

The project root (`../`) is mounted to `/home/mc4/MCCE4-Alpha` inside the container. Any changes made to the source code on your host machine will be immediately visible inside the container.

## Creating an Apptainer Image (.sif)

If you prefer to run MCCE4 using Apptainer (Singularity) directly on your host machine, you can convert the built Docker image into a `.sif` file.

1.  **Build the Docker image** (if not already done):
    ```bash
    cd Docker
    docker compose build
    ```

2.  **Build the SIF file** using Apptainer:
    ```bash
    apptainer build mcce4-alpha.sif docker-daemon://mcce4-alpha:latest
    ```

## Troubleshooting

If you encounter errors like `Failed to create user namespace` when running Apptainer:
1. Ensure you are starting the container with `docker compose up` (which applies the security config).
2. If running manually via `docker run`, you must include `--security-opt seccomp=unconfined --device /dev/fuse`.