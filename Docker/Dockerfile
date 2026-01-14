FROM continuumio/miniconda3

# Install system tools
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    build-essential \
    curl \
    git \
    iproute2 \
    net-tools \
    openssh-client \
    sudo \
    tree && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

RUN cd /tmp && \
    wget https://github.com/apptainer/apptainer/releases/download/v1.4.5/apptainer_1.4.5_amd64.deb && \
    apt-get update && apt-get install -y \
    ./apptainer_1.4.5_amd64.deb  && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Install Node via Conda
RUN conda install -c conda-forge nodejs

# Create user 'mc4', switch to the new user
RUN useradd -ms /bin/bash mc4
USER mc4
WORKDIR /home/mc4

# Optional: Ensure bash is used for RUN commands later (if you add any)
SHELL ["/bin/bash", "-l", "-c"]
