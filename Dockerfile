FROM ubuntu:22.04

# no interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# install packages
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    gfortran \
    pkg-config \
    wget \
    curl \
    git \
    vim \
    gdb \
    valgrind \
    && rm -rf /var/lib/apt/lists/*

RUN apt-get update && apt-get install -y \
    libopenmpi-dev \
    openmpi-bin \
    libgsl-dev \
    libfftw3-dev \
    libhdf5-dev \
    libhdf5-openmpi-dev \
    && rm -rf /var/lib/apt/lists/*

# Set up working directory
WORKDIR /workspace

# Set environment variables for MPI
ENV OMPI_ALLOW_RUN_AS_ROOT=1
ENV OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1

# Keep container running
CMD ["/bin/bash"]

