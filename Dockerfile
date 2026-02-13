# Use a standard Ubuntu base image with Python pre-installed
FROM python:3.10-slim-bullseye

ARG JULIA_RELEASE=1.10
ARG JULIA_VERSION=1.10.4

ENV DEBIAN_FRONTEND=noninteractive

# Install basic dependencies
RUN apt-get update && \
    apt-get install --yes --no-install-recommends \
    curl ca-certificates nano wget \
    build-essential software-properties-common \
    git cmake ninja-build zip unzip \
    libssl-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install Julia
RUN curl -s -L https://julialang-s3.julialang.org/bin/linux/x64/${JULIA_RELEASE}/julia-${JULIA_VERSION}-linux-x86_64.tar.gz | \
    tar -C /usr/local -x -z --strip-components=1 -f -

# Install additional required packages
RUN apt-get update && \
    apt-get install --yes --no-install-recommends \
    vim net-tools ffmpeg pkg-config && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install libcurl from source (keeping this for compatibility)
RUN wget https://curl.se/download/curl-7.81.0.tar.gz && \
    tar -xvf curl-7.81.0.tar.gz && cd curl-7.81.0 && \
    ./configure --with-openssl && make && make install && \
    cd / && rm -rf curl-7.81.0*

# Set up user environment
RUN mkdir -m 0777 /data
RUN mkdir -p $HOME/data

ENV JULIA_HISTORY=/data/logs/repl_history.jl

# Copy the data download script
COPY download_data.py /workspace/download_data.py
RUN python3 /workspace/download_data.py

# Setup thread management for Julia
COPY setup_threads.sh /etc/profile.d/

# Install Google Cloud SDK
RUN curl -sSL https://sdk.cloud.google.com | bash

# Install Julia packages
RUN julia -e 'using Pkg; Pkg.add(url="https://github.com/jakubMitura14/ImagePhantoms.jl.git")'
RUN julia -e 'using Pkg; Pkg.add(["Meshes","ImageFiltering","Accessors","UUIDs","JSON","HDF5","ImagePhantoms","MIRTjim","ImageGeoms","Sinograms", "PyCall", "FFTW", "LazyGrids", "Unitful", "Plots", "Revise"])'

# Set up Google Cloud credentials
ENV PATH $PATH:/root/google-cloud-sdk/bin
COPY infostrateg-b-f395071711c1.json /root/.devcontainer/infostrateg-b-f395071711c1.json
ENV GOOGLE_APPLICATION_CREDENTIALS="/root/.devcontainer/infostrateg-b-f395071711c1.json"
RUN gcloud auth activate-service-account --key-file="/root/.devcontainer/infostrateg-b-f395071711c1.json"

# Setup Weights & Biases
ENV WANDB_API_KEY=5a6c2e2caccb19a9c3cbfde388b61c7104eab632

# Copy application files
COPY generate-random_can_low_res.jl /root/.devcontainer/generate-random_can_low_res.jl
COPY get_geometry_main.jl /root/.devcontainer/get_geometry_main.jl
COPY CtFanArc_params.jl /root/.devcontainer/CtFanArc_params.jl

# Install Python packages (without CUDA dependencies)
RUN python3 -m pip install SimpleITK wandb adrt h5py scikit-image pydicom==2.4.4 
#pydicom-seg

# Install nii2dcm from source
RUN git clone https://github.com/tomaroberts/nii2dcm.git /tmp/nii2dcm && \
    cd /tmp/nii2dcm && \
    python3 -m pip install -r requirements.txt && \
    python3 -m pip install . && \
    cd / && \
    rm -rf /tmp/nii2dcm

# Ensure numpy compatibility
RUN python3 -m pip install numpy==1.23.2

# Copy entrypoint and organized code 
COPY entrypoint.sh /root/.devcontainer/entrypoint.sh
COPY in_docker_organized/ /root/.devcontainer/in_docker_organized/

# Make entrypoint executable
RUN chmod +x /root/.devcontainer/entrypoint.sh

WORKDIR /root/.devcontainer
ENTRYPOINT ["/bin/bash", "/root/.devcontainer/entrypoint.sh"]

#rozmycie i szum wyciagniecie tez do glownego ; ilosc do glownego
# w glownym sciezka do bucketa i tam lista jsonow