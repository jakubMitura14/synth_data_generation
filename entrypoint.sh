#!/bin/bash
set -e

echo "Starting Python script execution..."

python3 -m pip uninstall -y numpy && \
    python3 -m pip install numpy==1.23.2 nii2dcm

echo "Executing Python script..."
python3 /root/.devcontainer/in_docker_organized/coordinate_phantom_create.py

# Keep container running
echo "Script completed, keeping container alive..."
exec tail -f /dev/null


# julia in_docker_organized/main_create_phantom_can.jl 256x256x256 false true f21d088e-d5d3-4252-807f-efc036d3153c true true 0.1