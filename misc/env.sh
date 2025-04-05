#!/bin/bash

wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
chmod +x Miniforge3-Linux-x86_64.sh

echo "Installing Miniforge3..., follow the instructions on the screen"

./Miniforge3-Linux-x86_64.sh

echo "Miniforge3 installed successfully"

echo "Setting up the environment..."
conda create -n run3-vbsvvh -y
conda activate run3-vbsvvh

mamba install -y python==3.10.11
mamba install -y pip
mamba install root boost correctionlib mkl
pip install mplhep numpy pandas uproot

echo "Environment setup complete"