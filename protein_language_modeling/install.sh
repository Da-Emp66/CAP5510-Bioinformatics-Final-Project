#!/bin/bash

UV_PROJECT_ROOT=protein_language_modeling
DEV=""

if [ "$1" = "-dev" ] || [ "$1" = "--dev" ]; then
    DEV="true"
fi

init_conda () {
    ~/miniconda3/bin/conda init bash
    ~/miniconda3/bin/conda init zsh
    source ~/miniconda3/bin/activate
    conda init --all
}

install_miniconda () {
    # Try to initialize conda for the command check
    init_conda

    # If miniconda is not installed, install it
    if ! conda ; then
        mkdir -p ~/miniconda3
        wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
        bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
        rm ~/miniconda3/miniconda.sh
        init_conda
    fi
}

setup_miniconda_environment () {
    # Activate conda if it is not already activated
    init_conda

    # Install the conda environment if not already created,
    # activate the conda environment if it is already created
    if conda env list | grep "$CONDA_ENV" >/dev/null 2>/dev/null ; then
        conda activate "$CONDA_ENV"
    else
        conda create -n "$CONDA_ENV" python=$PYTHON_VERSION -y
        conda activate "$CONDA_ENV"
    fi
}

setup_uv () {
    # Install uv if not already installed
    UV_PATH=$(command -v uv)
    if [ -z "$UV_PATH" ]; then
        curl -LsSf https://astral.sh/uv/install.sh | sh
    fi

    # Activate uv for the current shell
    source $HOME/.local/bin/env

    # Install the project
    cd $UV_PROJECT_ROOT && \
        uv sync && \
        cd ..
}

source $UV_PROJECT_ROOT/.env

setup_uv

if [ "$DEV" = "true" ]; then
    code .
fi
