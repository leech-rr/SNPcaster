#!/bin/bash
# last update: 20241106 by t-yamagishi

# Create conda environments from conda env yaml files

# parameter check
if [ "$#" -lt 1 ]; then
    CMDNAME=`basename $0`
    echo "Usage: bash $CMDNAME CondaEnvDir"
    echo '  $1  Directory having conda yaml files'
    exit 1
fi

# directory check
env_dir=$1
if [ ! -d "$env_dir" ]; then
    echo "$env_dir is not a valid directory"
    exit 1
fi

# yaml file check
env_files=($(ls $env_dir/*.yml))
if [ ${#env_files[@]} -eq 0 ]; then
    echo "$env_dir has no yaml files"
    exit 1
fi

for env_yaml in ${env_files[@]}; do
    echo "Creating conda env name=${env_yaml}"
    mamba env create -f="${env_yaml}"
done

## Need to init the conda for bash
mamba init --system bash
