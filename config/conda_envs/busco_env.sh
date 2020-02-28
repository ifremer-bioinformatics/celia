#!/usr/bin/env bash
export PATH=$PATH:/appli/anaconda/latest/bin
source activate /appli/conda-env/bioinfo/busco-4.0.4
if [ ! -d $HOME/augustus_config/ ]
then
        mkdir -p $HOME/augustus_config/
        cp -r /appli/conda-env/bioinfo/busco-4.0.4/config/* $HOME/augustus_config/
fi
export AUGUSTUS_CONFIG_PATH=$HOME/augustus_config/
