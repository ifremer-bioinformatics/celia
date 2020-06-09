#!/usr/bin/env bash

BASEDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"

#activate nextflow environment
. ${BASEDIR}/config/conda_envs/nextflow_env.sh

#nextflow temp directory
export NXF_TEMP=$DATAWORK/nextflow/nxt-tmp

#run nextflow nextmb workflow ($1 is useful if you want to run resume)
nextflow -DXms=2G -DXmx=8G -trace nextflow.executor run CELIA.nf $1

#deactivate nextflow environment
. ${BASEDIR}/config/conda_envs/delenv.sh
