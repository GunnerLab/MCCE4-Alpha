#!/bin/bash

if [ $# -eq 0 ]; then
    echo "Usage: run <script.py> [args]"
    exit 1
fi

apptainer exec mcce4-alpha.sif "$@"