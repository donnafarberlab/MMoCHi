#!/bin/bash

if [ -d "./example_envs" ]; then
    rm -rf ./example_envs
fi

mkdir ./example_envs

for i in $(tox -l)
do
  if [ -d ./.tox/$i ]; then
    echo "Exporting environment for $i"
    conda activate ./.tox/$i
    conda list -e > ./example_envs/$i-deps.txt
    conda deactivate
  fi
done
