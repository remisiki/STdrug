#!/usr/bin/env bash

# Hyperparameters
dataName="hcc"
tissue="liver"
nclust=5

# File paths
projectDir=$(pwd)
outputDir="${projectDir}/output"
checkpointDir="${projectDir}/output/checkpoint"
dataDir="${projectDir}/data"

echo "Work in dir ${projectDir}, outputs will be saved at ${outputDir}"

# Spatial domain identification
source venv/bin/activate
echo "Using Python from $(which python3), $(python3 -V)"

python3 -u domain_matching/main.py \
  --data-name "${dataName}" \
  --nclust "${nclust}" \
  --output "${outputDir}"

# Drug repurposing
echo "Using R from $(which Rscript), $(Rscript --version)"

Rscript drug_score/Main.R \
  --data-name "${dataName}" \
  --tissue "${tissue}" \
  --output "${outputDir}" \
  --public-data-path "${dataDir}" \
  --cluster-output-path "${outputDir}" \
  --checkpoint-path "${checkpointDir}"
