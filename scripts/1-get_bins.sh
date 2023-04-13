#!/bin/bash -e

# Subset MAGs based on CheckM stats

# Directories
DIRIN=/nesi/project/ga02676/Waiwera_project/3.Binning/2.Final_bins
DIROUT=data/1.bins

mkdir -p $DIROUT

# Variables
COMP=70 # Completeness
CONT=5  # Contamination

# Subset CheckM results and copy MAGs to data/
awk -v comp="$COMP" -v cont="$CONT" -F "\t" '($3 > comp) && ($4 < cont) {print $6}' $DIRIN/checkm_data.txt \
  | xargs -I {} cp $DIRIN/{}.fna $DIROUT/
