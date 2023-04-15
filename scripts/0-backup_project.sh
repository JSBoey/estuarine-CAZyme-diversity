#!/bin/bash -e

TARGET=/nesi/project/ga02676/Waiwera_project/boey_work/estuarine-CAZyme-diversity

cp -v 2023-estuarine-cazyme-diversity.ipynb $TARGET
cp -v -r scripts/ $TARGET
cp -v -r docs/ $TARGET
cp -v -r results/ $TARGET

# Don't want to copy the interproscan stuff
rsync -av --exclude="interproscan*" bin $TARGET