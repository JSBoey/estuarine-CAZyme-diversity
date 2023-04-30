#!/bin/bash -e

TARGET=/nesi/project/ga02676/Waiwera_project/boey_work/estuarine-CAZyme-diversity

# Don't want to copy the interproscan stuff
rsync -av --exclude="interproscan*" bin $TARGET
rsync -av 2023-estuarine-cazyme-diversity.md $TARGET
rsync -av scripts $TARGET
rsync -av docs $TARGET
rsync -av results $TARGET
