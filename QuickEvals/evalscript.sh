#!/bin/zsh

#PBS -q N

echo Start job at `date` at `hostname`

export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
alias setupATLAS='source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh'

setupATLAS

localSetupROOT --rootVersion=5.34.32-HiggsComb-x86_64-slc6-gcc48-opt --skipConfirm

cd /mnt/Lustre/agrp/jonathas/btagntuple_analysis_framework

./run_framework -o QuickEvals/Parts/$OUTPUTNAME -r TaggerLists/$TAGGERFILE QuickEvals/FileLists/$TXTLIST
