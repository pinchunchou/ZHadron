#!/bin/sh
source /cvmfs/cms.cern.ch/cmsset_default.sh


CMSSWDir=/afs/cern.ch/user/p/pchou/CMSSW_12_5_2_patch1/src

cd $CMSSWDir

cmsenv

cd /afs/cern.ch/user/p/pchou/ZHadron

root -l -q -b "ZtrackAnalysis2.C($1)"
