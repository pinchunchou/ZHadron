#!/bin/sh
source /cvmfs/cms.cern.ch/cmsset_default.sh


CMSSWDir=/afs/cern.ch/user/p/pchou/CMSSW_12_5_2_patch1/src

cd $CMSSWDir

cmsenv

#root -l -q -b "/afs/cern.ch/user/p/pchou/ZHadron/Zhadron_loop_part.C(\"/eos/cms/store/group/phys_heavyions/chenyi/PbPb2018/Forest/DYJetsToLL_MLL-50_TuneCP5_HydjetDrumMB_5p02TeV-amcatnloFXFX-pythia8/HINPbPbSpring21MiniAOD-mva98_112X_upgrade2018_realistic_HI_v9-v1/MINIAODSIM/DYJetsToLL_MLL-50_TuneCP5_HydjetDrumMB_5p02TeV-amcatnloFXFX-pythia8/20221017_ZHadronInitialRunMCDY/221017_113040/000$1/\",\"/eos/cms/store/group/phys_heavyions_ops/pchou/MC/output_MC000$1_part$3.root\",$2,$3,50)"
root -l -q -b "/afs/cern.ch/user/p/pchou/ZHadron/Zhadron_loop_part.C+(\"/eos/cms/store/group/phys_heavyions/chenyi/PbPb2018/Forest/HIDoubleMuon/HIRun2018A-PbPb18_MiniAODv1-v1/MINIAOD/HIDoubleMuon/20221101_ZHadronDoubleMuonWithMuTree/221101_212941/000$1/\",\"/eos/cms/store/group/phys_heavyions_ops/pchou/doubleMu/output_doubleMu000$1_part$3.root\",$2,$3,50)"
#root -l -q -b "/afs/cern.ch/user/p/pchou/ZHadron/Zhadron_loop_part.C(\"/eos/cms/store/group/phys_heavyions/chenyi/PbPb2018/Forest/HIHardProbes/HIRun2018A-PbPb18_MiniAODv1-v1/MINIAOD/HIHardProbes/20221101_ZHadronHardProbesWithMuTree/221101_212923/000$1/\",\"/eos/cms/store/group/phys_heavyions_ops/pchou/HardProbes/output_HardProbes000$1_part$3.root\",$2,$3)"
#root -l -q -b "/afs/cern.ch/user/p/pchou/ZHadron/Zhadron_loop_part.C(\"/eos/cms/store/group/phys_heavyions/chenyi/PbPb2018/Forest/HIMinimumBias7/HIRun2018A-PbPb18_MiniAODv1-v1/MINIAOD/HIMinimumBias7/20221101_ZHadronMinimumBias7WithMuTree/221101_212932/000$1/\",\"/eos/cms/store/group/phys_heavyions_ops/pchou/MinBias7/output_MinBias7000$1_part$3.root\",$2,$3,2)"