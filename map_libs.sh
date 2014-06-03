#!/usr/bin/env bash

USE_LIBS="
 SUSYTools GoodRunsLists\
 ApplyJetResolutionSmearing ApplyJetCalibration CalibrationDataInterface\
 egammaAnalysisUtils\
 egammaEvent JetResolution JetUncertainties\
 MissingETUtility MuonEfficiencyCorrections MuonMomentumCorrections\
 PileupReweighting ReweightUtils\
 TauCorrUncert BCHCleaningTool PileupReweighting\
 TrigMuonEfficiency ElectronEfficiencyCorrection PATCore TileTripReader"

USE_HEADER="SUSYTools MissingETUtility egammaAnalysisUtils GoodRunsLists\
 CalibrationDataInterface TauCorrUncert PileupReweighting"


if [[ $1 == '-l' ]]
then
    lib_path=${ROOTCOREDIR}/lib/${ROOTCORECONFIG#/}
    echo -n -L${lib_path%/}" "
    echo -n -Wl,-rpath,${lib_path%/}" "
    for l in $USE_LIBS
    do
	echo -n -l${l}\ 
    done
elif [[ $1 == '-i' ]]
then
    for i in $USE_HEADER
    do 
	echo -n ${ROOTCOREDIR}/../$i\ 
    done
else 
    echo "ERROR: $0 needs flag -i or -l" >&2 
    false 
fi

