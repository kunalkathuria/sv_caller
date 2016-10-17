#!/bin/bash

while [[ $# -gt 1 ]]
do
key="$1"

case $key in
    -b|--BamFile)
    BAM="$2"
    shift # past argument
    ;;
    -bn|--BamFile_NS)
    BAM_NS="$2"
    shift # past argument
    ;;
    -r|--ReadThresh)
    READTHRESH="$2"
    shift # past argument
    ;;
    -m|--MatchRatioThresh)
    MATCHRATIO="$2"
    shift # past argument
    ;;
    -n|--NMatchRatioThresh)
    NMATCHRATIO="$2"
    shift # past argument
    ;;
    -np|--NMatchPctThresh)
    NMATCH_PCT="$2"
    shift # past argument
    ;;
    -ct|--CalcThresh)
    CALC_THRESH="$2"
    shift # past argument
    ;;
    -rd|--ReadDepthSignal)
    RDS="$2"
    shift # past argument
    ;;
    -s|--SamtoolsPath)
    SAM="$2"
    shift # past argument
    ;; 
    -sr|--MinClusterSize)
    SR_MARGIN="$2"
    shift # past argument
    ;;
    -bp|--bp_margin)
    BP_MARGIN="$2"
    shift # past argument
    ;;
    -sl|--slop)
    SLOP="$2"
    shift # past argument
    ;;
    -ro|--rec_overlap)
    RO="$2"
    shift # past argument
    ;;
    -vm|--var_margin)
    Variant_Margin="$2"
    shift # past argument
    ;;
    -rr|--refresh_rate)
    REF_RATE="$2"
    shift # past argument
    ;;
    -dt|--disj_thresh)
    DTHRESH="$2"
    shift # past argument
    ;;
    -sp|--splitter_path)
    S_PATH="$2"
    shift # past argument
    ;;
    -ss|--splitter_slop)
    SS="$2"
    shift # past argument
    ;;
    -risk|--splitter_risk)
    S_RISK="$2"
    shift # past argument
    ;;
    -srr|--splitter_refresh)
    SRR="$2"
    shift # past argument
    ;;		
    --default)
    DEFAULT=1
    ;;
    *)
    ;;
esac
shift # past argument or value
done

./runCode.sh $READTHRESH $MATCHRATIO $NMATCHRATIO $MIN_CS $NMATCH_PCT $CALC_THRESH $RDS $BAM $BAM_NS $SAM

# test case
# ./runAll.sh
