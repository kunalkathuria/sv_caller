#!/bin/bash

while [[ $# -gt 1 ]]
do
key="$1"

case $key in
    -b|--bamfile)
    BAM="$2"
    shift # past argument
    ;;
    -bn|--bam_ns)
    BAM_NS="$2"
    shift # past argument
    ;;
    -r|--read_thresh)
    READTHRESH="$2"
    shift # past argument
    ;;
    -m|--match_thresh)
    MATCHRATIO="$2"
    shift # past argument
    ;;
    -n|--n_match_thresh)
    NMATCHRATIO="$2"
    shift # past argument
    ;;
    -np|--n-pct_thresh)
    NMATCH_PCT="$2"
    shift # past argument
    ;;
    -ct|--calc_thresh)
    CALC_THRESH="$2"
    shift # past argument
    ;;
    -rd|--rd_signal)
    RDS="$2"
    shift # past argument
    ;;
    -s|--sam_path)
    SAM="$2"
    shift # past argument
    ;; 
    -sr|--split_margin)
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
