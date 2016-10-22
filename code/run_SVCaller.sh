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
    -np|--n_pct_thresh)
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
    -srr|--splitter_rate)
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
    -vm|--var_rate)
    VAR_RATE="$2"
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
    -mc|--min_cluster)
    MIN_CS="$2"
    shift # past argument
    ;;
    -sm|--sig_mult)
    SIG_MULT="$2"
    shift # past argument
    ;;
    -sb|--sig_bound)
    SIG_BOUND="$2"
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

./form_clusters.sh $BAM $BAM_NS $READ_THRESH $MATCHRATIO $NMATCHRATIO $CALC_THRESH $NMATCH_PCT $MIN_CS $SIG_MULT $BP_MARGIN $SIG_BOUND
./run_PE.sh $MIN_CS $VAR_RATE $SLOP $REF_RATE $DTHRESH # can replace $MIN_CS here with another value and run these 2 steps with new threshold rather than rerun whole cluster formation sequence. In this case, simply comment out previous line and run sv_caller.
./run_RD_SR.sh $S_PATH $S_RISK $SS $SR_MARGIN $RO $DTHRESH $RD
