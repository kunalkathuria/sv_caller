# sv_caller
Structural Variant Caller

Example call:

./run_SVCaller_d.sh (has all default values of parameters and files as noted/stored once in headers)

./run_SVCaller.sh [options]

Options:

    -b|--bamfile 		(position-sorted bam file)
    -bn|--bam_ns 		(name-sorted bam file)
    -r|--read_thresh 		(max number of discordant read mappings a single PE fragment mate can have so as not to be ignored completely, default = 20)
    -m|--match_thresh 		(min threshold of ratio of secondary alignment score/AS with AS of primary alignment needed for secondary alignment to be considered/used in set cover, def = 1)
    -n|--n_match_thresh 	(same as above for number of base pair matches, instead of AS, def = .95) 
    -np|--n_pct_thresh 		(min pct of ratio of n_bases of secondary alignment matching to n_bases of primary alignment matching needed for secondary alignment to be considered, default =  0)
    -ct|--calc_thresh 		(max entries of bamfile used to ascertain stats like coverage, mean insert length etc., def = 500,000)	
    -s|--sam_path 		(path to samtools)
    -bp|--bp_margin 		(in deciding breakpoint start, stop this is the minimum margin added to either side of breakpoint to account for sequencing errors etc. unless the calculated margin around breakpoints is anyway larger, def = 20 --> 10 added to either side)
    -sr|--splitter_rate		(to save time, only 1 split read mapping within these many bases of others is used to improve breakpoint precision, since only 1 is really required if coming from same junction/breakpoint, def = 250)
    -sl|--slop 			(in SV determination using PE mappings, this extra slop is used in calculating overlap as necessary for whatever reason, not needed ordinarily, def = 0)
    -rd|--rd_signal 		(whether or not to use RD signal, e.g. if CN file is not available set this to 0, def = 1)
    -ro|--rec_overlap		(min reciprocal overlap needed with CN region to boost/give more weight to intermediate PE-claimed SV)
    -vr|--var_rate		(this is the max distance between the breakpoint-start loc of position-sorted clusters for which overlap should be computed, while consolidating clusters into variants;  def = dynamic upto a max of 250)
    -rr|--refresh_rate	 	(too specific, broadly -- controls how often variant list cleaned up. See ClassifyVariants.py, def = 5)	
    -dt|--disj_thresh		(number of disjoint/unique fragment mappings required in variant set to pick set as claimed SV, def = 4)
    -sp|--splitter_path		(path to split-reads bam file)
    -ss|--splitter_slop		(extra slop in determining SR overlap with existing PE variants, shouldn't be needed ordinarily, def = 0)
    -risk|--splitter_risk	(too specific, broadly if this is set may gain both TPs and FPs, def = 0. See add_SR.py)

Results:

4 bedpe files are created in the standard 6-column format (chr1, start, stop, chr2, start, stop) except for insertions.bedpe which is in the standard (Chr_paste, start, stop, Chr_cut/Chr_copy, start, end) insertion format. The kind of insertion (cut vs copy vs inverted cut etc.) is also identified in the last column.

These are stored in the sv_caller/results/text directory. These make use of PE, SR alignments and RD signal. Intermediate bedpe files using only PE mappings are stored in the above directory as well in a folder titled 'pe_results'.
