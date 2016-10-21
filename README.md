# sv_caller
Structural Variant Caller

Example call:

./run_SVCaller_d.sh (has all default values of parameters and files as noted/stored once in headers)

./run_SVCaller.sh [options]

Options:

    -b|--bamfile 		(position-sorted bam file)
    -bn|--bam_ns 		(name-sorted bam file)
    -r|--read_thresh 		(max number of discordant read mappings a single PE fragment mate can have for it not to be ignored in uniqeness-based set cover, default = 20)
    -m|--match_thresh 		(min threshold of ratio of a given secondary alignment's score with AS of primary alignment needed for this secondary alignment to be considered/used in set cover, def = 1)
    -n|--n_match_thresh 	(same as above for number of base pair matches in alignment, instead of AS, def = .95) 
    -np|--n_pct_thresh 		(min ratio of matching bases in a given alignment to mean read length required for read to be used in analysis, default =  0)
    -ct|--calc_thresh 		(max entries of bamfile used to ascertain stats like coverage, mean insert length etc., def = 500,000)	
    -s|--sam_path 		(path to samtools)
    -bp|--bp_margin 		(in deciding breakpoint start, stop this is the minimum margin added to either side of breakpoint to account for sequencing errors etc. unless the calculated margin around breakpoints is anyway larger, def = 20 --> 10 added to either side)
    -sr|--splitter_rate		(to save time, only 1 split read mapping within these many bases of others is used to improve breakpoint precision, since only 1 is really required if coming from same junction/breakpoint, def = 250)
    -sl|--slop 			(in SV determination using PE mappings, this extra slop is used in calculating overlap as necessary for whatever reason, not needed ordinarily, def = 0)
    -rd|--rd_signal 		(whether or not to use read depth signal, e.g. if CN file is not available set this to 0, def = 1)
    -ro|--rec_overlap		(min reciprocal overlap needed with CN region to boost/give more weight to intermediate PE-claimed SV)
    -vr|--var_rate		(2 clusters are compared for overlap within this much max distance between breakpoint start locations, while consolidating clusters into variants; the lower this margin is, the faster the code will be;  def = dynamically calcuated with max threshold of 250)
    -rr|--refresh_rate	 	(too specific, broadly -- controls how often variant list cleaned up. See ClassifyVariants.py, def = 5)	
    -dt|--disj_thresh		(number of disjoint/unique fragment mappings required in variant set to pick set as claimed SV, def = 4)
    -sp|--splitter_path		(path to split-reads bam file)
    -ss|--splitter_slop		(extra slop in determining SR overlap with existing PE variants, shouldn't be needed ordinarily, def = 0)
    -risk|--splitter_risk	(too specific, broadly if this is set may gain both TPs and FPs, def = 0. See add_SR.py)
    -mc|--min_cluster		(the minimum number of reads supporting a cluster required for it to be considered in the analysis, def = 5) 	
    -sb|--sig_bound		(in forming breakpoint region, the value of mean_insert_length + sig_bound*sig_insert_length - mean_RDL - (outer dist b/w min and max location of cluster read) is added to the "exact" breakpoint given by end-point of read in direction of orientation. This is usually 3 for Poisson etc., def = 2)
    -sm||--sig_mult		(in forming clusters, a read's alignment location needs to fall within mean_IL + sig_mult*sigma_IL - 2*RDL of given cluster to be part of it as a necessary condition, def = 5)

Results:

4 bedpe files are created in the standard 6-column format (chr1, start, stop, chr2, start, stop) except for insertions.bedpe which is in the standard (Chr_paste, start, stop, Chr_cut/Chr_copy, start, end) insertion format. The kind of insertion (cut vs copy vs inverted cut etc.) is also identified in the last column.

These are stored in the sv_caller/results/text directory. These make use of PE, SR alignments and RD signal. Intermediate bedpe files using only PE mappings are stored in the above directory as well in a folder titled 'pe_results'.
