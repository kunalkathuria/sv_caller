# sv_caller
Structural Variant Caller

### SUMMARY

Tool that accepts BAM file of target as input and outputs 4 bedpe files containing deletions, insertions, inversions and tandem duplications repectively in standard bedpe format. The insertions are tagged as INS (copy-paste), INS_C (cut-paste), INS_I (inverted copy-paste), INS_C_I, INS_C_P (paste and cut locations/breakpoints are confirmed on same chromosome-- bp1 is indeed pasted location), INS_C_I_P. An Unknown.bedpe is also output containing unidentified variants.

sv_caller classifies discordant clusters into variants using PEM and then uses a uniqueness-based set-cover approach to pick those variants that are likely to be true. In between these 2 stages, it also uses split reads (to enhance breakpoint locations and add weight to existing variants) and read-depth signal (if desired, to disambiguate among purported variants).

### REQUIREMENTS

Unix-based OS with bash
python with pysam and other commonly used libraries
samtools
git
Best to have all executables (e.g. "samtools") on user path 

### INSTALLATION

Download sv_caller from GitHub and follow usage instructions below. 

### USAGE

Example call (from sv_caller/code):

./run_SVCaller.sh [options]

A good check to see if the tool is working properly is to give it as input the test BAM files in the sv_caller/data/bams/test folder and check if the resulting bedpe files found in the sv_caller/results/text folder match those contained in the sv_caller/results/test folder. So, first run the code with all default parameter values thus from the sv_caller/code folder:

./run_SVCaller.sh -a ../data/bams/test/test_bam.bam -b ../data/bams/test/test_bam.ns.bam -i "sam_path" -r ../data/bams/test/test_splitters.ns.bam

Then, make sure that the 4 files in sv_caller/results/test match the respective ones just created in sv_caller/results/text (deletions.bedpe etc.). If this is not the case, please recheck all the paths and checkout the master branch from GitHub again if necessary.

Options (also listed in command line call without arguments or with -h):

    -h                          (this menu)
    -a|--bamfile                (*REQUIRED: position-sorted bam file)
    -b|--bam_ns                 (*REQUIRED: name-sorted bam file)
    -c|--read_thresh            (max number of discordant read mappings a single PE fragment mate can have for it not to be ignored in uniqeness-based set cover, default = 20)
    -d|--match_thresh           (min threshold of ratio of a given secondary alignment's score with AS of primary alignment needed for this secondary alignment to be considered/used in set cover, def = 1)
    -e|--n_match_thresh         (same as above for number of base pair matches in alignment, instead of AS, def = .95)
    -f|--n_pct_thresh           (min ratio of matching bases in a given alignment to mean read length required for read to be used in analysis, default =  0)
    -g|--calc_thresh            (max entries of bamfile used to ascertain stats like coverage, mean insert length etc., def = 500,000)
    -i|--sam_path               (*REQUIRED: complete path to samtools command, including /samtools)
    -j|--bp_margin              (in deciding breakpoint start, stop this is the minimum margin added to either side of breakpoint to account for sequencing errors etc. unless the calculated margin around breakpoints is anyway larger, def = 20 --> 10 added to either side)
    -k|--splitter_rate          (to save time, only 1 split read mapping within these many bases of others is used to improve breakpoint precision, since only 1 is really required if coming from same junction/breakpoint, def = 250)
    -l|--slop                   (in SV determination using PE mappings, this extra slop is used in calculating overlap as necessary for whatever reason, not needed ordinarily, def = 0)
    -m|--rd_signal              (whether or not to use read depth signal, e.g. if CN file is not available set this to 0, def = 1)                                                                                                                    
    -n|--rec_overlap            (min reciprocal overlap needed with CN region to boost/give more weight to intermediate PE-claimed SV)
    -o|--var_rate               (2 clusters are compared for overlap within this much max distance between breakpoint start locations, while consolidating clusters into variants; the lower this margin is, the faster the code will be; value should be around high-percentile of breakpoint region width and so depends upon coverage, def = 250)
    -p|--refresh_rate           (too specific, broadly -- controls how often variant list cleaned up. See ClassifyVariants.py, def = 5)
    -q|--disj_thresh            (number of disjoint/unique fragment mappings required in variant set to pick set as claimed SV, def = 4)
    -r|--splitter_path          (*REQUIRED: path to split-reads bam file)
    -s|--splitter_slop          (extra slop in determining SR overlap with existing PE variants, shouldn't be needed ordinarily, def = 0)
    -t|--splitter_risk  (too specific, broadly if this is set may gain both TPs and FPs, def = 0. See add_SR.py)
    -u|--min_cluster            (the minimum number of reads supporting a cluster required for it to be considered in the analysis, def = 5)
    -v|--sig_bound              (in forming breakpoint region, the value of mean_insert_length + sig_bound*sig_insert_length - mean_RDL - (outer dist b/w min and max location of cluster read) is added to the "exact" breakpoint given by end-point of read in direction of orientation. This is usually 3 for Poisson etc., def = 2)
    -w|--sig_mult               (in forming clusters, a read's alignment location needs to fall within mean_IL + sig_mult*sigma_IL - 2*RDL of given cluster to be part of it as a necessary condition, def = 5)

#### RESULTS

4 bedpe files are created in the standard 6-column format (chr1, start, stop, chr2, start, stop) except for insertions.bedpe which is in the standard (Chr_paste, start, stop, Chr_cut/Chr_copy, start, end) insertion format. The kind of insertion (cut vs copy vs inverted cut etc.) is identified in the last column.

These files are stored in the sv_caller/results/text directory. Intermediate bedpe files using only PE mappings are stored in the above directory as well in a folder titled 'pe_results'.

### NOTES

1. Run-time depends on how many variants exist/are simulated. To speed up, consider raising the MIN_CS threshold supplied to run_PE.
2. Depending on proximity of SVs, reference copies etc., too much SLOP or BP_MARGIN can worsen results. The defaults are thus low, but feel free to change them and see if other values are better suited. In principle, higher values should not be needed.
3. There is a trade-off between the run-time and resulting benefit of run_RD_SR. It can be sped up by reducing the refresh margin.
4. For homozygous simulations, can raise MIN_CS to 25.
5. SIG_BOUND can be set to 1 or lower (breakpoint regions will now be narrower) but this can be compensated if too narrow by raising tolerance in one's comparison.
6. Var_margin should be set higher for low-coverage data. It should be set equal to the difference between max cluster length (usually mean + 3*sigma if Poisson etc.) and the length of the given breakpoint region. 
7. Min cluster size (MC) can be varied (raised from default of 5) to see if more variants are consequently claimed in the respective bedpe files. One can follow the approach of using different min CS's and picking the bedpe with max entries for each category as the final one (inversions.bedpe, deletions.bedpe etc.)
