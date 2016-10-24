# sv_caller
Structural Variant Caller

Example call:

./run_SVCaller.sh [options]

Options listed and explained in command line call without arguments or with -h.

Results:

4 bedpe files are created in the standard 6-column format (chr1, start, stop, chr2, start, stop) except for insertions.bedpe which is in the standard (Chr_paste, start, stop, Chr_cut/Chr_copy, start, end) insertion format. The kind of insertion (cut vs copy vs inverted cut etc.) is also identified in the last column.

These are stored in the sv_caller/results/text directory. These make use of PE, SR alignments and RD signal. Intermediate bedpe files using only PE mappings are stored in the above directory as well in a folder titled 'pe_results'.

Notes:

1. Run-time depends on how many variants exist/are simulated. To speed up, consider raising the MIN_CS threshold supplied to run_PE.
2. Depending on proximity of SVs, reference copies etc., too much SLOP or BP_MARGIN can worsen results. The defaults are thus low, but feel free to change them and see if other values are better suited. In principle, higher values should not be needed.
3. There is a trade-off between the run-time and resulting benefit of run_RD_SR. It can be sped up by reducing the refresh margin.
4. For homozygous simulations, can raise MIN_CS to 25.
5. SIG_BOUND can be set to 1 or lower (breakpoint regions will now be narrower) but this can be compensated if too narrow by raising tolerance in one's comparison.
6. Var_margin should be set higher for low-coverage data. It should be set equal to the difference between max cluster length (usually mean + 3*sigma if Poisson etc.) and the length of the given breakpoint region. 
7. Min cluster size (MC) can be varied (raised from default of 5) to see if more variants are consequently claimed in the respective bedpe files. One can follow the approach of using different min CS's and picking the bedpe with max entries for each category as the final one (inversions.bedpe, deletions.bedpe etc.)
