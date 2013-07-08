This note is about a simple implementation of the algorithm developed by Gene
Myers in 1986 (An O(ND) Difference Algorithm and Its Variations, Eugene W.
Myers, Algorithmica Vol 1, pp. 251--266). A PDF copy of the paper can be found
at http://neil.fraser.name/software/diff_match_patch/myers.pdf

While this paper has been cited more than 600 times according to Google
Scholar, I rarely read about this method in the recent literatures in
bioinformatics. This method is very efficent for sequence alignments of very
similar sequences.  It will be really useful for genome asssembly where the
sequence reads are typically very close to the assembly and alignment are used
for consensus rather than evolution comparison.  Unlikely the
Smith-Waterman-like dynamic programming approach for sequence alignment, this
method only concerns about editing distance rather than maximuming the
"alignment scores".  The greedy nature makes it really effiecient if the
difference between two sequences are small.  This notebook shows that is
acutally working quite well even the differences are big.  However, for using
it as a general purpose alignment tool for various bioinformatics tasks, some
pre- or post- process might be needed.  It will be interesting to do some study
to see how one can integrate such method into routine genome assembly workflow.


--Jason Chin, July 5, 2013
