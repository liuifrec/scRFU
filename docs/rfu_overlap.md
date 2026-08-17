# RFU Repertoire Overlap

`scrfu.tl.rfu_overlap()` accepts an existing sample-by-RFU matrix or constructs
one from row-level assignments when `sample_key` is supplied. Reusing a
`RFUPseudobulkResult` avoids silently changing pseudobulk semantics.

For presence sets A and B, Jaccard is `|A∩B|/|A∪B|`, Sørensen-Dice is
`2|A∩B|/(|A|+|B|)`, and overlap coefficient is
`|A∩B|/min(|A|,|B|)`. For non-negative abundance vectors x and y, cosine is
`x·y/(||x|| ||y||)`, weighted Jaccard is
`sum(min(x,y))/sum(max(x,y))`, and Bray-Curtis dissimilarity is
`sum(|x-y|)/sum(x+y)`; its similarity is one minus that value. Optional
Jensen-Shannon distance is the square root of the mean KL divergence to the
midpoint distribution, using natural logs.

The structured result labels the metric direction as `similarity` or
`distance`. A denominator-zero comparison is undefined and remains `NaN`; it
is never silently replaced by zero. Negative or non-finite matrices are
rejected, and `min_abundance` is applied before pairwise calculation.
