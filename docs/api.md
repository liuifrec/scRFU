# API Reference

## scrfu.tl.call_rfu

Parameters:

-   adata : AnnData
-   chain : str (default 'TRB')
-   airr_key : str (default 'airr')
-   out_key : str
-   rscript_path : path to RFU_call.R
-   atlas_path : path to km5000 centroids
-   extra_r\_args : optional list of extra arguments

Returns: - pandas.DataFrame with RFU assignments.
