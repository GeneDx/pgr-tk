.. pgr-tk documentation master file, created by
   sphinx-quickstart on Thu May 19 16:37:23 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to PGR-TK's documentation!
==================================


## Pangenome Research Toolkit

The PGR-TK is a python library for pangenome analysis. It includes features to (1) index human pangneome assembly and use the index for query (2) construct Minimizer Anchored Pangenome (MAP) Graph (3) Perform principal bundle decomposition for each haplotype to understand the repeat structures in pangenome.


PGR-TK implements the most computationally intensive algorithms with the Rust programming language and exposes an application interface in the Python programming language. This enables scripting and interactive analysis, e.g., with Jupyter Notebook. 

Here is a brief example using PGR-TK to generate a local MAP Graph of the MHC Class II locus::

    import pgrtk

    ref_db = pgrtk.AGCFile("/data/pgr-tk-HGRP-y1-evaluation-set-v0.agc") # lazy load an agc file of the reference without any SHIMMER index 

    sdb = pgrtk.SeqIndexDB()
    sdb.load_from_agc_index("/data/pgr-tk-HGRP-y1-evaluation-set-v0")

    ref_file_name, roi_chr, roi_b, roi_e = 'hg19_tagged.fa', "chr6_hg19", 32130918, 32959917
    padding = 0

    #get a segment of a reference
    roi_seq = ref_db.get_sub_seq(ref_file_name, roi_chr, roi_b-padding, roi_e+padding)


    # get the hits from the pangenome reference
    aln_range = pgrtk.query_sdb(sdb, roi_seq, merge_range_tol=200000)

    # collect the target sequences from the hits
    seq_list = []
    i = 0
    for k in list(aln_range.keys()):
        ctg_name, source, _ = seq_info[k]
        seq_id = k
        rgns = aln_range[k].copy()

        rgns = pgrtk.merge_regions(rgns,tol=1000)

        for rgn in rgns:
            b, e, length, orientation, aln = rgn

            seq =  sdb.get_sub_seq(source, ctg_name, b, e)
            if orientation == 1:
                seq = pgrtk.rc_byte_seq(seq)

            seq_list.append(("{}_{}_{}_{}".format(ctg_name, b, e, orientation), seq))
            
            i += 1

    # Create a shimmer index database with smaller window (denser shimmers)
    new_sdb = pgrtk.SeqIndexDB() 
    new_sdb.load_from_seq_list(seq_list, w=80, k=56, r=12, min_span=18)

    new_sdb.generate_mapg_gfa(0, "/results/HLA-ClassII.gfa")


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   pgrtk

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
