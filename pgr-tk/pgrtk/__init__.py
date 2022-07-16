# -*- coding: utf-8 -*-
"""This module is used to extract and compare sequences in a set of pan genome assemblies.
It includes a number modules for access the sequence data and query the sequence index.

Example
-------

This shows a simple example to query the pangenomics database::

    import pgrtk
  
    ## The AGCFile class is used to access the sequence data stored in a AGC file. 
    ref_db =pgrtk.AGCFile("hg19.agc")

    ## Load a pre-build index and sequence data from all humane genome assemblies 
    ## of the HPRC year one release. 
    sdb = pgrtk.SeqIndexDB()
    sdb.load_from_agc_index("HPRC-y1-rebuild-04252022")

    ## Extract a sequence from the hg19 AGC file.
    gene_seq = ref_db.get_sub_seq('hg19.fasta', 'chr6', 
                                  160952514, 161087407)

    ## find hits in the pangenomic reference
    alignment_ranges = pgrtk.query_sdb(sdb, gene_seq, 
                                         merge_range_tol=100000)

"""

import pgrtk
import numpy as np
from .pgrtk import *

__version__ = pgrtk.pgr_lib_version()

byte_rc_map = dict(zip([ord(c) for c in "ACGTNnacgt"],
                   [ord(c) for c in "TGCANntgca"]))


def rc_byte_seq(seq):
    """ Reverse complement a sequence as a list of bytes.

    Parameters
    ----------
    seq : list of bytes
        ascii code of the DNA sequence 

    Returns
    -------
    list of bytes 
        the list of bytes of the reverse complement DNA sequence

    """
    seq = [byte_rc_map[_] for _ in seq[::-1]]
    return seq


rc_map = dict(zip("ACGTNnactg", "TGCANntgca"))


def rc(seq):
    """ Reverse complement a sequence as a Python String.

    Parameters
    ----------
    seq : string
        a DNA sequence as a Python String 

    Returns
    -------
    string
        the reverse complement DNA sequence as a Python String

    """
    seq = "".join([rc_map[_] for _ in seq[::-1]])
    return seq


def string_to_u8(s):
    """ Convert a Python String to a list of bytes.

    Parameters
    ----------
    s : string
        a Python String of a DNA sequence

    Returns
    -------
    list of bytes
        a list of bytes representing the DNA sequence

    """
    return list(s.encode("utf-8"))


def u8_to_string(u8):
    """ Convert DNA sequene in a list of bytes to a Python String.

    Parameters
    ----------
    u8 : list of bytes
        a list of bytes representing the DNA sequence

    Returns
    -------
    string
        a Python String of a DNA sequence

    """
    return bytes(u8).decode("utf-8")


def query_sdb(seq_index_db, query_seq,
              gap_penality_factor=0.25,
              merge_range_tol=12,
              max_count=128,
              max_query_count=128,
              max_target_count=128,
              max_aln_span=8):
    """ Query a sequence index database for a query sequence. 

    Parameters
    ----------
    seq_index_db : SeqIndexDB object
        a sequence index database object

    query_seq : list of bytes
        a list of bytes representing the DNA sequence

    gap_penality_factor : float
        the gap penalty factor used in sparse dyanmic programming for finding the hits

    merge_range_tol : int
        a parameter used to merge the alignment ranges

    max_count : int
        only use the shimmer pairs that less than the ``max_count`` for sparse dynamic programming

    max_query_count : int
        only use the shimmer pairs that less than the ``max_count`` in the query sequence for sparse dynamic programming

    max_query_count : int
        only use the shimmer pairs that less than the ``max_count`` in the target sequence for sparse dynamic programming

    max_aln_span : int
        the size of span used in the sparse dynamic alignment for finding the hits


    Returns
    -------
    dict
        - a python dictionary with the key as the target sequence id and the value as a list of alignment ranges

        - each alignment ranges is a list of tuples, each tuple is (``start``, ``end``, ``length``, 
          ``orientation``, ``aln_records``)

        - the ``aln_records`` is a list of tuples of 
          (``target_sequence_id``, (``score``, ``list_of_the_hit_pairs``)), where 
          the ``list_of_the_hit_pairs`` is a list of tuples of 
          ((``query_start``, ``query_end``, ``query_orientation``), 
          (``target_start``, ``target_end``, ``target_orientation``))

    """
    r = seq_index_db.query_fragment_to_hps(
        query_seq,
        gap_penality_factor,
        max_count,
        max_query_count,
        max_target_count,
        max_aln_span)

    sid_to_alns = {}
    for (sid, alns) in r:
        aln_lens = []
        f_count = 0
        r_count = 0
        for s, aln in alns:
            if len(aln) > 2:
                aln_lens.append(len(aln))
                sid_to_alns.setdefault(sid, [])
                for hp in aln:
                    if hp[0][2] == hp[1][2]:
                        f_count += 1
                    else:
                        r_count += 1
                orientation = 0 if f_count > r_count else 1
                sid_to_alns[sid].append((aln, orientation))

    aln_range = {}
    for sid, alns in sid_to_alns.items():
        for aln, orientation in alns:
            target_coor = [(_[1][0], _[1][1]) for _ in aln]
            target_coor.sort()
            bgn = min(target_coor[0])
            end = max(target_coor[-1])
            aln_range.setdefault(sid, [])
            aln_range[sid].append((bgn, end, end-bgn, orientation, aln))

    if merge_range_tol > 0:
        for sid, rgns in aln_range.items():
            aln_range[sid] = merge_regions(
                rgns, tol=merge_range_tol)

    return aln_range


def map_intervals_in_sdb(seq_index_db, interval, query_seq,
              gap_penality_factor=0.001,
              merge_range_tol=100,
              max_count=32,
              max_query_count=32,
              max_target_count=32,
              max_aln_span=8):
    """
    TODO: Document
    
    """

    assert(len(interval) == 2)
    
    pos_map = seq_index_db.map_positions_in_seq(interval, query_seq, 0.001, 32, 32, 32, 100)
    
    seqid_to_positions = {}
    
    for res in pos_map:

        pos = res[0]
        sid, tpos, orientation = res[1]
        
        seqid_to_positions.setdefault(sid, {})
        seqid_to_positions[sid].setdefault(pos, [])
        seqid_to_positions[sid][pos].append((tpos, orientation))
    
    rtn = {}
    for sid in seqid_to_positions:
        #print(sid, seqid_to_positions[sid])
        if interval[0] in seqid_to_positions[sid] and interval[1] in seqid_to_positions[sid]:
            left_p = seqid_to_positions[sid][interval[0]]
            right_p = seqid_to_positions[sid][interval[1]]
            if len(left_p) != 1:
                continue
            if len(right_p) != 1:
                continue
            left_p, left_o = left_p[0]
            right_p, right_o = right_p[0]
            if left_o != right_o:
                continue
          
            rtn[sid] = (left_o, left_p, right_p)
    return rtn

def merge_regions(rgns, tol=1000):
    """ Take a list of ranges and merge them if two regions are within ``tol``.
    Parameters
    ----------
    rgns : list of tuples
        a list of tuples of (``start``, ``end``, ``length``, ``orientation``, ...)

    Returns
    -------
    list of tuples
        a list of tuples of (``start``, ``end``, ``length``, ``orientation``, ...)

    """
    # rgns is a list of (bgn, end, len, orientation)

    rgns.sort()
    frgns = [r for r in rgns if r[3] == 0]
    rrgns = [r for r in rgns if r[3] == 1]
    fwd_rgns = []
    last = None
    for r in frgns:
        r = list(r)
        if last is None:
            last = r[1]
            fwd_rgns.append(r)
            continue

        if r[1] < fwd_rgns[-1][1]:
            continue

        if r[0] - last < tol:  # merge
            fwd_rgns[-1][1] = r[1]
            fwd_rgns[-1][2] += r[2]
            fwd_rgns[-1][4] += r[4]
        else:
            fwd_rgns.append(r)
        last = fwd_rgns[-1][1]

    rev_rgns = []
    last = None
    for r in rrgns:
        r = list(r)
        if last is None:
            last = r[1]
            rev_rgns.append(r)
            continue

        if r[1] < rev_rgns[-1][1]:
            continue

        if r[0] - last < tol:  # merge
            rev_rgns[-1][1] = r[1]
            rev_rgns[-1][2] += r[2]
            rev_rgns[-1][4] += r[4]
        else:
            rev_rgns.append(r)

        last = rev_rgns[-1][1]
    return fwd_rgns + rev_rgns


def get_variant_calls(aln_segs, ref_bgn, ctg_bgn, rs0, cs0, strand):
    """ Generate a variant call internal representation from the alignment segments.

    Parameters
    ----------
    aln_segs : list of tuples
        -  a list of tuples of "alignment segments" generate by ``pgrtk.pgrtk.get_aln_segements()``
        -  the "alignment segments" are a list of ``(ref_loc: SeqLocus, tgt_loc: SeqLocus, align_type: AlnSegType)``. The data structures
           is defined as following Rust structs::

                pub struct SeqLocus {
                    pub id: u32,
                    pub bgn: u32,
                    pub len: u32,
                }

                pub enum AlnSegType {
                    Match,
                    Mismatch,
                    Insertion,
                    Deletion,
                    Unspecified,
                }

                pub struct AlnSegment {
                    pub ref_loc: SeqLocus,
                    pub tgt_loc: SeqLocus,
                    pub t: AlnSegType,
                }

    ref_bgn : int
        the reference sequence start position

    ctg_bgn : int
        the contig start position

    rs0 : string
        the reference sequence

    cs0 : string
        the contig sequence

    strand : int
        the contig strand


    Returns
    -------
    dict
        a dictionary mapping the key (ref_id, reference_position) to a set of variant 
        calls in the form of a dictionary mapping from (target_location, strand) to 
        a variant call record.
    """
    variant_calls = {}
    for s in aln_segs:
        ref_id = s.ref_loc[0]
        if s.t != ord('M'):
            if s.t == ord('X'):
                key = (ref_id, s.ref_loc[1]+ref_bgn+1)
                ref_bases = rs0[s.ref_loc[1]:s.ref_loc[1]+s.ref_loc[2]]
                alt_bases = cs0[s.tgt_loc[1]:s.tgt_loc[1]+s.tgt_loc[2]]

            if s.t == ord('I'):
                p0 = s.ref_loc[1]
                p1 = s.tgt_loc[1]

                while 1:
                    if rs0[p0-1] == cs0[p1+s.tgt_loc[2]-1] and rs0[p0-2] == cs0[p1-2]:
                        p0 = p0 - 1
                        p1 = p1 - 1
                    else:
                        break

                key = (ref_id, p0+ref_bgn)
                ref_bases = rs0[p0-1:p0+s.ref_loc[2]]
                alt_bases = cs0[p1-1:p1+s.tgt_loc[2]]

            if s.t == ord('D'):

                p0 = s.ref_loc[1]
                p1 = s.tgt_loc[1]

                while 1 and p0 > 0 and p1 > 0:
                    if rs0[p0+s.ref_loc[2]-1] == cs0[p1-1] and rs0[p0-2] == cs0[p1-2]:
                        p0 = p0 - 1
                        p1 = p1 - 1
                    else:
                        break

                key = (ref_id, p0+ref_bgn)
                ref_bases = rs0[p0-1:p0+s.ref_loc[2]]
                alt_bases = cs0[p1-1:p1+s.tgt_loc[2]]

            value = (chr(s.t), ref_bases, alt_bases,
                     (key[0], key[1], s.ref_loc[2]),
                     (s.tgt_loc[0], ctg_bgn+s.tgt_loc[1], s.tgt_loc[2]), strand)
            variant_calls.setdefault(key, {})
            variant_calls[key][(s.tgt_loc[0], strand)] = value

    return variant_calls


def output_variants_to_vcf_records(variant_calls, ref_name):
    """ Convert the variant calls to VCF records.

    Parameters
    ----------
    variant_calls : dict
        the variant calls generated by ``get_variant_calls()``

    ref_name : string
        reference sequence name

    Returns
    -------
    list
        list of VCF records
    """
    keys = sorted(list(variant_calls.keys()))
    vcf_recs = []
    for k in keys:
        v = variant_calls[k]
        ref_id = k[0]
        gt_set = set()
        gt = ["."]
        gt_idx = 0
        variants = []
        for kk in v:
            gt_idx += 1
            variants.append((gt_idx, v[kk][:3], ref_id))

        count = {}
        for v in variants:
            count.setdefault(v[0], [])
            count[v[0]].append(v[1])

        ht = [(str(x[2]), str(x[0])) for x in variants]
        ht.sort()
        ht = list(zip(*ht))

        for kk in sorted(count.keys()):
            ref_base = count[kk][0][1]
            alt_base = count[kk][0][2]
            # print(count)
            if ref_base == alt_base:
                continue

            vcf_recs.append((ref_name,  "{}".format(k[1]), ".", ref_base, alt_base,
                             "30", ".", ".", "GT:AD", "./1:0,1:"))

    return vcf_recs


def compute_graph_diffusion_entropy(gfa_fn, max_nodes = 6000):
    """ Give a GFA file name, compute an entropy by a simple diffusion model on the grap
        and generate the list of the final diffusion weight for each node
    
    Parameters
    ----------
    gfa_fn : string
        a gfa filename

    Returns
    -------
    tuple
        ``(entropy, list_of_diffusion_weight)``

        list_of_diffusion_weight = ``[(node_id, weight), ...]`` 
    """
    adj_list = {}

    with open(gfa_fn) as f:
        for r in f:
            r = r.strip().split("\t")
            if r[0] != "L":
                continue
            n1 = int(r[1])
            n2 = int(r[3])
            weight = None
            for f in r[6:]:
                f = f.split(":")
                if f[0] == "SC":
                    weight = int(f[2])
            if weight == None:
                weight = 1
            adj_list.setdefault( n1, [] )
            adj_list[ n1 ]. append((n2, weight))
            adj_list.setdefault( n2, [] )
            adj_list[ n2 ]. append((n1, weight))
    
    n_node = len(adj_list)
    if n_node > max_nodes:
        ## TODO: proper message to handle big graph
        return None

    adj_matrix = np.zeros( (n_node, n_node), dtype=np.float32 )
    for v, ws in adj_list.items():
        for w, weight in ws:
            adj_matrix[v][w] = weight
            
    n_adj_matrix = adj_matrix / np.sum(adj_matrix, axis=1)
    weights = []
    one = np.ones(n_node, dtype=np.float32)/n_node
    yy = one.copy()

    for i in range(n_node):
        yy = np.inner(n_adj_matrix, yy)

    entropy = -np.sum(yy * np.log2(yy))
    weight_list = list(enumerate(yy*n_node))
  
    return (entropy, weight_list)
