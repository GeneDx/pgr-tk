import pgrtk
import os, sys


def filter_aln(aln_segs):
    """
    ensure both target / query are strictly increasing
    """
 
    last_ts = aln_segs[0][1][0]
    last_te = aln_segs[0][1][1]

    last_qs = aln_segs[0][0][0]
    last_qe = aln_segs[0][0][1]

    
    rtn = [ ((last_ts, last_te), (last_qs, last_qe)) ]
    
    for seg in aln_segs:
      
        if seg[1][1] < seg[1][0]: continue
        if seg[0][-1] != seg[1][-1]: continue

        if seg[1][0] >= last_te: 
       
            last_ts = last_te
            last_te = seg[1][1]
 
            last_qs = last_qe
            last_qe = seg[0][1]
          
            if last_ts == last_te:
                continue
            
            
            rtn.append( ((last_ts, last_te), (last_qs, last_qe)) )
    return rtn
        

    
def filter_aln_rev(aln_segs):
    """
    ensure both target / query are strictly increasing
    """
    aln_segs = aln_segs.copy() 
    aln_segs.reverse()
    last_ts = aln_segs[0][1][0]
    last_te = aln_segs[0][1][1]

    last_qs = aln_segs[0][0][0]
    last_qe = aln_segs[0][0][1]

    
    rtn = [ ((last_ts, last_te), (last_qs, last_qe)) ]
    
    for seg in aln_segs:
      
        if seg[1][1] < seg[1][0]: continue
        if seg[0][-1] == seg[1][-1]: continue

        if seg[1][0] >= last_te: 
       
            last_ts = last_te
            last_te = seg[1][1]
 
            last_qe = last_qs
            last_qs = seg[0][0]
            
          
            if last_ts == last_te:
                continue
            
            
            rtn.append( ((last_ts, last_te), (last_qs, last_qe)) )
    return rtn 
    
def seq_align_to_sdb(seq_db, seq1):

    query_res = pgrtk.query_sdb(seq_db, seq1, 
                       merge_range_tol=0, 
                       gap_penalty_factor=0.001, 
                       max_query_count=1, 
                       max_target_count=1)
    
    _, kmer_size, _, _, _ = seq_db.get_shmmr_spec() 
    rtn = []

    for sid, alns in query_res.items():
        # print("#sid, hits:", sid, len(alns))
        
        ref_seq = seq_db.get_seq_by_id(sid)
        
        for aln in alns:
            ts, te, tl, orientation = aln[:-1]
            # print(ts, te, tl, orientation)
            aln = aln[-1]
            if orientation == 0 :
                filter_alignments = filter_aln(aln)
            else:
                filter_alignments = filter_aln_rev(aln)
            # print("# anchors: ", len(aln), len(filter_aln(aln)), len(filter_aln_rev(aln)))
            
            for seg in filter_alignments:
                
                last_ts, last_te = seg[0][:2]
                last_qs, last_qe = seg[1][:2]   

                last_ts -= kmer_size
                #last_te -= kmer_size
             
                s0str = pgrtk.u8_to_string(ref_seq[last_ts: last_te])
                if orientation == 0:
                    last_qs -= kmer_size
                    s1str =  pgrtk.u8_to_string(seq1[last_qs:last_qe])
                else:
                    last_qs -= kmer_size
                    s1str =  pgrtk.rc(pgrtk.u8_to_string(seq1[last_qs:last_qe]))
             
                if s0str[:16] != s1str[:16] or s0str[-16:] != s1str[-16:]:
                    print("XXXX1 {} :\n{}\n{}\n".format(orientation, s0str[:56],  s1str[:56]))
                    print("XXXX2 {} :\n{}\n{}\n".format(orientation, s0str[-56:],  s1str[-56:]))
                    diff = None
                elif min(len(s0str),len(s1str)) == 0 or abs(len(s0str)-len(s1str)) > 256:
                    diff = None
                else:
                    diff = pgrtk.get_variant_segments(s0str, s1str, max_wf_length=min(64, len(s0str), len(s1str)), max_diff_percent=1)

                if diff is not None:
                    if len(diff[0]) > 0:
                        for d in diff[0]:
                            rtn.append( ((sid, last_ts, last_te), (last_qs, last_qe), 
                                        (d[0] + last_ts, d[1] + last_qs, d[2], d[3], d[4]), orientation) )
                    else:
                        rtn.append( ((sid, last_ts, last_te), (last_qs, last_qe), 'ALL', orientation ) )
                elif diff is None:
                    rtn.append( ((sid, last_ts, last_te), (last_qs, last_qe), 'NULL', orientation ) )
    return rtn

def main(sdb_prefix, query_seq_fasta_path, out_prefix = "out"):
    target_sdb = pgrtk.SeqIndexDB()
    target_sdb.load_from_frg_index(sdb_prefix)
    query_sdb = pgrtk.SeqIndexDB()
    query_sdb.load_from_fastx(query_seq_fasta_path)

    target_sinfo = target_sdb.seq_info.copy() 
    sinfo = query_sdb.seq_info.copy()
    variant_file = open(out_prefix+".variants", "w")
    sv_candidate_file = open(out_prefix+".sv_candidate", "w")
    all_match_file = open(out_prefix+".all_match", "w")
    for sid in sinfo:
        ctg, src, length = sinfo[sid]
        query_seq = query_sdb.get_seq_by_id(sid)
        variants = seq_align_to_sdb(target_sdb, query_seq)
        for variant in variants: 
            t_sid, ts, te = variant[0]
            qs, qe = variant[1]
            t_ctg, _, _ = target_sinfo[t_sid]
            rec = variant[2]
            if rec in ['ALL', 'NULL']:
                print(t_ctg, ts, te, ctg, qs, qe, variant[2], variant[3], sep="\t", file=all_match_file)
            else:
                print(t_ctg, ts, te, ctg, qs, qe, rec[0], variant[3], sep="\t", file=all_match_file)
                print(t_ctg, rec[0], rec[2], rec[3], rec[4], ctg, sep="\t", file=variant_file)
            if rec == "NULL":
                print(t_ctg, variant[0][1], variant[0][2], ctg, variant[1][0], variant[1][1], sep="\t", file=sv_candidate_file)
    variant_file.close()
    sv_candidate_file.close()

            

    

if __name__ == "__main__":

    sdb_prefix = sys.argv[1]
    query_seq_fasta_path = sys.argv[2]
    prefix = sys.argv[3]
    main(sdb_prefix, query_seq_fasta_path, prefix)
