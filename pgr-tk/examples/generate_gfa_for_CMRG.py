import pgrlite
import os
import networkx as nx
from networkx.drawing import nx_pydot
from collections import Counter

def generate_gfa(cmrg_regions, gene_name, pg_db, out_dir):
    gene_seq = cmrg_regions[gene_name][1]
    aln_range0 = pgrlite.query_sdb(pg_db, gene_seq, merge_range_tol=len(gene_seq) * 0.25)
    print("The number of hits for {} is {}".format((gene_name), len(aln_range0)))
    count = 0
    for sid, rgns in aln_range0.items():
        count += len(rgns) 
    
    print("The total aligned regions {} is {}".format(gene_name, count))
    seq_info = pg_db.seq_info.copy()
    with open(os.path.join(out_dir, f"{gene_name}_hit.txt"), "w") as f:
        print("#source", "ctg", "len", "n_hit", sep="\t", file = f)
        for k in aln_range0:
            if len(aln_range0[k]) >= 1:
                ctg, src, len_ = seq_info[k]
                print(src, ctg, len_, len(aln_range0[k]), sep="\t", file = f)
    

    rgn_lengths = []
    with open(os.path.join(out_dir, f"{gene_name}_hit_range.txt"), "w") as f:
        print("#sourc", "ctg", "len", "t_rgn_start", "t_rgn_end", "t_rgn_len", sep="\t", file = f)
    
        for k in list(aln_range0.keys()):
            b, e = aln_range0[k][0][0:2]
            if e-b < len(gene_seq) * 0.25:
                continue
            ctg, src, len_ = seq_info[k]
            print(src, ctg, len_, b, e, e-b, sep="\t", file = f )
            rgn_lengths.append(e-b)
    
    with open(os.path.join(out_dir, f"{gene_name}_ht_copy_count.txt"), "w") as f:
        n_copy = {}
        for k in list(aln_range0.keys()):
            b, e = aln_range0[k][0][0:2]
            if e-b < len(gene_seq) * 0.25:
                continue
            n_copy[k] = len(aln_range0[k])
        copy_count = Counter(n_copy.values())
        for nc, nh in copy_count.items():
            print("{}\tnumber_of_copy: {}\tnumber_of_haplotype_contig: {}".format(gene_name, nc, nh), file = f)

    seq_list = []
    i = 0
    for k in list(aln_range0.keys()):
        ctg_name, source, _ = seq_info[k]
        seq_id = k
        rgns = aln_range0[k].copy()
        rgns = pgrlite.merge_regions(rgns, tol=int(len(gene_seq)*0.25))

        for rgn in rgns:
            b, e, length, orientation, aln = rgn
            if length < len(gene_seq)*0.25:
                continue
            seq = pg_db.get_sub_seq(source, ctg_name, b, e)
            if orientation == 1:
                seq = pgrlite.rc_byte_seq(seq)
            seq_list.append((i, "{}_{}_{}_{}".format(ctg_name, b, e, orientation), seq))
            i += 1

    with open( os.path.join(out_dir, f"{gene_name}.fa"), "w") as f:
        for sid, name, seq in seq_list:
            print(">{} {}".format(name, sid), file = f)
            print(pgrlite.u8_to_string(seq), file = f)

    new_sdb = pgrlite.SeqIndexDB() 
    new_sdb.load_from_seq_list(seq_list, w=48, k=48, r=1, min_span=24)
    new_sdb.generate_smp_gfa(0, os.path.join(out_dir, f"{gene_name}_48_48_1_24.gfa"))
    new_sdb.write_midx_to_text_file(os.path.join(out_dir, f"{gene_name}_48_48_1_24.midx"))

    links = new_sdb.get_smp_adj_list(0)
    link_count = Counter([(_[1],_[2]) for _ in links])
    G = nx.DiGraph()
    for sid, v, w in links:
        if sid == 0:
            continue
            
        penwidth = link_count[(v,w)] * 0.01
        weight = penwidth
        G.add_edge(tuple(v[:2]), tuple(w[:2]), weight= weight, penwidth=penwidth)

    nx_pydot.write_dot(G, os.path.join(out_dir, "{}_48_48_1_24.dot".format(gene_name)))
    nx.write_gexf(G, os.path.join(out_dir, "{}_48_48_1_24.gexf".format(gene_name)))


    new_sdb.load_from_seq_list(seq_list, w=48, k=48, r=8, min_span=24)
    new_sdb.generate_smp_gfa(0, os.path.join(out_dir, f"{gene_name}_48_48_8_24.gfa"))
    new_sdb.write_midx_to_text_file(os.path.join(out_dir, f"{gene_name}_48_48_8_24.midx"))
    
    links = new_sdb.get_smp_adj_list(0)
    link_count = Counter([(_[1],_[2]) for _ in links])
    G = nx.DiGraph()
    for sid, v, w in links:
        if sid == 0:
            continue
            
        penwidth = link_count[(v,w)] * 0.01
        weight = penwidth
        G.add_edge(tuple(v[:2]), tuple(w[:2]), weight= weight, penwidth=penwidth)
    
    nx_pydot.write_dot(G, os.path.join(out_dir, "{}_48_48_8_24.dot".format(gene_name)))
    nx.write_gexf(G, os.path.join(out_dir, "{}_48_48_8_24.gexf".format(gene_name)))

    new_sdb.generate_smp_gfa(0, os.path.join(out_dir, f"{gene_name}_48_48_8_24.gfa"))
    new_sdb.write_midx_to_text_file(os.path.join(out_dir, f"{gene_name}_48_48_8_24.midx"))


if __name__ == "__main__":
    ref_db =pgrlite.AGCFile("/data/HPRC-y1-rebuild-04252022/hg19.agc")
    pb_db = pgrlite.SeqIndexDB()
    pb_db.load_from_agc_index("/data/HPRC-y1-rebuild-04252022")
    CMRG_coordinates = {}
    padding = 20000
    with open("/data/HG002_GRCh37_CMRG_coordinates.bed") as f:
        for r in f:
            r = r.strip().split("\t")
            CMRG_coordinates[r[3]]=("chr{}".format(r[0]), int(r[1])-padding, int(r[2])+padding)

    CMRG_hg19_seq = {}
    for g, c in CMRG_coordinates.items():
        seq = ref_db.get_sub_seq('hg19.fasta', c[0], c[1], c[2])
        CMRG_hg19_seq[g] = (c, seq)

    for g_name in CMRG_hg19_seq:
        print("analyzing {}".format(g_name))
        generate_gfa(CMRG_hg19_seq, g_name, pb_db, "/scratch/GFA_files")
