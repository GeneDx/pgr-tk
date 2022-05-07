from .pgrlite import *

__version__ = pgr_lib_version()

byte_rc_map = dict(zip([ord(c) for c in "ACGTNnacgt"],
                   [ord(c) for c in "TGCANntgca"]))


def rc_byte_seq(seq):
    seq = [byte_rc_map[_] for _ in seq[::-1]]
    return seq


rc_map = dict(zip("ACGTNnactg", "TGCANntgca"))


def rc(seq):
    seq = "".join([dict(zip("ACGT", "TGCA"))[_] for _ in seq[::-1]])
    return seq


def string_to_u8(s):
    return list(s.encode("utf-8"))


def u8_to_string(u8):
    return bytes(u8).decode("utf-8")


class SeqShmmrIdxDB(object):

    def __init__(self, agc_prefix=None, fastx_file_spec=None, seq_list=None, backed_by="AGC"):
        if agc_prefix is not None and backed_by == "AGC":
            self.back_by = "AGC"
            self.agc_prefix = agc_prefix
            self.shmap = ShmmrFragMap()
            self.shmap.load_from_mdb(f"{agc_prefix}.mdb")
            self.agcfile = AGCFile(f"{agc_prefix}.agc")
            self.seq_index = {}
            with open(f"{agc_prefix}.midx") as f:
                for r in f:
                    # r = (id, seq_len, seq_id, seq_name)
                    r = r.strip().split("\t")
                    self.seq_index[int(r[0])] = int(r[1]), r[2], r[3]

        elif seq_list is not None and backed_by == "SeqList":
            # seq_vec is a tuple of (source,  list of (sid, fasta_seq_header, fasta_seq in u8), w, k, r, min_span)
            source, seq_l, w, k, r, min_span = seq_list
            self.back_by = "SeqVec"
            self.seq_list = seq_list
            self.shmap = ShmmrFragMap()
            self.shmap.load_from_seq_list(source, seq_l, w, k, r, min_span)
            self.seq_index = {}
            self.seqs = {}
            for sid, seq_len, name, source in self.shmap.seq_index:
                self.seq_index[sid] = seq_len, name, source

        elif fastx_file_spec is not None and backed_by == "Fastx":
            self.back_by = "Fastx"
            fastx_file_path, w, k, r, min_span = fastx_file_spec
            self.shmap = ShmmrFragMap()
            self.shmap.load_from_seq_fastx(fastx_file_path, w, k, r, min_span)
            self.seq_index = {}
            for sid, seq_len, name, source in self.shmap.seq_index:
                self.seq_index[sid] = seq_len, name, source
        else:
            print("Not implemented")

    def get_aln_ranges(self, query_seq, gap_penality_factor=0.25, merge_range_tol=0):
        r = self.shmap.query_fragment_to_hps(query_seq, gap_penality_factor)
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
                aln_range[sid] = self.merge_regions(
                    rgns, tol=merge_range_tol)

        return aln_range

    def merge_regions(self, rgns, tol=1000):
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
