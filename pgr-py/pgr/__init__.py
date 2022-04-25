from .pgr import *
import bisect


rc_map = dict(zip("ACGTNacgt", "TGCANTGCA"))


def rc(s):
    return "".join([rc_map[c] for c in s[::-1]])


class ContigMap(object):

    def __init__(self, ref_path, ctg_path):
        self.ref_path = ref_path
        self.ctg_path = ctg_path
        self.ref_db = None
        self.ctg_db = None
        self.shmmr_map = None
        self.shmmr_spec = (48,56,4)

    def load_seqs(self, build_shmmr_idx=True, shmmr_spec=None):
        self.ref_db = SeqDB(self.ref_path)
        self.ctg_db = SeqDB(self.ctg_path)
        self.ref_db.load_sequences()
        self.ctg_db.load_sequences()
        if shmmr_spec is not None:
            self.shmmr_spec = shmmr_spec 
        if build_shmmr_idx:
            self.ref_db.build_shmmrs_parallel(*self.shmmr_spec)
            self.ctg_db.build_shmmrs_parallel(*self.shmmr_spec)

    def build_shmmr_map(self, max_rpt_count=32):
        self.shmmr_map  = generate_shmmr_map(self.ref_db, self.ctg_db, max_rpt_count)
    
    def break_chains(self, chain, size):
        all_chains = []
        cur_chain = [chain[0]]
        
        for elm in chain[1:]:
            bgn = cur_chain[0][1]
            if elm[1] == cur_chain[-1][2] and elm[1] - bgn > size:
                all_chains.append(cur_chain)
                cur_chain = [elm]
            else:
                cur_chain.append(elm)
                
        if len(cur_chain) > 0:
            all_chains.append(cur_chain)
        
        return all_chains


    def generate_aln_chains(self):
        id_len_v = []
        for id_ in self.ref_db.get_all_ids():
            id_len_v.append([id_, self.ref_db.get_len_by_id(id_)])

        mr = map_seqs_with_db(id_len_v, self.shmmr_map)
        
        aln_chains = {}
        for m in mr:
            key = (m[0], m[3], m[6], m[7])
            aln_chains.setdefault(key, [])
            aln_chains[key].append(m)
        aln_chains2 = {}
        for k,v in aln_chains.items():
            v.sort(key=lambda x:x[1])
            aln_chains2[k] = self.break_chains(v, 16000) 
        return aln_chains2
    

    def aln_chains_to_variants(self, aln_chains):
        """
            going through each alignment chain and generate
            variant calls
        """
        ref_db = self.ref_db
        ctg_db = self.ctg_db

        variant_calls = {}
        aln_chunks = []
        aln_maps = {}
        cov = {}

        all_aln_chains = []
        for k, chain in aln_chains.items():
            sub_chain_id = 0
            for sub_chain in chain:
                new_key = (k[0], k[1], k[2], (k[3], sub_chain_id))
                all_aln_chains.append((new_key, sub_chain))
                sub_chain_id += 1
            

        for k, sub_chain in all_aln_chains:
            # print(k, len(aln_chains[k]))
            ref_id = k[0]
            ctg_id = k[1]
            strand = k[2]
            chain_id = k[3]
            first = sub_chain[0]
            last = sub_chain[-1]

            if strand == 0:
                ref_bgn = first[1]-55
                ctg_bgn = first[4]-55

                ref_end = last[2]
                ctg_end = last[5]
            else:
                ref_bgn = first[1]-55
                ctg_bgn = first[5]+1

                ref_end = last[2]
                ctg_end = last[4]-54
                ctg_bgn, ctg_end = ctg_end, ctg_bgn

            # print(ref_id, ref_bgn, ref_end, ctg_id, ctg_bgn, ctg_end, strand)

            rs0 = ref_db.get_subseq_by_id(ref_id, ref_bgn, ref_end)
            cs0 = ctg_db.get_subseq_by_id(ctg_id, ctg_bgn, ctg_end).upper()

            if strand == 1:
                cs0 = rc(cs0)
            assert(rs0[:56] == cs0[:56])

            aln_segs = get_aln_segements(ref_id, rs0, ctg_id, cs0)
            try: #work around the alinger's bug
                aln_map = get_aln_map(aln_segs, rs0, cs0)
            except:
                # TODO: record missing chunk
                continue
            aln_chunk = (ref_id, ref_bgn, ref_end, ctg_id,
                            ctg_bgn, ctg_end, strand, chain_id)
            
            
            if (ctg_id, ctg_bgn, strand, chain_id) not in aln_maps:
                aln_chunks.append(aln_chunk)
                aln_maps[(ctg_id, ctg_bgn, strand, chain_id)] = (ref_bgn, ctg_bgn, aln_map)
            else:
                (ref_bgn0, ctg_bgn0, _) = aln_maps[(ctg_id, ctg_bgn, strand, chain_id)]
                print((ctg_id, ctg_bgn, strand, chain_id), (ref_bgn0, ctg_bgn0), (ref_bgn, ctg_bgn))
                continue
            # print("".join([chr(c) for c in aln_map.ref_a_seq[:1000]]))
            # print("".join([chr(c) for c in aln_map.aln_seq[:1000]]))
            # print("".join([chr(c) for c in aln_map.tgt_a_seq[:1000]]))

            for s in aln_segs:
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
                    variant_calls[key][(s.tgt_loc[0], strand, chain_id)] = value

                    cov.setdefault(key, set())
                    cov[key].add((ctg_id, ctg_bgn, strand, chain_id))

        
        ## need a faster algorithm
        v_keys = list(variant_calls.keys())
        v_keys.sort()

        aln_chunks.sort()
        for r in aln_chunks:
            (ref_id, ref_bgn, ref_end, ctg_id,
                ctg_bgn, ctg_end, strand, chain_id) = r
            left = (ref_id, ref_bgn)
            right = (ref_id, ref_end)
            left_idx = bisect.bisect_right(v_keys, left)        
            right_idx = bisect.bisect_left(v_keys, right)
            
            for idx in range(left_idx, right_idx):
                k = v_keys[idx]
                assert(k[1] >= ref_bgn and k[1] < ref_end)
                cov.setdefault(k, set())        
                cov[k].add((ctg_id, ctg_bgn, strand, chain_id))
                (ref_bgn2, ctg_bgn2, aln_map) = aln_maps[(ctg_id, ctg_bgn, strand, chain_id)] 
                assert(ref_bgn==ref_bgn2)

        return aln_chunks, aln_maps, cov, variant_calls


    def output_variant_recs(self, variant_data, filename):
        """
            taking the alginemnt chains and variant calls from each chain
            and genrate merged output
        """
        ref_db = self.ref_db
        ctg_db = self.ctg_db
        aln_chunks, aln_maps, cov, variant_calls = variant_data

        file_handle = open(filename, "w")
        for r in aln_chunks:
            ref_id, ref_bgn, ref_end, ctg_id, ctg_bgn, ctg_end, strand, chain_id = r
            ref_name = ref_db.get_name_by_id(ref_id)
            ctg_name = ctg_db.get_name_by_id(ctg_id)
            print(
                "##", f"{ref_name}\t{ref_bgn}\t{ref_end}\t{ctg_name}\t{ctg_bgn}\t{ctg_end}\t{strand}\t{chain_id[0]}\t{chain_id[1]}", file=file_handle)

        keys = sorted(list(variant_calls.keys()))
        for k in keys:
            #if k not in cov: ## TODO, perhaps a bug
            #    continue
            v = variant_calls[k]
            #print(k, v, cov[k])
            ref_id = k[0]
            out = [ref_db.get_name_by_id(k[0]), "{}".format(k[1])]

            gt_set = set()
            gt = []
            gt_idx = 0
            out.append("{}".format(len(cov[k])))
            #print((ctg_id, strand), cov[k])
            for ctg_id, ctg_bgn, strand, chain_id in sorted(list(cov[k])):

                if (ctg_id, strand, chain_id) in v:
                    # print(len(v[ctg_id]))
                    vv = v[(ctg_id, strand, chain_id)]
                    p = vv[4][1]
                    strand = "+" if strand == 0 else "-"
                    ctg_name = ctg_db.get_name_by_id(ctg_id)

                    vtag = ":".join(vv[:3])
                    out.append(
                        f"{ctg_name}:{strand}:{chain_id[0]}:{chain_id[1]}:{p}:"+":".join(vv[:3]))
                    if vtag not in gt_set:
                        gt_set.add(vtag)
                        gt_idx += 1
                    gt.append("{}".format(gt_idx))

                elif (ctg_id, ctg_bgn, strand, chain_id) in cov[k]:
                    ref_bgn, ctg_bgn, aln_map = aln_maps[(
                        ctg_id, ctg_bgn, strand, chain_id)]
                    strand = "+" if strand == 0 else "-"
                
                    p = aln_map.pmap[k[1]-ref_bgn][1] + ctg_bgn
                    ctg_name = ctg_db.get_name_by_id(ctg_id)
                    ref_base = ref_db.get_subseq_by_id(ref_id, k[1]-1, k[1])
                    out.append(f"{ctg_name}:{strand}:{chain_id[0]}:{chain_id[1]}:{p}:" + f"M:{ref_base}:{ref_base}")
                    gt.append("0")

                else:
                    gt.append(".")

            out.append("|".join(gt))
            print("\t".join(out), file=file_handle)

def output_variant_recs_to_vcf(ctg_map, variant_data, filename):
    """
        taking the alginemnt chains and variant calls from each chain
        and genrate merged output
    """
    ref_db = ctg_map.ref_db
    ctg_db = ctg_map.ctg_db
    aln_chunks, aln_maps, cov, variant_calls = variant_data

    file_handle = open(filename, "w")
    """
    for r in aln_chunks:
        ref_id, ref_bgn, ref_end, ctg_id, ctg_bgn, ctg_end, strand, chain_id = r
        ref_name = ref_db.get_name_by_id(ref_id)
        ctg_name = ctg_db.get_name_by_id(ctg_id)
        print(
            "##", f"{ref_id}\t{ref_bgn}\t{ref_end}\t{ctg_name}\t{ctg_bgn}\t{ctg_end}\t{strand}\t{chain_id[0]}\t{chain_id[1]}", file=file_handle)
    """
    
    print(open("vcf_header.vcf").read().strip(), file=file_handle)
    
    keys = sorted(list(variant_calls.keys()))
    
    for k in keys:
      
        v = variant_calls[k]
        
        ref_id = k[0]

        out = [ref_db.get_name_by_id(ref_id), "{}".format(k[1]), "."]
        if out[0] == "2":
            print(out)
        gt_set = set()
        gt = []
        gt_idx = 0

        variants = []
        for ctg_id, ctg_bgn, strand, chain_id in sorted(list(cov[k])):

            if (ctg_id, strand, chain_id) in v:
                # print(len(v[ctg_id]))
                vv = v[(ctg_id, strand, chain_id)]
                p = vv[4][1]
                strand = "+" if strand == 0 else "-"
                ctg_name = ctg_db.get_name_by_id(ctg_id)

                vtag = tuple(vv[:3])
                
                if vtag not in gt_set:
                    gt_set.add(vtag)
                    gt_idx += 1
                    
                variants.append((gt_idx, vtag, ctg_id))
                gt.append("{}".format(gt_idx))

            elif (ctg_id, ctg_bgn, strand, chain_id) in cov[k]:
                ref_bgn, ctg_bgn, aln_map = aln_maps[(
                    ctg_id, ctg_bgn, strand, chain_id)]
                strand = "+" if strand == 0 else "-"

                p = aln_map.pmap[k[1]-ref_bgn][1] + ctg_bgn
                ctg_name = ctg_db.get_name_by_id(ctg_id)
                
                ref_base = ref_db.get_subseq_by_id(ref_id, k[1]-1, k[1])
                vtag = ("M", ref_base, ref_base)

                variants.append((0, vtag, ctg_id))
                gt.append("0")

            else:
                gt.append(".")
                
        variants = sorted(list(variants))
           
        if variants[-1][0] <= 1:
            count = {}
            for v in variants:
                count.setdefault(v[0], [])
                count[v[0]].append(v[1])
            if 0 in count and 1 in count:
                het = True
            else:
                het = False
                
            ht = [(str(x[2]), str(x[0])) for x in variants]
            ht.sort()
            ht = list(zip(*ht))
            if het == False:
                # print(count)
                ref_base=count[1][0][1]
                alt_base=count[1][0][2]
                #print(count)
                if ref_base == alt_base:
                    continue
                if 0 not in count and len(count[1]) == 1:
                    print("\t".join(out), ref_base, alt_base, 
                          "30", ".", ".", "GT:AD:CT:HT", "./1:0,1:" + ",".join(ht[0]) + ":" + ",".join(ht[1]), 
                          sep="\t", file=file_handle)
                else:
                    print("\t".join(out), ref_base, alt_base, 
                          "30", ".", ".", "GT:AD:CT:HT", "1/1:1,1:" + ",".join(ht[0]) + ":" + ",".join(ht[1]), 
                          sep="\t", file=file_handle)
            else:
                ref_base=count[1][0][1]
                alt_base=count[1][0][2]
                print("\t".join(out), ref_base, alt_base, 
                      "30", ".", ".", "GT:AD:CT:HT", "0/1:0,2:" + ",".join(ht[0]) + ":" + ",".join(ht[1]), 
                      sep="\t", file=file_handle)
        else:
            for v in sorted(list(set(variants))):
                ht = [(str(x[2]), str(x[0])) for x in [v]]
                ht = list(zip(*ht))
                #print(v, file=file_handle)
                ref_base=v[1][1]
                alt_base=v[1][2]
                if ref_base == alt_base:
                    continue
                print("\t".join(out), ref_base, alt_base, 
                      "30", ".", ".", "GT:AD:CT:HT", "./1:0,1:"+ ",".join(ht[0]) + ":" + ",".join(ht[1]), 
                      sep="\t", file=file_handle)