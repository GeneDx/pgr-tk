import pgrtk
sdb = pgrtk.SeqIndexDB()
sdb.load_from_fastx("test_seqs.fa")
sdb.write_frag_and_index_files("test_seqs_frag")
