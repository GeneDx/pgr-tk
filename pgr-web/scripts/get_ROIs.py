import gzip
import json


#gene2query = {}
gene2query = json.loads(open("ROIs_examples.json").read())
## we need the file of the coordinate of genes from https://s3.amazonaws.com/igv.org.genomes/hg38/ncbiRefSeq.sorted.txt.gz
with gzip.open("ncbiRefSeq.sorted.txt.gz") as f:
    for row in f:
        row = row.decode("utf-8")
        row = row.strip().split("\t")
        g = row[12]
        ch = row[2]
        if len(ch.split("_")) > 1:
            continue
        strand = row[3]
        bgn = int(row[4])
        end = int(row[5])
        if g not in gene2query:
            gene2query[g] = {
                "source": "hg38_tagged.fa",
                "ctg": f"{ch}_hg38",
                "bgn": bgn,
                "end": end,
                "padding": 10000,
                "merge_range_tol": 120000,
                "w": 48,
                "k": 56,
                "r": 1,
                "min_span": 12,
                "sketch": False,
                "min_cov": 2,
                "min_branch_size": 8,
                "bundle_length_cutoff": 500,
                "bundle_merge_distance": 10000
            }


print(json.dumps(gene2query)) 

