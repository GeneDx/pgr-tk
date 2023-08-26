import os
import json

if __name__ == "__main__":
    os.system("wget https://s3.amazonaws.com/igv.org.genomes/hg38/annotations/cytoBandIdeo.txt.gz")
    os.system("gunzip -f cytoBandIdeo.txt.gz")

    cypobands = {}
    with open("cytoBandIdeo.txt") as f:
        for row in f:
            row = row.strip().split("\t")
            cypobands.setdefault(row[0], [])
            cypobands[row[0]].append( (int(row[1]), int(row[2]), row[3], row[4]) )          
    
    out = open("cytoBandIdeo.json","w")
    json.dump({"cytobands": cypobands}, out)
    out.close()