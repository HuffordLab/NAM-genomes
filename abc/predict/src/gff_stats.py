


c_feat = "CDS"
with open("data/ref/Zm-B73-REFERENCE-NAM-5.0_Zm00001e.1.gff") as fl:
    for line in fl:
        if line[0] != "#":
            seqname, source, feature, start, end, score, strand, frame, attribute = line.strip().split()
            if feature == c_feat:
                print(f"{start}, {end}, {attribute}")
            
