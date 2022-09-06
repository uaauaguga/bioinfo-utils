#!/usr/bin/env python
import numpy as np
def main():
    sequences = {}
    sasas = {}
    with open("SASA.txt") as f:
        for line in f:
            assert line.startswith(">")
            seq_id = line[1:].strip()
            line = next(f)
            sequences[seq_id] = line.strip()
            line = next(f)
            sasa = np.array(line.strip().split(",")).astype(float)
            sasas[seq_id] = sasa
    paired_dict = {}
    for seq_id in sequences:
        paired_dict[seq_id] = np.full(len(sequences[seq_id]),False)
    # hydrog4079 hydrog ? ? GC C   1610 N3 ? ? ? 1_555 GC G   1651 N1 ? ? 2  C   1802 2  G   1843 1_555 ? ? ? ? ? ? WATSON-CRICK ?     ?
    pairs = [("A","U"),("U","A"),("C","G"),("G","C"),("G","U"),("U","G")]
    pairs = set(pairs)
    with open("6xyw.cif") as f:
        for line in f:
            line = line.strip()
            fields = line.strip().split()
            if not line.startswith("hydrog"):
                continue
            if len(fields) < 14:
                #print(line)
                continue
            c1, p1, c2, p2 = fields[5],fields[6],fields[13],fields[14] 
            p1, p2 = int(p1)-1, int(p2)-1
            
            seq_id = fields[4]
            if not  (c1 == sequences[seq_id][p1] and c2 == sequences[seq_id][p2]):
                #print(line)
                continue
                #print(c1,sequences[seq_id][p1],c2,sequences[seq_id][p2])
            if (c1,c2) in pairs:
                paired_dict[seq_id][p1] = True
                paired_dict[seq_id][p2] = True
    for seq_id in paired_dict:
        print(">" + seq_id)
        print(sequences[seq_id])
        print(",".join(paired_dict[seq_id].astype(int).astype(str)))
        print(",".join(np.round(sasas[seq_id],3).astype(str)))


if __name__ == "__main__":
    main()
