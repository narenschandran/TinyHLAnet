import re
import os

# Generates modified BLOSUM matrix from EMBOSS' BLOSUM

aa1lst = [
    "A",
    "C",
    "D",
    "E",
    "F",
    "G",
    "H",
    "I",
    "K",
    "L",
    "M",
    "N",
    "P",
    "Q",
    "R",
    "S",
    "T",
    "V",
    "W",
    "Y",
]


def line_proc(line):
    return re.sub("[ ]+", " ", line.strip()).split(" ")


def blosum_dictstr(f, name=None):
    fcon = open(f, "r")
    lines = [line_proc(line) for line in fcon.readlines() if line[0] != "#"]
    fcon.close()

    if name is None:
        name = os.path.basename(f)
    hd = lines[0]
    datlines = lines[1:]

    datlines = [datline for datline in datlines if datline[0] in aa1lst]

    res = {}
    for datline in datlines:
        aakey = datline[0]
        dct = {}
        for i in range(1, len(datline)):
            aaind = i - 1
            aa = hd[aaind]
            if aa in aa1lst:
                dct[aa] = int(datline[i])
        dct["X"] = 0
        res[aakey] = dct

    res["X"] = {aa: 0 for aa in hd if aa in aa1lst}
    res["X"]["X"] = 1
    k = list(res.keys())
    k.sort()

    resstr = []
    for i in range(len(k)):
        mainkey = k[i]
        valkeys = list(res[mainkey].keys())
        valkeys.sort()
        valstr = ["'%s': %d" % (valkey, res[mainkey][valkey]) for valkey in valkeys]
        resstr.append("'%s': {" % (mainkey) + ",".join(valstr) + "}")

    outstr = ",\n".join(resstr)
    outstr = "%s = {" % (name) + outstr + "}"
    return outstr


fs = ["EBLOSUM50", "EBLOSUM60", "EBLOSUM80"]

dctstrs = "\n".join([blosum_dictstr(f, os.path.basename(f)[1:]) for f in fs])

outfcon = open("BLOSUM.py", "w")
outfcon.writelines(dctstrs)
outfcon.close()
