# Radius using Cα as pivot calculated from a
# set of non-redundant pMHC-I  structures from PDB
pmhc_ca0 = {
    "A": 2.537,
    "C": 2.973,
    "D": 3.835,
    "E": 5.160,
    "F": 5.336,
    "G": 2.510,
    "H": 4.785,
    "I": 4.078,
    "K": 6.688,
    "L": 4.121,
    "M": 5.649,
    "N": 3.835,
    "P": 2.512,
    "Q": 5.099,
    "R": 7.577,
    "S": 2.597,
    "T": 2.722,
    "V": 2.683,
    "W": 6.793,
    "Y": 6.759,
}

vals = [v for v in pmhc_ca0.values()]
pmhc_ca = pmhc_ca0.copy()
pmhc_ca["X"] = sum(vals) / len(vals)
