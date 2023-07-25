

def aa_let_to_name(letter):
    letter_to_name = { "A" : "ALA", "R" : "ARG", "N" : "ASN", "D" : "ASP", \
      "C" : "CYS", "E" : "GLU", "Q" : "GLN", "G" : "GLY", \
      "H" : "HIS", "I" : "IIE", "L" : "LEU", "K" : "LYS", \
      "M" : "MET", "F" : "PHE", "P" : "PRO", "S" : "SER", \
      "T" : "THR", "W" : "TRP", "Y" : "TYR", "V" : "VAL" }

    return letter_to_name.get(letter)



