#! /usr/bin/python
import os, sys
sys.path.append('.')

import read_file

def parse_pdb_line(line):

    if (len(line) < 72):
        return("", 0, "", "", 0, 0.0, 0.0, 0.0, "", "",  "")
    try:
        field1 = line[:4]
        atom_num = int(read_file.rem_space(line[6:12]))
        atom_type = read_file.rem_space(line[12:17])
        mol_type = read_file.rem_space(line[17:21])
        mol_num = read_file.rem_space(line[24:27])
        xVal = float(read_file.rem_space(line[27:39]))
        yVal = float(read_file.rem_space(line[39:47]))
        zVal = float(read_file.rem_space(line[47:55]))
        colA = read_file.rem_space(line[55:61])
        colB = read_file.rem_space(line[61:67])
        label = ""
        if (len(line) >= 76):
            label = read_file.rem_space(line[67:77])
    except:
        print "ERROR _ line ", line
        return("ERROR", 0, "", "", 0, 0.0, 0.0, 0.0, "", "",  "")
    return ((field1, atom_num, atom_type, mol_type, mol_num, xVal, yVal, zVal, colA, colB, label))



