

import parse_pdb_class

def write_cos(aa_tupal, struct_num, mol_num):
    if (len(aa_tupal) != 0):
        myX = myY = myZ = 0.0
        for j in range(len(aa_tupal)):
            myX += aa_tupal[j][0]
            myY += aa_tupal[j][1]
            myZ += aa_tupal[j][2]
        myX = myX /len(aa_tupal)
        myY = myY /len(aa_tupal)
        myZ = myZ /len(aa_tupal)
        return (struct_num, mol_num, myX, myY, myZ)


def get_aa_com(pdbList, first_struct_num, chain):
    aa_com_list = []
    curr_struct = first_struct_num
    aa_tup = []
    prev_mol_num = 0
    for i in range(len(pdbList)):
        line = pdbList[i]
        myPdbLine = parse_pdb_class.Parse_pdb_line(line)
        if (myPdbLine.line_type == "END"):
            aa_com_list.append(write_cos(aa_tup, curr_struct, prev_mol_num))
            aa_tup = []
            prev_mon_num = 0
            curr_struct += 1
    
        if (myPdbLine.chain != chain):
            continue
        if (myPdbLine.line_type != "ATOM"):
            continue
        if (int(myPdbLine.mol_num) != prev_mol_num):
            if (len(aa_tup) != 0):
                aa_com_list.append(write_cos(aa_tup, curr_struct, prev_mol_num))
            prev_mol_type = myPdbLine.mol_type
            prev_mol_num = int(myPdbLine.mol_num)
            aa_tup = []
        atom_type = myPdbLine.atom_type
        if (atom_type[0] == "H"):
            continue
        if ((atom_type == "N") | (atom_type == "C") | (atom_type == "O") | \
          (atom_type == "NT")):
            continue
        if ((atom_type == "CA") & (myPdbLine.monomer_type != "GLY")):
            continue
        aa_tup.append( (myPdbLine.xVal, myPdbLine.yVal, myPdbLine.zVal) )
    

    return aa_com_list
