
import read_file

class Parse_pdb_line(object):
    def __init__(self, line):
        self.line_type = ""
        self.atom_num = 0
        self.atom_type = ""
        self.monomer_type = ""
        self.mol_type = ""
        self.mol_num = 0
        self.xVal = 0.0
        self.yVal = 0.0
        self.zVal = 0.0
        self.colA = ""
        self.colB = ""
        self.chain = ""

        if (len(line) < 6):
            self.line_type = line
            return

        self.line_type = read_file.rem_space(line[:6])
        if (self.line_type != "ATOM"):
            return

        if (len(line) < 72):
            return

        self.atom_num = int(read_file.rem_space(line[6:11]))
        self.atom_type = read_file.rem_space(line[12:16])
        element = read_file.comp_space(self.atom_type)
        self.element = element[0]
        self.monomer_type = read_file.rem_space(line[17:21])
        self.mol_type = read_file.rem_space(line[21:22])
        self.mol_num = read_file.rem_space(line[22:27])
        self.xVal = float(read_file.rem_space(line[30:38]))
        self.yVal = float(read_file.rem_space(line[38:46]))
        self.zVal = float(read_file.rem_space(line[46:54]))
        self.colA = read_file.rem_space(line[55:61])
        self.colB = read_file.rem_space(line[61:67])
        if (len(line) < 76):
            chain = ""
        else:
            self.chain = read_file.rem_space(line[67:77])

        return

    def add_spaces(self,inLine, num_spaces):
        for i in range(0,num_spaces):
            inLine += " "
        return inLine

    def write_pdb_line(self, fd):
        line = ""
        line = self.line_type
        line = self.add_spaces(line, 6-len(self.line_type))
        atom_num = str(self.atom_num)
        line = self.add_spaces(line, 11-(len(line)+len(atom_num)))
        line += atom_num
        atom_type = self.atom_type
        line = self.add_spaces(line, 16-(len(line)+len(atom_type)))
        line += atom_type
        line = self.add_spaces(line, 20 -(len(line) + len(self.monomer_type)))
        line += self.monomer_type
        line += " "
        line += self.mol_type
        mol_num = self.mol_num
        line = self.add_spaces(line, 26 -(len(line) + len(mol_num)))
        line += mol_num
        xVal = str(round(self.xVal,3))
        line = self.add_spaces(line, 38 -(len(line) + len(xVal)))
        line += xVal
        yVal = str(round(self.yVal,3))
        line = self.add_spaces(line, 46 -(len(line) + len(yVal)))
        line += yVal
        zVal = str(round(self.zVal,3))
        line = self.add_spaces(line, 54 -(len(line) + len(zVal)))
        line += zVal
        colA = self.colA
        line = self.add_spaces(line, 60 -(len(line) + len(colA)))
        line += colA
        colB = self.colB
        line = self.add_spaces(line, 66 -(len(line) + len(colB)))
        line += colB
        chain = self.chain
        line = self.add_spaces(line, 76 -(len(line) + len(chain)))
        line += chain
        line += "\n"
        fd.write(line)


    def print_class(self):
        print ("TYPE ", self.line_type)
        print ("ATOM NUM", self.atom_num)
        print ("ATOM TYPE", self.atom_type)
        print ("MONIMER TYPE", self.monomer_type)
        print ("MOL TYPE", self.mol_type)
        print ("MOL NUM", self.mol_num)
        print ("Xval", self.xVal)
        print ("Yval", self.yVal)
        print ("Zval", self.zVal)
        print ("Col A", self.colA)
        print ("Col B", self.colB)
        print ("Chain", self.chain)

