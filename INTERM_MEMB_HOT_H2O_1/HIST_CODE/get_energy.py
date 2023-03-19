#! /usr/bin/python3
#
#    This calculaes the energy from the dcd files created by namd.
#
import os,sys
MIN_ARGS=8
if (len(sys.argv) < MIN_ARGS):
    print ("NOT ENOUGH argv:\n   get_energy.py <tr> <rep> <dir> <template> <prefix>")
    exit()


for i in range(MIN_ARGS,len(sys.argv)):
    print ("add path", sys.argv[i])
    sys.path.append(sys.argv[i])

traNum = sys.argv[1]
firstRep = int(sys.argv[2])
lastRep = int(sys.argv[3])
dir = sys.argv[4]
template = sys.argv[5]
firstStr = sys.argv[6]
lastStr = sys.argv[7]


inFile = dir + "/RepTemp.dat"
import read_file


for rep in range(firstRep,lastRep+1):
    temps = read_file.read_file_list(inFile)
    print ("READ ", inFile, len(temps))
    temp = temps[rep]
    print ("TEMP is ", temp)

    filename = template
    energy1 = read_file.read_file_list(filename)
    outfile = "tmpE" + str(rep) + ".vmd"
    outFd = open(outfile, "w")
    for i in range(len(energy1)):
        line = energy1[i]
        off = line.find("TRAJ")
        if (off >= 0):
            outstr = " set rem " + traNum + "\n"
        elif (line.find("REP") >= 0):
            outstr = " set rep " + str(rep+1) + "\n"
        elif (line.find("TEMP") >= 0):
            outstr = " set t " + temp + "\n"
        elif (line.find("DIR") >= 0):
            outstr = ' set dir ' + dir + '\n'
        elif (line.find("START") >= 0):
            outstr = ' set firstStr ' + firstStr + '\n'
        elif (line.find("END") >= 0):
            outstr = ' set lastStr ' + lastStr + '\n'
        else:
           outstr = line + "\n"
        outFd.write(outstr)
        print ("WRITE ", outstr[:-1])
    outFd.close()

    cmd = "vmd -dispdev none -e " + outfile + " > en_" + template + "_" + str(rep+1) + ".out"
    print ("CMD ", traNum, rep, "(", cmd, ")")
    os.system(cmd)

