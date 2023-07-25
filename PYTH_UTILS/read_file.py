#! /usr/bin/python3



def Read_file(filename):
    return read_file(filename)
def read_file(filename):
    try:
        fp = open(filename,"r")
        data = fp.read()
        fp.close()
        return(data)
    except:
        return("")


def Read_file_list(filename):
    return read_file_list(filename)
def read_file_list(filename):
    try:
        fp = open(filename,"r")
        data = fp.read()
        fp.close()
        list = data.split("\n")
        outList = []
        for i in range(len(list)):
            if (len(list[i]) != 0):
                outList.append(list[i])
        return(outList)
    except:
        print("EMPTY or no file", filename)
        list = []
        return(list)

def read_file_list_rem_blank(filename):
    try:
        fp = open(filename,"r")
        data = fp.read()
        fp.close()
        list = data.split("\n")
        outList = []
        for i in range(len(list)):
            if (len(list[i]) > 1):
                outList.append(list[i])
        return(outList)
    except:
        print("EMPTY or no file", filename)
        list = []
        return(list)



def comp_space(line):
    if (len(line) > 0):
        oldLine = "HELLO"
        while (oldLine != line):
            oldLine = line
            line = line.replace('\t', " ")
            line = line.replace('\r', " ")
            line = line.replace('\n', " ")
        while (line[0] == ' '):
            line = line[1:]
            if (len(line) == 0):
                return ""
        while (line[-1] == ' '):
            line = line[:-1]
        oldLen = 99999
        newLen = len(line)
        while (oldLen != newLen):
            oldLen = newLen
            line = line.replace("  ", " ")
            newLen = len(line)
        while (line[0] == ' '):
            line = line[1:]
        while (line[-1] == ' '):
            line = line[:-1]
    else:
        print ("EMPTY LINE")
    return line

def rem_space(line):
    newLine = line.replace(" ","")
    return newLine

#
#   Change hash special characters to correct values
#   Remove leading and trailing quotes
def trans_special(line):
    if ((line[0] == '"') & (line[-1] == '"')):
        line = line[1:-1]
    if ((line[0] == "'") & (line[-1] == "'")):
        line = line[1:-1]
    done = False
    while (done == False):
        done = True
        for i in range(len(line)-1):
            print (i, "LINE ", line )
            if (line[i] == "\\"):
                if (line[i+1] == "t"):
                    line = line[:i] + "\t" + line[i+2:]
                    done = False
                    break
                if (line[i+1] == "s"):
                    line = line[:i] + " " + line[i+2:]
                    done = False
                    break

    return line
