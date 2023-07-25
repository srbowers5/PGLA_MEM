#! /usr/bin/python

#
#   In shell scripts don't forget to use \$ not just $
def make_name (inName, traj="", rep="", aa="", step=""):
    name = inName
    offset = 0
    while (offset >= 0):
        offset = name.find("$tr", offset)
        if (offset >= 0):
            name = name[:offset] + str(traj) + name[offset+3:]
    offset = 0
    while (offset >= 0):
        offset = name.find("$aa", offset)
        if (offset >= 0):
            name = name[:offset] + str(aa) + name[offset+3:]
    offset = 0
    while (offset >= 0):
        offset = name.find("$rep", offset)
        if (offset >= 0):
            name = name[:offset] + str(rep) + name[offset+4:]
       
    offset = 0
    while (offset >= 0):
        offset = name.find("$step", offset)
        if (offset >= 0):
            name = name[:offset] + str(step) + name[offset+5:]


    return name


#main
#inVal = "$trabc$trdef$trghi$tr"
#outVal = make_name(inVal, 3)

#print "IN  (", inVal, ")"
#print "OUT (", outVal, ")"
