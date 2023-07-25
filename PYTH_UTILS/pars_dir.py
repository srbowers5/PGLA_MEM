
import os, sys


def ParseDir(dirname, suffix, prefix="", min_len=0):
    suf_len = 0 - len(suffix)
    prefix_len = len(prefix)
    names = os.listdir(dirname)
    good_names = []
    path_list = []
    for i in range(len(names)):
        name_dir = dirname + "/" + names[i]
        file_size = os.stat(name_dir).st_size
        if (file_size < min_len):
            print "ZERO Len ", names[i]
            continue
        if (suf_len < 0):
            if (names[i][suf_len:] != suffix):
                continue
        if (prefix_len > 0):
            if (names[i][0:prefix_len:] == prefix):
                good_names.append(names[i])
        else:
            good_names.append(names[i])
    path_list = map (lambda x: dirname + "/"+x+"\n", good_names)
    return(good_names, path_list)

