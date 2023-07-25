#
#   This gets a PDB file given the PDB structure name and
#   directory in which to put the file.
#   If the file already exists, then it just returns 0.
#   If the file cannot be gotten, then the function returns -1.
#
#      find_pdb_file(pdb_struct, pdb_file_dir)
import os,sys
import urllib
import urllib2



Download_url="http://www.rcsb.org/pdb/files/"

def Load_pdb_struct(struct_name, out_name):
    url_string = Download_url + struct_name + ".pdb"
    try:
        resp1 = urllib2.urlopen(url_string);
        web_pg = resp1.read();
        if (len(web_pg) < 10):
            return(-1)
        fp  = open(out_name, "w")
        fp.write(web_pg)
        fp.close()
        return(0)
    except:
        return(-1)

def find_pdb_file(pdb_struct, pdb_file_dir):
    pdb_file = pdb_struct[0:4]
    out_pdb_name =  pdb_file_dir +  "/" + pdb_file + ".pdb"
    if ((os.path.isfile(out_pdb_name)) == False):
        ret_val = Load_pdb_struct(pdb_file, out_pdb_name)
        if (ret_val < 0):
            return(-1)
    if ((os.path.isfile(out_pdb_name)) == False):
        return(-1)
    return(0)


