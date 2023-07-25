
import read_file
import sys,os
import numpy

def gp_from_csv(in_file, out_file, gp_template, title, xlabel, ylabel, column_list):
    
    in_list = read_file.Read_file_list(in_file)
    j=0
    while (len(in_list[j]) <= 1):
        j += 1
    col_head = in_list[j].split(",")
    col_include = []
    for j in range (len(column_list)):
        found = 0
        for k in range (len(col_head)):
            off = col_head[k].find(column_list[j])
            if (off >= 0):
                col_include.append(k)
                found = 1
                break
        if (found == 0):
            print "COLUMN NOT FOUND", column_list[j]

    data_out_name = out_file + "_gp.txt"
    out_fp = open(data_out_name,"w")
    new_head = column_list[0].replace(" ", "_")
    out_string = new_head
    for i in range (1,len(column_list)):
        new_head = column_list[i].replace(" ", "_")
        out_string = out_string + " " + new_head
    out_string = out_string + "\n"
    out_fp.write(out_string)
    
    for i in range (1,len(in_list)):
        fields = in_list[i].split(",")
        if (len(fields) < len(column_list)):
            continue
        try:
            out_string = fields[col_include[0]]
            for k in range(1, len(col_include)):
                out_string = out_string + " " + fields[col_include[k]]
            out_string = out_string + "\n"
            out_fp.write(out_string)
        except:
            print "Problem with line", in_list[i]
    out_fp.close()

    templ_list = read_file.Read_file_list(gp_template)
    print "TEMPLATE LEN", len(templ_list)
    out_gp_name = out_file + ".gp"
    out_gif_name = out_file + ".gif"
    num_plots = len(column_list) - 1
    boxwidth = 0.6 / num_plots
    str_boxwidth = str("%.2f" % boxwidth)
    offset = boxwidth 
    out_gp_fp = open(out_gp_name,"w")
    for j in range (len(templ_list)):
        line = templ_list[j]
        line = line.replace("$TITLE", title)
        line = line.replace("$XLABEL", xlabel)
        line = line.replace("$YLABEL", ylabel)
        line = line.replace("$OUT_FILE", out_gif_name)
        line = line.replace("$BOXWIDTH", str_boxwidth)
        out_gp_fp.write((line+"\n"))

    off_array = numpy.zeros(num_plots, dtype=float)
    start_off = ((offset * num_plots) / (-2.0) )
    for j in range (num_plots):
        off_array[j] = start_off + (j * offset)
        print "OFF ARRAY", j, off_array[j]

    out_string = "plot '" + data_out_name + "' "
    for j in range(num_plots):
        if (off_array[j] > 0):
            out_string = out_string + "using ($1+" + str("%.2f" % off_array[j]) + '):"' + \
              column_list[j+1] + '" with boxes, ' + "'' "
        elif (off_array[j] < 0):
            out_string = out_string + "using ($1" + str("%.2f" % off_array[j]) + '):"' + \
              column_list[j+1] + '" with boxes, ' + "'' "
        else:
            out_string = out_string + "using ($1):" + '"' + \
              column_list[j+1] + '" with boxes, ' + "'' "
    out_string = out_string[:-5] + ";\n"
    out_gp_fp.write(out_string)
    out_string = "quit\n"
    out_gp_fp.write(out_string)

    out_gp_fp.close()
    cmd = "gnuplot " + out_gp_name
    os.system(cmd)
