
import read_file
#
#   Sum columns and get average for each column
#
#   Input:
#       List with lines of data
#   The output is 2 lists:
#      first list of sums
#      second list of averages
def sum_avg_col (data_list):
    colVal = []
    colAvg = []

    line = read_file.comp_space(data_list[0])
    fields = line.split(' ')
    for col in range(len(fields)):
        colVal.append(float(fields[col]))
    cnt = 1

    for row in range(1, len(data_list)):
        line = read_file.comp_space(data_list[row])
        fields = line.split(' ')
        for col in range(len(fields)):
            colVal[col] += float(fields[col])
        cnt += 1

    for i in range(len(colVal)):
        colAvg.append(colVal[i]/cnt)

    return colVal, colAvg

