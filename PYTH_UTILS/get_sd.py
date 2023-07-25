
#
#   sd, se, mean, num samples = get_sd_str( list_of_strings )
def get_sd(data_list):
    totVal = 0.0
    numVal = 0
    for i in range(len(data_list)):
        totVal += float(data_list[i])
        numVal += 1
    meanVal = totVal /numVal

    sqSum = 0.0
    for i in range(len(data_list)):
        totVal = (float(data_list[i]) - meanVal)
        sqSum += totVal ** 2

    sdVal = (sqSum / (numVal-1)) ** 0.5

    seVal = sdVal / (numVal ** 0.5)

    return sdVal, seVal, meanVal, numVal


#
#   sd, se, mean, num samples = get_sd_min(  )
def get_sd_min(data_list):
    totVal = 0.0
    numVal = 0
    minVal = 999999999999999999999999999999
    maxVal = -999999999999999999999999999999
    for i in range(len(data_list)):
        val = float(data_list[i])
        if (val > maxVal):
            maxVal = val
        if (val < minVal):
            minVal = val
        totVal += float(data_list[i])
        numVal += 1
    meanVal = totVal /numVal

    sqSum = 0.0
    for i in range(len(data_list)):
        totVal = (float(data_list[i]) - meanVal)
        sqSum += totVal ** 2

    sdVal = (sqSum / (numVal-1)) ** 0.5

    seVal = sdVal / (numVal ** 0.5)

    return sdVal, seVal, meanVal, numVal, minVal, maxVal


