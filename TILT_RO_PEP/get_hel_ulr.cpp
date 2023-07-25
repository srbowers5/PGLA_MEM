#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include <stdarg.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>
#include <arpa/inet.h>
#include <math.h>
//#include <cmath.h>
#include "get_hel_ulr.h"
#include "ParseLine.h"


int main(int argc, char * argv[]) {

    DIH_STRUCT dih_up;
    DIH_STRUCT dih_lo;
    FILE * outFpUp;
    FILE * outFpLo;
    char tmpLine[MAX_PARSE_LINE_LEN];
    char * parsedLine;
    size_t parsedLen;
    std::vector <char *> fields;
    int cnt;
    int nCnt;
    char * stringPtr;
    int i;
    int inHelCnt;

    char hel_up[20];
    char hel_lo[20];


//    arg[1] = num_structs
//    arg[2] = Input_file
//    arg[3] = Output_file

    parsedLine = &(tmpLine[0]);
    parsedLen = MAX_PARSE_LINE_LEN;
    ParseLine fileParsed(argv[1]);
    printf("OPEN %s\n", argv[1]);
    outFpUp = fopen(argv[2], "w");
    outFpLo = fopen(argv[3], "w");
    cnt = 0;
    while (nCnt = fileParsed.get_line_fields(&fields, &parsedLine, &parsedLen) > 0) {
        if (nCnt <= 0) {
            printf("DONE with reading\n");
            break;
        }
        if (fields.size() != NUM_FIELDS) {
            printf("ERROR RETURN fields size %ld %d struct %d\n", fields.size(), NUM_FIELDS, cnt);
            break;
        }

//   Starts at field[1]
        for (i=0; i<NUM_ANGLES; i++) {
            dih_up.phi_val[i] = atof(fields[i+1]);
        }
        for (i=0; i<NUM_ANGLES; i++) {
            dih_lo.phi_val[i] = atof(fields[i+1+(NUM_ANGLES)]);
        }
        for (i=0; i<NUM_ANGLES; i++) {
            dih_up.psi_val[i] = atof(fields[i+1+(NUM_ANGLES*2)]);
        }
        for (i=0; i<NUM_ANGLES; i++) {
            dih_lo.psi_val[i] = atof(fields[i+1+(NUM_ANGLES*3)]);
        }

        for (i=0; i<NUM_ANGLES; i++) {
            if ((dih_up.phi_val[i] > -87.0) && (dih_up.phi_val[i] < -27.0) &&
              (dih_up.psi_val[i] > -77.0) && (dih_up.psi_val[i] < -17.0)) {
                hel_up[i] = 'H';
            } else {
                hel_up[i] = '0';
            }
            if ((dih_lo.phi_val[i] > -87.0) && (dih_lo.phi_val[i] < -27.0) &&
              (dih_lo.psi_val[i] > -77.0) && (dih_lo.psi_val[i] < -17.0)) {
                hel_lo[i] = 'H';
            } else {
                hel_lo[i] = '0';
            }
        }
        inHelCnt = 0;
        for (i=0; i<NUM_ANGLES; i++) {
            if (hel_up[i] == 'H') {
                inHelCnt++;
            } else {
                if (inHelCnt == 1) {
                    hel_up[i-1] = '0';
                } else if (inHelCnt == 2) {
                    hel_up[i-1] = '0';
                    hel_up[i-2] = '0';
                }
                inHelCnt = 0;
            }
        }

        inHelCnt = 0;
        for (i=0; i<NUM_ANGLES; i++) {
            if (hel_lo[i] == 'H') {
                inHelCnt++;
            } else {
                if (inHelCnt == 1) {
                    hel_lo[i-1] = '0';
                } else if (inHelCnt == 2) {
                    hel_lo[i-1] = '0';
                    hel_lo[i-2] = '0';
                }
                inHelCnt = 0;
            }
        }
        hel_lo[19] = hel_up[19] = 0;

        fprintf(outFpUp, "%s\n", hel_up);
        fprintf(outFpLo, "%s\n", hel_lo);


//        printf("HEL %s    %s\n", hel_up, hel_lo);
        cnt++;
    }
    fclose(outFpUp);
    fclose(outFpLo);
    printf("DONE get_hel_url\n");
}
