#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include <stdarg.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>
#include <cmath>
#include <arpa/inet.h>
#include <math.h>
#include <set>
#include "dcdStruct.h"
#include "ReadPdb.h"
#include "ReadDcd.h"
#include "get_com_xsc.h"
#include "system_conf.h"
#include "PepCOM.h"
#include "get_water_region.h"


#define NUM_XY_REGION 3
#define MAX_Z_VAL  18.5
#define MIN_Z_VAL  -0.5
#define MAX_Z_BUCKET 18
#define MIN_Z_BUCKET 0
#define NUM_Z_BUCKS 19

int water_cnt[NUM_REG][NUM_Z_BUCK];

double get_dist_xy(float x1, float x2, float y1, float y2, float xSide) {
    double xDist;
    double yDist;
    double sqSum;
    double dist;
    double x1Cell, y1Cell, xSq, ySq, xSqCell, ySqCell, xCell, yCell;

    if (x2 > x1) {
        x1Cell = x1 + xSide;
    } else {
        x1Cell = x1 - xSide;
    }
    if (y2 > y1) {
        y1Cell = y1 + xSide;
    } else {
        y1Cell = y1 - xSide;
    }
    xDist = (double) (x1 - x2);
    xCell = (double) (x1Cell - x2);
    xSq = pow(xDist, 2);
    xSqCell = pow(xCell, 2);
    if (xSqCell < xSq) {
        xSq = xSqCell;
    }
    yDist = (double) (y1 - y2);
    yCell = (double) (y1Cell - y2);
    ySq = pow(yDist, 2);
    ySqCell = pow(yCell, 2);
    if (ySqCell < ySq) {
        ySq = ySqCell;
    }

    sqSum = xSq + ySq;
    dist = pow(sqSum, 0.5);
    return dist;

}




void get_pep_com(bool is_upper, float * xValCOM, float *yValCOM, float *zValCOM, ReadPdb * pdb, ReadDcd * dcd) {
    int atom_num;
    float xVal, yVal, zVal;
    int num_pep_atoms;

    xVal = 0.0;
    yVal = 0.0;
    zVal = 0.0;
    std::vector<int>::iterator first_it, end_it, it;
    if (is_upper == true) {
        first_it = std::begin(pdb->pep_up_heavy_atoms);
        end_it = std::end(pdb->pep_up_heavy_atoms);
    } else {
        first_it = std::begin(pdb->pep_low_heavy_atoms);
        end_it = std::end(pdb->pep_low_heavy_atoms);
    }
    for (it = first_it; it != end_it; ++it) {
        atom_num=*it;
        xVal += dcd->xValArray[0][atom_num-1];
        yVal += dcd->yValArray[0][atom_num-1];
        zVal += dcd->zValArray[0][atom_num-1];
    }
    num_pep_atoms = pdb->pep_low_heavy_atoms.size();
    *xValCOM = xVal / num_pep_atoms;
    *yValCOM = yVal / num_pep_atoms;
    *zValCOM = zVal / num_pep_atoms;
}


void cnt_water(bool isUp, ReadPdb * pdb, ReadDcd * dcd, float xSide, float xCOM, float yCOM) {

    std::vector<int>::iterator first_it, end_it, it;
    int atom_num;
    float xVal, yVal, zVal;
    double xyDist;
    int zBucket;

    first_it = std::begin(pdb->water_heavy_atoms);
    end_it = std::end(pdb->water_heavy_atoms);

    for (it = first_it; it != end_it; ++it) {
        atom_num=*it;
        zVal = dcd->zValArray[0][atom_num-1];
        if (isUp == false) {
            zVal = -zVal;
        }
        if ((zVal < MAX_Z_VAL) && (zVal > MIN_Z_VAL)) {
            xVal = dcd->xValArray[0][atom_num-1];
            yVal = dcd->yValArray[0][atom_num-1];
            zBucket = int(roundf(zVal));
            xyDist = get_dist_xy(xVal, xCOM, yVal, yCOM, xSide);
            if (xyDist <= PROXIMAL_DIST) {
                water_cnt[0][zBucket] += 1;
            } else if (xyDist <= FAR_DIST) {
                water_cnt[1][zBucket] += 1;
            } else {
                water_cnt[2][zBucket] += 1;
            }
        }
    }
}


            
int main(int argc, char * argv[]) {

    char dcd_file[MAX_FILENAME_LEN];
    char pdb_file[MAX_FILENAME_LEN];
    FILE * fp;
    FILE * fpUp;
    FILE * fpLo;
    bool isUp;

    int i, j, k;
    int aa;
    float xSide, ySide, zSide;

    int index = 0;

    float * distArray;
    float pep_xy_dist;
    int pepNum;
    ReadPdb::PEP_STRUCT_PTR pep_ptr;
//    arg[1] = Input_dir
//    arg[2] = tr
//    arg[3] = rep
//    arg[4] = prefix
//    arg[5] = first_str
//    arg[6] = last_str
//    arg[7] = CM list file
//    arg[8] = output_name

    for (i=0; i<10; i++) {
        printf("ARG %d = (%s)\n",i, argv[i]);
    }
    fflush(stdout);

    snprintf(pdb_file, MAX_FILENAME_LEN, "%s/spt.pdb", argv[1]);

    ReadPdb pdb(pdb_file);
    pdb.find_heavy_atoms(NUM_LIPID_LEAFLET, NUM_PEP_LEAFLET);
    printf("PDB %s\n", pdb_file);

    int first_step, last_step, step;
    first_step = atoi(argv[5]);
    last_step = atoi(argv[6]);
    ReadDcd * dcd;

    fp = fopen(argv[7],"w");
    fpUp = fopen(argv[8],"w");
    fpLo = fopen(argv[9],"w");
    printf("OPEN COM (%p) %s\n", fp, argv[7]);
    printf("OPEN UP (%p) %s\n", fpUp, argv[8]);
    printf("OPEN LO (%p) %s\n", fpLo, argv[9]);

    for (step=first_step; step < (last_step+1); step++) {
        for (i=0; i<NUM_Z_BUCK; i++) {
            for (j=0; j<NUM_REG; j++) {
                water_cnt[j][i] = 0;
            }
        }
        if ((step % 100) == 0) {
            printf("STEP %d\n", step);
        }
        snprintf(dcd_file, MAX_FILENAME_LEN, "%s/output/%s%s_%05d_%s.dcd", argv[1], argv[4], argv[2], step, argv[3]);
        dcd = new ReadDcd(dcd_file);
//        printf("READ DCD %s\n", dcd_file);
        fflush(stdout);

#ifdef AA_COM
        for (i=0; i<NUM_PEP; i++) {
            printf("GET COM %d\n", i);
            fflush(stdout);
            pepCOM[i] = new PepCOM(&pdb, dcd, i);
            printf("PEP COM %d %p\n", i, pepCOM[i]);
            fflush(stdout);
        }

        printf("DONE CENTER OF MASS %p\n", pepCOM[0]);
        fflush(stdout);
#endif /* AA_COM */


        xSide = dcd->boxArray[index]->ibox[0];
        ySide = dcd->boxArray[index]->ibox[2];
        zSide = dcd->boxArray[index]->ibox[5];

        float xUValCOM, yUValCOM, zUValCOM;
        float xLValCOM, yLValCOM, zLValCOM;

        get_pep_com(true, &xUValCOM, &yUValCOM, &zUValCOM, &pdb, dcd);
        get_pep_com(false, &xLValCOM, &yLValCOM, &zLValCOM, &pdb, dcd);
        fprintf(fp, "%f %f %f %f %f %f %f %f %f\n", xUValCOM, yUValCOM, zUValCOM, xLValCOM, yLValCOM, zLValCOM, xSide, ySide, zSide);

        isUp = true;
        cnt_water(isUp, &pdb, dcd, xSide, xUValCOM, yUValCOM);
        fprintf(fpUp, "%d", water_cnt[0][0]);
        for (j=0; j<NUM_REG; j++) {
            water_cnt[j][0] = 0;
        }
        for (i=1; i< NUM_Z_BUCK; i++) {
            fprintf(fpUp, " %d", water_cnt[0][i]);
            for (j=0; j<NUM_REG; j++) {
                water_cnt[j][i] = 0;
            }
        }
        fprintf(fpUp, "\n");

        isUp = false;
        cnt_water(isUp, &pdb, dcd, xSide, xUValCOM, yUValCOM);
        fprintf(fpLo, "%d", water_cnt[0][0]);
        for (j=0; j<NUM_REG; j++) {
            water_cnt[j][0] = 0;
        }
        for (i=1; i< NUM_Z_BUCK; i++) {
            fprintf(fpLo, " %d", water_cnt[0][i]);
            for (j=0; j<NUM_REG; j++) {
                water_cnt[j][i] = 0;
            }
        }
        fprintf(fpLo, "\n");

        delete dcd;
    }

    fclose(fpUp);
    fclose(fpLo);

    fclose(fp);

    printf("DONE get_water_cont - %s\n", argv[7]);

}

