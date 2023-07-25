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
#include "dcdStruct.h"
#include "ReadPdb.h"
#include "ReadDcd.h"
#include "get_com_xsc.h"

int Num_rem;



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


int main(int argc, char * argv[]) {

    char dcd_file[MAX_FILENAME_LEN];
    char pdb_file[MAX_FILENAME_LEN];
    FILE * fp;
    int i, j;
    float xSide, ySide, zSide;

    int index = 0;

    float * distArray;
    int num_dmpc_lipids;
    int num_dmpg_lipids;
    float pep_xy_dist;
//    arg[1] = Input_dir
//    arg[2] = tr
//    arg[3] = rep
//    arg[4] = prefix
//    arg[5] = first_str
//    arg[6] = last_str
//    arg[7] = CM list file
//    arg[8] = output_name

    snprintf(pdb_file, MAX_FILENAME_LEN, "%s/spt.pdb", argv[1]);


    ReadPdb pdb(pdb_file);
    pdb.find_heavy_atoms(NUM_LIPID_LEAFLET);
    printf("PDB %s\n", pdb_file);

    int first_step, last_step, step;
    first_step = atoi(argv[5]);
    last_step = atoi(argv[6]);
    ReadDcd * dcd;

    fp = fopen(argv[7],"w");
    for (step=first_step; step < (last_step+1); step++) {
        if ((step % 100) == 0) {
            printf("STEP %d\n", step);
        }
        snprintf(dcd_file, MAX_FILENAME_LEN, "%s/output/%s%s_%05d_%s.dcd", argv[1], argv[4], argv[2], step, argv[3]);
        dcd = new ReadDcd(dcd_file);


        xSide = dcd->boxArray[index]->ibox[0];
        ySide = dcd->boxArray[index]->ibox[2];
        zSide = dcd->boxArray[index]->ibox[5];

        float xUValCOM, yUValCOM, zUValCOM;
        float xLValCOM, yLValCOM, zLValCOM;

        get_pep_com(true, &xUValCOM, &yUValCOM, &zUValCOM, &pdb, dcd);
        get_pep_com(false, &xLValCOM, &yLValCOM, &zLValCOM, &pdb, dcd);
        fprintf(fp, "%f %f %f %f %f %f %f %f %f\n", xUValCOM, yUValCOM, zUValCOM, xLValCOM, yLValCOM, zLValCOM, xSide, ySide, zSide);
        delete dcd;
    }

    fclose(fp);
    printf("DONE get_com_xsc - %s\n", argv[7]);
}

