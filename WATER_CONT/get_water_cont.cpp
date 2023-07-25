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


#define NUM_X_BUCK 71
#define NUM_Y_BUCK 71
#define NUM_Z_BUCK 71
#define X_NEG_VALS   35
#define Y_NEG_VALS   35
#define Z_NEG_VALS   35
#define ATOM_CONT_DIST 4.5

std::vector <int> water_pos[NUM_X_BUCK][NUM_Y_BUCK][NUM_Z_BUCK];
int water_cont_aa_up[NUM_AA];
int water_cont_total_up[NUM_Z_BUCK];
int water_cont_aa_lo[NUM_AA];
int water_cont_total_lo[NUM_Z_BUCK];

double get_dist(float x1, float x2, float y1, float y2, float z1, float z2) {
    double xDist;
    double yDist;
    double zDist;
    double sqSum;
    double dist;

    xDist = (double) (x1 - x2);
    yDist = (double) (y1 - y2);
    zDist = (double) (z1 - z2);

    sqSum = pow(xDist, 2) + pow(yDist, 2) + pow(zDist, 2);
    dist = pow(sqSum, 0.5);
    return dist;

}



void get_water_array(ReadPdb * pdb, ReadDcd * dcd, float sideLen) {
    int atom_num;
    int zBuck, yBuck, xBuck;
    float xVal, yVal, zVal;
    float xValCenter, yValCenter, zValCenter;
    int xShift, yShift;

    std::vector<int>::iterator first_it, end_it, it;

    first_it = std::begin(pdb->water_heavy_atoms);
    end_it = std::end(pdb->water_heavy_atoms);

    for (it = first_it; it != end_it; ++it) {
        atom_num=*it;
        xValCenter = dcd->xValArray[0][atom_num-1];
        yValCenter = dcd->yValArray[0][atom_num-1];
        zVal = dcd->zValArray[0][atom_num-1];
        if (zVal < 0) {
            zBuck = (int(zVal - 0.5)) + Z_NEG_VALS;
        } else {
            zBuck = (int(zVal + 0.5)) + Z_NEG_VALS;
        }
        if ((zBuck > (NUM_Z_BUCK-1)) || (zBuck < 0)) {
            continue;
        }
        for (xShift=-1; xShift <=1; xShift++) {
            xVal = xValCenter + (sideLen * xShift);
            if (xVal < 0) {
                xBuck = (int(xVal - 0.5)) + X_NEG_VALS;
            } else {
                xBuck = (int(xVal + 0.5)) + X_NEG_VALS;
            }
            if ((xBuck > (NUM_X_BUCK-1)) || (xBuck < 0)) {
                continue;
            }

            for (yShift=-1; yShift <=1; yShift++) {
                yVal = yValCenter + (sideLen * yShift);

                if (yVal < 0) {
                    yBuck = (int(yVal - 0.5)) + Y_NEG_VALS;
                } else {
                    yBuck = (int(yVal + 0.5)) + Y_NEG_VALS;
                }
                if ((yBuck > (NUM_Y_BUCK-1)) || (yBuck < 0)) {
                    continue;
                }

                water_pos[xBuck][yBuck][zBuck].push_back(atom_num);
            }
        }
    }
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

void get_cont_vals(ReadPdb * pdb, ReadDcd * dcd, PepCOM * pepCOMUp_ptr, PepCOM * pepCOMLo_ptr) {

    int pepNum;
    ReadPdb::PEP_STRUCT_PTR pep_ptr;
    int aa;
    int atom_num, wat_atom_num;
    float xValPep, yValPep, zValPep;
    float xValWat, yValWat, zValWat;
    float zVal;
    int firstX, firstY, firstZ, lastX, lastY, lastZ;
    int xOff, yOff, zOff;
    float dist;
    std::vector<int>::iterator first_it, end_it, it;

    std::vector<int>::iterator wat_first_it, wat_end_it, wat_it;

    std::set<int> wat_cont_set_up;
    std::set<int> wat_aa_cont_set_up[NUM_AA];
    std::set<int> wat_cont_set_lo;
    std::set<int> wat_aa_cont_set_lo[NUM_AA];


    int bucketVal;

    wat_cont_set_up.erase(std::begin(wat_cont_set_up), std::end(wat_cont_set_up));
    wat_cont_set_lo.erase(std::begin(wat_cont_set_lo), std::end(wat_cont_set_lo));
    for (aa=0; aa<NUM_AA; aa++) {
        wat_aa_cont_set_up[aa].erase(std::begin(wat_aa_cont_set_up[aa]), std::end(wat_aa_cont_set_up[aa]));
        wat_aa_cont_set_lo[aa].erase(std::begin(wat_aa_cont_set_lo[aa]), std::end(wat_aa_cont_set_lo[aa]));
    }

    for (pepNum=0; pepNum<MAX_PEP; pepNum++) {
        pep_ptr = &(pdb->pep_struct[pepNum]);
        if ((pep_ptr->pep_heavy_atoms).size() == 0) {
//            printf("DO BREAK %d\n", pepNum);
            break;
        }
        for (aa=0; aa<NUM_AA; aa++) {
            first_it = std::begin(pep_ptr->pep_ha[aa]);
            end_it = std::end(pep_ptr->pep_ha[aa]);
            for (it = first_it; it != end_it; ++it) {
                atom_num=*it;
                xValPep = dcd->xValArray[0][atom_num-1];
                yValPep = dcd->yValArray[0][atom_num-1];
                zValPep = dcd->zValArray[0][atom_num-1];

                firstX = (int(xValPep) + X_NEG_VALS) - 4;
                lastX =  (int(xValPep) + X_NEG_VALS) + 4;
                firstY = (int(yValPep) + Y_NEG_VALS) - 4;
                lastY =  (int(yValPep) + Y_NEG_VALS) + 4;
                firstZ = (int(zValPep) + Z_NEG_VALS) - 4;
                lastZ =  (int(zValPep) + Z_NEG_VALS) + 4;
//                printf("PEP %d %f %d %d \n", pepNum, zValPep, firstZ, lastZ);
                if (firstX < 0) {
                    firstX = 0;
                } else if (lastX > (NUM_X_BUCK-1)) {
                    lastX = (NUM_X_BUCK-1);
                }
                if (firstY < 0) {
                    firstY = 0;
                } else if (lastY > (NUM_Y_BUCK-1)) {
                    lastY = (NUM_Y_BUCK-1);
                }
                if (firstZ < 0) {
                    firstZ = 0;
                } else if (firstZ > (NUM_Z_BUCK-1)) {
                    lastZ = (NUM_Z_BUCK-1);
                }
                for (xOff=firstX; xOff<=lastX; xOff++) {
                    for (yOff=firstY; yOff<=lastY; yOff++) {
                        for (zOff=firstZ; zOff<=lastZ; zOff++) {
                            wat_first_it = std::begin(water_pos[xOff][yOff][zOff]);
                            wat_end_it = std::end(water_pos[xOff][yOff][zOff]);
                            for (wat_it = wat_first_it; wat_it != wat_end_it; ++wat_it) {
                                wat_atom_num=*wat_it;
                                xValWat = dcd->xValArray[0][wat_atom_num-1];
                                yValWat = dcd->yValArray[0][wat_atom_num-1];
                                zValWat = dcd->zValArray[0][wat_atom_num-1];
                                dist = get_dist(xValPep, xValWat, yValPep, yValWat, zValPep, zValWat);
                                if (dist < ATOM_CONT_DIST) {
                                    if (pepNum < NUM_PEP_LEAFLET) {
                                        wat_cont_set_up.insert(wat_atom_num);
                                        wat_aa_cont_set_up[aa].insert(wat_atom_num);
                                    } else {
//                                        printf("ADD WAT LO CONT %d\n", wat_atom_num);
                                        wat_cont_set_lo.insert(wat_atom_num);
                                        wat_aa_cont_set_lo[aa].insert(wat_atom_num);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

//    std::set<int> wat_cont_set_up;
//    std::set<int> wat_aa_cont_set_up[NUM_AA];
//    std::set<int> wat_cont_set_lo;
//    std::set<int> wat_aa_cont_set_lo[NUM_AA];

    std::set<int>::iterator first_set_it, end_set_it, set_it;


    first_set_it = std::begin(wat_cont_set_up);
    end_set_it = std::end(wat_cont_set_up);
    for (set_it = first_set_it; set_it != end_set_it; ++set_it) {
        wat_atom_num = *set_it;
        zValWat = dcd->zValArray[0][wat_atom_num-1];
        zValWat += Z_NEG_VALS;
        bucketVal = int(zValWat);
        if (bucketVal < 0) {
            printf("Z VAL < MIN %d %f\n", bucketVal, dcd->zValArray[0][wat_atom_num-1]);
            bucketVal = 0;
        } else if (bucketVal  >= NUM_Z_BUCK) {
            printf("Z VAL > MAX %d %f\n", bucketVal, dcd->zValArray[0][wat_atom_num-1]);
            bucketVal = NUM_Z_BUCK-1;
        }

        water_cont_total_up[bucketVal]++;
        printf("ADD WAT UP %d %f\n", bucketVal, dcd->zValArray[0][wat_atom_num-1]);
    }

    first_set_it = std::begin(wat_cont_set_lo);
    end_set_it = std::end(wat_cont_set_lo);
    for (set_it = first_set_it; set_it != end_set_it; ++set_it) {
        wat_atom_num = *set_it;
        zValWat = dcd->zValArray[0][wat_atom_num-1];
        zValWat = -zValWat;
        zValWat += Z_NEG_VALS;
        bucketVal = int(zValWat);
        if (bucketVal < 0) {
            printf("Z VAL < MIN %d %f\n", bucketVal, dcd->zValArray[0][wat_atom_num-1]);
            bucketVal = 0;
        } else if (bucketVal  >= NUM_Z_BUCK) {
            printf("Z VAL > MAX %d %f\n", bucketVal, dcd->zValArray[0][wat_atom_num-1]);
            bucketVal = NUM_Z_BUCK-1;
        }
        water_cont_total_lo[bucketVal]++;
//        printf("ADD WAT LO %d %f\n", bucketVal, dcd->zValArray[0][wat_atom_num-1]);
    }

    for (aa=0; aa<NUM_AA; aa++) {
        water_cont_aa_up[aa] = 0;
        first_set_it = std::begin(wat_aa_cont_set_up[aa]);
        end_set_it = std::end(wat_aa_cont_set_up[aa]);
        for (set_it = first_set_it; set_it != end_set_it; ++set_it) {
            water_cont_aa_up[aa]++;
        }
    }


    for (aa=0; aa<NUM_AA; aa++) {
        water_cont_aa_lo[aa] = 0;
        first_set_it = std::begin(wat_aa_cont_set_lo[aa]);
        end_set_it = std::end(wat_aa_cont_set_lo[aa]);
        for (set_it = first_set_it; set_it != end_set_it; ++set_it) {
            water_cont_aa_lo[aa]++;
        }
    }

}

            
int main(int argc, char * argv[]) {

    char dcd_file[MAX_FILENAME_LEN];
    char pdb_file[MAX_FILENAME_LEN];
    char aa_cont_file[MAX_FILENAME_LEN];
    char pep_cont_file[MAX_FILENAME_LEN];
    FILE * fp;
    FILE * fpPepUp;
    FILE * fpPepLo;
    FILE * fpAa[NUM_AA];
    int i, j, k;
    int aa;
    float xSide, ySide, zSide;

    int index = 0;

    float * distArray;
    int num_dmpc_lipids;
    int num_dmpg_lipids;
    float pep_xy_dist;
    int pepNum;
    ReadPdb::PEP_STRUCT_PTR pep_ptr;
    PepCOM * pepCOM[NUM_PEP];
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
    pdb.find_heavy_atoms(NUM_LIPID_LEAFLET, NUM_PEP_LEAFLET);
    printf("PDB %s\n", pdb_file);

    int first_step, last_step, step;
    first_step = atoi(argv[5]);
    last_step = atoi(argv[6]);
    ReadDcd * dcd;

    printf("OPEN %s\n", argv[7]);
    fp = fopen(argv[7],"w");

    snprintf(pep_cont_file, MAX_FILENAME_LEN, "%s_up.dat", argv[8]);
    fpPepUp = fopen(pep_cont_file,"w");
    snprintf(pep_cont_file, MAX_FILENAME_LEN, "%s_lo.dat", argv[8]);
    fpPepLo = fopen(pep_cont_file,"w");
    for (aa=0; aa<NUM_AA; aa++) {
        snprintf(aa_cont_file, MAX_FILENAME_LEN, "%s%d.dat", argv[9], (aa+1));
        fpAa[aa] = fopen(aa_cont_file,"w");
    }
    for (step=first_step; step < (last_step+1); step++) {
        for (i=0; i<NUM_Z_BUCK; i++) {
            water_cont_total_up[i] = 0;
            water_cont_total_lo[i] = 0;
            for (j=0; j<NUM_Y_BUCK; j++) {
                for (k=0; k<NUM_X_BUCK; k++) {
                    water_pos[k][j][i].clear();
                }
            }
        }
        for (aa=0; aa<NUM_AA; aa++) {
             water_cont_aa_up[aa] = 0;
             water_cont_aa_lo[aa] = 0;
        }
        if ((step % 100) == 0) {
            printf("STEP %d\n", step);
        }
        snprintf(dcd_file, MAX_FILENAME_LEN, "%s/output/%s%s_%05d_%s.dcd", argv[1], argv[4], argv[2], step, argv[3]);
        dcd = new ReadDcd(dcd_file);

        for (i=0; i<NUM_PEP; i++) {
            printf("GET COM %d\n", i);
            fflush(stdout);
            pepCOM[i] = new PepCOM(&pdb, dcd, i);
            printf("PEP COM %d %p\n", i, pepCOM[i]);
            fflush(stdout);
        }

        printf("DONE CENTER OF MASS %p\n", pepCOM[0]);
        fflush(stdout);
        xSide = dcd->boxArray[index]->ibox[0];
        ySide = dcd->boxArray[index]->ibox[2];
        zSide = dcd->boxArray[index]->ibox[5];

        float xUValCOM, yUValCOM, zUValCOM;
        float xLValCOM, yLValCOM, zLValCOM;

        get_pep_com(true, &xUValCOM, &yUValCOM, &zUValCOM, &pdb, dcd);
        get_pep_com(false, &xLValCOM, &yLValCOM, &zLValCOM, &pdb, dcd);
        fprintf(fp, "%f %f %f %f %f %f %f %f %f\n", xUValCOM, yUValCOM, zUValCOM, xLValCOM, yLValCOM, zLValCOM, xSide, ySide, zSide);

        get_water_array(&pdb, dcd, xSide);

/*   ASSUMES one pep per leaflet. Change for DIMER */
        get_cont_vals(&pdb, dcd, pepCOM[0], pepCOM[1]);

/*
*   Write out the results for the DCD and zero the values for the next DCD.
*/
        fprintf(fpPepUp, "%d", water_cont_total_up[0]);
        water_cont_total_up[0] = 0;
        for (i=1; i< NUM_Z_BUCK; i++) {
            fprintf(fpPepUp, " %d", water_cont_total_up[i]);
            water_cont_total_up[i] = 0;
        }
        fprintf(fpPepUp, "\n");

        for (aa=0; aa<NUM_AA; aa++) {
            fprintf(fpAa[aa], "%d %d\n", water_cont_aa_up[aa], water_cont_aa_lo[aa]);
        }


        fprintf(fpPepLo, "%d", water_cont_total_lo[0]);
        water_cont_total_lo[0] = 0;
        for (i=1; i< NUM_Z_BUCK; i++) {
            fprintf(fpPepLo, " %d", water_cont_total_lo[i]);
            water_cont_total_lo[i] = 0;
        }
        fprintf(fpPepLo, "\n");

        delete dcd;
    }

    fclose(fpPepUp);
    fclose(fpPepLo);
    for (aa=0; aa<NUM_AA; aa++) {
        fclose(fpAa[aa]);
    }

    fclose(fp);

    printf("DONE get_water_cont - %s\n", argv[7]);

}

