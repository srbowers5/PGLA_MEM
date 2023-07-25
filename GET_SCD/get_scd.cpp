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
#include "get_scd.h"

int Num_rem;

float uCell_xy;
float uCell_z;
float avg_uCell_xy;
float avg_uCell_z;
float side_xy;
float side_z;
int num_steps;


double chain_scd[NUM_CHAINS][NUM_C_TAIL];
int chain_scd_cnt[NUM_CHAINS][NUM_C_TAIL];

double step_chain_scd[NUM_CHAINS][NUM_C_TAIL];
int step_chain_scd_cnt[NUM_CHAINS][NUM_C_TAIL];

double chain_scd_r[NUM_CHAINS][NUM_C_TAIL][MAX_XY_DIST];
int chain_scd_cnt_r[NUM_CHAINS][NUM_C_TAIL][MAX_XY_DIST];
double pc_chain_scd_r[NUM_CHAINS][NUM_C_TAIL][MAX_XY_DIST];
int pc_chain_scd_cnt_r[NUM_CHAINS][NUM_C_TAIL][MAX_XY_DIST];
double pg_chain_scd_r[NUM_CHAINS][NUM_C_TAIL][MAX_XY_DIST];
int pg_chain_scd_cnt_r[NUM_CHAINS][NUM_C_TAIL][MAX_XY_DIST];

double scd_r[NUM_CHAINS][MAX_XY_DIST];
int scd_r_cnt[NUM_CHAINS][MAX_XY_DIST];
double pc_scd_r[NUM_CHAINS][MAX_XY_DIST];
int pc_scd_r_cnt[NUM_CHAINS][MAX_XY_DIST];
double pg_scd_r[NUM_CHAINS][MAX_XY_DIST];
int pg_scd_r_cnt[NUM_CHAINS][MAX_XY_DIST];



double chain_scd_reg[NUM_CHAINS][NUM_C_TAIL][NUM_REG];
int chain_scd_cnt_reg[NUM_CHAINS][NUM_C_TAIL][NUM_REG];

double chain_scd_reg_tot[NUM_CHAINS][NUM_REG];
int chain_scd_cnt_reg_tot[NUM_CHAINS][NUM_REG];


double pc_chain_scd[NUM_CHAINS][NUM_C_TAIL];
int pc_chain_scd_cnt[NUM_CHAINS][NUM_C_TAIL];

double pc_step_chain_scd[NUM_CHAINS][NUM_C_TAIL];
int pc_step_chain_scd_cnt[NUM_CHAINS][NUM_C_TAIL];


double pc_chain_scd_reg[NUM_CHAINS][NUM_C_TAIL][NUM_REG];
int pc_chain_scd_cnt_reg[NUM_CHAINS][NUM_C_TAIL][NUM_REG];

double pc_chain_scd_reg_tot[NUM_CHAINS][NUM_REG];
int pc_chain_scd_cnt_reg_tot[NUM_CHAINS][NUM_REG];
double pc_chain_scd_reg_tot_step[NUM_CHAINS][NUM_REG];
int pc_chain_scd_cnt_reg_tot_step[NUM_CHAINS][NUM_REG];


double pg_chain_scd[NUM_CHAINS][NUM_C_TAIL];
int pg_chain_scd_cnt[NUM_CHAINS][NUM_C_TAIL];


double pg_step_chain_scd[NUM_CHAINS][NUM_C_TAIL];
int pg_step_chain_scd_cnt[NUM_CHAINS][NUM_C_TAIL];


double pg_step_chain_scd_near[NUM_CHAINS][NUM_C_TAIL];
int pg_step_chain_scd_near_cnt[NUM_CHAINS][NUM_C_TAIL];
double pg_step_chain_scd_far[NUM_CHAINS][NUM_C_TAIL];
int pg_step_chain_scd_far_cnt[NUM_CHAINS][NUM_C_TAIL];

double pc_step_chain_scd_near[NUM_CHAINS][NUM_C_TAIL];
int pc_step_chain_scd_near_cnt[NUM_CHAINS][NUM_C_TAIL];
double pc_step_chain_scd_far[NUM_CHAINS][NUM_C_TAIL];
int pc_step_chain_scd_far_cnt[NUM_CHAINS][NUM_C_TAIL];

double pg_chain_scd_reg[NUM_CHAINS][NUM_C_TAIL][NUM_REG];
int pg_chain_scd_cnt_reg[NUM_CHAINS][NUM_C_TAIL][NUM_REG];

double pg_chain_scd_reg_tot[NUM_CHAINS][NUM_REG];
int pg_chain_scd_cnt_reg_tot[NUM_CHAINS][NUM_REG];
double pg_chain_scd_reg_tot_step[NUM_CHAINS][NUM_REG];
int pg_chain_scd_cnt_reg_tot_step[NUM_CHAINS][NUM_REG];


int scd_cnt_r[NUM_C_TAIL][MAX_XY_DIST];
int pc_scd_cnt_r[NUM_C_TAIL][MAX_XY_DIST];
int pg_scd_cnt_r[NUM_C_TAIL][MAX_XY_DIST];


double scd_r_OPP[NUM_CHAINS][MAX_XY_DIST];
int scd_cnt_r_OPP[NUM_CHAINS][MAX_XY_DIST];
double pc_scd_r_OPP[NUM_CHAINS][MAX_XY_DIST];
int pc_scd_cnt_r_OPP[NUM_CHAINS][MAX_XY_DIST];
double pg_scd_r_OPP[NUM_CHAINS][MAX_XY_DIST];
int pg_scd_cnt_r_OPP[NUM_CHAINS][MAX_XY_DIST];









//
//   SCD = ((3 * cosSq) - 1) / 2;
//
//   cos = x / ( (x**2 + y**2) **0.5)
//   cos**2 = z**2 / (x**2 + y**2 + z**2)

double get_one_scd(float xC, float yC, float zC, float xH, float yH, float zH) {
    float xDiffSq, yDiffSq, zDiffSq;
    float distSq;
    float cosSq, scd;

    xDiffSq = pow( (xC - xH), 2);
    yDiffSq = pow( (yC - yH), 2);
    zDiffSq = pow( (zC - zH), 2);
    distSq = xDiffSq + yDiffSq + zDiffSq;
    cosSq = zDiffSq / distSq;

    scd = ((3 * cosSq) - 1)/ 2;

    return scd;
}


float get_min_xy_dist(float pepX, float pepY, float lipX, float lipY, float uCellXY) {

    float xDiff, yDiff;
    float uCellHalf;
    float distSq;
    float myDist;

    uCellHalf = uCellXY/2;

    xDiff = pepX - lipX;
    yDiff = pepY - lipY;

//    printf("GET MIN %f %f %f\n", xDiff, yDiff, uCellXY);
    if (xDiff > uCellHalf) {
        xDiff = xDiff - uCellXY;
    } else if (xDiff < -uCellHalf) {
        xDiff = xDiff + uCellXY;
    }

    if (yDiff > uCellHalf) {
        yDiff = yDiff - uCellXY;
    } else if (yDiff < -uCellHalf) {
        yDiff = yDiff + uCellXY;
    }
   // printf("GET MOD MIN %f %f %f\n", xDiff, yDiff, uCellXY);

    distSq = (pow(xDiff,2) + pow(yDiff,2));
    myDist = pow(distSq, 0.5);

   // printf("GET MOD MIN %f %f %f %f %f\n", xDiff, yDiff, uCellXY, distSq, myDist);

    return myDist;
}

void write_scd_vals_OPP(ReadPdb::LIPID_C_TAIL_PTR lip_c_ptr, double scdValue, int  lipReg, int lipDistInt) {

//    printf("IN write_scd_vals_OPP %f %d\n", scdValue, lipDistInt);
    if (lipDistInt < MAX_XY_DIST) {
        if (lip_c_ptr->mon_type == ReadPdb::IS_DMPC) {
            pc_scd_r_OPP[lip_c_ptr->chain_num-1][lipDistInt] += scdValue;
            pc_scd_cnt_r_OPP[lip_c_ptr->chain_num-1][lipDistInt]++;
        } else if (lip_c_ptr->mon_type == ReadPdb::IS_DMPG) {
            pg_scd_r_OPP[lip_c_ptr->chain_num-1][lipDistInt] += scdValue;
            pg_scd_cnt_r_OPP[lip_c_ptr->chain_num-1][lipDistInt]++;
        }
        scd_r_OPP[lip_c_ptr->chain_num-1][lipDistInt] += scdValue;
        scd_cnt_r_OPP[lip_c_ptr->chain_num-1][lipDistInt]++;
    }

}
void write_scd_vals(ReadPdb::LIPID_C_TAIL_PTR lip_c_ptr, double scdValue, int  lipReg, int lipDistInt) {

    fflush(stdout);
    if (lip_c_ptr->mon_type == ReadPdb::IS_DMPC) {
        pc_chain_scd[lip_c_ptr->chain_num-1][lip_c_ptr->chain_pos-1] += scdValue;
        pc_chain_scd_cnt[lip_c_ptr->chain_num-1][lip_c_ptr->chain_pos-1]++;
        pc_step_chain_scd[lip_c_ptr->chain_num-1][lip_c_ptr->chain_pos-1] += scdValue;
        pc_step_chain_scd_cnt[lip_c_ptr->chain_num-1][lip_c_ptr->chain_pos-1]++;
        if (lipReg == 0) {
            pc_step_chain_scd_near[lip_c_ptr->chain_num-1][lip_c_ptr->chain_pos-1] += scdValue;
            pc_step_chain_scd_near_cnt[lip_c_ptr->chain_num-1][lip_c_ptr->chain_pos-1]++;
        } else if (lipReg == 2) {
            pc_step_chain_scd_far[lip_c_ptr->chain_num-1][lip_c_ptr->chain_pos-1] += scdValue;
            pc_step_chain_scd_far_cnt[lip_c_ptr->chain_num-1][lip_c_ptr->chain_pos-1]++;
        }

        if (lipDistInt < MAX_XY_DIST) {
            pc_chain_scd_r[lip_c_ptr->chain_num-1][lip_c_ptr->chain_pos-1][lipDistInt] += scdValue;
            pc_chain_scd_cnt_r[lip_c_ptr->chain_num-1][lip_c_ptr->chain_pos-1][lipDistInt]++;
            pc_scd_r[lip_c_ptr->chain_num-1][lipDistInt] += scdValue;
            pc_scd_cnt_r[lip_c_ptr->chain_num-1][lipDistInt]++;
        }
        pc_chain_scd_reg[lip_c_ptr->chain_num-1][lip_c_ptr->chain_pos-1][lipReg] += scdValue;
        pc_chain_scd_cnt_reg[lip_c_ptr->chain_num-1][lip_c_ptr->chain_pos-1][lipReg]++;
        pc_chain_scd_reg_tot[lip_c_ptr->chain_num-1][lipReg] += scdValue;
        pc_chain_scd_cnt_reg_tot[lip_c_ptr->chain_num-1][lipReg]++;
        pc_chain_scd_reg_tot_step[lip_c_ptr->chain_num-1][lipReg] += scdValue;
        pc_chain_scd_cnt_reg_tot_step[lip_c_ptr->chain_num-1][lipReg]++;

    } else if (lip_c_ptr->mon_type == ReadPdb::IS_DMPG) {
        pg_chain_scd[lip_c_ptr->chain_num-1][lip_c_ptr->chain_pos-1] += scdValue;
        pg_chain_scd_cnt[lip_c_ptr->chain_num-1][lip_c_ptr->chain_pos-1]++;
        pg_step_chain_scd[lip_c_ptr->chain_num-1][lip_c_ptr->chain_pos-1] += scdValue;
        pg_step_chain_scd_cnt[lip_c_ptr->chain_num-1][lip_c_ptr->chain_pos-1]++;
        if (lipReg == 0) {
            pg_step_chain_scd_near[lip_c_ptr->chain_num-1][lip_c_ptr->chain_pos-1] += scdValue;
            pg_step_chain_scd_near_cnt[lip_c_ptr->chain_num-1][lip_c_ptr->chain_pos-1]++;
        } else if (lipReg == 2) {
            pg_step_chain_scd_far[lip_c_ptr->chain_num-1][lip_c_ptr->chain_pos-1] += scdValue;
            pg_step_chain_scd_far_cnt[lip_c_ptr->chain_num-1][lip_c_ptr->chain_pos-1]++;
        }
        if (lipDistInt < MAX_XY_DIST) {
            pg_chain_scd_r[lip_c_ptr->chain_num-1][lip_c_ptr->chain_pos-1][lipDistInt] += scdValue;
            pg_chain_scd_cnt_r[lip_c_ptr->chain_num-1][lip_c_ptr->chain_pos-1][lipDistInt]++;
            pg_scd_r[lip_c_ptr->chain_num-1][lipDistInt] += scdValue;
            pg_scd_cnt_r[lip_c_ptr->chain_num-1][lipDistInt]++;
        }
        pg_chain_scd_reg[lip_c_ptr->chain_num-1][lip_c_ptr->chain_pos-1][lipReg] += scdValue;
        pg_chain_scd_cnt_reg[lip_c_ptr->chain_num-1][lip_c_ptr->chain_pos-1][lipReg]++;
        pg_chain_scd_reg_tot[lip_c_ptr->chain_num-1][lipReg] += scdValue;
        pg_chain_scd_cnt_reg_tot[lip_c_ptr->chain_num-1][lipReg]++;
        pg_chain_scd_reg_tot_step[lip_c_ptr->chain_num-1][lipReg] += scdValue;
        pg_chain_scd_cnt_reg_tot_step[lip_c_ptr->chain_num-1][lipReg]++;
    }
    chain_scd[lip_c_ptr->chain_num-1][lip_c_ptr->chain_pos-1] += scdValue;
    chain_scd_cnt[lip_c_ptr->chain_num-1][lip_c_ptr->chain_pos-1]++;
    step_chain_scd[lip_c_ptr->chain_num-1][lip_c_ptr->chain_pos-1] += scdValue;
    step_chain_scd_cnt[lip_c_ptr->chain_num-1][lip_c_ptr->chain_pos-1]++;

    if (lipDistInt < MAX_XY_DIST) {
        chain_scd_r[lip_c_ptr->chain_num-1][lip_c_ptr->chain_pos-1][lipDistInt] += scdValue;
        chain_scd_cnt_r[lip_c_ptr->chain_num-1][lip_c_ptr->chain_pos-1][lipDistInt]++;
        scd_r[lip_c_ptr->chain_num-1][lipDistInt] += scdValue;
        scd_cnt_r[lip_c_ptr->chain_num-1][lipDistInt]++;
    }
    chain_scd_reg[lip_c_ptr->chain_num-1][lip_c_ptr->chain_pos-1][lipReg] += scdValue;
    chain_scd_cnt_reg[lip_c_ptr->chain_num-1][lip_c_ptr->chain_pos-1][lipReg]++;
    chain_scd_reg_tot[lip_c_ptr->chain_num-1][lipReg] += scdValue;
    chain_scd_cnt_reg_tot[lip_c_ptr->chain_num-1][lipReg]++;
}


void get_scd_leaf(bool isUpper, bool isOpp, float pepX, float pepY, float pepZ, ReadPdb * pdb, ReadDcd * dcd) {

    ReadPdb::LIPID_C_TAIL_PTR lip_c_ptr;
    float xC, yC, zC, xH, yH, zH;
    double scdVal;
    float lipDist;
    int lipDistInt;
    float uCellXY;
    int  lipReg;
    int  i,j;

    uCellXY = dcd->uCell[0];
    std::vector<ReadPdb::LIPID_C_TAIL_PTR>::iterator first_it, end_it, it;
    if (isUpper) {
        first_it = std::begin(pdb->lipid_tail_c_h_up);
        end_it = std::end(pdb->lipid_tail_c_h_up);
    } else {
        first_it = std::begin(pdb->lipid_tail_c_h_low);
        end_it = std::end(pdb->lipid_tail_c_h_low);
    }
   // printf("Leaf LOOP\n"); 
    fflush(stdout);
    int cnt = 0;
    for (it = first_it; it != end_it; ++it) {
        cnt++;
        lip_c_ptr = *it;
        if (lip_c_ptr->hxAtom_num != 0) {
            xC = dcd->xValArray[0][lip_c_ptr->cAtom_num-1];
            yC = dcd->yValArray[0][lip_c_ptr->cAtom_num-1];
            zC = dcd->zValArray[0][lip_c_ptr->cAtom_num-1];

            xH = dcd->xValArray[0][lip_c_ptr->hxAtom_num-1];
            yH = dcd->yValArray[0][lip_c_ptr->hxAtom_num-1];
            zH = dcd->zValArray[0][lip_c_ptr->hxAtom_num-1];

            lipDist = get_min_xy_dist(pepX, pepY, xC, yC, uCellXY);
            lipDistInt = int(lipDist);
            if (lipDist < PROXIMAL_DIST) {
                lipReg = 0;
            } else if (lipDist < FAR_DIST) {
                lipReg = 1;
            } else {
                lipReg = 2;
            }
            scdVal = get_one_scd(xC, yC, zC, xH, yH, zH);
            if (isOpp == false) {
                write_scd_vals(lip_c_ptr, scdVal, lipReg, lipDistInt);
            } else {
                write_scd_vals_OPP(lip_c_ptr, scdVal, lipReg, lipDistInt);
            }

            if (lip_c_ptr->hyAtom_num != 0) {
                xH = dcd->xValArray[0][lip_c_ptr->hyAtom_num-1];
                yH = dcd->yValArray[0][lip_c_ptr->hyAtom_num-1];
                zH = dcd->zValArray[0][lip_c_ptr->hyAtom_num-1];
                scdVal = get_one_scd(xC, yC, zC, xH, yH, zH);
                if (isOpp == false) {
                    write_scd_vals(lip_c_ptr, scdVal, lipReg, lipDistInt);
                } else {
 //printf("Write SCD OPP\n"); 
 fflush(stdout);
                    write_scd_vals_OPP(lip_c_ptr, scdVal, lipReg, lipDistInt);
                }

                if (lip_c_ptr->hzAtom_num != 0) {
                    xH = dcd->xValArray[0][lip_c_ptr->hzAtom_num-1];
                    yH = dcd->yValArray[0][lip_c_ptr->hzAtom_num-1];
                    zH = dcd->zValArray[0][lip_c_ptr->hzAtom_num-1];
   // printf("GET3 one SCD\n"); 
    fflush(stdout);
                    scdVal = get_one_scd(xC, yC, zC, xH, yH, zH);
                    if (isOpp == false) {
   // printf("Write SCD\n"); 
    fflush(stdout);
                        write_scd_vals(lip_c_ptr, scdVal, lipReg, lipDistInt);
                    } else {
   // printf("Write SCD OPP\n"); 
    fflush(stdout);
                        write_scd_vals_OPP(lip_c_ptr, scdVal, lipReg, lipDistInt);
                    }
                }
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
    if (num_pep_atoms == 0) {
        *xValCOM = 0.0;
        *yValCOM = 0.0;
        *zValCOM = 0.0;
    } else {
        *xValCOM = xVal / num_pep_atoms;
        *yValCOM = yVal / num_pep_atoms;
        *zValCOM = zVal / num_pep_atoms;
    }
}

int main(int argc, char * argv[]) {

    char dcd_file[MAX_FILENAME_LEN];
    char pdb_file[MAX_FILENAME_LEN];
    char cm_file[MAX_FILENAME_LEN];
    char outFile[MAX_FILENAME_LEN];
    FILE *fpOut;
    FILE *fpScdStep;
    FILE *fpPGScdStep;
    FILE *fpPCScdStep;
    FILE *fpPGScdRegStep;
    FILE *fpPCScdRegStep;
    FILE *fpPGScdNearStep;
    FILE *fpPGScdFarStep;
    FILE *fpPCScdNearStep;
    FILE *fpPCScdFarStep;
    int i, j, k;
    float uCell[3];
    ReadPdb::LIPID_C_TAIL_PTR lip_c_ptr;
    float xC, yC, zC, xH, yH, zH;
    double scdVal;


//    arg[1] = Input_dir
//    arg[2] = tr
//    arg[3] = rep
//    arg[4] = prefix
//    arg[5] = first_str
//    arg[6] = last_str

    snprintf(pdb_file, MAX_FILENAME_LEN, "%s/spt.pdb", argv[1]);


    ReadPdb pdb(pdb_file);
    pdb.find_lipid_tails(NUM_LIPID_LEAFLET);
    printf("PDB %s\n", pdb_file);
    pdb.find_heavy_atoms(NUM_LIPID_LEAFLET);


    int first_step, last_step, step;
    first_step = atoi(argv[5]);
    last_step = atoi(argv[6]);
    ReadDcd * dcd;


    for (i=0; i< NUM_CHAINS; i++) {
        for (j=0; j< NUM_C_TAIL; j++) {
            chain_scd[i][j] = 0.0;
            chain_scd_cnt[i][j] = 0;
            pc_chain_scd[i][j] = 0.0;
            pc_chain_scd_cnt[i][j] = 0;
            pg_chain_scd[i][j] = 0.0;
            pg_chain_scd_cnt[i][j] = 0;

            for (k=0; k<MAX_XY_DIST; k++) {
                scd_r[i][k] = 0.0;
                scd_cnt_r[i][k] = 0;
                pc_scd_r[i][k] = 0.0;
                pc_scd_cnt_r[i][k] = 0;
                pg_scd_r[i][k] = 0.0;
                pg_scd_cnt_r[i][k] = 0;

                chain_scd_r[i][j][k] = 0.0;
                chain_scd_cnt_r[i][j][k] = 0;
                pc_chain_scd_r[i][j][k] = 0.0;
                pc_chain_scd_cnt_r[i][j][k] = 0;
                pg_chain_scd_r[i][j][k] = 0.0;
                pg_chain_scd_cnt_r[i][j][k] = 0;

                scd_r_OPP[i][k] = 0.0;
                scd_cnt_r_OPP[i][k] = 0;
                pc_scd_r_OPP[i][k] = 0.0;
                pc_scd_cnt_r_OPP[i][k] = 0;
                pg_scd_r_OPP[i][k] = 0.0;
                pg_scd_cnt_r_OPP[i][k] = 0;
            }

            for (k=0; k<NUM_REG; k++) {
                chain_scd_reg[i][j][k] = 0.0;
                chain_scd_cnt_reg[i][j][k] = 0;
                pc_chain_scd_reg[i][j][k] = 0.0;
                pc_chain_scd_cnt_reg[i][j][k] = 0;
                pg_chain_scd_reg[i][j][k] = 0.0;
                pg_chain_scd_cnt_reg[i][j][k] = 0;

                chain_scd_reg_tot[i][k] = 0.0;
                chain_scd_cnt_reg_tot[i][k] = 0;
                pc_chain_scd_reg_tot[i][k] = 0.0;
                pc_chain_scd_cnt_reg_tot[i][k] = 0;
                pc_chain_scd_reg_tot_step[i][k] = 0.0;
                pc_chain_scd_cnt_reg_tot_step[i][k] = 0;
                pg_chain_scd_reg_tot[i][k] = 0.0;
                pg_chain_scd_cnt_reg_tot[i][k] = 0;
                pg_chain_scd_reg_tot_step[i][k] = 0.0;
                pg_chain_scd_cnt_reg_tot_step[i][k] = 0;
            }
        }
    }
    snprintf(outFile, MAX_FILENAME_LEN, "%s/SCD_DIR/pg_scd_step_%s_%s_%s-%s.dat", argv[1], argv[2], argv[3], argv[5], argv[6]);
    fpPGScdStep = fopen(outFile,"w");
    snprintf(outFile, MAX_FILENAME_LEN, "%s/SCD_DIR/pc_scd_step_%s_%s_%s-%s.dat", argv[1], argv[2], argv[3], argv[5], argv[6]);
    fpPCScdStep = fopen(outFile,"w");
    snprintf(outFile, MAX_FILENAME_LEN, "%s/SCD_DIR/scd_step_%s_%s_%s-%s.dat", argv[1], argv[2], argv[3], argv[5], argv[6]);
    fpScdStep = fopen(outFile,"w");
    snprintf(outFile, MAX_FILENAME_LEN, "%s/SCD_DIR/pg_scd_reg_step_%s_%s_%s-%s.dat", argv[1], argv[2], argv[3], argv[5], argv[6]);
    fpPGScdRegStep = fopen(outFile,"w");
    snprintf(outFile, MAX_FILENAME_LEN, "%s/SCD_DIR/pc_scd_reg_step_%s_%s_%s-%s.dat", argv[1], argv[2], argv[3], argv[5], argv[6]);
    fpPCScdRegStep = fopen(outFile,"w");

    snprintf(outFile, MAX_FILENAME_LEN, "%s/SCD_DIR/pc_scd_sn2_near_step_%s_%s_%s-%s.dat", argv[1], argv[2], argv[3], argv[5], argv[6]);
    fpPCScdNearStep = fopen(outFile,"w");
    snprintf(outFile, MAX_FILENAME_LEN, "%s/SCD_DIR/pc_scd_sn2_far_step_%s_%s_%s-%s.dat", argv[1], argv[2], argv[3], argv[5], argv[6]);
    fpPCScdFarStep = fopen(outFile,"w");
    snprintf(outFile, MAX_FILENAME_LEN, "%s/SCD_DIR/pg_scd_sn2_near_step_%s_%s_%s-%s.dat", argv[1], argv[2], argv[3], argv[5], argv[6]);
    fpPGScdNearStep = fopen(outFile,"w");
    snprintf(outFile, MAX_FILENAME_LEN, "%s/SCD_DIR/pg_scd_sn2_far_step_%s_%s_%s-%s.dat", argv[1], argv[2], argv[3], argv[5], argv[6]);
    fpPGScdFarStep = fopen(outFile,"w");

   // printf("START LOOP\n");
    fflush(stdout);
    for (step=first_step; step < (last_step+1); step++) {
        snprintf(dcd_file, MAX_FILENAME_LEN, "%s/output/%s%s_%05d_%s.dcd", argv[1], argv[4], argv[2], step, argv[3]);
        dcd = new ReadDcd(dcd_file);


        dcd->uCell[0] = dcd->boxArray[0]->ibox[0];
        dcd->uCell[1] = dcd->boxArray[0]->ibox[2];
        dcd->uCell[2] = dcd->boxArray[0]->ibox[5];
        uCell_xy += dcd->uCell[0];
        uCell_z += dcd->uCell[2];

        float xUValCOM, yUValCOM, zUValCOM;
        float xLValCOM, yLValCOM, zLValCOM;


        get_pep_com(true, &xUValCOM, &yUValCOM, &zUValCOM, &pdb, dcd);
        get_pep_com(false, &xLValCOM, &yLValCOM, &zLValCOM, &pdb, dcd);


        for (i=0; i<NUM_CHAINS; i++) {
            for (j=0; j<NUM_C_TAIL; j++) {
                step_chain_scd[i][j] = 0.0;
                step_chain_scd_cnt[i][j] = 0;
                pg_step_chain_scd[i][j] = 0.0;
                pg_step_chain_scd_cnt[i][j] = 0;
                pc_step_chain_scd[i][j] = 0.0;
                pc_step_chain_scd_cnt[i][j] = 0;
                pg_step_chain_scd_near[i][j] = 0.0;
                pg_step_chain_scd_near_cnt[i][j] = 0;
                pg_step_chain_scd_far[i][j] = 0.0;
                pg_step_chain_scd_far_cnt[i][j] = 0;
                pc_step_chain_scd_near[i][j] = 0.0;
                pc_step_chain_scd_near_cnt[i][j] = 0;
                pc_step_chain_scd_far[i][j] = 0.0;
                pc_step_chain_scd_far_cnt[i][j] = 0;
            }
        }
        for (i=0; i<NUM_CHAINS; i++) {
            for (j=0; j<NUM_REG; j++) {
                pg_chain_scd_reg_tot_step[i][j] = 0.0;
                pc_chain_scd_reg_tot_step[i][j] = 0.0;
                pg_chain_scd_cnt_reg_tot_step[i][j] = 0;
                pc_chain_scd_cnt_reg_tot_step[i][j] = 0;
            }
        }
//  Upper leaf
//  Upper leaf OPP
//  Lower leaf
//  Lower leaf OPP
        get_scd_leaf(true, false, xUValCOM, yUValCOM, zUValCOM, &pdb, dcd);
        get_scd_leaf(true, true, xLValCOM, yLValCOM, zLValCOM, &pdb, dcd);
        get_scd_leaf(false, false, xLValCOM, yLValCOM, zLValCOM, &pdb, dcd);
        get_scd_leaf(false, true, xUValCOM, yUValCOM, zUValCOM, &pdb, dcd);
        for (i=0; i<NUM_CHAINS; i++) {
            for (j=1; j<NUM_C_TAIL; j++) {
                fprintf(fpScdStep, "%.5f ", step_chain_scd[i][j]/step_chain_scd_cnt[i][j]);
                fprintf(fpPGScdStep, "%.5f ", pg_step_chain_scd[i][j]/pg_step_chain_scd_cnt[i][j]);
                fprintf(fpPCScdStep, "%.5f ", pc_step_chain_scd[i][j]/pc_step_chain_scd_cnt[i][j]);
                fprintf(fpPGScdNearStep, "%.5f ", pg_step_chain_scd_near[i][j]/ pg_step_chain_scd_near_cnt[i][j]);
                fprintf(fpPCScdNearStep, "%.5f ", pc_step_chain_scd_near[i][j]/ pc_step_chain_scd_near_cnt[i][j]);
                fprintf(fpPGScdFarStep, "%.5f ", pg_step_chain_scd_far[i][j]/ pg_step_chain_scd_far_cnt[i][j]);
                fprintf(fpPCScdFarStep, "%.5f ", pc_step_chain_scd_far[i][j]/ pc_step_chain_scd_far_cnt[i][j]);
            }
        }
        fprintf(fpScdStep, "\n");
        fprintf(fpPGScdStep, "\n");
        fprintf(fpPCScdStep, "\n");
        fprintf(fpPGScdNearStep, "\n");
        fprintf(fpPGScdFarStep, "\n");
        fprintf(fpPCScdNearStep, "\n");
        fprintf(fpPCScdFarStep, "\n");
        for (i=0; i<NUM_CHAINS; i++) {
            for (j=0; j<NUM_REG; j++) {
                if (pg_chain_scd_cnt_reg_tot_step[i][j]) {
                    fprintf(fpPGScdRegStep, "%.5f ", pg_chain_scd_reg_tot_step[i][j]/pg_chain_scd_cnt_reg_tot_step[i][j]);
                } else {
                    fprintf(fpPGScdRegStep, "%.5f ", pg_chain_scd_reg_tot_step[i][j+1]/pg_chain_scd_cnt_reg_tot_step[i][j+1]);
                }
                if (pc_chain_scd_cnt_reg_tot_step[i][j]) {
                    fprintf(fpPCScdRegStep, "%.5f ", pc_chain_scd_reg_tot_step[i][j]/pc_chain_scd_cnt_reg_tot_step[i][j]);
                } else {
                    fprintf(fpPCScdRegStep, "%.5f ", pc_chain_scd_reg_tot_step[i][j+1]/pc_chain_scd_cnt_reg_tot_step[i][j+1]);
                }
            }
        }
        fprintf(fpPGScdRegStep, "\n");
        fprintf(fpPCScdRegStep, "\n");

        delete dcd;
    }
    fclose(fpScdStep);
    fclose(fpPGScdStep);
    fclose(fpPCScdStep);
    fclose(fpPGScdRegStep);
    fclose(fpPCScdRegStep);
    fclose(fpPGScdNearStep);
    fclose(fpPGScdFarStep);
    fclose(fpPCScdNearStep);
    fclose(fpPCScdFarStep);

   // printf("Write SCD_\n");
    fflush(stdout);
    snprintf(outFile, MAX_FILENAME_LEN, "%s/SCD_DIR/scd_%s_%s_%s-%s.dat", argv[1], argv[2], argv[3], argv[5], argv[6]);
    fpOut = fopen(outFile,"w");
    for (i=1; i<NUM_C_TAIL; i++) {
        fprintf(fpOut, "%d %f %f %f %f %f %f\n", i+1, chain_scd[0][i]/chain_scd_cnt[0][i], chain_scd[1][i]/chain_scd_cnt[1][i],
          pc_chain_scd[0][i]/pc_chain_scd_cnt[0][i], pc_chain_scd[1][i]/pc_chain_scd_cnt[1][i],
          pg_chain_scd[0][i]/pg_chain_scd_cnt[0][i], pg_chain_scd[1][i]/pg_chain_scd_cnt[1][i]);
    }
    fclose(fpOut);
    printf("Wrote %s\n", outFile);

    float reg_tot[2][3];
    float pc_reg_tot[2][3];
    float pg_reg_tot[2][3];
    int reg_tot_cnt[2][3];
    int pc_reg_tot_cnt[2][3];
    int pg_reg_tot_cnt[2][3];
    float reg_tot_all[3];
    float pc_reg_tot_all[3];
    float pg_reg_tot_all[3];
    int reg_tot_cnt_all[3];
    int pc_reg_tot_cnt_all[3];
    int pg_reg_tot_cnt_all[3];
    for (i=0; i<2; i++) {
        for (j=0; j<3; j++) {
            reg_tot[i][j] = 0.0;
            reg_tot_cnt[i][j] = 0;
            pc_reg_tot[i][j] = 0.0;
            pc_reg_tot_cnt[i][j] = 0;
            pg_reg_tot[i][j] = 0.0;
            pg_reg_tot_cnt[i][j] = 0;
        }
    }
    snprintf(outFile, MAX_FILENAME_LEN, "%s/SCD_DIR/scd_reg_%s_%s_%s-%s.dat", argv[1], argv[2], argv[3], argv[5], argv[6]);
    fpOut = fopen(outFile,"w");
    for (i=1; i<NUM_C_TAIL; i++) {
        fprintf(fpOut,"%d", i+1);
        for (j=0; j<NUM_REG; j++) {
            fprintf(fpOut, " %f %f %f %f %f %f", chain_scd_reg[0][i][j]/chain_scd_cnt_reg[0][i][j], chain_scd_reg[1][i][j]/chain_scd_cnt_reg[1][i][j],
              pc_chain_scd_reg[0][i][j]/pc_chain_scd_cnt_reg[0][i][j], pc_chain_scd_reg[1][i][j]/pc_chain_scd_cnt_reg[1][i][j],
              pg_chain_scd_reg[0][i][j]/pg_chain_scd_cnt_reg[0][i][j], pg_chain_scd_reg[1][i][j]/pg_chain_scd_cnt_reg[1][i][j]);
            reg_tot[0][j] += chain_scd_reg[0][i][j];
            reg_tot[1][j] += chain_scd_reg[1][i][j];
            reg_tot_all[j] += (chain_scd_reg[0][i][j] + chain_scd_reg[1][i][j]);

            reg_tot_cnt[0][j] += chain_scd_cnt_reg[0][i][j];
            reg_tot_cnt[1][j] += chain_scd_cnt_reg[1][i][j];
            reg_tot_cnt_all[j] += (chain_scd_cnt_reg[0][i][j] + chain_scd_cnt_reg[1][i][j]);

            pc_reg_tot[0][j] += pc_chain_scd_reg[0][i][j];
            pc_reg_tot[1][j] += pc_chain_scd_reg[1][i][j];
            pc_reg_tot_all[j] += (pc_chain_scd_reg[0][i][j] + pc_chain_scd_reg[1][i][j]);

            pc_reg_tot_cnt[0][j] += pc_chain_scd_cnt_reg[0][i][j];
            pc_reg_tot_cnt[1][j] += pc_chain_scd_cnt_reg[1][i][j];
            pc_reg_tot_cnt_all[j] += (pc_chain_scd_cnt_reg[0][i][j] + pc_chain_scd_cnt_reg[1][i][j]);

            pg_reg_tot[0][j] += pg_chain_scd_reg[0][i][j];
            pg_reg_tot[1][j] += pg_chain_scd_reg[1][i][j];
            pg_reg_tot_all[j] += (pg_chain_scd_reg[0][i][j] + pg_chain_scd_reg[1][i][j]);

            pg_reg_tot_cnt[0][j] += pg_chain_scd_cnt_reg[0][i][j];
            pg_reg_tot_cnt[1][j] += pg_chain_scd_cnt_reg[1][i][j];
            pg_reg_tot_cnt_all[j] += (pg_chain_scd_cnt_reg[0][i][j] + pg_chain_scd_cnt_reg[1][i][j]);

        }
        fprintf(fpOut,"\n");
    }
    fclose(fpOut);
    printf("Wrote %s\n", outFile);

    snprintf(outFile, MAX_FILENAME_LEN, "%s/SCD_DIR/scd_reg_tot_%s_%s_%s-%s.dat", argv[1], argv[2], argv[3], argv[5], argv[6]);
    fpOut = fopen(outFile,"w");
    for (j=0; j<NUM_REG; j++) {
        fprintf(fpOut, "%d %f %f %f %f %f  %f %f %f %f %f  %f %f %f %f %f\n", j, 
          chain_scd_reg_tot[0][j]/chain_scd_cnt_reg_tot[0][j], 
          pc_chain_scd_reg_tot[0][j]/pc_chain_scd_cnt_reg_tot[0][j], 
          pg_chain_scd_reg_tot[0][j]/pg_chain_scd_cnt_reg_tot[0][j], 
          chain_scd_reg_tot[1][j]/chain_scd_cnt_reg_tot[1][j],
          pc_chain_scd_reg_tot[1][j]/pc_chain_scd_cnt_reg_tot[1][j],
          pg_chain_scd_reg_tot[1][j]/pg_chain_scd_cnt_reg_tot[1][j],
          reg_tot[0][j]/reg_tot_cnt[0][j],
          pc_reg_tot[0][j]/pc_reg_tot_cnt[0][j],
          pg_reg_tot[0][j]/pg_reg_tot_cnt[0][j],
          reg_tot[1][j]/reg_tot_cnt[1][j],
          pc_reg_tot[1][j]/pc_reg_tot_cnt[1][j],
          pg_reg_tot[1][j]/pg_reg_tot_cnt[1][j],
          reg_tot_all[j]/reg_tot_cnt_all[j],
          pc_reg_tot_all[j]/pc_reg_tot_cnt_all[j],
          pg_reg_tot_all[j]/pg_reg_tot_cnt_all[j]);
    }
    fclose(fpOut);
    printf("Wrote %s\n", outFile);




    snprintf(outFile, MAX_FILENAME_LEN, "%s/SCD_DIR/scd_r_sn1_%s_%s_%s-%s.dat", argv[1], argv[2], argv[3], argv[5], argv[6]);
    fpOut = fopen(outFile,"w");
    for (i=1; i<NUM_C_TAIL; i++) {
        for (j=0; j<MAX_XY_DIST; j++) {
           // printf("%d %d %f %d\n", i, j, chain_scd_r[i][j],chain_scd_cnt_r[i][j]);
            fprintf(fpOut, "%f ", chain_scd_r[1][i][j]/chain_scd_cnt_r[1][i][j]);
        }
        fprintf(fpOut,"\n");
    }
    fclose(fpOut);
    printf("Wrote %s\n", outFile);


    snprintf(outFile, MAX_FILENAME_LEN, "%s/SCD_DIR/scd_r_sn2_%s_%s_%s-%s.dat", argv[1], argv[2], argv[3], argv[5], argv[6]);
    fpOut = fopen(outFile,"w");
    for (i=1; i<NUM_C_TAIL; i++) {
        for (j=0; j<MAX_XY_DIST; j++) {
           // printf("%d %d %f %d\n", i, j, chain_scd_r[i][j],chain_scd_cnt_r[i][j]);
            fprintf(fpOut, "%f ", chain_scd_r[0][i][j]/chain_scd_cnt_r[0][i][j]);
        }
        fprintf(fpOut,"\n");
    }
    fclose(fpOut);
    printf("Wrote %s\n", outFile);

    for (i=0; i<NUM_CHAINS; i++) {
       for (j=0; j<MAX_XY_DIST; j++) {
           scd_r[i][j] = 0.0;
           scd_r_cnt[i][j] = 0;
           pc_scd_r[i][j] = 0.0;
           pc_scd_r_cnt[i][j] = 0;
           pg_scd_r[i][j] = 0.0;
           pg_scd_r_cnt[i][j] = 0;
           for (k=0; k<NUM_C_TAIL; k++) {
               scd_r[i][j] += chain_scd_r[i][k][j];
               scd_r_cnt[i][j] += chain_scd_cnt_r[i][k][j];
               pc_scd_r[i][j] += pc_chain_scd_r[i][k][j];
               pc_scd_r_cnt[i][j] += pc_chain_scd_cnt_r[i][k][j];
               pg_scd_r[i][j] += pg_chain_scd_r[i][k][j];
               pg_scd_r_cnt[i][j] += pg_chain_scd_cnt_r[i][k][j];
           }
       }
   }



    snprintf(outFile, MAX_FILENAME_LEN, "%s/SCD_DIR/pc_scd_r_sn2_%s_%s_%s-%s.dat", 
      argv[1], argv[2], argv[3], argv[5], argv[6]);
    fpOut = fopen(outFile,"w");
    for (j=0; j<MAX_XY_DIST; j++) {
        fprintf(fpOut, "%f %f\n", (j+0.5), pc_scd_r[0][j]/pc_scd_r_cnt[0][j]);
    }
    fclose(fpOut);
    printf("Wrote %s\n", outFile);


    snprintf(outFile, MAX_FILENAME_LEN, "%s/SCD_DIR/pg_scd_r_sn2_%s_%s_%s-%s.dat",
      argv[1], argv[2], argv[3], argv[5], argv[6]);
    fpOut = fopen(outFile,"w");
    for (j=0; j<MAX_XY_DIST; j++) {
        fprintf(fpOut, "%f %f\n", (j+0.5), pg_scd_r[0][j]/pg_scd_r_cnt[0][j]);
    }
    fclose(fpOut);
    printf("Wrote %s\n", outFile);





    snprintf(outFile, MAX_FILENAME_LEN, "%s/SCD_DIR/pc_scd_r_sn2_OPP_%s_%s_%s-%s.dat",
      argv[1], argv[2], argv[3], argv[5], argv[6]);
    fpOut = fopen(outFile,"w");
    for (j=0; j<MAX_XY_DIST; j++) {
        fprintf(fpOut, "%f %f\n", (j+0.5), pc_scd_r_OPP[0][j]/pc_scd_cnt_r_OPP[0][j]);
    }
    fclose(fpOut);
    printf("Wrote %s\n", outFile);


    snprintf(outFile, MAX_FILENAME_LEN, "%s/SCD_DIR/pg_scd_r_sn2_OPP_%s_%s_%s-%s.dat",
      argv[1], argv[2], argv[3], argv[5], argv[6]);
    fpOut = fopen(outFile,"w");
    for (j=0; j<MAX_XY_DIST; j++) {
        fprintf(fpOut, "%f %f\n", (j+0.5), pg_scd_r_OPP[0][j]/pg_scd_cnt_r_OPP[0][j]);
    }
    fclose(fpOut);
    printf("Wrote %s\n", outFile);




    printf("DONE get_scd\n");
}
