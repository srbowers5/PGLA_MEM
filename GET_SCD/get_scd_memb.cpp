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

double chain_scd_all[NUM_CHAINS];
int chain_scd_all_cnt[NUM_CHAINS];
double pc_chain_scd_all[NUM_CHAINS];
int pc_chain_scd_all_cnt[NUM_CHAINS];
double pg_chain_scd_all[NUM_CHAINS];
int pg_chain_scd_all_cnt[NUM_CHAINS];

double chain_scd[NUM_CHAINS][NUM_C_TAIL];
int chain_scd_cnt[NUM_CHAINS][NUM_C_TAIL];

double step_chain_scd[NUM_CHAINS][NUM_C_TAIL];
int step_chain_scd_cnt[NUM_CHAINS][NUM_C_TAIL];

double chain_scd_r[NUM_C_TAIL][MAX_XY_DIST];
int chain_scd_cnt_r[NUM_C_TAIL][MAX_XY_DIST];

double chain_scd_reg[NUM_CHAINS][NUM_C_TAIL][NUM_REG];
int chain_scd_cnt_reg[NUM_CHAINS][NUM_C_TAIL][NUM_REG];

double chain_scd_reg_tot[NUM_CHAINS][NUM_REG];
int chain_scd_cnt_reg_tot[NUM_CHAINS][NUM_REG];


double pc_chain_scd[NUM_CHAINS][NUM_C_TAIL];
int pc_chain_scd_cnt[NUM_CHAINS][NUM_C_TAIL];

double pc_step_chain_scd[NUM_CHAINS][NUM_C_TAIL];
int pc_step_chain_scd_cnt[NUM_CHAINS][NUM_C_TAIL];

double pc_chain_scd_r[NUM_C_TAIL][MAX_XY_DIST];
int pc_chain_scd_cnt_r[NUM_C_TAIL][MAX_XY_DIST];

double pc_chain_scd_reg[NUM_CHAINS][NUM_C_TAIL][NUM_REG];
int pc_chain_scd_cnt_reg[NUM_CHAINS][NUM_C_TAIL][NUM_REG];

double pc_chain_scd_reg_tot[NUM_CHAINS][NUM_REG];
int pc_chain_scd_cnt_reg_tot[NUM_CHAINS][NUM_REG];


double pg_chain_scd[NUM_CHAINS][NUM_C_TAIL];
int pg_chain_scd_cnt[NUM_CHAINS][NUM_C_TAIL];

double pg_step_chain_scd[NUM_CHAINS][NUM_C_TAIL];
int pg_step_chain_scd_cnt[NUM_CHAINS][NUM_C_TAIL];

double pg_chain_scd_r[NUM_C_TAIL][MAX_XY_DIST];
int pg_chain_scd_cnt_r[NUM_C_TAIL][MAX_XY_DIST];

double pg_chain_scd_reg[NUM_CHAINS][NUM_C_TAIL][NUM_REG];
int pg_chain_scd_cnt_reg[NUM_CHAINS][NUM_C_TAIL][NUM_REG];

double pg_chain_scd_reg_tot[NUM_CHAINS][NUM_REG];
int pg_chain_scd_cnt_reg_tot[NUM_CHAINS][NUM_REG];











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

void write_scd_vals(ReadPdb::LIPID_C_TAIL_PTR lip_c_ptr, double scdValue, int  lipReg, int lipDistInt) {

    if (lip_c_ptr->mon_type == ReadPdb::IS_DMPC) {
        pc_chain_scd_all[lip_c_ptr->chain_num-1] += scdValue;
        pc_chain_scd_all_cnt[lip_c_ptr->chain_num-1]++;
        pc_chain_scd[lip_c_ptr->chain_num-1][lip_c_ptr->chain_pos-1] += scdValue;
        pc_chain_scd_cnt[lip_c_ptr->chain_num-1][lip_c_ptr->chain_pos-1]++;
        pc_step_chain_scd[lip_c_ptr->chain_num-1][lip_c_ptr->chain_pos-1] += scdValue;
        pc_step_chain_scd_cnt[lip_c_ptr->chain_num-1][lip_c_ptr->chain_pos-1]++;

        if (lipDistInt < MAX_XY_DIST) {
            pc_chain_scd_r[lip_c_ptr->chain_pos-1][lipDistInt] += scdValue;
            pc_chain_scd_cnt_r[lip_c_ptr->chain_pos-1][lipDistInt]++;
        }
        pc_chain_scd_reg[lip_c_ptr->chain_num-1][lip_c_ptr->chain_pos-1][lipReg] += scdValue;
        pc_chain_scd_cnt_reg[lip_c_ptr->chain_num-1][lip_c_ptr->chain_pos-1][lipReg]++;
        pc_chain_scd_reg_tot[lip_c_ptr->chain_num-1][lipReg] += scdValue;
        pc_chain_scd_cnt_reg_tot[lip_c_ptr->chain_num-1][lipReg]++;
//        if (lip_c_ptr->chain_num == 1) {
//            printf("PC %f %f %d %f ", scdValue, pg_chain_scd_all[lip_c_ptr->chain_num-1], 
//              pg_chain_scd_all_cnt[lip_c_ptr->chain_num-1], 
//              pg_chain_scd_all[lip_c_ptr->chain_num-1]/pg_chain_scd_all_cnt[lip_c_ptr->chain_num-1]);
//        }

    } else if (lip_c_ptr->mon_type == ReadPdb::IS_DMPG) {
        pg_chain_scd_all[lip_c_ptr->chain_num-1] += scdValue;
        pg_chain_scd_all_cnt[lip_c_ptr->chain_num-1]++;
        pg_chain_scd[lip_c_ptr->chain_num-1][lip_c_ptr->chain_pos-1] += scdValue;
        pg_chain_scd_cnt[lip_c_ptr->chain_num-1][lip_c_ptr->chain_pos-1]++;
        pg_step_chain_scd[lip_c_ptr->chain_num-1][lip_c_ptr->chain_pos-1] += scdValue;
        pg_step_chain_scd_cnt[lip_c_ptr->chain_num-1][lip_c_ptr->chain_pos-1]++;
        if (lipDistInt < MAX_XY_DIST) {
            pg_chain_scd_r[lip_c_ptr->chain_pos-1][lipDistInt] += scdValue;
            pg_chain_scd_cnt_r[lip_c_ptr->chain_pos-1][lipDistInt]++;
        }
        pg_chain_scd_reg[lip_c_ptr->chain_num-1][lip_c_ptr->chain_pos-1][lipReg] += scdValue;
        pg_chain_scd_cnt_reg[lip_c_ptr->chain_num-1][lip_c_ptr->chain_pos-1][lipReg]++;
        pg_chain_scd_reg_tot[lip_c_ptr->chain_num-1][lipReg] += scdValue;
        pg_chain_scd_cnt_reg_tot[lip_c_ptr->chain_num-1][lipReg]++;
//        if (lip_c_ptr->chain_num == 1) {
 //           printf("PG %f %f %d %f ", scdValue, pg_chain_scd_all[lip_c_ptr->chain_num-1], 
  //            pg_chain_scd_all_cnt[lip_c_ptr->chain_num-1], 
   //           pg_chain_scd_all[lip_c_ptr->chain_num-1]);
//        }
    } else {
        printf("ERROR invalid lipid type %d\n", lip_c_ptr->mon_type);
    }
    chain_scd_all[lip_c_ptr->chain_num-1] += scdValue;
    chain_scd_all_cnt[lip_c_ptr->chain_num-1]++;
//    if (lip_c_ptr->chain_num == 1) {
//        printf("Tot %f %f %f %d %f\n", chain_scd_all[lip_c_ptr->chain_num-1], 
//          (pc_chain_scd_all[lip_c_ptr->chain_num-1] + pg_chain_scd_all[lip_c_ptr->chain_num-1]),
//          chain_scd_all[lip_c_ptr->chain_num-1] - (pc_chain_scd_all[lip_c_ptr->chain_num-1] + pg_chain_scd_all[lip_c_ptr->chain_num-1]),
//          chain_scd_all_cnt[lip_c_ptr->chain_num-1],
//          chain_scd_all[lip_c_ptr->chain_num-1]/chain_scd_all_cnt[lip_c_ptr->chain_num-1]); 
 //   }

    chain_scd[lip_c_ptr->chain_num-1][lip_c_ptr->chain_pos-1] += scdValue;
    chain_scd_cnt[lip_c_ptr->chain_num-1][lip_c_ptr->chain_pos-1]++;
    step_chain_scd[lip_c_ptr->chain_num-1][lip_c_ptr->chain_pos-1] += scdValue;
    step_chain_scd_cnt[lip_c_ptr->chain_num-1][lip_c_ptr->chain_pos-1]++;

    if (lipDistInt < MAX_XY_DIST) {
        chain_scd_r[lip_c_ptr->chain_pos-1][lipDistInt] += scdValue;
        chain_scd_cnt_r[lip_c_ptr->chain_pos-1][lipDistInt]++;
    }
    chain_scd_reg[lip_c_ptr->chain_num-1][lip_c_ptr->chain_pos-1][lipReg] += scdValue;
    chain_scd_cnt_reg[lip_c_ptr->chain_num-1][lip_c_ptr->chain_pos-1][lipReg]++;
    chain_scd_reg_tot[lip_c_ptr->chain_num-1][lipReg] += scdValue;
    chain_scd_cnt_reg_tot[lip_c_ptr->chain_num-1][lipReg]++;
}


void get_scd_leaf(bool isUpper, float pepX, float pepY, float pepZ, ReadPdb * pdb, ReadDcd * dcd, int frame) {

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
   
    for (it = first_it; it != end_it; ++it) {
        lip_c_ptr = *it;
        if (lip_c_ptr->hxAtom_num != 0) {
            xC = dcd->xValArray[frame][lip_c_ptr->cAtom_num-1];
            yC = dcd->yValArray[frame][lip_c_ptr->cAtom_num-1];
            zC = dcd->zValArray[frame][lip_c_ptr->cAtom_num-1];

            xH = dcd->xValArray[frame][lip_c_ptr->hxAtom_num-1];
            yH = dcd->yValArray[frame][lip_c_ptr->hxAtom_num-1];
            zH = dcd->zValArray[frame][lip_c_ptr->hxAtom_num-1];

//            printf("ATOM VAL frame %d (%f, %f, %f) (%f,%f %f)\n", frame, xC, yC, zC, xH, yH, zH);
            lipDist = get_min_xy_dist(pepX, pepY, xC, yC, uCellXY);
            lipDistInt = int(lipDist);
            if (lipDist < PROXIMAL_DIST) {
                lipReg = 0;
            } else if (lipDist < FAR_DIST) {
                lipReg = 1;
            } else {
                lipReg = 2;
            }

           // printf("(%f, %f), (%f, %f), %f, %f, %d, %d\n", pepX, pepY, xC, yC, uCellXY, lipDist, lipDistInt, lipReg);

            scdVal = get_one_scd(xC, yC, zC, xH, yH, zH);
            write_scd_vals(lip_c_ptr, scdVal, lipReg, lipDistInt);

            if (lip_c_ptr->hyAtom_num != 0) {
                xH = dcd->xValArray[frame][lip_c_ptr->hyAtom_num-1];
                yH = dcd->yValArray[frame][lip_c_ptr->hyAtom_num-1];
                zH = dcd->zValArray[frame][lip_c_ptr->hyAtom_num-1];
                scdVal = get_one_scd(xC, yC, zC, xH, yH, zH);
                write_scd_vals(lip_c_ptr, scdVal, lipReg, lipDistInt);

                if (lip_c_ptr->hzAtom_num != 0) {
                    xH = dcd->xValArray[frame][lip_c_ptr->hzAtom_num-1];
                    yH = dcd->yValArray[frame][lip_c_ptr->hzAtom_num-1];
                    zH = dcd->zValArray[frame][lip_c_ptr->hzAtom_num-1];
                    scdVal = get_one_scd(xC, yC, zC, xH, yH, zH);
                    write_scd_vals(lip_c_ptr, scdVal, lipReg, lipDistInt);
                }
            }
        }
    }
}



void get_pep_com(bool is_upper, float * xValCOM, float *yValCOM, float *zValCOM, ReadPdb * pdb, ReadDcd * dcd, int frame) {
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
        xVal += dcd->xValArray[frame][atom_num-1];
        yVal += dcd->yValArray[frame][atom_num-1];
        zVal += dcd->zValArray[frame][atom_num-1];
//        printf("PEP VAL frame %d (%f, %f, %f)\n", frame, xVal, yVal, zVal);
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
//    arg[7] = pdbFile


    snprintf(pdb_file, MAX_FILENAME_LEN, "%s", argv[7]);
    printf("PDB %s\n", pdb_file);

//    snprintf(pdb_file, MAX_FILENAME_LEN, "%s/spt.pdb", argv[1]);
    ReadPdb pdb(pdb_file);
    pdb.find_lipid_tails(NUM_LIPID_LEAFLET);
    printf("PDB %s\n", pdb_file);
    pdb.find_heavy_atoms(NUM_LIPID_LEAFLET);


    int first_step, last_step, step;
    first_step = atoi(argv[5]);
    last_step = atoi(argv[6]);
    ReadDcd * dcd;


    for (i=0; i< NUM_CHAINS; i++) {
        chain_scd_all[i] = 0.0;
        pc_chain_scd_all[i] = 0.0;
        pg_chain_scd_all[i] = 0.0;
        chain_scd_all_cnt[i] = 0;
        pc_chain_scd_all_cnt[i] = 0;
        pg_chain_scd_all_cnt[i] = 0;
        for (j=0; j< NUM_C_TAIL; j++) {
            chain_scd[i][j] = 0.0;
            chain_scd_cnt[i][j] = 0;
            pc_chain_scd[i][j] = 0.0;
            pc_chain_scd_cnt[i][j] = 0;
            pg_chain_scd[i][j] = 0.0;
            pg_chain_scd_cnt[i][j] = 0;

            for (k=0; k<MAX_XY_DIST; k++) {
                chain_scd_r[j][k] = 0.0;
                chain_scd_cnt_r[j][k] = 0;
                pc_chain_scd_r[j][k] = 0.0;
                pc_chain_scd_cnt_r[j][k] = 0;
                pg_chain_scd_r[j][k] = 0.0;
                pg_chain_scd_cnt_r[j][k] = 0;
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
                pg_chain_scd_reg_tot[i][k] = 0.0;
                pg_chain_scd_cnt_reg_tot[i][k] = 0;
            }
        }
    }
// One DCD file
    printf("READ %s/output/%s%s_%s.dcd", argv[1], argv[4], argv[2], argv[3]);
    snprintf(dcd_file, MAX_FILENAME_LEN, "%s/output/%s%s_%s.dcd", argv[1], argv[4], argv[2], argv[3]);
      dcd = new ReadDcd(dcd_file);

    if (dcd->nStr < last_step) {
        last_step = dcd->nStr;
    }
    printf("Frames %d to %d total NUM =     %d\n", first_step, last_step, (last_step + 1 - first_step));

    for (int frame=first_step-1; frame<last_step; frame++) {


        dcd->uCell[0] = dcd->boxArray[frame]->ibox[0];
        dcd->uCell[1] = dcd->boxArray[frame]->ibox[2];
        dcd->uCell[2] = dcd->boxArray[frame]->ibox[5];
        uCell_xy += dcd->uCell[0];
        uCell_z += dcd->uCell[2];

        float xUValCOM, yUValCOM, zUValCOM;
        float xLValCOM, yLValCOM, zLValCOM;


        get_pep_com(true, &xUValCOM, &yUValCOM, &zUValCOM, &pdb, dcd, frame);
        get_pep_com(false, &xLValCOM, &yLValCOM, &zLValCOM, &pdb, dcd, frame);


        for (i=0; i<NUM_CHAINS; i++) {
            for (j=0; j<NUM_C_TAIL; j++) {
                step_chain_scd[i][j] = 0.0;
                step_chain_scd_cnt[i][j] = 0;
            }
        }

        get_scd_leaf(true, xUValCOM, yUValCOM, zUValCOM, &pdb, dcd, frame);
        get_scd_leaf(false, xLValCOM, yLValCOM, zLValCOM, &pdb, dcd, frame);
    }
    delete dcd;

    snprintf(outFile, MAX_FILENAME_LEN, "%s/scd_all_%s_%s_%s-%s.dat", argv[1], argv[2], argv[3], argv[5], argv[6]);
    fpOut = fopen(outFile,"w");
//    printf("VAL %f %f %f , CNT %d %d %d\n", chain_scd_all[0]/chain_scd_all_cnt[0], 
//          pc_chain_scd_all[0]/pc_chain_scd_all_cnt[0], 
//          pg_chain_scd_all[0]/pg_chain_scd_all_cnt[0],
//          chain_scd_all_cnt[0], pc_chain_scd_all_cnt[0], pg_chain_scd_all_cnt[0]);
          
    fprintf(fpOut, "0 %f %f %f %f %f %f\n", chain_scd_all[0]/chain_scd_all_cnt[0], chain_scd_all[1]/chain_scd_all_cnt[1],
          pc_chain_scd_all[0]/pc_chain_scd_all_cnt[0], pc_chain_scd_all[1]/pc_chain_scd_all_cnt[1],
          pg_chain_scd_all[0]/pg_chain_scd_all_cnt[0], pg_chain_scd_all[1]/pg_chain_scd_all_cnt[1]);
    fclose(fpOut);
    printf("Wrote %s\n", outFile);


    snprintf(outFile, MAX_FILENAME_LEN, "%s/scd_%s_%s_%s-%s.dat", argv[1], argv[2], argv[3], argv[5], argv[6]);
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
    snprintf(outFile, MAX_FILENAME_LEN, "%s/scd_reg_%s_%s_%s-%s.dat", argv[1], argv[2], argv[3], argv[5], argv[6]);
    fpOut = fopen(outFile,"w");
    for (i=1; i<NUM_C_TAIL; i++) {
        fprintf(fpOut,"%d", i+1);
        for (j=0; j<NUM_REG; j++) {
//            printf("SCD REG 0 cnt %d %d %d\n", i, j, chain_scd_cnt_reg[0][i][j]);
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

    snprintf(outFile, MAX_FILENAME_LEN, "%s/scd_reg_tot_%s_%s_%s-%s.dat", argv[1], argv[2], argv[3], argv[5], argv[6]);
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




    snprintf(outFile, MAX_FILENAME_LEN, "%s/scd_r_%s_%s_%s-%s.dat", argv[1], argv[2], argv[3], argv[5], argv[6]);
    fpOut = fopen(outFile,"w");
    for (i=1; i<NUM_C_TAIL; i++) {
        for (j=0; j<MAX_XY_DIST; j++) {
           // printf("%d %d %f %d\n", i, j, chain_scd_r[i][j],chain_scd_cnt_r[i][j]);
            fprintf(fpOut, "%f ", chain_scd_r[i][j]/chain_scd_cnt_r[i][j]);
        }
        fprintf(fpOut,"\n");
    }
    fclose(fpOut);
    printf("Wrote %s\n", outFile);

    printf("DONE get_scd\n");
}
