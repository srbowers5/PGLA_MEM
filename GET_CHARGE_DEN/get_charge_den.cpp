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
#include "ReadPsf.h"
#include "ReadCM.h"
#include "get_charge_den.h"

extern void get_areas(float xVal, float yVal, float * a1Val, float * a2Val, float * a3Val);
void get_Efield ();

#define LIP_CHO         0
#define LIP_P           1
#define LIP_GLY         2
#define LIP_TAIL        3
#define NUM_LIPID_SEC   4

int Num_rem;

#define E_FIELD_SIZE 30
#define MAX_FORCE_DIST   39
#define KVAL  8988000000   //  N⋅A2⋅C−2
#define C_CONVERT 1.60217646e-19 
#define Msq_to_Asq   1.0e+20


float tot_charge_reg[NUM_REG_R][NUM_REG_Z];
float tot_wat_charge_reg[NUM_REG_R][NUM_REG_Z];

float tot_PC_tail_charge[NUM_LIPID_SEC];
float tot_PG_tail_charge[NUM_LIPID_SEC];

float charge_reg[NUM_REG_R][NUM_REG_Z];
float wat_charge_reg[NUM_REG_R][NUM_REG_Z];
float lip_charge_reg[NUM_REG_R][NUM_REG_Z];
float pep_charge_reg[NUM_REG_R][NUM_REG_Z];
float ion_charge_reg[NUM_REG_R][NUM_REG_Z];
float chargeVol[NUM_REG_R][NUM_REG_Z];
float chargeArea[NUM_REG_R];

float PG_charge_reg_sec[NUM_REG_R][NUM_REG_Z][NUM_LIPID_SEC];
float PC_charge_reg_sec[NUM_REG_R][NUM_REG_Z][NUM_LIPID_SEC];

float PG_charge_reg_sec_step[NUM_REG_R][NUM_REG_Z][NUM_LIPID_SEC];
float PC_charge_reg_sec_step[NUM_REG_R][NUM_REG_Z][NUM_LIPID_SEC];

float PG_charge_reg_z_sec[NUM_REG_R][NUM_Z_BUCKETS][NUM_LIPID_SEC];
float PC_charge_reg_z_sec[NUM_REG_R][NUM_Z_BUCKETS][NUM_LIPID_SEC];

float net_E_Z_field[E_FIELD_SIZE][E_FIELD_SIZE];
float net_E_R_field[E_FIELD_SIZE][E_FIELD_SIZE];

float denAreas[UCELL_MAX_X][UCELL_MAX_Y][NUM_AREAS];
float uCell_charge[UCELL_MAX_X][UCELL_MAX_Y][UCELL_MAX_Z];
float uCell_OPP_charge[UCELL_MAX_X][UCELL_MAX_Y][UCELL_MAX_Z];

float tot_charge;
float lipid_tot_charge;
float ion_tot_charge;
float wat_tot_charge;
float pep_tot_charge;

float net_charge_dist[NUM_BUCKETS][NUM_Z_BUCKETS];
float neg_charge_dist[NUM_BUCKETS][NUM_Z_BUCKETS];
float pos_charge_dist[NUM_BUCKETS][NUM_Z_BUCKETS];


float pep_pos_charge_dist[NUM_BUCKETS][NUM_Z_BUCKETS];
float pep_neutral_charge_dist[NUM_BUCKETS][NUM_Z_BUCKETS];
float pep_side_charge_dist[NUM_BUCKETS][NUM_Z_BUCKETS];
float pep_bb_charge_dist[NUM_BUCKETS][NUM_Z_BUCKETS];
float lipid_net_charge_dist[NUM_BUCKETS][NUM_Z_BUCKETS];
float lipid_pos_charge_dist[NUM_BUCKETS][NUM_Z_BUCKETS];
float lipid_neg_charge_dist[NUM_BUCKETS][NUM_Z_BUCKETS];
float lipid_net_OPP_charge_dist[NUM_BUCKETS][NUM_Z_BUCKETS];
float lipid_pos_OPP_charge_dist[NUM_BUCKETS][NUM_Z_BUCKETS];
float lipid_neg_OPP_charge_dist[NUM_BUCKETS][NUM_Z_BUCKETS];
float pep_net_charge_dist[NUM_BUCKETS][NUM_Z_BUCKETS];
float pep_net_OPP_charge_dist[NUM_BUCKETS][NUM_Z_BUCKETS];
float ion_net_charge_dist[NUM_BUCKETS][NUM_Z_BUCKETS];
float ion_net_OPP_charge_dist[NUM_BUCKETS][NUM_Z_BUCKETS];


float wat_net_charge_dist[NUM_BUCKETS][NUM_Z_BUCKETS];
float wat_net_charge_den[NUM_BUCKETS][NUM_Z_BUCKETS];


float lipid_net_charge_den[NUM_BUCKETS][NUM_Z_BUCKETS];
float pep_pos_charge_den[NUM_BUCKETS][NUM_Z_BUCKETS];
float pep_neutral_charge_den[NUM_BUCKETS][NUM_Z_BUCKETS];

float pep_side_charge_den[NUM_BUCKETS][NUM_Z_BUCKETS];
float pep_bb_charge_den[NUM_BUCKETS][NUM_Z_BUCKETS];
float lipid_pos_charge_den[NUM_BUCKETS][NUM_Z_BUCKETS];
float lipid_neg_charge_den[NUM_BUCKETS][NUM_Z_BUCKETS];
float lipid_net_OPP_charge_den[NUM_BUCKETS][NUM_Z_BUCKETS];
float lipid_pos_OPP_charge_den[NUM_BUCKETS][NUM_Z_BUCKETS];
float lipid_neg_OPP_charge_den[NUM_BUCKETS][NUM_Z_BUCKETS];
float pep_net_charge_den[NUM_BUCKETS][NUM_Z_BUCKETS];
float pep_net_OPP_charge_den[NUM_BUCKETS][NUM_Z_BUCKETS];
float ion_net_charge_den[NUM_BUCKETS][NUM_Z_BUCKETS];
float ion_net_OPP_charge_den[NUM_BUCKETS][NUM_Z_BUCKETS];


int net_OPP_charge_dist[NUM_BUCKETS][NUM_Z_BUCKETS];
int neg_OPP_charge_dist[NUM_BUCKETS][NUM_Z_BUCKETS];
int pos_OPP_charge_dist[NUM_BUCKETS][NUM_Z_BUCKETS];

float net_charge_den[NUM_BUCKETS][NUM_Z_BUCKETS];
float net_OPP_charge_den[NUM_BUCKETS][NUM_Z_BUCKETS];

float uCell_xy;
float uCell_z;
float avg_uCell_xy;
float avg_uCell_z;
float side_xy;
float side_z;
int num_steps;


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


int pepZmin = 1000;
int pepZmax = -1000;
void add_2d_hist(bool in_upper, bool is_dmpc, bool is_dmpg, bool is_ion, bool is_pep, bool is_wat, bool is_side, bool is_charged_aa, float xyDist, float zVal, float charge, int atom_num, int lipid_sec) {

    int histBucket;
    int histZBucket;
    float xyVal;
    int rReg;
    int zReg;
    float sum_reg;
    float sum_diff;

    if ((zVal > MAX_Z_DIST) | (zVal < MIN_Z_DIST)) {
         return;
    }

    if (xyDist > MAX_XY_DIST) {
        return;
    }


    xyVal = (xyDist * BUCKETS_PER_DIST) + 0.00001;
    histBucket = (int) xyVal;
    if (histBucket >= NUM_BUCKETS) {
        histBucket = NUM_BUCKETS-1;
    }
    histZBucket = (int) (zVal - MIN_Z_DIST);
    if ((histZBucket >= NUM_Z_BUCKETS) || (histZBucket < 0))  {
        return;
    }
    if (is_dmpc) {
        lipid_net_charge_dist[histBucket][histZBucket] += charge;
        lipid_tot_charge += charge;
    } else if (is_dmpg) {
        lipid_net_charge_dist[histBucket][histZBucket] += charge;
        lipid_tot_charge += charge;
    } else if (is_ion) {
        ion_net_charge_dist[histBucket][histZBucket] += charge;
        ion_tot_charge += charge;
    } else if (is_pep) {
        pep_net_charge_dist[histBucket][histZBucket] += charge;
        pep_tot_charge += charge;
        if (is_side) {
            pep_side_charge_dist[histBucket][histZBucket] += charge;
        } else {
            pep_bb_charge_dist[histBucket][histZBucket] += charge;
        }
      
        if (is_charged_aa) {
            pep_pos_charge_dist[histBucket][histZBucket] += charge;
        } else {
           pep_neutral_charge_dist[histBucket][histZBucket] += charge;
        }
    } else if (is_wat) {
        wat_net_charge_dist[histBucket][histZBucket] += charge;
        wat_tot_charge += charge;
    }
    net_charge_dist[histBucket][histZBucket] += charge;
    tot_charge += charge;

    if (xyDist < PROXIMAL_DIST) {
        rReg = 0;
    } else if (xyDist < DISTANT_DIST) {
        rReg = 1;
    } else {
        rReg = 2;
    }

    if (zVal < 0) {    /* ONLY USE Reference leaflet */
        return;
    } else if (zVal < HYDROPHOBIC_Z_DIST) {
        zReg = 0;
    } else if (zVal < PHOS_Z_DIST) {
        zReg = 1;
    } else if (zVal < INTERFACE_Z_DIST) {
        zReg = 2;
    } else {
        zReg = 3;
    }

    charge_reg[rReg][zReg] += charge;
    tot_charge_reg[rReg][zReg] += charge;
    if (is_wat) {
        wat_charge_reg[rReg][zReg] += charge;
        tot_wat_charge_reg[rReg][zReg] += charge;
    } else if (is_dmpc) {
        PC_charge_reg_sec[rReg][zReg][lipid_sec] += charge;
        PC_charge_reg_sec_step[rReg][zReg][lipid_sec] += charge;
        lip_charge_reg[rReg][zReg] += charge;
        PC_charge_reg_z_sec[rReg][histZBucket][lipid_sec] += charge;
        tot_PC_tail_charge[lipid_sec] += charge;
    } else if (is_dmpg) {
        lip_charge_reg[rReg][zReg] += charge;
        PG_charge_reg_sec[rReg][zReg][lipid_sec] += charge;
        PG_charge_reg_sec_step[rReg][zReg][lipid_sec] += charge;
        PG_charge_reg_z_sec[rReg][histZBucket][lipid_sec] += charge;
        tot_PG_tail_charge[lipid_sec] += charge;
    } else if (is_pep) {
        pep_charge_reg[rReg][zReg] += charge;
    } else if (is_ion) {
        ion_charge_reg[rReg][zReg] += charge;
    } else {
        printf("ERROR - not type\n");
    }
    sum_reg = wat_charge_reg[rReg][zReg] + pep_charge_reg[rReg][zReg] + lip_charge_reg[rReg][zReg] + ion_charge_reg[rReg][zReg];
    sum_diff = abs(sum_reg - charge_reg[rReg][zReg]);
    if (sum_diff > 0.001) {
        printf("ERROR - (%f != %f), %d, %d %f\n", sum_reg, charge_reg[rReg][zReg], rReg, zReg, charge);
    } else {
//        printf("%f, %d, %d %f\n", sum_reg, rReg, zReg, charge);
    }
}

void get_xy_dists(float x1Diff, float x2Diff, float y1Diff, float y2Diff, float max_xyDist, float *d1,float *d2,float *d3,float *d4) {

    if (x1Diff < max_xyDist) {
        if (y1Diff < max_xyDist) {
            *d1 = pow((pow(x1Diff,2) + pow(y1Diff,2)),0.5);
        }
        if (y2Diff < max_xyDist) {
            *d2 = pow((pow(x1Diff,2) + pow(y2Diff,2)),0.5);
        }
    }
    if (x2Diff < max_xyDist) {
        if (y1Diff < max_xyDist) {
            *d3 = pow((pow(x2Diff,2) + pow(y1Diff,2)),0.5);
        }
        if (y2Diff < max_xyDist) {
            *d4 = pow((pow(x2Diff,2) + pow(y2Diff,2)),0.5);
        }
    }
}

void get_charge_hist(bool is_upper, float xCOM, float yCOM, ReadPdb * pdb, ReadDcd * dcd, ReadPsf * psf) {

    std::vector<ReadPsf::PSF_ATOM_PTR>::iterator first_it, end_it, it;
    int atom_num;
    float xVal, yVal, zVal;
    float xCell, x1Diff, y1Diff, z1Diff, x2Diff, y2Diff;
    float xyDists[4];
    int histBucket;
    int zBucket;
    int i, j, k;
    int mol_type;
    float charge;
    ReadPsf::PSF_ATOM_PTR atom_ptr;
    bool is_dmpc, is_dmpg, is_ion, is_pep, is_wat, is_side, charged_aa;
    int prev_lipid_num;
    int lipid_num;
    int lipid_atom_num;
    int lipid_sec;
    

    prev_lipid_num = 0;

    first_it = std::begin(psf->atoms);
    end_it = std::end(psf->atoms);

    for (it = first_it; it != end_it; ++it) {
        atom_ptr=*it;
        atom_num = atom_ptr->atom_num;

        mol_type = atom_ptr->mol_type_val; 
        if (atom_ptr->atom_type == SIDE_MON) {
            is_side = true;
        } else {
            is_side = false;
        }
        charged_aa = 0;
        if (mol_type == DMPC_MOL) {
            is_dmpc = true;
            is_dmpg = false;
            is_ion =  false;
            is_pep =  false;
            is_wat =  false;
            lipid_num = atom_ptr->mon_num;
            if (lipid_num != prev_lipid_num) {
                prev_lipid_num = lipid_num;
                lipid_atom_num = 1;
            } else {
                lipid_atom_num++;
            }
            if (lipid_atom_num <= 19) {
                lipid_sec = LIP_CHO;
            } else if (lipid_atom_num <= 24) {
                lipid_sec = LIP_P;
            } else if ((lipid_atom_num <= 30) || ((lipid_atom_num >=36) && (lipid_atom_num <= 39))) {
                lipid_sec = LIP_GLY;
            } else {
                 lipid_sec = LIP_TAIL;
            }
        } else if (mol_type == DMPG_MOL) { 
            is_dmpc = false;
            is_dmpg = true;
            is_ion =  false;
            is_pep =  false;
            is_wat =  false;
            lipid_num = atom_ptr->mon_num;
            if (lipid_num != prev_lipid_num) {
                prev_lipid_num = lipid_num;
                lipid_atom_num = 1;
            } else {
                lipid_atom_num++;
            }
            if (lipid_atom_num <= 12) {
                lipid_sec = LIP_CHO;
            } else if (lipid_atom_num <= 17) {
                lipid_sec = LIP_P;
            } else if ((lipid_atom_num <= 23) || ((lipid_atom_num >=29) && (lipid_atom_num <= 32))) {
                lipid_sec = LIP_GLY;
            } else {
                 lipid_sec = LIP_TAIL;
            }
        } else if ((mol_type == SOD_MOL) || (mol_type == CLA_MOL)) {
            is_dmpc = false;
            is_dmpg = false;
            is_ion =  true;
            is_pep =  false;
            is_wat =  false;
            lipid_sec = 0;
        } else if ((mol_type == PGL1_MOL) || (mol_type == PGL2_MOL)) {

            if ((atom_ptr->mon_num == 1) ||
              (atom_ptr->mon_num == 5) ||
              (atom_ptr->mon_num == 12) ||
              (atom_ptr->mon_num == 15) ||
              (atom_ptr->mon_num == 19)) {
                charged_aa = true;
            } else {
                charged_aa = false;
            }

            is_dmpg = false;
            is_dmpc = false;
            is_ion =  false;
            is_pep =  true;
            is_wat =  false;
            lipid_sec = 0;
        } else {   /* else water */
            is_dmpg = false;
            is_dmpc = false;
            is_ion =  false;
            is_pep =  false;
            is_wat =  true;
            lipid_sec = 0;
        }
        charge = atom_ptr->atom_charge;

        xCell = dcd->uCell[0];
        xVal = dcd->xValArray[0][atom_num-1];
        yVal = dcd->yValArray[0][atom_num-1];
        zVal = dcd->zValArray[0][atom_num-1];
        x1Diff = abs(xCOM - xVal);
        y1Diff = abs(yCOM - yVal);
        x2Diff = abs(xCell - x1Diff);
        y2Diff = abs(xCell - y1Diff);

        if (is_upper == false) {
            zVal = -zVal;
        }

        xyDists[0] =  xyDists[1] =  xyDists[2] =  xyDists[3] = MAX_XY_DIST;
        get_xy_dists(x1Diff, x2Diff, y1Diff, y2Diff, MAX_XY_DIST, &xyDists[0], &xyDists[1], &xyDists[2], &xyDists[3]);

        for (i=0; i<4;i++) {
             if (xyDists[i] < MAX_XY_DIST) {
                add_2d_hist(is_upper,  is_dmpc, is_dmpg, is_ion, is_pep, is_wat, is_side, charged_aa, xyDists[i], zVal, charge, atom_num, lipid_sec);

             }
        }
    }
}



int main(int argc, char * argv[]) {

    char dcd_file[MAX_FILENAME_LEN];
    char pdb_file[MAX_FILENAME_LEN];
    char psf_file[MAX_FILENAME_LEN];
    char cm_file[MAX_FILENAME_LEN];
    int i, j, k;
    float uCell[3];


//    arg[1] = Input_dir
//    arg[2] = tr
//    arg[3] = rep
//    arg[4] = prefix
//    arg[5] = first_str
//    arg[6] = last_str
//    arg[7] = input_dir

    for (i=0; i< NUM_LIPID_SEC; i++) {
        tot_PG_tail_charge[i] = 0.0;
        tot_PC_tail_charge[i] = 0.0;
    }


    for (i=0; i<NUM_REG_R; i++) {
        for (j=0; j<NUM_REG_Z; j++) {
            for (k=0; k<NUM_LIPID_SEC; k++) {
                PG_charge_reg_sec[i][j][k] = 0.0;
                PC_charge_reg_sec[i][j][k] = 0.0;
            }
        }
    }
    for (i=0; i<NUM_REG_R; i++) {
        for (j=0; j<NUM_REG_Z; j++) {
            for (k=0; k<NUM_Z_BUCKETS; k++) {
                PG_charge_reg_z_sec[i][j][k] = 0.0;
                PC_charge_reg_z_sec[i][j][k] = 0.0;
            }
        }
    }




    fflush(stdout);
    snprintf(pdb_file, MAX_FILENAME_LEN, "%s/spt.pdb", argv[7]);
    snprintf(psf_file, MAX_FILENAME_LEN, "%s/mem22_pep.psf", argv[7]);


    ReadPdb pdb(pdb_file);
    ReadPsf psf(psf_file);
    pdb.find_heavy_atoms(NUM_LIPID_LEAFLET);
    pdb.find_charged_atoms(NUM_LIPID_LEAFLET);


    int first_step, last_step, step;
    first_step = atoi(argv[5]);
    last_step = atoi(argv[6]);
    ReadDcd * dcd;

    tot_charge = 0;
    lipid_tot_charge = 0;
    ion_tot_charge = 0;
    pep_tot_charge = 0;
    wat_tot_charge = 0;

    for (i=0; i<NUM_BUCKETS; i++) {
        for (j=0; j<NUM_Z_BUCKETS; j++) {
            net_charge_dist[i][j] = 0;
            pos_charge_dist[i][j] = 0;
            neg_charge_dist[i][j] = 0;
            net_OPP_charge_dist[i][j] = 0;
            pos_OPP_charge_dist[i][j] = 0;
            neg_OPP_charge_dist[i][j] = 0;

            lipid_net_charge_dist[i][j] = 0;
            pep_side_charge_dist[i][j] = 0;
            pep_bb_charge_dist[i][j] = 0;
            pep_pos_charge_dist[i][j] = 0;
            pep_neutral_charge_dist[i][j] = 0;
            lipid_pos_charge_dist[i][j] = 0;
            lipid_neg_charge_dist[i][j] = 0;
            lipid_net_OPP_charge_dist[i][j] = 0;
            lipid_pos_OPP_charge_dist[i][j] = 0;
            lipid_neg_OPP_charge_dist[i][j] = 0;

            wat_net_charge_dist[i][j] = 0;

            pep_net_charge_dist[i][j] = 0;
            pep_net_OPP_charge_dist[i][j] = 0;
            ion_net_charge_dist[i][j] = 0;
            ion_net_OPP_charge_dist[i][j] = 0;
        }
    }



    uCell_xy = 0.0;
    uCell_z = 0.0;


    FILE *fpStepOut;
    char outFile[MAX_FILENAME_LEN];
    snprintf(outFile, MAX_FILENAME_LEN, "%s/charge_reg_%s_%s_%s-%s.dat", argv[1], argv[2], argv[3], argv[5], argv[6]);
    fpStepOut = fopen(outFile,"w");

    FILE *fpWatStepOut;
    snprintf(outFile, MAX_FILENAME_LEN, "%s/charge_reg_wat_%s_%s_%s-%s.dat", argv[1], argv[2], argv[3], argv[5], argv[6]);
    fpWatStepOut = fopen(outFile,"w");

    FILE *fpLipStepOut;
    snprintf(outFile, MAX_FILENAME_LEN, "%s/charge_reg_lip_%s_%s_%s-%s.dat", argv[1], argv[2], argv[3], argv[5], argv[6]);
    fpLipStepOut = fopen(outFile,"w");

    FILE *fpPepStepOut;
    snprintf(outFile, MAX_FILENAME_LEN, "%s/charge_reg_pep_%s_%s_%s-%s.dat", argv[1], argv[2], argv[3], argv[5], argv[6]);
    fpPepStepOut = fopen(outFile,"w");

    FILE *fpIonStepOut;
    snprintf(outFile, MAX_FILENAME_LEN, "%s/charge_reg_ion_%s_%s_%s-%s.dat", argv[1], argv[2], argv[3], argv[5], argv[6]);
    fpIonStepOut = fopen(outFile,"w");

    FILE * fpPGRegStep[NUM_LIPID_SEC];
    FILE * fpPCRegStep[NUM_LIPID_SEC];
    snprintf(outFile, MAX_FILENAME_LEN, "%s/charge_reg_PG_cho_step_%s_%s_%s-%s.dat", argv[1], argv[2], argv[3], argv[5], argv[6]);
    fpPGRegStep[0] = fopen(outFile,"w");
    snprintf(outFile, MAX_FILENAME_LEN, "%s/charge_reg_PG_pho_step_%s_%s_%s-%s.dat", argv[1], argv[2], argv[3], argv[5], argv[6]);
    fpPGRegStep[1] = fopen(outFile,"w");
    snprintf(outFile, MAX_FILENAME_LEN, "%s/charge_reg_PG_gly_step_%s_%s_%s-%s.dat", argv[1], argv[2], argv[3], argv[5], argv[6]);
    fpPGRegStep[2] = fopen(outFile,"w");
    snprintf(outFile, MAX_FILENAME_LEN, "%s/charge_reg_PG_tail_step_%s_%s_%s-%s.dat", argv[1], argv[2], argv[3], argv[5], argv[6]);
    fpPGRegStep[3] = fopen(outFile,"w");

    snprintf(outFile, MAX_FILENAME_LEN, "%s/charge_reg_PC_cho_step_%s_%s_%s-%s.dat", argv[1], argv[2], argv[3], argv[5], argv[6]);
    fpPCRegStep[0] = fopen(outFile,"w");
    snprintf(outFile, MAX_FILENAME_LEN, "%s/charge_reg_PC_pho_step_%s_%s_%s-%s.dat", argv[1], argv[2], argv[3], argv[5], argv[6]);
    fpPCRegStep[1] = fopen(outFile,"w");
    snprintf(outFile, MAX_FILENAME_LEN, "%s/charge_reg_PC_gly_step_%s_%s_%s-%s.dat", argv[1], argv[2], argv[3], argv[5], argv[6]);
    fpPCRegStep[2] = fopen(outFile,"w");
    snprintf(outFile, MAX_FILENAME_LEN, "%s/charge_reg_PC_tail_step_%s_%s_%s-%s.dat", argv[1], argv[2], argv[3], argv[5], argv[6]);
    fpPCRegStep[3] = fopen(outFile,"w");


    chargeArea[0] = (pow(PROXIMAL_DIST ,2) * M_PI);
    chargeArea[1] = (pow(DISTANT_DIST ,2) * M_PI) - chargeArea[0];
    chargeArea[2] = (pow(MAX_XY_DIST ,2) * M_PI) - (chargeArea[0] + chargeArea[1]);


    chargeVol[0][0] = (HYDROPHOBIC_Z_DIST * pow(PROXIMAL_DIST ,2) * M_PI);
    chargeVol[1][0] = (HYDROPHOBIC_Z_DIST * pow(DISTANT_DIST ,2) * M_PI) - chargeVol[0][0];
    chargeVol[2][0] = (HYDROPHOBIC_Z_DIST * pow(MAX_XY_DIST ,2) * M_PI) - (chargeVol[0][0] + chargeVol[1][0]);

    chargeVol[0][1] = ((PHOS_Z_DIST - HYDROPHOBIC_Z_DIST) * pow(PROXIMAL_DIST ,2) * M_PI);
    chargeVol[1][1] = ((PHOS_Z_DIST - HYDROPHOBIC_Z_DIST) * pow(DISTANT_DIST ,2) * M_PI) - chargeVol[0][1];
    chargeVol[2][1] = ((PHOS_Z_DIST - HYDROPHOBIC_Z_DIST) * pow(MAX_XY_DIST ,2) * M_PI) - (chargeVol[0][1] + chargeVol[1][1]);

    chargeVol[0][2] = ((INTERFACE_Z_DIST - PHOS_Z_DIST) * pow(PROXIMAL_DIST ,2) * M_PI);
    chargeVol[1][2] = ((INTERFACE_Z_DIST - PHOS_Z_DIST) * pow(DISTANT_DIST ,2) * M_PI)  - chargeVol[0][2];
    chargeVol[2][2] = ((INTERFACE_Z_DIST - PHOS_Z_DIST) * pow(MAX_XY_DIST ,2) * M_PI) - (chargeVol[0][2] + chargeVol[1][2]);

    chargeVol[0][3] = ((BULK_Z_DIST - INTERFACE_Z_DIST) * pow(PROXIMAL_DIST ,2) * M_PI);
    chargeVol[1][3] = ((BULK_Z_DIST - INTERFACE_Z_DIST) * pow(DISTANT_DIST ,2) * M_PI)  - chargeVol[0][3];
    chargeVol[2][3] = ((BULK_Z_DIST - INTERFACE_Z_DIST) * pow(MAX_XY_DIST ,2) * M_PI) - (chargeVol[0][3] + chargeVol[1][3]);


    for (i=0; i<NUM_REG_R; i++) {
        for (j=0; j<NUM_REG_Z; j++) {
            tot_charge_reg[i][j] = 0.0;
            tot_wat_charge_reg[i][j] = 0.0;
        }
    }

    num_steps = (last_step+1) - first_step;
    for (step=first_step; step < (last_step+1); step++) {
        snprintf(dcd_file, MAX_FILENAME_LEN, "%s/output/%s%s_%05d_%s.dcd", argv[7], argv[4], argv[2], step, argv[3]);
        printf("READ DCD %s\n", dcd_file);
        fflush(stdout);
        dcd = new ReadDcd(dcd_file);

        for (i=0; i<NUM_REG_R; i++) {
            for (j=0; j<NUM_REG_Z; j++) {
                for (k=0; k<NUM_LIPID_SEC; k++) {
                    PG_charge_reg_sec_step[i][j][k] = 0.0;
                    PC_charge_reg_sec_step[i][j][k] = 0.0;
                }
            }
        }

        for (i=0; i<NUM_REG_R; i++) {
            for (j=0; j<NUM_REG_Z; j++) {
                charge_reg[i][j] = 0;
                wat_charge_reg[i][j] = 0;
                lip_charge_reg[i][j] = 0;
                pep_charge_reg[i][j] = 0;
                ion_charge_reg[i][j] = 0;
            }
        }

        dcd->uCell[0] = dcd->boxArray[0]->ibox[0];
        dcd->uCell[1] = dcd->boxArray[0]->ibox[2];
        dcd->uCell[2] = dcd->boxArray[0]->ibox[5];
        uCell_xy += dcd->uCell[0];
        uCell_z += dcd->uCell[2];

        float xUValCOM, yUValCOM, zUValCOM;
        float xLValCOM, yLValCOM, zLValCOM;


        get_pep_com(true, &xUValCOM, &yUValCOM, &zUValCOM, &pdb, dcd);
        get_pep_com(false, &xLValCOM, &yLValCOM, &zLValCOM, &pdb, dcd);

        printf("UPPER \n");
        fflush(stdout);
        get_charge_hist(true, xUValCOM, yUValCOM, &pdb, dcd, &psf);
        printf("LOWER \n");
        get_charge_hist(false, xLValCOM, yLValCOM, &pdb, dcd, &psf);


        for (i=0; i<NUM_REG_R; i++) {
            for (j=0; j<NUM_REG_Z; j++) {
                for (k=0; k<NUM_LIPID_SEC; k++) {
                    PG_charge_reg_sec_step[i][j][k] = PG_charge_reg_sec_step[i][j][k] / chargeVol[i][j];
                    PC_charge_reg_sec_step[i][j][k] = PC_charge_reg_sec_step[i][j][k] / chargeVol[i][j];
                }
            }
        }


       for (i=0; i<NUM_REG_R; i++) {
            for (j=0; j<NUM_REG_Z; j++) {
                for (k=0; k<NUM_LIPID_SEC; k++) {
                    fprintf(fpPGRegStep[k], "%f ", PG_charge_reg_sec_step[i][j][k]);
                    fprintf(fpPCRegStep[k], "%f ", PC_charge_reg_sec_step[i][j][k]);
                }
            }
        }
        for (k=0; k<NUM_LIPID_SEC; k++) {
            fprintf(fpPGRegStep[k], "\n");
            fprintf(fpPCRegStep[k], "\n");
        }


        for (i=0; i<NUM_REG_R; i++) {
            for (j=0; j<NUM_REG_Z; j++) {
                charge_reg[i][j] = charge_reg[i][j] / chargeVol[i][j];
                wat_charge_reg[i][j] = wat_charge_reg[i][j] / chargeVol[i][j];
                lip_charge_reg[i][j] = lip_charge_reg[i][j] / chargeVol[i][j];
                pep_charge_reg[i][j] = pep_charge_reg[i][j] / chargeVol[i][j];
                ion_charge_reg[i][j] = ion_charge_reg[i][j] / chargeVol[i][j];
            }
        }
        for (i=0; i<NUM_REG_R; i++) {
            for (j=0; j<NUM_REG_Z; j++) {
                fprintf(fpStepOut, "%f ", charge_reg[i][j]);
                fprintf(fpWatStepOut, "%f ", wat_charge_reg[i][j]);
                fprintf(fpLipStepOut, "%f ", lip_charge_reg[i][j]);
                fprintf(fpPepStepOut, "%f ", pep_charge_reg[i][j]);
                fprintf(fpIonStepOut, "%f ", ion_charge_reg[i][j]);
            }
        }

        fprintf(fpStepOut, "\n");
        fprintf(fpWatStepOut, "\n");
        fprintf(fpLipStepOut, "\n");
        fprintf(fpPepStepOut, "\n");
        fprintf(fpIonStepOut, "\n");
        delete dcd;
    }
    fclose(fpStepOut);
    fclose(fpWatStepOut);
    fclose(fpLipStepOut);
    fclose(fpPepStepOut);
    fclose(fpIonStepOut);

    for (k=0; k<NUM_LIPID_SEC; k++) {
        fclose(fpPGRegStep[k]);
        fclose(fpPCRegStep[k]);
    }
    for (i=0; i<NUM_REG_R; i++) {
        for (j=0; j<NUM_REG_Z; j++) {
            tot_charge_reg[i][j] = tot_charge_reg[i][j] / (chargeVol[i][j] * num_steps);
        }
    }
    for (i=0; i<NUM_REG_R; i++) {
        printf("CHARGE Rad Region %d   %f %f %f\n", i, tot_charge_reg[i][0], tot_charge_reg[i][1], tot_charge_reg[i][2]);
    }

    for (i=0; i<NUM_REG_R; i++) {
        for (j=0; j<NUM_REG_Z; j++) {
            tot_wat_charge_reg[i][j] = tot_wat_charge_reg[i][j] / (chargeVol[i][j] * num_steps);
        }
    }
    for (i=0; i<NUM_REG_R; i++) {
        printf("CHARGE Water Rad Region %d   %f %f %f\n", i, tot_wat_charge_reg[i][0], tot_wat_charge_reg[i][1], tot_wat_charge_reg[i][2]);
    }


    float tot_vol;
    tot_vol = (pow(((float)MAX_XY_DIST), 2) * M_PI) * (MAX_Z_DIST - MIN_Z_DIST);
    tot_charge = tot_charge/ (num_steps * tot_vol); 
    printf("Tot Charge %f\n", tot_charge);




//    printf("DONE STEPS\n");
    num_steps = (last_step+1) - first_step;
    avg_uCell_xy = uCell_xy/num_steps;
    avg_uCell_z = uCell_z/num_steps;
    side_xy = avg_uCell_xy/2.0;
    side_z = avg_uCell_z/2.0;


    float curr_dist;
    float currRad;
    float prevRad;
    float myArea;
    curr_dist = FIRST_DIST;
    prevRad = 0.0;
    for (i=0; i<NUM_BUCKETS; i++) {
        currRad = (i+1) * BUCKET_DIST;
        myArea =  (pow(currRad,2) - pow(prevRad, 2)) * M_PI;
        prevRad = currRad;

        for (j=0; j<NUM_Z_BUCKETS; j++) {
            net_charge_den[i][j] = (net_charge_dist[i][j]/myArea)/num_steps;


            lipid_net_charge_den[i][j] = (lipid_net_charge_dist[i][j]/myArea)/num_steps;
            pep_side_charge_den[i][j] = (pep_side_charge_dist[i][j]/myArea)/num_steps;
            pep_bb_charge_den[i][j] = (pep_bb_charge_dist[i][j]/myArea)/num_steps;
            pep_pos_charge_den[i][j] = (pep_pos_charge_dist[i][j]/myArea)/num_steps;
            pep_neutral_charge_den[i][j] = (pep_neutral_charge_dist[i][j]/myArea)/num_steps;

            lipid_pos_charge_den[i][j] =  (lipid_pos_charge_dist[i][j]/myArea)/num_steps;

            lipid_neg_charge_den[i][j] = (lipid_neg_charge_dist[i][j]/myArea)/num_steps;

            pep_net_charge_den[i][j] = (pep_net_charge_dist[i][j]/myArea)/num_steps;
            ion_net_charge_den[i][j] = (ion_net_charge_dist[i][j]/myArea)/num_steps;
            wat_net_charge_den[i][j] = (wat_net_charge_dist[i][j]/myArea)/num_steps;
        }
    }



//
//    Write charge density hist file
    FILE *fpOut;
    snprintf(outFile, MAX_FILENAME_LEN, "%s/charge_den_r_%s_%s_%s-%s.dat", argv[1], argv[2], argv[3], argv[5], argv[6]);
    fpOut = fopen(outFile,"w");
    for (i=0; i<NUM_BUCKETS; i++) {
        for (j=0; j<NUM_Z_BUCKETS; j++) {
            fprintf(fpOut, "%f ", net_charge_den[i][j]);
        }
        fprintf(fpOut, "\n");
    }
    fclose(fpOut);

//
//    Write lipid charge density hist file
    snprintf(outFile, MAX_FILENAME_LEN, "%s/lipid_charge_den_r_%s_%s_%s-%s.dat", argv[1], argv[2], argv[3], argv[5], argv[6]);
    fpOut = fopen(outFile,"w");
    for (i=0; i<NUM_BUCKETS; i++) {
        for (j=0; j<NUM_Z_BUCKETS; j++) {
            fprintf(fpOut, "%f ", lipid_net_charge_den[i][j]);
        }
        fprintf(fpOut, "\n");
    }
    fclose(fpOut);


//
//    Write pep side charge density hist file
    snprintf(outFile, MAX_FILENAME_LEN, "%s/pep_side_charge_den_r_%s_%s_%s-%s.dat", argv[1], argv[2], argv[3], argv[5], argv[6]);
    fpOut = fopen(outFile,"w");
    for (i=0; i<NUM_BUCKETS; i++) {
        for (j=0; j<NUM_Z_BUCKETS; j++) {
            fprintf(fpOut, "%f ", pep_side_charge_den[i][j]);
        }
        fprintf(fpOut, "\n");
    }
    fclose(fpOut);

 //
//    Write pep bb charge charge density hist file
    snprintf(outFile, MAX_FILENAME_LEN, "%s/pep_bb_charge_den_r_%s_%s_%s-%s.dat", argv[1], argv[2], argv[3], argv[5], argv[6]);
    fpOut = fopen(outFile,"w");
    for (i=0; i<NUM_BUCKETS; i++) {
        for (j=0; j<NUM_Z_BUCKETS; j++) {
            fprintf(fpOut, "%f ", pep_bb_charge_den[i][j]);
        }
        fprintf(fpOut, "\n");
    }
    fclose(fpOut);

//
//    Write pep charged AA charge charge density hist file
    snprintf(outFile, MAX_FILENAME_LEN, "%s/pep_pos_aa_charge_den_r_%s_%s_%s-%s.dat", argv[1], argv[2], argv[3], argv[5], argv[6]);
    fpOut = fopen(outFile,"w");
    for (i=0; i<NUM_BUCKETS; i++) {
        for (j=0; j<NUM_Z_BUCKETS; j++) {
            fprintf(fpOut, "%f ", pep_pos_charge_den[i][j]);
        }
        fprintf(fpOut, "\n");
    }
    fclose(fpOut);

//
//    Write pep neutral AA charge charge density hist file
    snprintf(outFile, MAX_FILENAME_LEN, "%s/pep_neutral_aa_charge_den_r_%s_%s_%s-%s.dat", argv[1], argv[2], argv[3], argv[5], argv[6]);
    fpOut = fopen(outFile,"w");
    for (i=0; i<NUM_BUCKETS; i++) {
        for (j=0; j<NUM_Z_BUCKETS; j++) {
            fprintf(fpOut, "%f ", pep_neutral_charge_den[i][j]);
        }
        fprintf(fpOut, "\n");
    }
    fclose(fpOut);


//
//    Write lipid charge density hist file
    snprintf(outFile, MAX_FILENAME_LEN, "%s/wat_charge_den_r_%s_%s_%s-%s.dat", argv[1], argv[2], argv[3], argv[5], argv[6]);
    fpOut = fopen(outFile,"w");
    for (i=0; i<NUM_BUCKETS; i++) {
        for (j=0; j<NUM_Z_BUCKETS; j++) {
            fprintf(fpOut, "%f ", wat_net_charge_den[i][j]);
        }
        fprintf(fpOut, "\n");
    }
    fclose(fpOut);



    snprintf(outFile, MAX_FILENAME_LEN, "%s/pep_charge_den_r_%s_%s_%s-%s.dat", argv[1], argv[2], argv[3], argv[5], argv[6]);
    fpOut = fopen(outFile,"w");
    for (i=0; i<NUM_BUCKETS; i++) {
        for (j=0; j<NUM_Z_BUCKETS; j++) {
            fprintf(fpOut, "%f ", pep_net_charge_den[i][j]);
        }
        fprintf(fpOut, "\n");
    }
    fclose(fpOut);


    snprintf(outFile, MAX_FILENAME_LEN, "%s/ion_charge_den_r_%s_%s_%s-%s.dat", argv[1], argv[2], argv[3], argv[5], argv[6]);
    fpOut = fopen(outFile,"w");
    for (i=0; i<NUM_BUCKETS; i++) {
        for (j=0; j<NUM_Z_BUCKETS; j++) {
            fprintf(fpOut, "%f ", ion_net_charge_den[i][j]);
        }
        fprintf(fpOut, "\n");
    }
    fclose(fpOut);

    for (i=0; i<NUM_REG_R; i++) {
        for (j=0; j<NUM_REG_Z; j++) {
            for (k=0; k<NUM_LIPID_SEC; k++) {
                PG_charge_reg_sec[i][j][k] = PG_charge_reg_sec[i][j][k] / (chargeVol[i][j] * num_steps);
                PC_charge_reg_sec[i][j][k] = PC_charge_reg_sec[i][j][k] / (chargeVol[i][j] * num_steps);
            }
        }
    }



    FILE * fpOutCho;
    FILE * fpOutPho;
    FILE * fpOutGly;
    FILE * fpOutTail;

    snprintf(outFile, MAX_FILENAME_LEN, "%s/PC_cho_z_charge_den_r_%s_%s_%s-%s.dat", argv[1], argv[2], argv[3], argv[5], argv[6]);
    fpOutCho = fopen(outFile,"w");
    snprintf(outFile, MAX_FILENAME_LEN, "%s/PC_pho_z_charge_den_r_%s_%s_%s-%s.dat", argv[1], argv[2], argv[3], argv[5], argv[6]);
    fpOutPho = fopen(outFile,"w");
    snprintf(outFile, MAX_FILENAME_LEN, "%s/PC_gly_z_charge_den_r_%s_%s_%s-%s.dat", argv[1], argv[2], argv[3], argv[5], argv[6]);
    fpOutGly = fopen(outFile,"w");
    snprintf(outFile, MAX_FILENAME_LEN, "%s/PC_tail_z_charge_den_r_%s_%s_%s-%s.dat", argv[1], argv[2], argv[3], argv[5], argv[6]);
    fpOutTail = fopen(outFile,"w");
    for (i=0; i<NUM_REG_R; i++) {
        for (j=0; j<NUM_Z_BUCKETS; j++) {
            fprintf(fpOutCho, "%f ", PC_charge_reg_z_sec[i][j][0]/(chargeArea[i] * num_steps));
            fprintf(fpOutPho, "%f ", PC_charge_reg_z_sec[i][j][1]/(chargeArea[i] * num_steps));
            fprintf(fpOutGly, "%f ", PC_charge_reg_z_sec[i][j][2]/(chargeArea[i] * num_steps));
            fprintf(fpOutTail, "%f ", PC_charge_reg_z_sec[i][j][3]/(chargeArea[i] * num_steps));
        }
        fprintf(fpOutCho, "\n");
        fprintf(fpOutPho, "\n");
        fprintf(fpOutGly, "\n");
        fprintf(fpOutTail, "\n");
    }
    fclose(fpOutCho);
    fclose(fpOutPho);
    fclose(fpOutGly);
    fclose(fpOutTail);

    snprintf(outFile, MAX_FILENAME_LEN, "%s/PG_cho_z_charge_den_r_%s_%s_%s-%s.dat", argv[1], argv[2], argv[3], argv[5], argv[6]);
    fpOutCho = fopen(outFile,"w");
    snprintf(outFile, MAX_FILENAME_LEN, "%s/PG_pho_z_charge_den_r_%s_%s_%s-%s.dat", argv[1], argv[2], argv[3], argv[5], argv[6]);
    fpOutPho = fopen(outFile,"w");
    snprintf(outFile, MAX_FILENAME_LEN, "%s/PG_gly_z_charge_den_r_%s_%s_%s-%s.dat", argv[1], argv[2], argv[3], argv[5], argv[6]);
    fpOutGly = fopen(outFile,"w");
    snprintf(outFile, MAX_FILENAME_LEN, "%s/PG_tail_z_charge_den_r_%s_%s_%s-%s.dat", argv[1], argv[2], argv[3], argv[5], argv[6]);
    fpOutTail = fopen(outFile,"w");
    for (i=0; i<NUM_REG_R; i++) { 
        for (j=0; j<NUM_Z_BUCKETS; j++) {
            fprintf(fpOutCho, "%f ", PG_charge_reg_z_sec[i][j][0]/(chargeArea[i] * num_steps));
            fprintf(fpOutPho, "%f ", PG_charge_reg_z_sec[i][j][1]/(chargeArea[i] * num_steps));
            fprintf(fpOutGly, "%f ", PG_charge_reg_z_sec[i][j][2]/(chargeArea[i] * num_steps));
            fprintf(fpOutTail, "%f ", PG_charge_reg_z_sec[i][j][3]/(chargeArea[i] * num_steps));
        }
        fprintf(fpOutCho, "\n");
        fprintf(fpOutPho, "\n");
        fprintf(fpOutGly, "\n");
        fprintf(fpOutTail, "\n");
    }
    fclose(fpOutCho);
    fclose(fpOutPho);
    fclose(fpOutGly);
    fclose(fpOutTail);



    snprintf(outFile, MAX_FILENAME_LEN, "%s/PC_cho_charge_den_r_%s_%s_%s-%s.dat", argv[1], argv[2], argv[3], argv[5], argv[6]);
    fpOutCho = fopen(outFile,"w");
    snprintf(outFile, MAX_FILENAME_LEN, "%s/PC_pho_charge_den_r_%s_%s_%s-%s.dat", argv[1], argv[2], argv[3], argv[5], argv[6]);
    fpOutPho = fopen(outFile,"w");
    snprintf(outFile, MAX_FILENAME_LEN, "%s/PC_gly_charge_den_r_%s_%s_%s-%s.dat", argv[1], argv[2], argv[3], argv[5], argv[6]);
    fpOutGly = fopen(outFile,"w");
    snprintf(outFile, MAX_FILENAME_LEN, "%s/PC_tail_charge_den_r_%s_%s_%s-%s.dat", argv[1], argv[2], argv[3], argv[5], argv[6]);
    fpOutTail = fopen(outFile,"w");
    for (i=0; i<NUM_REG_R; i++) {
        for (j=0; j<NUM_REG_Z; j++) {
            fprintf(fpOutCho, "%f ", PC_charge_reg_sec[i][j][0]);
            fprintf(fpOutPho, "%f ", PC_charge_reg_sec[i][j][1]);
            fprintf(fpOutGly, "%f ", PC_charge_reg_sec[i][j][2]);
            fprintf(fpOutTail, "%f ", PC_charge_reg_sec[i][j][3]);
        }
        fprintf(fpOutCho, "\n");
        fprintf(fpOutPho, "\n");
        fprintf(fpOutGly, "\n");
        fprintf(fpOutTail, "\n");
    }
    fclose(fpOutCho);
    fclose(fpOutPho);
    fclose(fpOutGly);
    fclose(fpOutTail);

   snprintf(outFile, MAX_FILENAME_LEN, "%s/PG_cho_charge_den_r_%s_%s_%s-%s.dat", argv[1], argv[2], argv[3], argv[5], argv[6]);
    fpOutCho = fopen(outFile,"w");
    snprintf(outFile, MAX_FILENAME_LEN, "%s/PG_pho_charge_den_r_%s_%s_%s-%s.dat", argv[1], argv[2], argv[3], argv[5], argv[6]);
    fpOutPho = fopen(outFile,"w");
    snprintf(outFile, MAX_FILENAME_LEN, "%s/PG_gly_charge_den_r_%s_%s_%s-%s.dat", argv[1], argv[2], argv[3], argv[5], argv[6]);
    fpOutGly = fopen(outFile,"w");
    snprintf(outFile, MAX_FILENAME_LEN, "%s/PG_tail_charge_den_r_%s_%s_%s-%s.dat", argv[1], argv[2], argv[3], argv[5], argv[6]);
    fpOutTail = fopen(outFile,"w");
    for (i=0; i<NUM_REG_R; i++) {
        for (j=0; j<NUM_REG_Z; j++) {
            fprintf(fpOutCho, "%f ", PG_charge_reg_sec[i][j][0]);
            fprintf(fpOutPho, "%f ", PG_charge_reg_sec[i][j][1]);
            fprintf(fpOutGly, "%f ", PG_charge_reg_sec[i][j][2]);
            fprintf(fpOutTail, "%f ", PG_charge_reg_sec[i][j][3]);
        }
        fprintf(fpOutCho, "\n");
        fprintf(fpOutPho, "\n");
        fprintf(fpOutGly, "\n");
        fprintf(fpOutTail, "\n");
    }
    fclose(fpOutCho);
    fclose(fpOutPho);
    fclose(fpOutGly);
    fclose(fpOutTail);

    printf("PC charge %f %f %f %f\n", 
      tot_PC_tail_charge[0],
      tot_PC_tail_charge[1],
      tot_PC_tail_charge[2],
      tot_PC_tail_charge[3]);

    printf("PG charge %f %f %f %f\n", 
      tot_PG_tail_charge[0],
      tot_PG_tail_charge[1],
      tot_PG_tail_charge[2],
      tot_PG_tail_charge[3]);

}

