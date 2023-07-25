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


int Num_rem;

#define E_FIELD_SIZE 30
#define MAX_FORCE_DIST   39
#define KVAL  8988000000   //  N⋅A2⋅C−2
#define C_CONVERT 1.60217646e-19 
#define Msq_to_Asq   1.0e+20

#define NUM_Z_VALS 30

float tot_charge_Z[NUM_Z_VALS];
float tot_charge_dist_Z[NUM_Z_VALS];

float lip_charge_Z[NUM_Z_VALS];
float lip_charge_dist_Z[NUM_Z_VALS];

float ion_charge_Z[NUM_Z_VALS];
float ion_charge_dist_Z[NUM_Z_VALS];

float wat_charge_Z[NUM_Z_VALS];
float wat_charge_dist_Z[NUM_Z_VALS];

float pep_charge_Z[NUM_Z_VALS];
float pep_charge_dist_Z[NUM_Z_VALS];


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
void add_2d_hist(bool in_upper, bool is_dmpc, bool is_dmpg, bool is_ion, bool is_pep, bool is_wat, bool is_side, bool is_charged_aa, float xyDist, float zVal, float charge, int atom_num) {

    int histBucket;
    int histZBucket;
    float xyVal;
    int rReg;
    int zReg;

    if (xyDist > MAX_XY_DIST) {
        return;
    }

    histZBucket = (int) zVal;
    if ((histZBucket >= NUM_Z_VALS) || (histZBucket < 0))  {
        return;
    }

    tot_charge_Z[histZBucket] += charge;
    if (xyDist > DISTANT_DIST) {
        tot_charge_dist_Z[histZBucket] += charge;
    }

    if ((is_dmpc) || (is_dmpg)) {
        lip_charge_Z[histZBucket] += charge;
        if (xyDist > DISTANT_DIST) {
            lip_charge_dist_Z[histZBucket] += charge;
        }

    } else if (is_ion) {
        ion_charge_Z[histZBucket] += charge;
        if (xyDist > DISTANT_DIST) {
            ion_charge_dist_Z[histZBucket] += charge;
        }
    } else if (is_pep) {
        pep_charge_Z[histZBucket] += charge;
        if (xyDist > DISTANT_DIST) {
            pep_charge_dist_Z[histZBucket] += charge;
        }
    } else if (is_wat) {
        wat_charge_Z[histZBucket] += charge;
        if (xyDist > DISTANT_DIST) {
            wat_charge_dist_Z[histZBucket] += charge;
        }
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
    int i;
    int mol_type;
    float charge;
    ReadPsf::PSF_ATOM_PTR atom_ptr;
    bool is_dmpc, is_dmpg, is_ion, is_pep, is_wat, is_side, charged_aa;
    


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
        } else if (mol_type == DMPG_MOL) { 
            is_dmpc = true;
            is_dmpg = false;
            is_ion =  false;
            is_pep =  false;
            is_wat =  false;
        } else if ((mol_type == SOD_MOL) || (mol_type == CLA_MOL)) {
            is_dmpc = false;
            is_dmpg = false;
            is_ion =  true;
            is_pep =  false;
            is_wat =  false;
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
        } else {   /* else water */
            is_dmpg = false;
            is_dmpc = false;
            is_ion =  false;
            is_pep =  false;
            is_wat =  true;
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
                add_2d_hist(is_upper,  is_dmpc, is_dmpg, is_ion, is_pep, is_wat, is_side, charged_aa, xyDists[i], zVal, charge, atom_num);

             }
        }
    }
}



int main(int argc, char * argv[]) {

    char dcd_file[MAX_FILENAME_LEN];
    char pdb_file[MAX_FILENAME_LEN];
    char psf_file[MAX_FILENAME_LEN];
    char cm_file[MAX_FILENAME_LEN];
    int i, j;
    float uCell[3];
    float chargeVolDist, chargeVolTot;
    FILE * fpOut;


//    arg[1] = Input_dir
//    arg[2] = tr
//    arg[3] = rep
//    arg[4] = prefix
//    arg[5] = first_str
//    arg[6] = last_str
//    arg[7] = input_dir
//    arg[8] = outFile

    fflush(stdout);
    snprintf(pdb_file, MAX_FILENAME_LEN, "%s/spt.pdb", argv[1]);
    snprintf(psf_file, MAX_FILENAME_LEN, "%s/mem22_pep.psf", argv[1]);


    ReadPdb pdb(pdb_file);
    ReadPsf psf(psf_file);
    pdb.find_heavy_atoms(NUM_LIPID_LEAFLET);
    pdb.find_charged_atoms(NUM_LIPID_LEAFLET);


    int first_step, last_step, step;
    first_step = atoi(argv[5]);
    last_step = atoi(argv[6]);
    ReadDcd * dcd;

    chargeVolDist = (pow(MAX_XY_DIST ,2) * M_PI) - (pow(DISTANT_DIST ,2) * M_PI);
    chargeVolTot = pow(MAX_XY_DIST ,2) * M_PI;

    for (i=0; i<NUM_Z_VALS; i++) {
        tot_charge_Z[i] = 0;
        tot_charge_dist_Z[i] = 0;
        lip_charge_Z[i] = 0;
        lip_charge_dist_Z[i] = 0;
        ion_charge_Z[i] = 0;
        ion_charge_dist_Z[i] = 0;
        wat_charge_Z[i] = 0;
        wat_charge_dist_Z[i] = 0;
        pep_charge_Z[i] = 0;
        pep_charge_dist_Z[i] = 0;
    }



    num_steps = (last_step+1) - first_step;
    for (step=first_step; step < (last_step+1); step++) {
        snprintf(dcd_file, MAX_FILENAME_LEN, "%s/output/%s%s_%05d_%s.dcd", argv[1], argv[4], argv[2], step, argv[3]);
//        printf("READ DCD %s\n", dcd_file);
        dcd = new ReadDcd(dcd_file);

        dcd->uCell[0] = dcd->boxArray[0]->ibox[0];
        dcd->uCell[1] = dcd->boxArray[0]->ibox[2];
        dcd->uCell[2] = dcd->boxArray[0]->ibox[5];
        uCell_xy = dcd->uCell[0];
        uCell_z = dcd->uCell[2];

        float xUValCOM, yUValCOM, zUValCOM;
        float xLValCOM, yLValCOM, zLValCOM;


        get_pep_com(true, &xUValCOM, &yUValCOM, &zUValCOM, &pdb, dcd);
        get_pep_com(false, &xLValCOM, &yLValCOM, &zLValCOM, &pdb, dcd);

        get_charge_hist(true, xUValCOM, yUValCOM, &pdb, dcd, &psf);
        get_charge_hist(false, xLValCOM, yLValCOM, &pdb, dcd, &psf);

        delete dcd;
    }
    num_steps = (last_step+1) - first_step;

    float val;
    fpOut = fopen(argv[7],"w");
    for (i=0; i<NUM_Z_VALS; i++) {
        val = i + 0.5;
        fprintf(fpOut, "%f ",val);
        fprintf(fpOut, "%2.8f ", (tot_charge_Z[i] / (num_steps * chargeVolTot)));
        fprintf(fpOut, "%2.8f ", (tot_charge_dist_Z[i] / (num_steps * chargeVolDist)));

        fprintf(fpOut, "%2.8f ", (lip_charge_Z[i] / (num_steps * chargeVolTot)));
        fprintf(fpOut, "%2.8f ", (lip_charge_dist_Z[i] / (num_steps * chargeVolDist)));

        fprintf(fpOut, "%2.8f ", (ion_charge_Z[i] / (num_steps * chargeVolTot)));
        fprintf(fpOut, "%2.8f ", (ion_charge_dist_Z[i] / (num_steps * chargeVolDist)));

        fprintf(fpOut, "%2.8f ", (wat_charge_Z[i] / (num_steps * chargeVolTot)));
        fprintf(fpOut, "%2.8f ", (wat_charge_dist_Z[i] / (num_steps * chargeVolDist)));

        fprintf(fpOut, "%2.8f ", (pep_charge_Z[i] / (num_steps * chargeVolTot)));
        fprintf(fpOut, "%2.8f\n", (pep_charge_dist_Z[i] / (num_steps * chargeVolDist)));

    }
    fclose(fpOut);

}


