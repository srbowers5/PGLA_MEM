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
#include "ReadPsf.h"

#define LIP_CHO         0
#define LIP_P           1
#define LIP_GLY         2
#define LIP_TAIL        3
#define NUM_LIPID_SEC   4

float tot_PG_charge[NUM_LIPID_SEC];
float tot_PC_charge[NUM_LIPID_SEC];


void get_charge_group(ReadPsf * psf) {

    std::vector<ReadPsf::PSF_ATOM_PTR>::iterator first_it, end_it, it;
    int atom_num;
    int mol_type;
    float charge;
    ReadPsf::PSF_ATOM_PTR atom_ptr;
    int prev_lipid_num;
    int lipid_num;
    int lipid_atom_num;
    int lipid_sec;
    int foundDMPC, foundDMPG;
    

    prev_lipid_num = 0;
    foundDMPC = 0;
    foundDMPG = 0;

    first_it = std::begin(psf->atoms);
    end_it = std::end(psf->atoms);

    for (it = first_it; it != end_it; ++it) {
        atom_ptr=*it;
        atom_num = atom_ptr->atom_num;

        mol_type = atom_ptr->mol_type_val; 
        if ((mol_type == DMPC_MOL)  && (foundDMPC < 2)) {
            lipid_num = atom_ptr->mon_num;
            if (lipid_num != prev_lipid_num) {
                foundDMPC++;
                if (foundDMPC == 2) {   /* Done with molecule */
                    continue;
                } else {                /* Start new molecule */
                    prev_lipid_num = lipid_num;
                    lipid_atom_num = 1;
                }
            } else {                    /* Continue molecule */
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
            tot_PC_charge[lipid_sec] += atom_ptr->atom_charge;
        } else if ((mol_type == DMPG_MOL) && (foundDMPG < 2)) { 

            lipid_num = atom_ptr->mon_num;
            if (lipid_num != prev_lipid_num) {
                foundDMPG++;
                if (foundDMPG == 2) {   /* Done with molecule */
                    continue;
                } else {                /* Start new molecule */
                    prev_lipid_num = lipid_num;
                    lipid_atom_num = 1;
                }
            } else {                    /* Continue molecule */
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
            tot_PG_charge[lipid_sec] += atom_ptr->atom_charge;
        }
    }
}



int main(int argc, char * argv[]) {

    int i;
    FILE * fpOut;



//    arg[1] = Input_psf_file
//    arg[2] = Output file

    for (i=0; i< NUM_LIPID_SEC; i++) {
        tot_PG_charge[i] = 0.0;
        tot_PC_charge[i] = 0.0;
    }


    ReadPsf psf(argv[1]);

    get_charge_group(&psf);

    fpOut = fopen(argv[2],"w");

    fprintf(fpOut, "PC_Choline %f\n", tot_PC_charge[LIP_CHO]);
    fprintf(fpOut, "PC_Gly %f\n", tot_PC_charge[LIP_GLY]);
    fprintf(fpOut, "PC_Pho %f\n", tot_PC_charge[LIP_P]);
    fprintf(fpOut, "PC_Tail %f\n", tot_PC_charge[LIP_TAIL]);


    fprintf(fpOut, "PG_Terminal_GLy %f\n", tot_PG_charge[LIP_CHO]);
    fprintf(fpOut, "PG_Gly %f\n", tot_PG_charge[LIP_GLY]);
    fprintf(fpOut, "PG_Pho %f\n", tot_PG_charge[LIP_P]);
    fprintf(fpOut, "PG_Tail %f\n", tot_PG_charge[LIP_TAIL]);
}

