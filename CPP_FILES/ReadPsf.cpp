#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include <stdarg.h>
#include <iostream>
#include <fstream>
#include <arpa/inet.h>
#include <vector>
#include "ReadPsf.h"

#define NUM_LIPID_LEAFLET 49

int ReadPsf::parse_psf_line(const char * line, int prev_atom_num) {

    PSF_ATOM_PTR atom_ptr;
    char tmp_atom_string[6];
    char tmp_mon_type_string[5];
    char tmp_aa_type_string[5];
    char tmp_charge_string[10];
    char tmp_mol_type_string[10];
    char tmp_mol_num_string[30];
    char tmp_mon_num_string[5];
    int atom_num;
    int mon_num;
    int mol_num;
    double atom_charge;

    std::string memb_string = "MEMB";
    std::string dmpc_string = "DMPC";
    std::string dmpg_string = "DMPG";
    std::string pgl1_string = "PGL1";
    std::string pgl2_string = "PGL2";
    std::string pgl3_string = "PGL3";
    std::string pgl4_string = "PGL4";
    std::string wat_string  = "WAT ";
    std::string sod_string  = "SOD ";
    std::string cla_string  = "CLA ";

    const char * memb_str = memb_string.c_str();
    const char * dmpc_str = dmpc_string.c_str();
    const char * dmpg_str = dmpg_string.c_str();
    const char * pgl1_str = pgl1_string.c_str();
    const char * pgl2_str = pgl2_string.c_str();
    const char * pgl3_str = pgl3_string.c_str();
    const char * pgl4_str = pgl4_string.c_str();
    const char * wat_str =  wat_string.c_str();
    const char * sod_str =  sod_string.c_str();
    const char * cla_str =  cla_string.c_str();



    if (strlen(line) < sizeof(MIN_PSF_ATOM_SIZE)) {
        // Empty line continue.
        return -1;
    }
    memcpy(tmp_atom_string, &(line[3]), 5);
    tmp_atom_string[5] = 0;
    atom_num = atoi(tmp_atom_string);
    if (atom_num != (prev_atom_num + 1)) {
        printf("EXIT ATOM NUM %d %d\n", atom_num, prev_atom_num);
        return 0;
    }
    memcpy(tmp_mol_num_string, &(line[14]), 2);
    tmp_mol_num_string[2] = 0;
    mol_num = atoi(tmp_mol_num_string);


    memcpy(tmp_charge_string, &(line[35]), 9);
    tmp_charge_string[9] = 0;
    atom_charge = atof(tmp_charge_string);

    memcpy(tmp_mol_type_string, &(line[9]), 4);
    tmp_mol_type_string[4] = 0;
    if (strcmp(tmp_mol_type_string, memb_str) == 0) {
        memcpy(tmp_mol_type_string, &(line[19]), 4);
        tmp_mol_type_string[4] = 0;
    }


    memcpy(tmp_mon_num_string, &(line[14]), 4);
    tmp_mon_num_string[4] = 0;
    mon_num = atof(tmp_mon_num_string);

    memcpy(tmp_aa_type_string, &(line[19]), 4);
    tmp_aa_type_string[4] = 0;


    atom_ptr = (PSF_ATOM_PTR) malloc(sizeof(PSF_ATOM));
    atom_ptr->atom_num = atom_num;
    atom_ptr->mon_num = mon_num;
    atom_ptr->atom_charge = atom_charge;
    memcpy(atom_ptr->aa_type, tmp_aa_type_string,5);
    
    memcpy(atom_ptr->mol_type, tmp_mol_type_string,5);

    if (memcmp(dmpc_str, atom_ptr->mol_type, 4) == 0) {
        atom_ptr->mol_type_val = DMPC_MOL;
        if (mol_num <= NUM_LIPID_LEAFLET) {
            atom_ptr->leaf = LEAF_UP;
        } else {
            atom_ptr->leaf = LEAF_LO;
        }
    } else if (memcmp(dmpg_str, atom_ptr->mol_type, 4) == 0) {
        atom_ptr->mol_type_val = DMPG_MOL;
        if (mol_num <= NUM_LIPID_LEAFLET) {
            atom_ptr->leaf = LEAF_UP;
        } else {
            atom_ptr->leaf = LEAF_LO;
        }
    } else if (memcmp(pgl1_str, atom_ptr->mol_type, 4) == 0) {
        atom_ptr->mol_type_val = PGL1_MOL;
        atom_ptr->leaf = LEAF_UP;
    } else if (memcmp(pgl2_str, atom_ptr->mol_type, 4) == 0) {
        atom_ptr->mol_type_val = PGL2_MOL;
        atom_ptr->leaf = LEAF_UP;
    } else if (memcmp(pgl3_str, atom_ptr->mol_type, 4) == 0) {
        atom_ptr->mol_type_val = PGL3_MOL;
        atom_ptr->leaf = LEAF_LO;
    } else if (memcmp(pgl4_str, atom_ptr->mol_type, 4) == 0) {
        atom_ptr->mol_type_val = PGL4_MOL;
        atom_ptr->leaf = LEAF_LO;
    } else if (memcmp(wat_str, atom_ptr->mol_type, 4) == 0) {
        atom_ptr->mol_type_val = WAT_MOL;
        atom_ptr->leaf = LEAF_UN;
    } else if (memcmp(sod_str, atom_ptr->mol_type, 4) == 0) {
        atom_ptr->mol_type_val = SOD_MOL;
        atom_ptr->leaf = LEAF_UN;
    } else if (memcmp(cla_str, atom_ptr->mol_type, 4) == 0) {
        atom_ptr->mol_type_val = CLA_MOL;
        atom_ptr->leaf = LEAF_UN;
    }

    memcpy(tmp_mon_type_string, &(line[24]), 4);
    tmp_mon_type_string[4] = 0;
    memcpy(atom_ptr->monomer_type, tmp_mon_type_string,5);
    if ((atom_ptr->mol_type_val  == PGL1_MOL) || 
      (atom_ptr->mol_type_val  == PGL2_MOL) ||
      (atom_ptr->mol_type_val  == PGL3_MOL) ||
      (atom_ptr->mol_type_val  == PGL4_MOL)) {
        if ((strncmp(tmp_mon_type_string, "N ", 2) == 0) ||
          (strncmp(tmp_mon_type_string, "HT", 2) == 0) ||
          (strncmp(tmp_mon_type_string, "HA", 2) == 0) ||
          (strncmp(tmp_mon_type_string, "HN", 2) == 0) ||
          (strncmp(tmp_mon_type_string, "CA ", 3) == 0) ||
          (strncmp(tmp_mon_type_string, "C ", 2) == 0) ||
          (strncmp(tmp_mon_type_string, "O ", 2) == 0)) {
            atom_ptr->atom_type = BB_MON;
        } else {
            atom_ptr->atom_type = SIDE_MON;
        }
    } else {
        atom_ptr->atom_type = NOT_PEP_MON;
    }
    
//    printf("MOL %d SIDE = %d AA TYPE %s MON_TYPE %s AA NUM %d %s\n", atom_ptr->atom_num, 
//      atom_ptr->atom_type, atom_ptr->aa_type, atom_ptr->monomer_type, atom_ptr->mon_num, tmp_mon_num_string);




    atoms.push_back(atom_ptr);

    return 1;
}

