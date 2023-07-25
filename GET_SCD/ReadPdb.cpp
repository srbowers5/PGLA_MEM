#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include <stdarg.h>
#include <iostream>
#include <fstream>
#include <arpa/inet.h>
#include <vector>
#include "ReadPdb.h"

void rem_lead_space(char *outString, char *inString, int max_len) {
     int idx = 0;
     while  (inString[idx] == ' ') {
         idx++;
     }
     strncpy(outString, &(inString[idx]), max_len);


}

unsigned char get_first_char(char * inString) {
    int i = 0;
    unsigned char next_char;
    next_char = inString[i];
    while (next_char == ' ') {
        i++;
        next_char = inString[i];
    }
    return next_char;
}


int ReadPdb::find_heavy_atoms(int lipids_per_leaf) {
    PDB_ATOM_PTR atom_ptr;
    char tmp_field[128];
    char chain[128];
    char lipid_type[128];
    unsigned char atom_type;
    char atom_type_string[128];

    int cnt = 0;
    for (std::vector<PDB_ATOM_PTR>::iterator it = std::begin(atoms); it != std::end(atoms); ++it) {
//        atom_ptr=atoms[i];
        atom_ptr=*it;
        cnt++;
        rem_lead_space(chain, atom_ptr->chain, 128);
        if ( (strcmp(chain, "MEMB")) == 0) {
            atom_type = get_first_char(atom_ptr->atom_type);
            if (atom_type != 'H') {
                rem_lead_space(lipid_type, atom_ptr->monomer_type, 128);
                if (atom_ptr->mol_num <= lipids_per_leaf) {
                    if ( (strcmp(lipid_type, "DMPC")) == 0) {
                        dmpc_up_heavy_atoms.push_back(atom_ptr->atom_num);
                        if (atom_type == 'P') {
                            dmpc_up_P_atoms.push_back(atom_ptr->atom_num);
                        }
                    } else if ( (strcmp(lipid_type, "DMPG")) == 0) {
                        dmpg_up_heavy_atoms.push_back(atom_ptr->atom_num);
                        if (atom_type == 'P') {
                            dmpg_up_P_atoms.push_back(atom_ptr->atom_num);
                        }
                    } else if ( (strcmp(lipid_type, "DMCE")) == 0) {
                        dmpe_up_heavy_atoms.push_back(atom_ptr->atom_num);
                        if (atom_type == 'P') {
                            dmpe_up_P_atoms.push_back(atom_ptr->atom_num);
                        }
                    }
                    lip_up_heavy_atoms.push_back(atom_ptr->atom_num);
                    if (atom_type == 'P') {
                        lip_up_Phos.push_back(atom_ptr->atom_num);
                    }
                } else {
                    if ( (strcmp(lipid_type, "DMPC")) == 0) {
                        dmpc_low_heavy_atoms.push_back(atom_ptr->atom_num);
                        if (atom_type == 'P') {
                            dmpc_low_P_atoms.push_back(atom_ptr->atom_num);
                        }
                    } else if ( (strcmp(lipid_type, "DMPG")) == 0) {
                        dmpg_low_heavy_atoms.push_back(atom_ptr->atom_num);
                        if (atom_type == 'P') {
                            dmpg_low_P_atoms.push_back(atom_ptr->atom_num);
                        }
                    } else if ( (strcmp(lipid_type, "DMPE")) == 0) {
                        dmpe_low_heavy_atoms.push_back(atom_ptr->atom_num);
                        if (atom_type == 'P') {
                            dmpe_low_P_atoms.push_back(atom_ptr->atom_num);
                        }
                    }
                    lip_low_heavy_atoms.push_back(atom_ptr->atom_num);
                    if (atom_type == 'P') {
                        lip_low_Phos.push_back(atom_ptr->atom_num);
                    }
                }
            }
        } else if ( (strncmp(chain, "PGL1", 4)) == 0) {
            atom_type = get_first_char(atom_ptr->atom_type);
            if (atom_type != 'H') {
                pep_up_heavy_atoms.push_back(atom_ptr->atom_num);
            }
            rem_lead_space(atom_type_string, atom_ptr->atom_type, 128);
            if (strncmp(atom_type_string, "CA", 2) == 0) {
                pep_up_CA.push_back(atom_ptr->atom_num);
            }
        } else if ( (strncmp(chain, "PGL2", 4)) == 0) {
            atom_type = get_first_char(atom_ptr->atom_type);
            if (atom_type != 'H') {
                pep_low_heavy_atoms.push_back(atom_ptr->atom_num);
            }
            rem_lead_space(atom_type_string, atom_ptr->atom_type, 128);
            if (strncmp(atom_type_string, "CA", 2) == 0) {
                pep_low_CA.push_back(atom_ptr->atom_num);
            }
        } else if ( (strncmp(chain, "WAT", 3)) == 0) {
            atom_type = get_first_char(atom_ptr->atom_type);
            if (atom_type == 'O') {
                water_heavy_atoms.push_back(atom_ptr->atom_num);
            }
        }
    }
//    printf("NUMBER OF HA %ul %ul %ul %ul\n", dmpc_up_heavy_atoms.size(), dmpc_low_heavy_atoms.size(),dmpg_up_heavy_atoms.size(), dmpg_low_heavy_atoms.size());
//    printf("NUMBER OF LIP P %d %d NUM PGLA CA %d %d NUM WAT %d\n", 
//      lip_up_Phos.size(), lip_low_Phos.size(),
//      pep_up_CA.size(), pep_low_CA.size(),
//      water_heavy_atoms.size());
}

int ReadPdb::find_lipid_tails(int lipids_per_leaf) {
    PDB_ATOM_PTR atom_ptr;
    char tmp_field[128];
    char chain[128];
    char lipid_type[128];
    unsigned char atom_type;
    char atom_type_string[128];
    char hType;
    int currC, currHx, currHy, currHz;
    LIPID_C_TAIL_PTR lipid_tail_ptr;
    int len_string;

    lipid_tail_ptr = 0;
    currC = currHx = currHy = currHz = 0;
    int cnt = 0;
    int mol_num = 0;
    for (std::vector<PDB_ATOM_PTR>::iterator it = std::begin(atoms); it != std::end(atoms); ++it) {
//        atom_ptr=atoms[i];
        atom_ptr=*it;
        cnt++;
        rem_lead_space(chain, atom_ptr->chain, 128);

        if ( (strcmp(chain, "MEMB")) == 0) {
            rem_lead_space(atom_type_string, atom_ptr->atom_type, 128);
            if (atom_type_string[0] == 'C') {
                if (lipid_tail_ptr != 0) {
                    if (mol_num <= lipids_per_leaf) {
                        lipid_tail_c_h_up.push_back(lipid_tail_ptr);
                    } else {
                        lipid_tail_c_h_low.push_back(lipid_tail_ptr);
                    }
                    lipid_tail_ptr = 0;
                }
                mol_num = atom_ptr->mol_num;
                rem_lead_space(lipid_type, atom_ptr->monomer_type, 128);

                if ((atom_type_string[1] == '2') && (atom_type_string[2] != 0) && (atom_type_string[2] != ' ')) {
                    lipid_tail_ptr = (LIPID_C_TAIL_PTR) malloc(sizeof(LIPID_C_TAIL));
                    lipid_tail_ptr->chain_num = 1;
                    lipid_tail_ptr->chain_pos =  atoi( &(atom_type_string[2]));
                    lipid_tail_ptr->cAtom_num = atom_ptr->atom_num;
                    lipid_tail_ptr->hxAtom_num = 0;
                    lipid_tail_ptr->hyAtom_num = 0;
                    lipid_tail_ptr->hzAtom_num = 0;
                    if ( (strcmp(lipid_type, "DMPC")) == 0) {
                        lipid_tail_ptr->mon_type = IS_DMPC;
                    } else if ( (strcmp(lipid_type, "DMPG")) == 0) {
                        lipid_tail_ptr->mon_type = IS_DMPG;
                    }
                }
                if ((atom_type_string[1] == '3') && (atom_type_string[2] != 0) && (atom_type_string[2] != ' ')) {
                    lipid_tail_ptr = (LIPID_C_TAIL_PTR) malloc(sizeof(LIPID_C_TAIL));
                    lipid_tail_ptr->chain_num = 2;
                    lipid_tail_ptr->chain_pos =  atoi( &(atom_type_string[2]));
                    lipid_tail_ptr->cAtom_num = atom_ptr->atom_num;
                    lipid_tail_ptr->hxAtom_num = 0;
                    lipid_tail_ptr->hyAtom_num = 0;
                    lipid_tail_ptr->hzAtom_num = 0;
                    if ( (strcmp(lipid_type, "DMPC")) == 0) {
                        lipid_tail_ptr->mon_type = IS_DMPC;
                    } else if ( (strcmp(lipid_type, "DMPG")) == 0) {
                        lipid_tail_ptr->mon_type = IS_DMPG;
                    }
                }

            } else if ((lipid_tail_ptr) && (atom_type_string[0] == 'H')) {
                len_string = strlen(atom_type_string);
                hType = atom_type_string[len_string-1];
                if (hType == 'R') {
                    lipid_tail_ptr->hxAtom_num = atom_ptr->atom_num;
                    if (lipid_tail_ptr->chain_num != 1) {
                        printf("ERROR - bad chain numR %s %d (%c) \n", atom_type_string, lipid_tail_ptr->chain_num, hType);
                    }
                } else if (hType == 'S') {
                    lipid_tail_ptr->hyAtom_num = atom_ptr->atom_num;
                    if (lipid_tail_ptr->chain_num != 1) {
                        printf("ERROR - bad chain numS %s %d (%c)\n", atom_type_string, lipid_tail_ptr->chain_num, hType);
                    }
                } else if (hType == 'T') {
                    lipid_tail_ptr->hzAtom_num = atom_ptr->atom_num;
                    if (lipid_tail_ptr->chain_num != 1) {
                        printf("ERROR - bad chain numT %s %d\n", atom_type_string, lipid_tail_ptr->chain_num);
                    }
                } else if (hType == 'X') {
                    lipid_tail_ptr->hxAtom_num = atom_ptr->atom_num;
                    if (lipid_tail_ptr->chain_num != 2) {
                        printf("ERROR - bad chain numX  %s %d\n", atom_type_string, lipid_tail_ptr->chain_num);
                    }
                } else if (hType == 'Y') {
                    lipid_tail_ptr->hyAtom_num = atom_ptr->atom_num;
                    if (lipid_tail_ptr->chain_num != 2) {
                        printf("ERROR - bad chain numY %s %d\n", atom_type_string, lipid_tail_ptr->chain_num);
                    }
                } else if (hType == 'Z') {
                    lipid_tail_ptr->hzAtom_num = atom_ptr->atom_num;
                    if (lipid_tail_ptr->chain_num != 2) {
                        printf("ERROR - bad chain numZ %s %d\n", atom_type_string, lipid_tail_ptr->chain_num);
                    }
                }
                atom_type_string[len_string-1] = 0;
                if (lipid_tail_ptr->chain_pos != atoi( &(atom_type_string[1]))) {
                    printf("ERROR - bad chain pos %s\n", atom_type_string);
                }
            }
        }
    }
    if (lipid_tail_ptr != 0) {
        if (mol_num <= lipids_per_leaf) {
            lipid_tail_c_h_up.push_back(lipid_tail_ptr);
        } else {
            lipid_tail_c_h_low.push_back(lipid_tail_ptr);
        }
        lipid_tail_ptr = 0;
    }
}

int ReadPdb::parse_pdb_line(const char * line) {
    ATOM_LINE_PTR atom_line_ptr;
    BOX_LINE_PTR box_line_ptr;
    PDB_ATOM_PTR atom_ptr;
    char tmp_string[32];
    int retVal;

    if (strlen(line) < sizeof(TYPE_LEN)) {
        // Empty line continue.
        return 1;
    }
//    printf("LINE %s\n", line);
    retVal = 1;
    if (strncmp(ATOM_TYPE, line, TYPE_LEN) == 0) {
        if (strlen(line) < sizeof(ATOM_LINE)) {
            printf("LINE to short of ATOM  %d %s\n", (int) strlen(line), line);
            fflush(stdout);
            exit(98);
        }
        atom_line_ptr = (ATOM_LINE_PTR) line;
        atom_ptr = (PDB_ATOM_PTR) malloc(sizeof(PDB_ATOM));

        memcpy(tmp_string, atom_line_ptr->atom_num, 5);
        tmp_string[5] = 0;
        atom_ptr->atom_num = atoi(tmp_string);

        memcpy(tmp_string, atom_line_ptr->mol_num, 5);
        tmp_string[5] = 0;
        atom_ptr->mol_num = atoi(tmp_string);

        memcpy(tmp_string, atom_line_ptr->xVal, 8);
        tmp_string[8] = 0;
        atom_ptr->xVal = atof(tmp_string);

        memcpy(tmp_string, atom_line_ptr->yVal, 8);
        tmp_string[8] = 0;
        atom_ptr->yVal = atof(tmp_string);

        memcpy(tmp_string, atom_line_ptr->zVal, 8);
        tmp_string[8] = 0;
        atom_ptr->zVal = atof(tmp_string);

        memcpy(atom_ptr->atom_type, atom_line_ptr->atom_type, 5);
        atom_ptr->atom_type[5] = 0;
//        printf("ATOM_TYPE %s\n", atom_ptr->atom_type);

        memcpy(atom_ptr->monomer_type, atom_line_ptr->monomer_type, 5);
        atom_ptr->monomer_type[5] = 0;
//        printf("MON_TYPE %s\n", atom_ptr->monomer_type);

        atom_ptr->mol_type = atom_line_ptr->mol_type;
//        printf("ATOM_TYPE %s\n", atom_ptr->atom_type);

        memcpy(atom_ptr->chain, atom_line_ptr->chain, 10);
        atom_ptr->chain[10] = 0;
//        printf("CHAIN %s\n", atom_ptr->chain);

        memcpy(atom_ptr->colA, atom_line_ptr->colA, 6);
        atom_ptr->colA[6] = 0;
        memcpy(atom_ptr->colB, atom_line_ptr->colB, 6);
        atom_ptr->colB[6] = 0;
        atoms.push_back(atom_ptr);
        retVal = 1;

    } else if (strncmp(BOX_TYPE, line, TYPE_LEN) == 0) {
        if (strlen(line) < sizeof(BOX_LINE)) {
            printf("LINE to short of BOX  %d %s\n", (int) strlen(line), line);
            exit(98);
        }
        box_line_ptr = (BOX_LINE_PTR) line;
        memcpy(tmp_string, box_line_ptr->xLen, 8);
        tmp_string[8] = 0;
        box.box[0] = atof(tmp_string);

        memcpy(tmp_string, box_line_ptr->yLen, 8);
        tmp_string[8] = 0;
        box.box[1] = atof(tmp_string);

        memcpy(tmp_string, box_line_ptr->zLen, 8);
        tmp_string[8] = 0;
        box.box[2] = atof(tmp_string);
        retVal = 1;
    } else if (strncmp(END_TYPE, line, TYPE_LEN) == 0) {
//        printf("FOUND END (%s)\n", line);
        retVal = 0; // Done with struct.
    }
    return retVal;
}

void ReadPdb::print_pdb() {
    int i;
    printf("unit cell %f %f %f\n", box.box[0], box.box[1], box.box[2]);
    for (i=0; i<atoms.size(); i++) {
        printf("%d %d %s %s %c %s %s %s %f %f %f\n", 
          atoms[i]->atom_num, atoms[i]->mol_num, atoms[i]->atom_type,
          atoms[i]->monomer_type, atoms[i]->mol_type, atoms[i]->chain,
          atoms[i]->colA, atoms[i]->colB, atoms[i]->xVal, atoms[i]->yVal,
          atoms[i]->zVal);
    }
    fflush(stdout);
}


