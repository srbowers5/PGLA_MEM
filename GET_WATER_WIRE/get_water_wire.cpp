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
#include "system_conf.h"
#include "get_water_wire.h"

//    arg[1] = Input_dir
//    arg[2] = Input_pdb
//    arg[3] = tr
//    arg[4] = rep
//    arg[5] = prefix
//    arg[6] = first_str
//    arg[7] = last_str
//    arg[8] = output_name


std::vector<WATER_MOL_PTR> water_vec;

HBOND bond_array[MAX_HBONDS];
int cnt_hbonds;

int clust_array[MAX_CLUSTS*2][MAX_HBONDS*2];
int num_clusts;
int cnt_clusts[MAX_CLUSTS*2];
float clust_min[MAX_CLUSTS*2];
float clust_max[MAX_CLUSTS*2];
float clust_diff[MAX_CLUSTS*2];

float clust_min_up[MAX_CLUSTS*2];
float clust_min_lo[MAX_CLUSTS*2];
float clust_max_uo[MAX_CLUSTS*2];
float clust_max_lo[MAX_CLUSTS*2];

/*
*   Function: get_water_mol()
*
*   Reads DCD structure of picks out the water molecules where the O atom
*   has abs(Z) < INSERT_DEPTH, and saves coordinates in vector.
*/
void get_water_mol(ReadPdb * pdb, ReadDcd * dcd) {
    std::vector<int>::iterator first_it, end_it, it;
    int water_cnt;
    int atom_num;
    float zVal;
    WATER_MOL_PTR water_ptr;

    first_it = std::begin(pdb->water_heavy_atoms);
    end_it = std::end(pdb->water_heavy_atoms);

    water_cnt = 0;
    for (it = first_it; it != end_it; ++it) {
        atom_num=*it;
        zVal = dcd->zValArray[0][atom_num-1];
        if ((zVal > -INSERT_DEPTH) && (zVal < INSERT_DEPTH)) {
            water_ptr = (WATER_MOL_PTR) malloc(sizeof(WATER_MOL));
            water_ptr->ox_mol_coord[0] = dcd->xValArray[0][atom_num-1];
            water_ptr->ox_mol_coord[1] = dcd->yValArray[0][atom_num-1];
            water_ptr->ox_mol_coord[2] = zVal;

            atom_num += 1;
            water_ptr->h1_mol_coord[0] = dcd->xValArray[0][atom_num-1];
            water_ptr->h1_mol_coord[1] = dcd->yValArray[0][atom_num-1];
            water_ptr->h1_mol_coord[2] = dcd->zValArray[0][atom_num-1];

            atom_num += 1;
            water_ptr->h2_mol_coord[0] = dcd->xValArray[0][atom_num-1];
            water_ptr->h2_mol_coord[1] = dcd->yValArray[0][atom_num-1];
            water_ptr->h2_mol_coord[2] = dcd->zValArray[0][atom_num-1];
            water_ptr->mol_index = water_cnt;
            water_cnt += 1;
            water_vec.push_back(water_ptr);
        }
    }
}

float get_dist(float x1,float y1,float z1, float x2, float y2, float z2) {
    float valx, valy, valz, val;

//    printf("IN %f %f %f,  %f %f %f\n", x1, y1, z1, x2, y2, z2);
    valx = pow((x1-x2), 2);
    valy = pow((y1-y2), 2);
    valz = pow((z1-z2), 2);

    val = pow((valx + valy + valz), 0.5);
//    printf (" VALS %f %f %f  = %f\n", valx, valy, valz, val);
    return val;
}

float get_angle( float ax, float ay, float az, float bx, float by, float bz, float cx, float cy, float cz)
{
    float abx, aby, abz;
    float bcx, bcy, bcz;
    float ablen, bclen;
    float vec1x, vec2x, vec1y, vec2y, vec1z, vec2z;
    float val, rads, degrees;
    abx = bx - ax;
    aby = by - ay;
    abz = bz - az;

    bcx = bx - cx;
    bcy = by - cy;
    bcz = bz - cz;

    ablen = pow( (pow(abx,2) + pow(aby,2) + pow(abz,2)), 0.5);
    bclen = pow( (pow(bcx,2) + pow(bcy,2) + pow(bcz,2)), 0.5);

    vec1x = abx / ablen;
    vec1y = aby / ablen;
    vec1z = abz / ablen;

    vec2x = bcx / bclen;
    vec2y = bcy / bclen;
    vec2z = bcz / bclen;

    val = ((vec1x *vec2x) + (vec1y *vec2y) + (vec1z *vec2z));
    rads = acos(val);
    degrees = (rads * 180) / M_PI;

    return degrees;
}

int add_hbond(int m1, int m2, int hVal) {

    if (cnt_hbonds > MAX_HBONDS) {
        printf("ERROR - to many HBONDS %d %d\n", cnt_hbonds, MAX_HBONDS);
        return -1;
    }

    bond_array[cnt_hbonds].mol1 = m1;
    bond_array[cnt_hbonds].mol2 = m2;
    bond_array[cnt_hbonds].hVal = hVal;
    cnt_hbonds += 1;
    return 0;

}



int find_clusts() {

    int free_hbonds;
    int curr_clust;
    int off;
    int offMin;
    int mol1Val;
    int clust_off;
    int curr_mol;

    offMin = 0;
    num_clusts = 0;
    free_hbonds = cnt_hbonds;
    offMin = 0;
    for (curr_clust=0; curr_clust<MAX_CLUSTS; curr_clust++) {
        mol1Val = bond_array[offMin].mol1;
        while (mol1Val < 0) {
            offMin += 1;
            mol1Val = bond_array[offMin].mol1;
        }
        clust_array[curr_clust][0] = bond_array[offMin].mol1;
        clust_array[curr_clust][1] = bond_array[offMin].mol2;
        cnt_clusts[curr_clust] = 2;


        bond_array[offMin].mol1 = -1;
        free_hbonds--;
        offMin++;

        clust_off = 0;
        while (clust_off < cnt_clusts[curr_clust]) {
            curr_mol = clust_array[curr_clust][clust_off];
            for (off=offMin ; off<cnt_hbonds; off++) {
                if (bond_array[off].mol1 < 0) {
                    continue;
                }
                if (bond_array[off].mol1 == curr_mol) {
                    clust_array[curr_clust][cnt_clusts[curr_clust]] = bond_array[off].mol2;
                    cnt_clusts[curr_clust]++;
                    bond_array[off].mol1 = -1;
                    free_hbonds--;
                } else if (bond_array[off].mol2 == curr_mol) {
                    clust_array[curr_clust][cnt_clusts[curr_clust]] = bond_array[off].mol1;
                    cnt_clusts[curr_clust]++;
                    bond_array[off].mol1 = -1;
                    free_hbonds--;
                }
            }
            clust_off++;
        }
        num_clusts++;

        if (free_hbonds == 0) {
            break;
        }
    }
    return 0;
}


int find_hbonds() {
    int vector_len;
    float oxx, oxy, oxz;
    float oxx2, oxy2, oxz2;
    float h1x, h1y, h1z;
    float h2x, h2y, h2z;
    float angle_deg;
    int i,j;
    float dist;
    int retCode;

    vector_len = water_vec.size();
    cnt_hbonds = 0;
    for (i=0; i<vector_len; i++) {
        oxx = water_vec[i]->ox_mol_coord[0];
        oxy = water_vec[i]->ox_mol_coord[1];
        oxz = water_vec[i]->ox_mol_coord[2];
        for (j=0; j<vector_len; j++) {
            if (i==j) {
                continue;
            }
            h1x = water_vec[j]->h1_mol_coord[0];
            h1y = water_vec[j]->h1_mol_coord[1];
            h1z = water_vec[j]->h1_mol_coord[2];

            h2x = water_vec[j]->h2_mol_coord[0];
            h2y = water_vec[j]->h2_mol_coord[1];
            h2z = water_vec[j]->h2_mol_coord[2];

            dist = get_dist(oxx, oxy, oxz, h1x, h1y, h1z);
            if (dist < MAX_HBOND_DIST) {
                oxx2 = water_vec[j]->ox_mol_coord[0];
                oxy2 = water_vec[j]->ox_mol_coord[1];
                oxz2 = water_vec[j]->ox_mol_coord[2];
                angle_deg = get_angle(oxx, oxy, oxz, h1x, h1y, h1z, oxx2, oxy2, oxz2);
                if (angle_deg < MAX_HBOND_ANGLE) {
                    add_hbond(i,j,1);
                }
            }

           dist = get_dist(oxx, oxy, oxz, h2x, h2y, h2z);
            if (dist < MAX_HBOND_DIST) {
                oxx2 = water_vec[j]->ox_mol_coord[0];
                oxy2 = water_vec[j]->ox_mol_coord[1];
                oxz2 = water_vec[j]->ox_mol_coord[2];
                angle_deg = get_angle(oxx, oxy, oxz, h2x, h2y, h2z, oxx2, oxy2, oxz2);
                if (angle_deg < MAX_HBOND_ANGLE) {
                    retCode = add_hbond(i,j,2);
                    if (retCode < 0) {
                        printf("FAILED HBONDS\n");
                        return -1;
                    }
                }
            }
        }
    }
    return 0;
}



void get_clust_min_max() {
    int i, j;
    float minZ, maxZ;
    float minZup, maxZup;
    float minZlo, maxZlo;
    float val;
    WATER_MOL_PTR water_ptr;

    for (i=0; i<num_clusts; i++) {
        minZ = 9999.0;
        maxZ = -9999.0;
        for (j=0; j<cnt_clusts[i]; j++) {
            water_ptr = water_vec[clust_array[i][j]];
            val = water_ptr->ox_mol_coord[2];
            if (val < minZ) {
                minZ = val;
            }
            if (val > maxZ) {
                maxZ = val;
            }
        }
        clust_diff[i] = maxZ - minZ;
        clust_min[i] = minZ;
        clust_max[i] = maxZ;
    }
}

void write_long_clust(FILE * fp, int curr_step) {
    float max_diff;
    float max_diffHi;
    float max_diffLo;
    int max_clust;
    int max_clustHi;
    int max_clustLo;
    bool clustUp;
    int i;

    max_diff = 0.0;
    max_clust = 0;
    max_diffHi = 0.0;
    max_clustHi = 0;
    max_diffLo = 0.0;
    max_clustLo = 0;
    for (i=0; i<num_clusts; i++) {
        if (clust_min[i] > 0) {
            clustUp = true;
        } else if (clust_max[i] < 0)   {
            clustUp = false;
        } else if (clust_max[i] > -clust_min[i]) {
            clustUp = true;
        } else {
            clustUp = false;
        }

        if (clust_diff[i] > max_diff) {
            max_diff = clust_diff[i];
            max_clust = i;
        }
        if (clustUp == true) {
            if (clust_diff[i] > max_diffHi) {
                max_diffHi = clust_diff[i];
                max_clustHi = i;
            }
        } else {
            if (clust_diff[i] > max_diffLo) {
                max_diffLo = clust_diff[i];
                max_clustLo = i;
            }
        }
    }
    fprintf(fp, "%d %f %f %f %f %f %f %f %f %f\n", curr_step, 
      clust_diff[max_clust], clust_min[max_clust], clust_max[max_clust],
      clust_diff[max_clustLo], clust_min[max_clustLo], clust_max[max_clustLo],
      clust_diff[max_clustHi], clust_min[max_clustHi], clust_max[max_clustHi]);
}



int main(int argc, char * argv[]) {

    char dcd_file[MAX_FILENAME_LEN];
    char pdb_file[MAX_FILENAME_LEN];
    char water_wire_file[MAX_FILENAME_LEN];
    WATER_MOL_PTR water_ptr;
    FILE * out_fp;
    int retCode;
    int i,j;

    int first_step, last_step, step;
    ReadDcd * dcd;


    ReadPdb::PEP_STRUCT_PTR pep_ptr;


    snprintf(pdb_file, MAX_FILENAME_LEN, "%s/%s", argv[1], argv[2]);
    printf("PDB %s\n", pdb_file);

    ReadPdb pdb(pdb_file);
    pdb.find_heavy_atoms(NUM_LIPID_LEAFLET, NUM_PEP_LEAFLET);


    first_step = atoi(argv[6]);
    last_step = atoi(argv[7]);

    snprintf(water_wire_file, MAX_FILENAME_LEN, "%s/%s", argv[1], argv[8]);
    printf("OutFile %s\n", water_wire_file);
    out_fp = fopen(water_wire_file,"w");


    for (step=first_step; step < (last_step+1); step++) {
        snprintf(dcd_file, MAX_FILENAME_LEN, "%s/output/%s%s_%05d_%s.dcd", argv[1], argv[5], argv[3], step, argv[4]);
//        printf("DCD %s\n", dcd_file);
        dcd = new ReadDcd(dcd_file);
        get_water_mol(&pdb, dcd);

//        std::vector<WATER_MOL_PTR>::iterator first_it, end_it, it;
//        first_it = std::begin(water_vec);
//        end_it = std::end(water_vec);

        retCode = find_hbonds();


        int tot_mol_num;
//        printf("DO FIND CLUSTS \n");
        fflush(stdin);
        find_clusts();
        get_clust_min_max();
/****
        tot_mol_num = 0;
        printf("NUM_CLUSTS %d\n", num_clusts);
        for (i=0; i<num_clusts; i++) {
            printf("CLUST %d size %d %f %f %f\n", (i+1), cnt_clusts[i], clust_diff[i], clust_min[i], clust_max[i]);
            tot_mol_num += cnt_clusts[i];
            for (j=0; j<cnt_clusts[i]; j++) {
                printf("   %d\n", clust_array[i][j]);
            }
        }
        printf("GET %d clusts %d mol %d %d hbonds\n", num_clusts, tot_mol_num, cnt_hbonds, (tot_mol_num - num_clusts) );
*****/
        write_long_clust(out_fp, step);

        while (!water_vec.empty())
        {
            water_ptr = water_vec.back();
            water_vec.pop_back();
            free(water_ptr);
        }
        delete dcd;
        if (retCode < 0) {
            printf("ERROR\n");
        }
    }
    fclose(out_fp);

}
