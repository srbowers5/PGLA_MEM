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
#include "system_conf.h"
#include "dcdStruct.h"
#include "ReadPdb.h"
#include "ReadDcd.h"
#include "ReadCM.h"
#include "get_lipid_den.h"

int Num_rem;

int dmpc_OPP_p_dist[NUM_BUCKETS];
int dmpg_OPP_p_dist[NUM_BUCKETS];
int dmpc_p_dist[NUM_BUCKETS];
int dmpg_p_dist[NUM_BUCKETS];
int dmpc_p_up_dist[NUM_BUCKETS];
int dmpg_p_up_dist[NUM_BUCKETS];
int dmpc_p_lo_dist[NUM_BUCKETS];
int dmpg_p_lo_dist[NUM_BUCKETS];
int pep_dist[NUM_BUCKETS];
int pep_dist_fixed[NUM_BUCKETS];
int pep_dist_bound[NUM_BUCKETS];
int pep_OPP_dist[NUM_BUCKETS];
int pep_OPP_dist_fixed[NUM_BUCKETS];
int pep_OPP_dist_bound[NUM_BUCKETS];

int dmpc_OPP_den_reg[NUM_REG];
int dmpg_OPP_den_reg[NUM_REG];
int dmpc_SAM_den_reg[NUM_REG];
int dmpg_SAM_den_reg[NUM_REG];
int dmpc_SAM_vol_fix_den_reg[NUM_REG];
int dmpc_SAM_vol_ref_den_reg[NUM_REG];
int dmpc_OPP_vol_fix_den_reg[NUM_REG];
int dmpc_OPP_vol_ref_den_reg[NUM_REG];
int dmpg_SAM_vol_fix_den_reg[NUM_REG];
int dmpg_SAM_vol_ref_den_reg[NUM_REG];
int dmpg_OPP_vol_fix_den_reg[NUM_REG];
int dmpg_OPP_vol_ref_den_reg[NUM_REG];

int lipid_cont_dist_hist[NUM_BUCKETS];
int dmpc_dist_hist[NUM_BUCKETS];
int dmpg_dist_hist[NUM_BUCKETS];

int sod_ion_dist_hist[NUM_BUCKETS];
int cla_ion_dist_hist[NUM_BUCKETS];

int dmpc_same_dist_hist[NUM_BUCKETS];
int dmpc_same_dist_hist_fixed[NUM_BUCKETS];
int dmpc_same_dist_hist_bound[NUM_BUCKETS];


int dmpc_opp_dist_hist[NUM_BUCKETS];
int dmpc_opp_dist_hist_fixed[NUM_BUCKETS];
int dmpc_opp_dist_hist_bound[NUM_BUCKETS];

int dmpg_same_dist_hist[NUM_BUCKETS];
int dmpg_same_dist_hist_fixed[NUM_BUCKETS];
int dmpg_same_dist_hist_bound[NUM_BUCKETS];

int dmpg_opp_dist_hist[NUM_BUCKETS];
int dmpg_opp_dist_hist_fixed[NUM_BUCKETS];
int dmpg_opp_dist_hist_bound[NUM_BUCKETS];


int dmpc_SAM_xy_reg[NUM_BUCKETS][4];
int dmpg_SAM_xy_reg[NUM_BUCKETS][4];
int dmpc_OPP_xy_reg[NUM_BUCKETS][4];
int dmpg_OPP_xy_reg[NUM_BUCKETS][4];


int cnt_dmpc_SAM_xy[4];
int cnt_dmpc_OPP_xy[4];
int cnt_dmpg_SAM_xy[4];
int cnt_dmpg_OPP_xy[4];

int dmpc_2d_dist_hist[NUM_BUCKETS][NUM_Z_BUCKETS];
int dmpg_2d_dist_hist[NUM_BUCKETS][NUM_Z_BUCKETS];
int lipid_2d_dist_hist[NUM_BUCKETS][NUM_Z_BUCKETS];
int water_2d_dist_hist[NUM_BUCKETS][NUM_Z_BUCKETS];

int peptide_dist_hist[NUM_BUCKETS];
int peptide_2d_dist_hist[NUM_BUCKETS][NUM_Z_BUCKETS];

float ref_bound[NUM_BUCKETS];
float opp_bound[NUM_BUCKETS];

float    near_fix_vol;
float    near_ref_vol;
float    near_fix_opp_vol;
float    near_ref_opp_vol;

float    prox_fix_vol;
float    prox_ref_vol;
float    prox_fix_opp_vol;
float    prox_ref_opp_vol;

float    far_fix_vol;
float    far_ref_vol;
float    far_fix_opp_vol;
float    far_ref_opp_vol;

float pc_region_vol_den_step[NUM_REG], pg_region_vol_den_step[NUM_REG], pc_region_area_den_step[NUM_REG], pg_region_area_den_step[NUM_REG];


void get_xy_dists(float x1Diff, float x2Diff, float y1Diff, float y2Diff, float max_xyDist, float *d1,float *d2,float *d3,float *d4);


void read_bound() {
    int i;
    float dist;
    float tmp_ref_bound[NUM_BUCKETS];
    float tmp_opp_bound[NUM_BUCKETS];
    FILE * fp;
    float myArea;
    float junk1, junk2;
  
    fp = fopen ("lipid_bound_r_all.dat", "r");
    for (i=0; i< NUM_BUCKETS; i++) {
        fscanf(fp, "%f %f %f %f %f", &dist, &(tmp_opp_bound[i]), &(tmp_ref_bound[i]), &junk1, &junk2);
        ref_bound[i] = abs(tmp_ref_bound[i]);
        opp_bound[i] = abs(tmp_opp_bound[i]);
    }

    fclose(fp);
    printf("REF DISTENCES ");
    for (i=0; i<NUM_BUCKETS; i++) {
        printf("%f ", ref_bound[i]);
    }
    printf("\n");

    printf("OPP DISTENCES ");
    for (i=0; i<NUM_BUCKETS; i++) {
        printf("%f ", opp_bound[i]);
    }
    printf("\n");

    near_fix_vol = 0.0;
    near_fix_vol = 0.0;
    near_ref_opp_vol = 0.0;
    near_ref_opp_vol = 0.0;
    for (i=0; i< MAX_Z_DIST; i++) {
        myArea =  (pow((i+1),2) - pow(i, 2)) * M_PI;
        if (i < PROXIMAL_DIST) {
            near_fix_vol += myArea * FIXED_LIPID_BOUND;
            near_ref_vol += myArea * ref_bound[i];
            near_fix_opp_vol += myArea * FIXED_LIPID_BOUND;
            near_ref_opp_vol += myArea * opp_bound[i];
        } else if (i<FAR_DIST) {
            prox_fix_vol += myArea * FIXED_LIPID_BOUND;
            prox_ref_vol += myArea * ref_bound[i];
            prox_fix_opp_vol += myArea * FIXED_LIPID_BOUND;
            prox_ref_opp_vol += myArea * opp_bound[i];
        } else if (i<MAX_Z_DIST) {
            far_fix_vol += myArea * FIXED_LIPID_BOUND;
            far_ref_vol += myArea * ref_bound[i];
            far_fix_opp_vol += myArea * FIXED_LIPID_BOUND;
            far_ref_opp_vol += myArea * opp_bound[i];
        }
    }

}


void add_cont_dist_hist(float * distArray, unsigned char * lip_cont_array, int size_array, int * hist_array) {
    int i;
    float val;
    int histBucket = 0;
    for (i=0; i<size_array; i++) {
        if (lip_cont_array[i] != 0) {
////            printf("LIPID %d contact DIST=%f\n", i, distArray[i]);
            val = (distArray[i] * BUCKETS_PER_DIST) + 0.00001;
            histBucket = (int) val;
            if (histBucket >= NUM_BUCKETS) {
//                histBucket = NUM_BUCKETS-1;
                continue;
            }
            hist_array[histBucket] += 1;
//            printf("ADD COM DIST %d val %d\n", histBucket, hist_array[histBucket]);
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


void add_2d_ion_hist(mol_defines mol_type, bool isUpper, float pepXVal, float pepYVal, float pepZVal, ReadPdb * pdb, ReadDcd * dcd ) {

    float xCell;
    float x1Diff, y1Diff, z1Diff, x2Diff, y2Diff;
    float xyDist, xyVal;
    float xyDists[4];
    int histBucket;
    std::vector<int>::iterator first_it, end_it, it;
    int atom_num;
    float xVal, yVal, zVal;

    if (mol_type == IS_CLA_ION) {
        first_it = std::begin(pdb->cl_ion);
        end_it = std::end(pdb->cl_ion);

    } else if (mol_type == IS_SOD_ION) {
        first_it = std::begin(pdb->sod_ion);
        end_it = std::end(pdb->sod_ion);
    } else {
        return;
    }
    for (it = first_it; it != end_it; ++it) {
        atom_num=*it;
        xVal = dcd->xValArray[0][atom_num-1];
        yVal = dcd->yValArray[0][atom_num-1];
        zVal = dcd->zValArray[0][atom_num-1];
        if ((zVal < 0) && (isUpper == true)) {
            continue;
        }
        if ((zVal > 0) && (isUpper == false)) {
            continue;
        }

        xCell = dcd->uCell[0];
        x1Diff = abs(xVal - pepXVal);
        y1Diff = abs(yVal - pepYVal);
        z1Diff = abs(zVal - pepZVal);

        x2Diff = abs(xCell - x1Diff);
        y2Diff = abs(xCell - y1Diff);
        xyDists[0] =  xyDists[1] =  xyDists[2] =  xyDists[3] = MAX_XY_DIST;
        get_xy_dists(x1Diff, x2Diff, y1Diff, y2Diff, MAX_XY_DIST, &xyDists[0], &xyDists[1], &xyDists[2], &xyDists[3]);
        for (int i=0; i<4; i++) {
            xyDist = xyDists[i];
            if (xyDist >= MAX_XY_DIST) {
                continue;
            }
            xyVal = (xyDist * BUCKETS_PER_DIST) + 0.00001;
            histBucket = (int) xyVal;
            if (histBucket >= NUM_BUCKETS) {
                continue;
            }
            if (mol_type == IS_CLA_ION) {
                cla_ion_dist_hist[histBucket]++;
            } else {
                sod_ion_dist_hist[histBucket]++;
            }
        }
    }
}


void add_2d_hist(mol_defines mol_type, bool isUpper, float xyDist, float zVal, float pep_xy_dist) {

    int histBucket;
    int histZBucket;
    int reg_index;
    int pep_xy_index;
    float xyVal;
    float zValAbs;
    if ((zVal > MAX_Z_DIST) | (zVal < MIN_Z_DIST)) {
         return;
    }

    if (xyDist > MAX_XY_DIST) {
        return;
    }

    zValAbs = abs(zVal);
/**** OLD JUNK
    if (xyDist > FAR_DIST) {
        reg_index = 3;
    } else if (xyDist > DISTANT_DIST) {
        reg_index = 2;
    } else if (xyDist > PROXIMAL_DIST) {
        reg_index = 1;
    } else {
        reg_index = 0;
    }
**   END OLD JUNK  *****/

    if (xyDist < NEAR_END) {
        reg_index = 0;
    } else if (xyDist < PROXIMAL_END) {
        reg_index = 1;
    } else if (xyDist < FAR_END) {
        reg_index = 2;
    } else {
        reg_index = 3;
    }


    pep_xy_index = 0;
    if (pep_xy_dist > 10.0) {
        if (pep_xy_dist < 20.0) {
            pep_xy_index = 1;
        } else if (pep_xy_dist < 30) {
            pep_xy_index = 2;
        }
        else {
            pep_xy_index = 3;
        }
    }
//    printf("XY DIST %f index %d\n", pep_xy_dist, pep_xy_index);

    if (isUpper == false) {
        zVal = -zVal;
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
//    if (isUpper) {
//        printf("Upper Z VAL = %f bucket %d XY %f (%d)\n", zVal, histZBucket, xyDist, histBucket);
//    } else {
//        printf("Lower Z VAL = %f bucket %d XY %f (%d)\n", zVal, histZBucket, xyDist, histBucket);
//    }
    if (mol_type == IS_WATER) {
        water_2d_dist_hist[histBucket][histZBucket]++;
    } else if (mol_type == IS_DMPC) {
        dmpc_2d_dist_hist[histBucket][histZBucket]++;
        dmpc_dist_hist[histBucket]++;
        if (zVal > 0) {
            dmpc_same_dist_hist[histBucket]++;
            if (zValAbs < FIXED_LIPID_BOUND) {
                dmpc_same_dist_hist_fixed[histBucket]++;
                dmpc_SAM_vol_fix_den_reg[reg_index]++;
                pc_region_vol_den_step[reg_index]++;
            }
            if (zValAbs < ref_bound[histBucket]) {
                dmpc_same_dist_hist_bound[histBucket]++;
                dmpc_SAM_vol_ref_den_reg[reg_index]++;
            }
            dmpc_SAM_den_reg[reg_index]++;
            dmpc_SAM_xy_reg[histBucket][pep_xy_index]++;
            cnt_dmpc_SAM_xy[pep_xy_index]++;
        } else {
            dmpc_opp_dist_hist[histBucket]++;
            if (zValAbs < FIXED_LIPID_BOUND) {
                dmpc_opp_dist_hist_fixed[histBucket]++;
                dmpc_OPP_vol_fix_den_reg[reg_index]++;
            }
            if (zValAbs < opp_bound[histBucket]) {
                dmpc_opp_dist_hist_bound[histBucket]++;
                dmpc_OPP_vol_ref_den_reg[reg_index]++;
            }
            dmpc_OPP_den_reg[reg_index]++;
            dmpc_OPP_xy_reg[histBucket][pep_xy_index]++;
            cnt_dmpc_OPP_xy[pep_xy_index]++;
        }
        lipid_2d_dist_hist[histBucket][histZBucket]++;


    } else if (mol_type == IS_DMPG) {
        dmpg_2d_dist_hist[histBucket][histZBucket]++;
        dmpg_dist_hist[histBucket]++;
        if (zVal > 0) {
            dmpg_same_dist_hist[histBucket]++;
            if (zValAbs < FIXED_LIPID_BOUND) {
                dmpg_same_dist_hist_fixed[histBucket]++;
                dmpg_SAM_vol_fix_den_reg[reg_index]++;
                pg_region_vol_den_step[reg_index]++;
            }
            if (zValAbs < ref_bound[histBucket]) {
                dmpg_same_dist_hist_bound[histBucket]++;
                dmpg_SAM_vol_ref_den_reg[reg_index]++;
            }
            dmpg_SAM_den_reg[reg_index]++;
            dmpg_SAM_xy_reg[histBucket][pep_xy_index]++;
            cnt_dmpg_SAM_xy[pep_xy_index]++;
        } else {
            dmpg_opp_dist_hist[histBucket]++;
            if (zValAbs < FIXED_LIPID_BOUND) {
                dmpg_opp_dist_hist_fixed[histBucket]++;
                dmpg_OPP_vol_fix_den_reg[reg_index]++;
            }
            if (zValAbs < ref_bound[histBucket]) {
                dmpg_opp_dist_hist_bound[histBucket]++;
                dmpg_OPP_vol_ref_den_reg[reg_index]++;
            }
            dmpg_OPP_den_reg[reg_index]++;
            dmpg_OPP_xy_reg[histBucket][pep_xy_index]++;
            cnt_dmpg_OPP_xy[pep_xy_index]++;
        }
        lipid_2d_dist_hist[histBucket][histZBucket]++;
    } else if (mol_type == IS_DMPE) {
        printf("DMPE not Now supported. Add code HERE\n");
        lipid_2d_dist_hist[histBucket][histZBucket]++;
    } else {
        printf("OTHER\n");
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


void add_peptide_dist_hist(float x1, float y1, float z1, float x2, float y2, float z2, ReadDcd * dcd) {
    float xCell, x1Diff, y1Diff, z1Diff, x2Diff, y2Diff;
    float xyDists[4];
    int histBucket;
    int zBucket;
    int i;

    xCell = dcd->uCell[0];
    x1Diff = abs(x1 - x2);
    y1Diff = abs(y1 - y2);
    z1Diff = abs(z1 - z2);

    zBucket = int(z1Diff + 0.00001);
    if (zBucket >= NUM_Z_BUCKETS) {
        zBucket = NUM_Z_BUCKETS;
    }

    x2Diff = abs(xCell - x1Diff);
    y2Diff = abs(xCell - y1Diff);
    xyDists[0] =  xyDists[1] =  xyDists[2] =  xyDists[3] = MAX_XY_DIST;
    get_xy_dists(x1Diff, x2Diff, y1Diff, y2Diff, MAX_XY_DIST, &xyDists[0], &xyDists[1], &xyDists[2], &xyDists[3]);
    for (i=0; i<4;i++) {
        if (xyDists[i] < MAX_XY_DIST) {
            histBucket = int (xyDists[i] + 0.00001);
            if (histBucket < NUM_BUCKETS) {
                peptide_dist_hist[histBucket]++;
                peptide_2d_dist_hist[histBucket][zBucket]++;
            }
        }
    }
}

float get_peptide_dist(float x1, float y1, float z1, float x2, float y2, float z2, ReadDcd * dcd) {
   float xCell, x1Diff, y1Diff, z1Diff, x2Diff, y2Diff;
    float xyDists[4];
    float min_dist;
    int histBucket;
    int zBucket;
    int i;

    xCell = dcd->uCell[0];
    x1Diff = abs(x1 - x2);
    y1Diff = abs(y1 - y2);
    z1Diff = abs(z1 - z2);

    zBucket = int(z1Diff + 0.00001);
    if (zBucket >= NUM_Z_BUCKETS) {
        zBucket = NUM_Z_BUCKETS;
    }

    x2Diff = abs(xCell - x1Diff);
    y2Diff = abs(xCell - y1Diff);
    xyDists[0] =  xyDists[1] =  xyDists[2] =  xyDists[3] = MAX_XY_DIST;
    get_xy_dists(x1Diff, x2Diff, y1Diff, y2Diff, MAX_XY_DIST, &xyDists[0], &xyDists[1], &xyDists[2], &xyDists[3]);
    min_dist = 9999.0;
    for (i=0; i<4;i++) {
        if (xyDists[i] < min_dist) {
            min_dist = xyDists[i];
        }
    }
    return min_dist;
}



void get_pep_array(bool is_upper, bool pep_upper, float xCOM, float yCOM, float zCOM, ReadPdb * pdb, ReadDcd * dcd) {
    float xCell;
    int atom_num;
    float xVal, yVal, zVal;
    float zValAbs;
    float x1Diff, y1Diff, z1Diff;
    float x2Diff, y2Diff, z2Diff;
    float xDiff, yDiff, zDiff;
    float xyDist;
    float xyDists[4];
    mol_defines molVal;
    int i;
    int histBucket;
    std::vector<int>::iterator first_it, end_it, it;


    if (is_upper) {
        first_it = std::begin(pdb->pep_up_heavy_atoms);
        end_it = std::end(pdb->pep_up_heavy_atoms);
    } else {
        first_it = std::begin(pdb->pep_low_heavy_atoms);
        end_it = std::end(pdb->pep_low_heavy_atoms);
    }

    for (it = first_it; it != end_it; ++it) {
        atom_num=*it;
        xCell = dcd->uCell[0];
        xVal = dcd->xValArray[0][atom_num-1];
        yVal = dcd->yValArray[0][atom_num-1];
        zVal = dcd->zValArray[0][atom_num-1];
        x1Diff = abs(xCOM - xVal);
        y1Diff = abs(yCOM - yVal);
        zValAbs = abs(zVal);

        xyDist = pow((pow(x1Diff,2) + pow(y1Diff,2)),0.5);
        if (is_upper == pep_upper) {
            histBucket = int (xyDist + 0.00001);
            if (histBucket < NUM_BUCKETS) {
                if (is_upper == pep_upper) {
                    pep_dist[histBucket]++;
                    if (zValAbs < FIXED_LIPID_BOUND) {
                        pep_dist_fixed[histBucket]++;
                    }
                    if (zValAbs < ref_bound[histBucket]) {
                        pep_dist_bound[histBucket]++;
                    }
                }
            }
        } else {
            x2Diff = abs(xCell - x1Diff);
            y2Diff = abs(xCell - y1Diff);
            xyDists[0] =  xyDists[1] =  xyDists[2] =  xyDists[3] = MAX_XY_DIST;
            get_xy_dists(x1Diff, x2Diff, y1Diff, y2Diff, MAX_XY_DIST, &xyDists[0], &xyDists[1], &xyDists[2], &xyDists[3]);
            for (i=0; i<4;i++) {
                if (xyDists[i] < MAX_XY_DIST) {
                    histBucket = int (xyDists[i] + 0.00001);
                    if (histBucket < NUM_BUCKETS) {
                        pep_OPP_dist[histBucket]++;
                    }
                    if (zValAbs < FIXED_LIPID_BOUND) {
                        pep_OPP_dist_fixed[histBucket]++;
                    }
                    if (zValAbs < opp_bound[histBucket]) {
                        pep_OPP_dist_bound[histBucket]++;
                    }
                }
            }
        }
    }
}
void get_wat_array(bool pep_upper, float xCOM, float yCOM, float zCOM, ReadPdb * pdb, ReadDcd * dcd, float pep_xy) {

    float xCell;
    int atom_num;
    float xVal, yVal, zVal;
    float x1Diff, y1Diff, z1Diff;
    float x2Diff, y2Diff, z2Diff;
    float xDiff, yDiff, zDiff;
    float xyDist;
    float xyDists[4];
    mol_defines molVal;
    int i;


//    printf ("UCELL %f %f %f\n", dcd->uCell[0], dcd->uCell[1], dcd->uCell[2]);

    std::vector<int>::iterator first_it, end_it, it;
    first_it = std::begin(pdb->water_heavy_atoms);
    end_it = std::end(pdb->water_heavy_atoms);
    molVal = IS_WATER;

    for (it = first_it; it != end_it; ++it) {
        atom_num=*it;
        xCell = dcd->uCell[0];
        xVal = dcd->xValArray[0][atom_num-1];
        yVal = dcd->yValArray[0][atom_num-1];
        zVal = dcd->zValArray[0][atom_num-1];
        x1Diff = abs(xCOM - xVal);
        y1Diff = abs(yCOM - yVal);
        z1Diff = abs(zCOM - zVal);
        x2Diff = abs(xCell - x1Diff);
        y2Diff = abs(xCell - y1Diff);

        xyDists[0] =  xyDists[1] =  xyDists[2] =  xyDists[3] = MAX_XY_DIST;
        get_xy_dists(x1Diff, x2Diff, y1Diff, y2Diff, MAX_XY_DIST, &xyDists[0], &xyDists[1], &xyDists[2], &xyDists[3]);
        for (i=0; i<4;i++) {
             if (xyDists[i] < MAX_XY_DIST) {
                add_2d_hist(molVal, pep_upper, xyDists[i], zVal, pep_xy);
             }
        }
    }
}

//
//  Get R distance for lipids
//

void get_r_array(bool is_upper, bool pep_upper, std::string mol_type, float xCOM, float yCOM, float zCOM, ReadPdb * pdb, ReadDcd * dcd, float pep_xy) {


    float xCell;
    int atom_num;
    float xVal, yVal, zVal;
    float x1Diff, y1Diff, z1Diff;
    float x2Diff, y2Diff, z2Diff;
    float xDiff, yDiff, zDiff;
    float xyDist;
    float xyDists[4];
    mol_defines molVal;
    int i;

//    printf ("UCELL %f %f %f\n", dcd->uCell[0], dcd->uCell[1], dcd->uCell[2]);

    int atom_cnt;
    atom_cnt = 0;

    std::vector<int>::iterator first_it, end_it, it;
    if (is_upper == true) {
        if (mol_type.compare("DMPC") == 0) {
            first_it = std::begin(pdb->dmpc_up_heavy_atoms);
            end_it = std::end(pdb->dmpc_up_heavy_atoms);
            molVal = IS_DMPC;
        } else if (mol_type.compare("DMPG") == 0) {
            first_it = std::begin(pdb->dmpg_up_heavy_atoms);
            end_it = std::end(pdb->dmpg_up_heavy_atoms);
            molVal = IS_DMPG;
        } else if (mol_type.compare("DMPE") == 0) {
            first_it = std::begin(pdb->dmpe_up_heavy_atoms);
            end_it = std::end(pdb->dmpe_up_heavy_atoms);
            molVal = IS_DMPE;
        }
    } else {
        if (mol_type.compare("DMPC") == 0) {
            first_it = std::begin(pdb->dmpc_low_heavy_atoms);
            end_it = std::end(pdb->dmpc_low_heavy_atoms);
            molVal = IS_DMPC;
        } else if (mol_type.compare("DMPG") == 0) {
            first_it = std::begin(pdb->dmpg_low_heavy_atoms);
            end_it = std::end(pdb->dmpg_low_heavy_atoms);
            molVal = IS_DMPG;
        } else if (mol_type.compare("DMPE") == 0) {
            first_it = std::begin(pdb->dmpe_low_heavy_atoms);
            end_it = std::end(pdb->dmpe_low_heavy_atoms);
            molVal = IS_DMPE;
        }
    }
    for (it = first_it; it != end_it; ++it) {
        atom_num=*it;
        xCell = dcd->uCell[0];
        xVal = dcd->xValArray[0][atom_num-1];
        yVal = dcd->yValArray[0][atom_num-1];
        zVal = dcd->zValArray[0][atom_num-1];
        x1Diff = abs(xCOM - xVal);
        y1Diff = abs(yCOM - yVal);
        z1Diff = abs(zCOM - zVal);
        x2Diff = abs(xCell - x1Diff);
        y2Diff = abs(xCell - y1Diff);

        xyDists[0] =  xyDists[1] =  xyDists[2] =  xyDists[3] = MAX_XY_DIST;
        get_xy_dists(x1Diff, x2Diff, y1Diff, y2Diff, MAX_XY_DIST, &xyDists[0], &xyDists[1], &xyDists[2], &xyDists[3]);
        int cnt;
        atom_cnt++;
        cnt = 0;
        for (i=0; i<4;i++) {
             if (xyDists[i] < MAX_XY_DIST) {
                cnt++;
                add_2d_hist(molVal, pep_upper, xyDists[i], zVal, pep_xy);
             }
        }
    }
}
//
//     Now get the distance to each Phos atom from COM.
//
//

void get_p_array(bool is_upper, bool pep_upper, std::string mol_type, float xCOM, float yCOM, float zCOM, unsigned char * contArray, ReadPdb * pdb, ReadDcd * dcd) {


    float xCell;
    int atom_num;
    int mol_num;
    float xVal, yVal, zVal;
    float x1Diff, y1Diff, z1Diff;
    float x2Diff, y2Diff, z2Diff;
    float xDiff, yDiff, zDiff;
    float xyDist;
    float xyDists[4];
    mol_defines molVal;
    int i;
    int histBucket = 0;
    std::vector<int>::iterator first_it, end_it, it;
//    printf("Start get_p_array is_up=%d  pep_up=%d\n", is_upper, pep_upper);

    if (is_upper == true) {
        if (mol_type.compare("DMPC") == 0) {
            first_it = std::begin(pdb->dmpc_up_P_atoms);
            end_it = std::end(pdb->dmpc_up_P_atoms);
            molVal = IS_DMPC;
        } else if (mol_type.compare("DMPG") == 0) {
            first_it = std::begin(pdb->dmpg_up_P_atoms);
            end_it = std::end(pdb->dmpg_up_P_atoms);
            molVal = IS_DMPG;
        } else if (mol_type.compare("DMPE") == 0) {
            first_it = std::begin(pdb->dmpe_up_P_atoms);
            end_it = std::end(pdb->dmpe_up_P_atoms);
            molVal = IS_DMPE;
        }
    } else {
        if (mol_type.compare("DMPC") == 0) {
            first_it = std::begin(pdb->dmpc_low_P_atoms);
            end_it = std::end(pdb->dmpc_low_P_atoms);
            molVal = IS_DMPC;
        } else if (mol_type.compare("DMPG") == 0) {
            first_it = std::begin(pdb->dmpg_low_P_atoms);
            end_it = std::end(pdb->dmpg_low_P_atoms);
            molVal = IS_DMPG;
//            printf("IS DMPG \n");
        } else if (mol_type.compare("DMPE") == 0) {
            first_it = std::begin(pdb->dmpe_low_P_atoms);
            end_it = std::end(pdb->dmpe_low_P_atoms);
            molVal = IS_DMPE;
        }
    }

//    printf("get_p_array\n");
    for (it = first_it; it != end_it; ++it) {
//        printf("READ it\n");
        atom_num=*it;
        xCell = dcd->uCell[0];
        mol_num = pdb->atoms[atom_num-1]->mol_num;
//        if (molVal == IS_DMPG) {
//            printf("DMPG upper=%d,   atom=%d,    mol=%d\n", is_upper, atom_num-1, pdb->atoms[atom_num-1]->mol_num);
//        } else if (molVal == IS_DMPC) {
//            printf("DMPC upper=%d,   atom=%d,    mol=%d\n", is_upper, atom_num-1, pdb->atoms[atom_num-1]->mol_num);
//        }
        xVal = dcd->xValArray[0][atom_num-1];
        yVal = dcd->yValArray[0][atom_num-1];
        zVal = dcd->zValArray[0][atom_num-1];
        x1Diff = abs(xCOM - xVal);
        y1Diff = abs(yCOM - yVal);
        z1Diff = abs(zCOM - zVal);
        x2Diff = abs(xCell - x1Diff);
        y2Diff = abs(xCell - y1Diff);
//        printf("COM %f %f   PHOS %f %f    DISTS %f %f\n", xCOM, yCOM, xVal, yVal, x1Diff, y1Diff);
        xyDists[0] =  xyDists[1] =  xyDists[2] =  xyDists[3] = MAX_XY_DIST;
//        printf("XVAL %f %f YVAL %f %f\n", x1Diff, x2Diff, y1Diff, y2Diff);
        get_xy_dists(x1Diff, x2Diff, y1Diff, y2Diff, MAX_XY_DIST, &xyDists[0], &xyDists[1], &xyDists[2], &xyDists[3]);
        for (i=0; i<4;i++) {
             if (xyDists[i] < MAX_XY_DIST) {
                 if (is_upper == pep_upper) {
                     if (xyDists[i] < NEAR_END) {
                         if (molVal == IS_DMPC) {
                              pc_region_area_den_step[0]++;
                         } else {
                             pg_region_area_den_step[0]++;
                         }
                     } else if (xyDists[i] < PROXIMAL_END) {
                         if (molVal == IS_DMPC) {
                              pc_region_area_den_step[1]++;
                         } else {
                              pg_region_area_den_step[1]++;
                         }

                     } else if (xyDists[i] < FAR_END) {
                         if (molVal == IS_DMPC) {
                              pc_region_area_den_step[2]++;
                         } else {
                              pg_region_area_den_step[2]++;
                         }
                     }
                }

                histBucket = int (xyDists[i] + 0.00001);
                if (histBucket < NUM_BUCKETS) {
                    if (is_upper == pep_upper) {
                        if (molVal == IS_DMPC) {
                            dmpc_p_dist[histBucket]++;
                        } else {
                            dmpg_p_dist[histBucket]++;
                        }
                        if (is_upper) {
                            if (molVal == IS_DMPC) {
                                dmpc_p_up_dist[histBucket]++;
                            } else {
                                dmpg_p_up_dist[histBucket]++;
                            }
                        } else {
                            if (molVal == IS_DMPC) {
                                dmpc_p_lo_dist[histBucket]++;
                            } else {
                                dmpg_p_lo_dist[histBucket]++;
                            }
                        }
 
                    } else {
                        if (molVal == IS_DMPC) {
                            dmpc_OPP_p_dist[histBucket]++;
                        } else {
                            dmpg_OPP_p_dist[histBucket]++;
                        }
                    }
                    if (contArray[mol_num] != 0) {
                        lipid_cont_dist_hist[histBucket]++;
                    }
                }
            }
        }
    }
}

    


int main(int argc, char * argv[]) {

    char dcd_file[MAX_FILENAME_LEN];
    char pdb_file[MAX_FILENAME_LEN];
    char cm_file[MAX_FILENAME_LEN];
    char * outFile;
    int i, j;
    float uCell[3];
    float nearArea, proxArea, farArea;



#ifdef OLD
    std::vector <int> lip_up_heavy_atom_off;
    std::vector <int> lip_low_heavy_atom_off;
    std::vector <int> pep_up_heavy_atom_off;
    std::vector <int> pep_low_heavy_atom_off;
#endif

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

    read_bound();

    snprintf(pdb_file, MAX_FILENAME_LEN, "%s/spt.pdb", argv[1]);


    ReadPdb pdb(pdb_file);
    pdb.find_heavy_atoms(NUM_LIPID_LEAFLET, NUM_PEP_LEAFLET);
    printf("PDB %s\n", pdb_file);

    ReadCM cm(argv[7]);
    printf("CM %s\n", argv[7]);

    int first_step, last_step, step;
    first_step = atoi(argv[5]);
    last_step = atoi(argv[6]);
    std::vector <float> up_r_dist;
    std::vector <float> low_r_dist;
    ReadDcd * dcd;

    for (i=0; i<NUM_REG; i++) {
        dmpc_OPP_den_reg[i] = 0;
        dmpg_OPP_den_reg[i] = 0;
        dmpc_SAM_den_reg[i] = 0;
        dmpg_SAM_den_reg[i] = 0;
        dmpc_SAM_vol_fix_den_reg[i] = 0;
        dmpc_SAM_vol_ref_den_reg[i] = 0;
        dmpc_OPP_vol_fix_den_reg[i] = 0;
        dmpc_OPP_vol_ref_den_reg[i] = 0;
        dmpg_SAM_vol_fix_den_reg[i] = 0;
        dmpg_SAM_vol_ref_den_reg[i] = 0;
        dmpg_OPP_vol_fix_den_reg[i] = 0;
        dmpg_OPP_vol_ref_den_reg[i] = 0;
    }


    for (i=0; i<NUM_BUCKETS; i++) {
        dmpc_dist_hist[i] = 0;
        dmpc_same_dist_hist[i] = 0;
        dmpc_same_dist_hist_fixed[i] = 0;
        dmpc_same_dist_hist_bound[i] = 0;

        dmpc_opp_dist_hist[i] = 0;
        dmpc_opp_dist_hist_fixed[i] = 0;
        dmpc_opp_dist_hist_bound[i] = 0;

        dmpg_dist_hist[i] = 0;
        dmpg_same_dist_hist[i] = 0;
        dmpg_same_dist_hist_fixed[i] = 0;
        dmpg_same_dist_hist_bound[i] = 0;
        dmpg_opp_dist_hist[i] = 0;
        dmpg_opp_dist_hist_fixed[i] = 0;
        dmpg_opp_dist_hist_bound[i] = 0;
        lipid_cont_dist_hist[i] = 0;
        dmpc_OPP_p_dist[i] = 0;
        dmpg_OPP_p_dist[i] = 0;
        dmpc_p_dist[i] = 0;
        dmpg_p_dist[i] = 0;
        dmpc_p_up_dist[i] = 0;
        dmpg_p_up_dist[i] = 0;
        dmpc_p_lo_dist[i] = 0;
        dmpg_p_lo_dist[i] = 0;
        pep_dist[i] = 0;
        pep_dist_fixed[i] = 0;
        pep_dist_bound[i] = 0;
        pep_OPP_dist[i] = 0;
        pep_OPP_dist_fixed[i] = 0;
        pep_OPP_dist_bound[i] = 0;
        for (j=0; j<NUM_Z_BUCKETS; j++) {
            dmpc_2d_dist_hist[i][j] = 0;
            dmpg_2d_dist_hist[i][j] = 0;
            lipid_2d_dist_hist[i][j] = 0;
            water_2d_dist_hist[i][j] = 0;
        }
        for (j=0; j<4; j++) {
            dmpc_SAM_xy_reg[i][j] = 0;
            dmpg_SAM_xy_reg[i][j] = 0;
            dmpc_OPP_xy_reg[i][j] = 0;
            dmpg_OPP_xy_reg[i][j] = 0;
        }
    }
    for (j=0; j<4; j++) {
        cnt_dmpc_SAM_xy[j] = 0;
        cnt_dmpc_OPP_xy[j] = 0;
        cnt_dmpg_SAM_xy[j] = 0;
        cnt_dmpg_OPP_xy[j] = 0;
    }

    printf("Start loop %d\n", step);
    fflush(stdout);



    char step_reg_file[MAX_FILENAME_LEN];
    FILE *fp_reg_step;
    snprintf(step_reg_file, MAX_FILENAME_LEN, "%s/reg_den_step_%s_%s_%s-%s.dat", argv[1], argv[2], argv[3], argv[5], argv[6]);
    fp_reg_step = fopen(step_reg_file,"w");


    nearArea =  pow(NEAR_END,2) * M_PI;
    proxArea =  (pow(PROXIMAL_END,2) - pow(NEAR_END,2)) * M_PI;
    farArea =  (pow(FAR_END,2) - pow(PROXIMAL_END,2)) * M_PI;


    for (step=first_step; step < (last_step+1); step++) {
        snprintf(dcd_file, MAX_FILENAME_LEN, "%s/output/%s%s_%05d_%s.dcd", argv[1], argv[4], argv[2], step, argv[3]);
        if ((step % 100) == 0) {
            printf("STEP %d - DCD %s \n", step, dcd_file);
        }
        dcd = new ReadDcd(dcd_file);


        for (i=0; i<NUM_REG; i++) {
             pc_region_vol_den_step[i] = 0.0;
             pg_region_vol_den_step[i] = 0.0;
             pc_region_area_den_step[i] = 0.0;
             pg_region_area_den_step[i] = 0.0;
        }


        dcd->uCell[0] = dcd->boxArray[index]->ibox[0];
        dcd->uCell[1] = dcd->boxArray[index]->ibox[2];
        dcd->uCell[2] = dcd->boxArray[index]->ibox[5];

        float xUValCOM, yUValCOM, zUValCOM;
        float xLValCOM, yLValCOM, zLValCOM;
        unsigned char * cm_ptr;


        cm_ptr = &(cm.contacts_ptr[step-first_step]->num_lip_cont[0]);
        get_pep_com(true, &xUValCOM, &yUValCOM, &zUValCOM, &pdb, dcd);
        get_pep_com(false, &xLValCOM, &yLValCOM, &zLValCOM, &pdb, dcd);


        pep_xy_dist =   get_peptide_dist(xUValCOM, yUValCOM, zUValCOM, xLValCOM, yLValCOM, zLValCOM, dcd);
//        printf("XY DIST = %f\n", pep_xy_dist);

        add_peptide_dist_hist(xUValCOM, yUValCOM, zUValCOM, xLValCOM, yLValCOM, zLValCOM, dcd);

        get_pep_array(true, true, xUValCOM, yUValCOM, zUValCOM, &pdb, dcd);
        get_pep_array(true, false, xLValCOM, yLValCOM, zLValCOM, &pdb, dcd);
        get_pep_array(false, true, xUValCOM, yUValCOM, zUValCOM, &pdb, dcd);
        get_pep_array(false, false, xLValCOM, yLValCOM, zLValCOM, &pdb, dcd);

        get_wat_array(true, xUValCOM, yUValCOM, zUValCOM, &pdb, dcd, pep_xy_dist);
        get_wat_array(false, xLValCOM, yLValCOM, zLValCOM, &pdb, dcd, pep_xy_dist);

        add_2d_ion_hist(IS_SOD_ION, true, xUValCOM, yUValCOM, zUValCOM, &pdb, dcd);
        add_2d_ion_hist(IS_SOD_ION, false, xLValCOM, yLValCOM, zLValCOM, &pdb, dcd);
        add_2d_ion_hist(IS_CLA_ION, true, xUValCOM, yUValCOM, zUValCOM, &pdb, dcd);
        add_2d_ion_hist(IS_CLA_ION, false, xLValCOM, yLValCOM, zLValCOM, &pdb, dcd);

        num_dmpc_lipids = pdb.dmpc_up_heavy_atoms.size();
        get_r_array(false, true, "DMPC", xUValCOM, yUValCOM, zUValCOM, &pdb, dcd, pep_xy_dist);
        get_r_array(true, true, "DMPC", xUValCOM, yUValCOM, zUValCOM,  &pdb, dcd, pep_xy_dist);
        get_p_array(true, true, "DMPC", xUValCOM, yUValCOM, zUValCOM, cm_ptr, &pdb, dcd);
        get_p_array(false, true, "DMPC", xUValCOM, yUValCOM, zUValCOM, cm_ptr, &pdb, dcd);
        
        num_dmpc_lipids = pdb.dmpc_low_heavy_atoms.size();
        get_r_array(true, false, "DMPC", xLValCOM, yLValCOM, zLValCOM, &pdb, dcd, pep_xy_dist);
        get_r_array(false, false, "DMPC", xLValCOM, yLValCOM, zLValCOM, &pdb, dcd, pep_xy_dist);
        get_p_array(true, false, "DMPC", xLValCOM, yLValCOM, zLValCOM, cm_ptr, &pdb, dcd);
        get_p_array(false, false, "DMPC", xLValCOM, yLValCOM, zLValCOM, cm_ptr, &pdb, dcd);

        num_dmpg_lipids = pdb.dmpg_up_heavy_atoms.size();
        get_r_array(false, true, "DMPG", xUValCOM, yUValCOM, zUValCOM, &pdb, dcd, pep_xy_dist);
        get_r_array(true, true, "DMPG", xUValCOM, yUValCOM, zUValCOM,  &pdb, dcd, pep_xy_dist);
        get_p_array(true, true, "DMPG", xUValCOM, yUValCOM, zUValCOM, cm_ptr, &pdb, dcd);
        get_p_array(false, true, "DMPG", xUValCOM, yUValCOM, zUValCOM, cm_ptr, &pdb, dcd);
        
        num_dmpg_lipids = pdb.dmpg_low_heavy_atoms.size();
        get_r_array(true, false, "DMPG", xLValCOM, yLValCOM, zLValCOM, &pdb, dcd, pep_xy_dist);
        get_r_array(false, false, "DMPG", xLValCOM, yLValCOM, zLValCOM, &pdb, dcd, pep_xy_dist);
        get_p_array(true, false, "DMPG", xLValCOM, yLValCOM, zLValCOM, cm_ptr, &pdb, dcd);
        get_p_array(false, false, "DMPG", xLValCOM, yLValCOM, zLValCOM, cm_ptr, &pdb, dcd);
        delete dcd;


        fprintf(fp_reg_step, "%f %f %f %f %f %f %f %f %f %f %f %f\n",    
             pc_region_vol_den_step[0]/(nearArea*FIXED_LIPID_BOUND*2), pg_region_vol_den_step[0]/(nearArea*FIXED_LIPID_BOUND*2), 
             pc_region_area_den_step[0]/(nearArea*2), pg_region_area_den_step[0]/(nearArea*2),
             pc_region_vol_den_step[1]/(proxArea*FIXED_LIPID_BOUND*2), pg_region_vol_den_step[1]/(proxArea*FIXED_LIPID_BOUND*2), 
             pc_region_area_den_step[1]/(proxArea*2), pg_region_area_den_step[1]/(proxArea*2),
             pc_region_vol_den_step[2]/(farArea*FIXED_LIPID_BOUND*2), pg_region_vol_den_step[2]/(farArea*FIXED_LIPID_BOUND*2), 
             pc_region_area_den_step[2]/(farArea*2), pg_region_area_den_step[2]/(farArea*2));

    }

//
//   Write the SOD density hist file
    FILE *fp_ion_den;
    float curr_dist;
    float currRad;
    float prevRad;
    float myArea;
    int num_steps;

    char ion_den_file[MAX_FILENAME_LEN];
    snprintf(ion_den_file, MAX_FILENAME_LEN, "%s/ion_den_r_%s_%s_%s-%s.dat", argv[1], argv[2], argv[3], argv[5], argv[6]);
    fp_ion_den = fopen(ion_den_file,"w");


    num_steps = (last_step+1) - first_step;
    curr_dist = FIRST_DIST;
    prevRad = 0.0;
    for (i=0; i<NUM_BUCKETS; i++) {
        currRad = (i+1) * BUCKET_DIST;
        myArea =  (pow(currRad,2) - pow(prevRad, 2)) * M_PI;
        prevRad = currRad;
        fprintf(fp_ion_den, "%f %f %f\n", curr_dist,
          (sod_ion_dist_hist[i]/ (myArea * num_steps * 2 * FIXED_LIPID_BOUND)), (cla_ion_dist_hist[i]/ (myArea * num_steps * 2 * FIXED_LIPID_BOUND)));
    }
    fclose(fp_ion_den);


//
//    Write lipid density hist file
    FILE *fp_lip_den;
    char lip_den_file[MAX_FILENAME_LEN];
    snprintf(lip_den_file, MAX_FILENAME_LEN, "%s/lip_den_r_%s_%s_%s-%s.dat", argv[1], argv[2], argv[3], argv[5], argv[6]);
    fp_lip_den = fopen(lip_den_file,"w");
    num_steps = (last_step+1) - first_step;
    curr_dist = FIRST_DIST;
    prevRad = 0.0;
    for (i=0; i<NUM_BUCKETS; i++) {
        currRad = (i+1) * BUCKET_DIST;
        myArea =  (pow(currRad,2) - pow(prevRad, 2)) * M_PI;
        prevRad = currRad;
        fprintf(fp_lip_den, "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", 
          curr_dist, 
//          dmpc_dist_hist[i], dmpc_same_dist_hist[i], dmpc_opp_dist_hist[i], 
//          dmpg_dist_hist[i], dmpg_same_dist_hist[i], dmpg_opp_dist_hist[i],

          dmpc_dist_hist[i]/(myArea * num_steps * 2), 
          dmpc_same_dist_hist[i]/(myArea * num_steps * 2), 
          dmpc_same_dist_hist_fixed[i]/(myArea * num_steps * 2* FIXED_LIPID_BOUND), 
          dmpc_same_dist_hist_bound[i]/(myArea * num_steps * 2 * ref_bound[i]), 

          dmpc_opp_dist_hist[i]/(myArea * num_steps * 2), 
          dmpc_opp_dist_hist_fixed[i]/(myArea * num_steps * 2* FIXED_LIPID_BOUND), 
          dmpc_opp_dist_hist_bound[i]/(myArea * num_steps * 2 * opp_bound[i]), 

          dmpg_dist_hist[i]/(myArea * num_steps * 2),
          dmpg_same_dist_hist[i]/(myArea * num_steps * 2),
          dmpg_same_dist_hist_fixed[i]/(myArea * num_steps * 2* FIXED_LIPID_BOUND), 
          dmpg_same_dist_hist_bound[i]/(myArea * num_steps * 2 * ref_bound[i]), 
          dmpg_opp_dist_hist[i]/(myArea * num_steps * 2),
          dmpg_opp_dist_hist_fixed[i]/(myArea * num_steps * 2* FIXED_LIPID_BOUND), 
          dmpg_opp_dist_hist_bound[i]/(myArea * num_steps * 2 * opp_bound[i]), 
          ((float) dmpg_same_dist_hist[i])/dmpc_same_dist_hist[i],
          ((float) dmpg_opp_dist_hist[i])/dmpc_opp_dist_hist[i]);
        curr_dist += BUCKET_DIST;

    }
    fclose(fp_lip_den);

//
//    Write 2D lipid density hist file
    FILE * fp_out;
    char fileName[MAX_FILENAME_LEN];

//
// Write the water file
    snprintf(fileName, MAX_FILENAME_LEN, "%s/water_den_r_%s_%s_%s-%s.dat", argv[1], argv[2], argv[3], argv[5], argv[6]);
    fp_out = fopen(fileName,"w");
    prevRad = NUM_BUCKETS*BUCKET_DIST;
    for (i=NUM_BUCKETS-1; i>=0; i--) {
        currRad = i * BUCKET_DIST;
        myArea =  (pow(prevRad,2) - pow(currRad, 2)) * M_PI;
        prevRad = currRad;
        for (j=0; j<NUM_Z_BUCKETS; j++) {
            fprintf(fp_out, "%f ", (water_2d_dist_hist[i][j]/myArea)/(num_steps*2));
        }
        fprintf(fp_out, "\n");
    }

    prevRad = 0.0;
    for (i=0; i<NUM_BUCKETS; i++) {
        currRad = (i+1) * BUCKET_DIST;
        myArea =  (pow(currRad,2) - pow(prevRad, 2)) * M_PI;
        prevRad = currRad;
        for (j=0; j<NUM_Z_BUCKETS; j++) {
            fprintf(fp_out, "%f ", (water_2d_dist_hist[i][j]/myArea)/(num_steps*2));
        }
        fprintf(fp_out, "\n");
    }
    fclose(fp_out);


//
// Write the DMPC file
    snprintf(fileName, MAX_FILENAME_LEN, "%s/dmpc_den_r_%s_%s_%s-%s.dat", argv[1], argv[2], argv[3], argv[5], argv[6]);
    fp_out = fopen(fileName,"w");
    prevRad = NUM_BUCKETS*BUCKET_DIST;
    for (i=NUM_BUCKETS-1; i>=0; i--) {
        currRad = i * BUCKET_DIST;
        myArea =  (pow(prevRad,2) - pow(currRad, 2)) * M_PI;
        prevRad = currRad;
        for (j=0; j<NUM_Z_BUCKETS; j++) {
            fprintf(fp_out, "%f ", (dmpc_2d_dist_hist[i][j]/myArea)/(num_steps*2));
        }
        fprintf(fp_out, "\n");
    }

    prevRad = 0.0;
    for (i=0; i<NUM_BUCKETS; i++) {
        currRad = (i+1) * BUCKET_DIST;
        myArea =  (pow(currRad,2) - pow(prevRad, 2)) * M_PI;
        prevRad = currRad;
        for (j=0; j<NUM_Z_BUCKETS; j++) {
            fprintf(fp_out, "%f ", (dmpc_2d_dist_hist[i][j]/myArea)/(num_steps*2));
        }
        fprintf(fp_out, "\n");
    }
    fclose(fp_out);
//
// Write DMPG file.
//
    snprintf(fileName, MAX_FILENAME_LEN, "%s/dmpg_den_r_%s_%s_%s-%s.dat", argv[1], argv[2], argv[3], argv[5], argv[6]);
    fp_out = fopen(fileName,"w");
    prevRad = NUM_BUCKETS*BUCKET_DIST;
    for (i=NUM_BUCKETS-1; i>=0; i--) {
        currRad = i * BUCKET_DIST;
        myArea =  (pow(prevRad,2) - pow(currRad, 2)) * M_PI;
        prevRad = currRad;
        for (j=0; j<NUM_Z_BUCKETS; j++) {
            fprintf(fp_out, "%f ", (dmpg_2d_dist_hist[i][j]/myArea)/(num_steps*2));
        }
        fprintf(fp_out, "\n");
    }
    prevRad = 0.0;
    for (i=0; i<NUM_BUCKETS; i++) {
        currRad = (i+1) * BUCKET_DIST;
        myArea =  (pow(currRad,2) - pow(prevRad, 2)) * M_PI;
        prevRad = currRad;
        for (j=0; j<NUM_Z_BUCKETS; j++) {
            fprintf(fp_out, "%f ", (dmpg_2d_dist_hist[i][j]/myArea)/(num_steps*2));
        }
        fprintf(fp_out, "\n");
    }
    fclose(fp_out);
//
//   Write All lipids file.
//
    snprintf(fileName, MAX_FILENAME_LEN, "%s/all_den_r_%s_%s_%s-%s.dat", argv[1], argv[2], argv[3], argv[5], argv[6]);
    fp_out = fopen(fileName,"w");

    prevRad = NUM_BUCKETS*BUCKET_DIST;
    for (i=NUM_BUCKETS-1; i>=0; i--) {
        currRad = i * BUCKET_DIST;
        myArea =  (pow(prevRad,2) - pow(currRad, 2)) * M_PI;
        prevRad = currRad;
        for (j=0; j<NUM_Z_BUCKETS; j++) {
            fprintf(fp_out, "%f ", (lipid_2d_dist_hist[i][j]/myArea)/(num_steps*2));
        }
        fprintf(fp_out, "\n");
    }
    prevRad = 0.0;
    for (i=0; i<NUM_BUCKETS; i++) {
        currRad = (i+1) * BUCKET_DIST;
        myArea =  (pow(currRad,2) - pow(prevRad, 2)) * M_PI;
        prevRad = currRad;
        for (j=0; j<NUM_Z_BUCKETS; j++) {
            fprintf(fp_out, "%f ", (lipid_2d_dist_hist[i][j]/myArea)/(num_steps*2));
        }
        fprintf(fp_out, "\n");
    }
    fclose(fp_out);
//
//  Get the lipid boundry

    snprintf(fileName, MAX_FILENAME_LEN, "%s/lipid_bound_r_%s_%s_%s-%s.dat", argv[1], argv[2], argv[3], argv[5], argv[6]);
    float bucket1, bucket2;
    float bucket1_den, bucket2_den;
    float boundry_den;
    int den_cnt = 0;
    boundry_den = 0.0;
    
    fp_out = fopen(fileName,"w");
    for (i=0; i<NUM_BUCKETS; i++) {
        bucket1 = 0;
        bucket2 = 0;

        for (j=0; j<NUM_Z_BUCKETS; j++) {
            if (lipid_2d_dist_hist[i][j] > water_2d_dist_hist[i][j]) {
                bucket1 = (float) j;
                myArea =  (pow(i+1,2) - pow(i, 2)) * M_PI;
                bucket1_den = (lipid_2d_dist_hist[i][j]/myArea) / (num_steps * 2);
                break;
            }
        }

        for (j=NUM_Z_BUCKETS-1; j>=0; j--) {
            if (lipid_2d_dist_hist[i][j] > water_2d_dist_hist[i][j]) {
                bucket2 = (float) j;
                bucket2_den = (lipid_2d_dist_hist[i][j]/myArea) / (num_steps * 2)  ;
                break;
            }
        }

        fprintf(fp_out, "%4.1f %4.1f %4.1f %4.4f %4.4f\n", (((float) i)*BUCKET_DIST) + FIRST_DIST,
          MIN_Z_DIST + (((float)bucket1) * BUCKET_Z_DIST),
          MIN_Z_DIST + (((float)bucket2) * BUCKET_Z_DIST), bucket1_den, bucket2_den);


        if (i >= (NUM_BUCKETS-8)) {
            boundry_den += (bucket1_den + bucket2_den);
            den_cnt += 2;
        }
    }
    boundry_den = boundry_den / den_cnt;
    fclose(fp_out);

    snprintf(fileName, MAX_FILENAME_LEN, "%s/bound_den_%s_%s_%s-%s.dat", argv[1], argv[2], argv[3], argv[5], argv[6]);
    fp_out = fopen(fileName,"w");
    fprintf(fp_out, "%4.4f Density at boundry\n", boundry_den);
    fclose(fp_out);


    snprintf(fileName, MAX_FILENAME_LEN, "%s/lipid_den_bound_r_%s_%s_%s-%s.dat", argv[1], argv[2], argv[3], argv[5], argv[6]);
    
    fp_out = fopen(fileName,"w");
    for (i=0; i<NUM_BUCKETS; i++) {
        bucket1 = 0;
        bucket2 = 0;
        myArea =  (pow(i+1,2) - pow(i, 2)) * M_PI;

        for (j=0; j<NUM_Z_BUCKETS; j++) {
            bucket1_den = (lipid_2d_dist_hist[i][j]/myArea) / (num_steps * 2);
            if (bucket1_den > boundry_den) {
                bucket1 = (float) j;
                break;
            }
        }
        for (j=NUM_Z_BUCKETS-1; j>=0; j--) {
            bucket2_den = (lipid_2d_dist_hist[i][j]/myArea) / (num_steps * 2);
            if (bucket2_den > boundry_den) {
                bucket2 = (float) j;
                break;
            }
        }
        fprintf(fp_out, "%4.1f %4.1f %4.1f\n", (((float) i)*BUCKET_DIST) + FIRST_DIST,
          MIN_Z_DIST + (((float)bucket1) * BUCKET_Z_DIST), MIN_Z_DIST + (((float)bucket2) * BUCKET_Z_DIST));
    }
    fclose(fp_out);


//
//    Write lipid density contact file
//
    FILE *fp_lip_cont;
    char lip_cont_file[MAX_FILENAME_LEN];
    snprintf(lip_cont_file, MAX_FILENAME_LEN, "%s/lip_cont_r_%s_%s_%s-%s.dat", argv[1], argv[2], argv[3], argv[5], argv[6]);
    fp_lip_cont = fopen(lip_cont_file,"w");
    curr_dist = FIRST_DIST;
    prevRad = 0.0;
    for (i=0; i<NUM_BUCKETS; i++) {
        currRad = (i+1) * BUCKET_DIST;
        myArea =  (pow(currRad,2) - pow(prevRad, 2)) * M_PI;
        prevRad = currRad;
        fprintf(fp_lip_cont, "%f %d %f\n", curr_dist, lipid_cont_dist_hist[i], lipid_cont_dist_hist[i]/(myArea*num_steps*2));
        curr_dist += BUCKET_DIST;
    }
    fclose(fp_lip_cont);

//
//    Write P density file
//
    char p_file[MAX_FILENAME_LEN];
    snprintf(p_file, MAX_FILENAME_LEN, "%s/Phos_r_%s_%s_%s-%s.dat", argv[1], argv[2], argv[3], argv[5], argv[6]);
    fp_out = fopen(p_file,"w");
    curr_dist = FIRST_DIST;
    prevRad = 0.0;
    for (i=0; i<NUM_BUCKETS; i++) {
        currRad = (i+1) * BUCKET_DIST;
        myArea =  (pow(currRad,2) - pow(prevRad, 2)) * M_PI;
        prevRad = currRad;

        fprintf(fp_out, "%f %d %f %d %f %d %f %d %f\n", curr_dist, 
          dmpg_p_dist[i], dmpg_p_dist[i]/(myArea*num_steps*2),
          dmpg_OPP_p_dist[i], dmpg_OPP_p_dist[i]/(myArea*num_steps*2),
          dmpc_p_dist[i], dmpc_p_dist[i]/(myArea*num_steps*2),
          dmpc_OPP_p_dist[i], dmpc_OPP_p_dist[i]/(myArea*num_steps*2));
        curr_dist += BUCKET_DIST;
    }
    fclose(fp_out);


//
//    Write Region P density file
//
    float GNearVal, GProxVal, GDistVal, GFarVal;
    float GNearValUP, GProxValUP, GDistValUP, GFarValUP;
    float GNearValLO, GProxValLO, GDistValLO, GFarValLO;

    float CNearVal, CProxVal, CDistVal, CFarVal;
    float CNearValUP, CProxValUP, CDistValUP, CFarValUP;
    float CNearValLO, CProxValLO, CDistValLO, CFarValLO;
    float totNearVal, totProxVal, totDistVal, totFarVal;
    float totNearValUP, totProxValUP, totDistValUP, totFarValUP;
    float totNearValLO, totProxValLO, totDistValLO, totFarValLO;
    float gVal, cVal, totVal;
    float gValUP, cValUP, totValUP;
    float gValLO, cValLO, totValLO;
    snprintf(p_file, MAX_FILENAME_LEN, "%s/Phos_region_%s_%s_%s-%s.dat", argv[1], argv[2], argv[3], argv[5], argv[6]);
    fp_out = fopen(p_file,"w");
    gVal = 0.0;
    cVal = 0.0;
    totVal = 0.0;
    gValUP = 0.0;
    cValUP = 0.0;
    totValUP = 0.0;
    gValLO = 0.0;
    cValLO = 0.0;
    totValLO = 0.0;
    myArea =  pow(NEAR_END,2) * M_PI;
    for (i=0; i<NEAR_END; i++) {
        gVal += dmpg_p_dist[i];
        cVal += dmpc_p_dist[i];
        totVal += dmpc_p_dist[i] + dmpg_p_dist[i];
        gValUP += dmpg_p_up_dist[i];
        cValUP += dmpc_p_up_dist[i];
        totValUP += dmpc_p_up_dist[i] + dmpg_p_up_dist[i];
        gValLO += dmpg_p_lo_dist[i];
        cValLO += dmpc_p_lo_dist[i];
        totValLO += dmpc_p_lo_dist[i] + dmpg_p_lo_dist[i];
    }
    GNearVal = (gVal / myArea)/ (num_steps * 2);
    CNearVal = (cVal / myArea)/ (num_steps * 2);
    totNearVal = (totVal / myArea)/ (num_steps * 2);
    GNearValUP = (gValUP / myArea)/ (num_steps);
    CNearValUP = (cValUP / myArea)/ (num_steps);
    totNearValUP = (totValUP / myArea)/ (num_steps);
    GNearValLO = (gValLO / myArea)/ (num_steps);
    CNearValLO = (cValLO / myArea)/ (num_steps);
    totNearValLO = (totValLO / myArea)/ (num_steps);

    gVal = 0.0;
    cVal = 0.0;
    totVal = 0.0;
    gValUP = 0.0;
    cValUP = 0.0;
    totValUP = 0.0;
    gValLO = 0.0;
    cValLO = 0.0;
    totValLO = 0.0;
    myArea =  (pow(PROXIMAL_END,2) - pow(NEAR_END,2)) * M_PI;
    for (i=NEAR_END; i<PROXIMAL_END; i++) {
        gVal += dmpg_p_dist[i];
        cVal += dmpc_p_dist[i];
        totVal += dmpc_p_dist[i] + dmpg_p_dist[i];
        gValUP += dmpg_p_up_dist[i];
        cValUP += dmpc_p_up_dist[i];
        totValUP += dmpc_p_up_dist[i] + dmpg_p_up_dist[i];
        gValLO += dmpg_p_lo_dist[i];
        cValLO += dmpc_p_lo_dist[i];
        totValLO += dmpc_p_lo_dist[i] + dmpg_p_lo_dist[i];
    }
    GProxVal = (gVal / myArea)/ (num_steps * 2);
    CProxVal = (cVal / myArea)/ (num_steps * 2);
    totProxVal = (totVal / myArea)/ (num_steps * 2);
    GProxValUP = (gValUP / myArea)/ (num_steps);
    CProxValUP = (cValUP / myArea)/ (num_steps);
    totProxValUP = (totValUP / myArea)/ (num_steps);
    GProxValLO = (gValLO / myArea)/ (num_steps);
    CProxValLO = (cValLO / myArea)/ (num_steps);
    totProxValLO = (totValLO / myArea)/ (num_steps);

//    gVal = 0.0;
//    cVal = 0.0;
//    totVal = 0.0;
//    myArea =  (pow(DIST_END,2) - pow(PROXIMAL_END,2)) * M_PI;
//    for (i=PROXIMAL_END; i<DIST_END; i++) {
//        gVal += dmpg_p_dist[i];
//        cVal += dmpc_p_dist[i];
//        totVal += dmpc_p_dist[i] + dmpg_p_dist[i];
//    }
//    GDistVal = (gVal / myArea)/ (num_steps * 2);
//    CDistVal = (cVal / myArea)/ (num_steps * 2);
//    totDistVal = (totVal / myArea)/ (num_steps * 2);

    gVal = 0.0;
    cVal = 0.0;
    totVal = 0.0;
    gValUP = 0.0;
    cValUP = 0.0;
    totValUP = 0.0;
    gValLO = 0.0;
    cValLO = 0.0;
    totValLO = 0.0;
    myArea =  (pow(FAR_END,2) - pow(PROXIMAL_END,2)) * M_PI;
    for (i=DIST_END; i<FAR_END; i++) {
        gVal += dmpg_p_dist[i];
        cVal += dmpc_p_dist[i];
        totVal += dmpc_p_dist[i] + dmpg_p_dist[i];
        gValUP += dmpg_p_up_dist[i];
        cValUP += dmpc_p_up_dist[i];
        totValUP += dmpc_p_up_dist[i] + dmpg_p_up_dist[i];
        gValLO += dmpg_p_lo_dist[i];
        cValLO += dmpc_p_lo_dist[i];
        totValLO += dmpc_p_lo_dist[i] + dmpg_p_lo_dist[i];
    }
    GFarVal = (gVal / myArea)/ (num_steps * 2);
    CFarVal = (cVal / myArea)/ (num_steps * 2);
    totFarVal = (totVal / myArea)/ (num_steps * 2);
    GFarValUP = (gValUP / myArea)/ (num_steps);
    CFarValUP = (cValUP / myArea)/ (num_steps);
    totFarValUP = (totValUP / myArea)/ (num_steps);
    GFarValLO = (gValLO / myArea)/ (num_steps);
    CFarValLO = (cValLO / myArea)/ (num_steps);
    totFarValLO = (totValLO / myArea)/ (num_steps);

    fprintf(fp_out, "DMPG Near %4.5f Prox %4.5f, Dist %4.5f, Far %4.5f\n", GNearVal, GProxVal, GDistVal, GFarVal);
    fprintf(fp_out, "DMPC Near %4.5f Prox %4.5f, Dist %4.5f, Far %4.5f\n", CNearVal, CProxVal, CDistVal, CFarVal);
    fprintf(fp_out, "Total Near %4.5f Prox %4.5f, Dist %4.5f, Far %4.5f\n", totNearVal, totProxVal, totDistVal, totFarVal);

    fprintf(fp_out, "UP DMPG Near %4.5f Prox %4.5f, Dist %4.5f, Far %4.5f\n", GNearValUP, GProxValUP, GDistValUP, GFarValUP);
    fprintf(fp_out, "UP DMPC Near %4.5f Prox %4.5f, Dist %4.5f, Far %4.5f\n", CNearValUP, CProxValUP, CDistValUP, CFarValUP);
    fprintf(fp_out, "UP Total Near %4.5f Prox %4.5f, Dist %4.5f, Far %4.5f\n", totNearValUP, totProxValUP, totDistValUP, totFarValUP);

    fprintf(fp_out, "LO DMPG Near %4.5f Prox %4.5f, Dist %4.5f, Far %4.5f\n", GNearValLO, GProxValLO, GDistValLO, GFarValLO);
    fprintf(fp_out, "LO DMPC Near %4.5f Prox %4.5f, Dist %4.5f, Far %4.5f\n", CNearValLO, CProxValLO, CDistValLO, CFarValLO);
    fprintf(fp_out, "LO Total Near %4.5f Prox %4.5f, Dist %4.5f, Far %4.5f\n", totNearValLO, totProxValLO, totDistValLO, totFarValLO);


    fclose(fp_out);

//
//   Write PEP density file
//
    char pep_file[MAX_FILENAME_LEN];
    snprintf(pep_file, MAX_FILENAME_LEN, "%s/Pep_r_%s_%s_%s-%s.dat", argv[1], argv[2], argv[3], argv[5], argv[6]);
    fp_out = fopen(pep_file,"w");
    curr_dist = FIRST_DIST;
    prevRad = 0.0;
    for (i=0; i<NUM_BUCKETS; i++) {
        currRad = (i+1) * BUCKET_DIST;
        myArea =  (pow(currRad,2) - pow(prevRad, 2)) * M_PI;
        prevRad = currRad;

        fprintf(fp_out, "%f %d %f %f %f %d %f %f %f \n", curr_dist,
          pep_dist[i], pep_dist[i]/(myArea*num_steps*2),
          pep_dist_fixed[i]/(myArea*num_steps*2*FIXED_LIPID_BOUND),
          pep_dist_bound[i]/(myArea*num_steps*2*ref_bound[i]),
          pep_OPP_dist[i], pep_OPP_dist[i]/(myArea*num_steps*2),
          pep_OPP_dist_fixed[i]/(myArea*num_steps*2*FIXED_LIPID_BOUND),
          pep_OPP_dist_bound[i]/(myArea*num_steps*2*opp_bound[i]));
        curr_dist += BUCKET_DIST;
    }
    fclose(fp_out);
//
//  Write region density file
    char reg_file[MAX_FILENAME_LEN];
    snprintf(reg_file, MAX_FILENAME_LEN, "%s/lip_reg_%s_%s_%s-%s.dat", argv[1], argv[2], argv[3], argv[5], argv[6]);
    fp_out = fopen(reg_file,"w");

// Near region

    prevRad = NEAR_DIST;
    currRad = PROXIMAL_DIST;
    myArea =  pow(currRad,2) * M_PI;
    fprintf(fp_out, "Near %d %f %d %f %d %f %d %f %d %f %d %f %d %f %d %f %d %f %d %f %d %f %d %f\n", 
      dmpc_SAM_den_reg[0], dmpc_SAM_den_reg[0] / (myArea * num_steps*2),
      dmpg_SAM_den_reg[0], dmpg_SAM_den_reg[0] / (myArea * num_steps*2),
      dmpc_OPP_den_reg[0], dmpc_OPP_den_reg[0] / (myArea * num_steps*2),
      dmpg_OPP_den_reg[0], dmpg_OPP_den_reg[0] / (myArea * num_steps*2),
      dmpc_SAM_vol_fix_den_reg[0], dmpc_SAM_vol_fix_den_reg[0] / (myArea * num_steps*2*FIXED_LIPID_BOUND),
      dmpg_SAM_vol_fix_den_reg[0], dmpg_SAM_vol_fix_den_reg[0] / (myArea * num_steps*2*FIXED_LIPID_BOUND),
      dmpc_OPP_vol_fix_den_reg[0], dmpc_OPP_vol_fix_den_reg[0] / (myArea * num_steps*2*FIXED_LIPID_BOUND),
      dmpg_OPP_vol_fix_den_reg[0], dmpg_OPP_vol_fix_den_reg[0] / (myArea * num_steps*2*FIXED_LIPID_BOUND),

      dmpc_SAM_vol_fix_den_reg[0], dmpc_SAM_vol_fix_den_reg[0] / (near_ref_vol * num_steps*2),
      dmpg_SAM_vol_fix_den_reg[0], dmpg_SAM_vol_fix_den_reg[0] / (near_ref_vol * num_steps*2),
      dmpc_OPP_vol_fix_den_reg[0], dmpc_OPP_vol_fix_den_reg[0] / (near_ref_opp_vol * num_steps*2),
      dmpg_OPP_vol_fix_den_reg[0], dmpg_OPP_vol_fix_den_reg[0] / (near_ref_opp_vol * num_steps*2));


    currRad = FAR_DIST;
    prevRad = PROXIMAL_DIST;
    myArea =  (pow(currRad,2) - pow(prevRad, 2)) * M_PI;
    fprintf(fp_out, "Proximal %d %f %d %f %d %f %d %f %d %f %d %f %d %f %d %f %d %f %d %f %d %f %d %f\n", 
      dmpc_SAM_den_reg[1], dmpc_SAM_den_reg[1] / (myArea * num_steps*2),
      dmpg_SAM_den_reg[1], dmpg_SAM_den_reg[1] / (myArea * num_steps*2),
      dmpc_OPP_den_reg[1], dmpc_OPP_den_reg[1] / (myArea * num_steps*2),
      dmpg_OPP_den_reg[1], dmpg_OPP_den_reg[1] / (myArea * num_steps*2),
      dmpc_SAM_vol_fix_den_reg[1], dmpc_SAM_vol_fix_den_reg[1] / (myArea * num_steps*2*FIXED_LIPID_BOUND),
      dmpg_SAM_vol_fix_den_reg[1], dmpg_SAM_vol_fix_den_reg[1] / (myArea * num_steps*2*FIXED_LIPID_BOUND),
      dmpc_OPP_vol_fix_den_reg[1], dmpc_OPP_vol_fix_den_reg[1] / (myArea * num_steps*2*FIXED_LIPID_BOUND),
      dmpg_OPP_vol_fix_den_reg[1], dmpg_OPP_vol_fix_den_reg[1] / (myArea * num_steps*2*FIXED_LIPID_BOUND),

      dmpc_SAM_vol_fix_den_reg[1], dmpc_SAM_vol_fix_den_reg[1] / (prox_ref_vol * num_steps*2),
      dmpg_SAM_vol_fix_den_reg[1], dmpg_SAM_vol_fix_den_reg[1] / (prox_ref_vol * num_steps*2),
      dmpc_OPP_vol_fix_den_reg[1], dmpc_OPP_vol_fix_den_reg[1] / (prox_ref_opp_vol * num_steps*2),
      dmpg_OPP_vol_fix_den_reg[1], dmpg_OPP_vol_fix_den_reg[1] / (prox_ref_opp_vol * num_steps*2));


//    currRad = FAR_DIST;
//    prevRad = DISTANT_DIST;
//    myArea =  (pow(currRad,2) - pow(prevRad, 2)) * M_PI;
//    fprintf(fp_out, "Distant %d %f %d %f %d %f %d %f %d %f %d %f %d %f %d %f\n", 
//      dmpc_SAM_den_reg[2], dmpc_SAM_den_reg[2] / (myArea * num_steps*2),
//      dmpg_SAM_den_reg[2], dmpg_SAM_den_reg[2] / (myArea * num_steps*2),
//      dmpc_OPP_den_reg[2], dmpc_OPP_den_reg[2] / (myArea * num_steps*2),
 //     dmpg_OPP_den_reg[2], dmpg_OPP_den_reg[2] / (myArea * num_steps*2));

    prevRad = FAR_DIST;
    currRad = MAX_XY_DIST;
    myArea =  (pow(currRad,2) - pow(prevRad, 2)) * M_PI;
    fprintf(fp_out, "Far %d %f %d %f %d %f %d %f %d %f %d %f %d %f %d %f %d %f %d %f %d %f %d %f\n", 
      dmpc_SAM_den_reg[3], dmpc_SAM_den_reg[3] / (myArea * num_steps*2),
      dmpg_SAM_den_reg[3], dmpg_SAM_den_reg[3] / (myArea * num_steps*2),
      dmpc_OPP_den_reg[3], dmpc_OPP_den_reg[3] / (myArea * num_steps*2),
      dmpg_OPP_den_reg[3], dmpg_OPP_den_reg[3] / (myArea * num_steps*2),
      dmpc_SAM_vol_fix_den_reg[3], dmpc_SAM_vol_fix_den_reg[3] / (myArea * num_steps*2*FIXED_LIPID_BOUND),
      dmpg_SAM_vol_fix_den_reg[3], dmpg_SAM_vol_fix_den_reg[3] / (myArea * num_steps*2*FIXED_LIPID_BOUND),
      dmpc_OPP_vol_fix_den_reg[3], dmpc_OPP_vol_fix_den_reg[3] / (myArea * num_steps*2*FIXED_LIPID_BOUND),
      dmpg_OPP_vol_fix_den_reg[3], dmpg_OPP_vol_fix_den_reg[3] / (myArea * num_steps*2*FIXED_LIPID_BOUND),
      dmpc_SAM_vol_fix_den_reg[3], dmpc_SAM_vol_fix_den_reg[3] / (far_ref_vol * num_steps*2),
      dmpg_SAM_vol_fix_den_reg[3], dmpg_SAM_vol_fix_den_reg[3] / (far_ref_vol * num_steps*2),
      dmpc_OPP_vol_fix_den_reg[3], dmpc_OPP_vol_fix_den_reg[3] / (far_ref_opp_vol * num_steps*2),
      dmpg_OPP_vol_fix_den_reg[3], dmpg_OPP_vol_fix_den_reg[3] / (far_ref_opp_vol * num_steps*2));

    fclose(fp_out);
//
// Write the histogram of distances between XY of COM of peptides.

    char pep_hist_file[MAX_FILENAME_LEN];
    snprintf(pep_hist_file, MAX_FILENAME_LEN, "%s/Pep_dist_hist_%s_%s_%s-%s.dat", argv[1], argv[2], argv[3], argv[5], argv[6]);
    fp_out = fopen(pep_hist_file,"w");

    curr_dist = FIRST_DIST;
    prevRad = 0.0;
    for (i=0; i<NUM_BUCKETS; i++) {
        currRad = (i+1) * BUCKET_DIST;
        myArea =  (pow(currRad,2) - pow(prevRad, 2)) * M_PI;
        prevRad = currRad;

        fprintf(fp_out, "%f %d %f %f\n", curr_dist,
          peptide_dist_hist[i], peptide_dist_hist[i]/(myArea*num_steps*2),
          peptide_dist_hist[i]/(1.0*num_steps*2));
        curr_dist += BUCKET_DIST;
    }
    fclose(fp_out);

//
// Write the 2D histogram of distances between XY of COM of peptides, and the Z difference between COM.

    char pep_2D_hist_file[MAX_FILENAME_LEN];
    snprintf(pep_2D_hist_file, MAX_FILENAME_LEN, "%s/Pep_2d_dist_hist_%s_%s_%s-%s.dat", argv[1], argv[2], argv[3], argv[5], argv[6]);
    fp_out = fopen(pep_2D_hist_file,"w");

    curr_dist = FIRST_DIST;
    prevRad = 0.0;
    for (i=0; i<NUM_BUCKETS; i++) {
        currRad = (i+1) * BUCKET_DIST;
        myArea =  (pow(currRad,2) - pow(prevRad, 2)) * M_PI;
        prevRad = currRad;
        for (j=0; j<NUM_Z_BUCKETS; j++) {
            fprintf(fp_out, "%f ", peptide_2d_dist_hist[i][j]/(myArea*num_steps*2));
        }
        fprintf(fp_out, "\n");
    }
    fclose(fp_out);

//
//  Val of DEN for each PEP xy Dist PC SAME side
    char xy_dist_file[MAX_FILENAME_LEN];
    snprintf(xy_dist_file, MAX_FILENAME_LEN, "%s/PC_SAM_pepxy_hist_%s_%s_%s-%s.dat", argv[1], argv[2], argv[3], argv[5], argv[6]);
    fp_out = fopen(xy_dist_file,"w");

    for (i=0; i<3; i++) {
        curr_dist = 1;
        prevRad = 0.0;
        for (j=0; j<NUM_BUCKETS; j++) {
            currRad = (j+1) * BUCKET_DIST;
            myArea =  (pow(currRad,2) - pow(prevRad, 2)) * M_PI;
            prevRad = currRad;
            fprintf(fp_out, "%f ", dmpc_SAM_xy_reg[j][i]/((float)cnt_dmpc_SAM_xy[i]));
        }
        fprintf(fp_out, "\n");
    }
    fclose(fp_out);
//
//  Val DEN for each PEP xy Dist PC OPP side
    snprintf(xy_dist_file, MAX_FILENAME_LEN, "%s/PC_OPP_pepxy_hist_%s_%s_%s-%s.dat", argv[1], argv[2], argv[3], argv[5], argv[6]);
    fp_out = fopen(xy_dist_file,"w");

    for (i=0; i<3; i++) {
        curr_dist = FIRST_DIST;
        prevRad = 0.0;
        for (j=0; j<NUM_BUCKETS; j++) {
            currRad = (j+1) * BUCKET_DIST;
            myArea =  (pow(currRad,2) - pow(prevRad, 2)) * M_PI;
            prevRad = currRad;
            fprintf(fp_out, "%f ", dmpc_OPP_xy_reg[j][i]/((float)cnt_dmpc_OPP_xy[i]));
        }
        fprintf(fp_out, "\n");
    }
    fclose(fp_out);
//
//
//  Val of DEN for each PEP xy Dist PG SAME side
    snprintf(xy_dist_file, MAX_FILENAME_LEN, "%s/PG_SAM_pepxy_hist_%s_%s_%s-%s.dat", argv[1], argv[2], argv[3], argv[5], argv[6]);
    fp_out = fopen(xy_dist_file,"w");

    for (i=0; i<3; i++) {
        curr_dist = FIRST_DIST;
        prevRad = 0.0;
        for (j=0; j<NUM_BUCKETS; j++) {
            currRad = (j+1) * BUCKET_DIST;
            myArea =  (pow(currRad,2) - pow(prevRad, 2)) * M_PI;
            prevRad = currRad;
            fprintf(fp_out, "%f ", dmpg_SAM_xy_reg[j][i]/((float)cnt_dmpg_SAM_xy[i]));
        }
        fprintf(fp_out, "\n");
    }
    fclose(fp_out);
//
//  Val DEN for each PEP xy Dist PG OPP side
    snprintf(xy_dist_file, MAX_FILENAME_LEN, "%s/PG_OPP_pepxy_hist_%s_%s_%s-%s.dat", argv[1], argv[2], argv[3], argv[5], argv[6]);
    fp_out = fopen(xy_dist_file,"w");

    for (i=0; i<3; i++) {
        curr_dist = FIRST_DIST;
        prevRad = 0.0;
        for (j=0; j<NUM_BUCKETS; j++) {
            currRad = (j+1) * BUCKET_DIST;
            myArea =  (pow(currRad,2) - pow(prevRad, 2)) * M_PI;
            prevRad = currRad;
            fprintf(fp_out, "%f ", dmpg_OPP_xy_reg[j][i]/((float)cnt_dmpg_OPP_xy[i]));
        }
        fprintf(fp_out, "\n");
    }
    fclose(fp_out);

//
//  Val of DEN for each PEP xy Dist PC SAME side
    snprintf(xy_dist_file, MAX_FILENAME_LEN, "%s/PC_SAM_pepxy_hist_den_%s_%s_%s-%s.dat", argv[1], argv[2], argv[3], argv[5], argv[6]);
    fp_out = fopen(xy_dist_file,"w");

    for (i=0; i<3; i++) {
        curr_dist = 1;
        prevRad = 0.0;
        for (j=0; j<NUM_BUCKETS; j++) {
            currRad = (j+1) * BUCKET_DIST * 10;
            myArea =  (pow(currRad,2) - pow(prevRad, 2)) * M_PI;
            prevRad = currRad;
            fprintf(fp_out, "%f ", dmpc_SAM_xy_reg[j][i]/(myArea*num_steps*2));
        }
        fprintf(fp_out, "\n");
    }
    fclose(fp_out);
//
//  Val DEN for each PEP xy Dist PC OPP side
    snprintf(xy_dist_file, MAX_FILENAME_LEN, "%s/PC_OPP_pepxy_hist_den_%s_%s_%s-%s.dat", argv[1], argv[2], argv[3], argv[5], argv[6]);
    fp_out = fopen(xy_dist_file,"w");

    for (i=0; i<3; i++) {
        curr_dist = FIRST_DIST;
        prevRad = 0.0;
        for (j=0; j<NUM_BUCKETS; j++) {
            currRad = (j+1) * BUCKET_DIST * 10;
            myArea =  (pow(currRad,2) - pow(prevRad, 2)) * M_PI;
            prevRad = currRad;
            fprintf(fp_out, "%f ", dmpc_OPP_xy_reg[j][i]/(myArea*num_steps*2));
        }
        fprintf(fp_out, "\n");
    }
    fclose(fp_out);
//
//
//  Val of DEN for each PEP xy Dist PG SAME side
    snprintf(xy_dist_file, MAX_FILENAME_LEN, "%s/PG_SAM_pepxy_hist_den_%s_%s_%s-%s.dat", argv[1], argv[2], argv[3], argv[5], argv[6]);
    fp_out = fopen(xy_dist_file,"w");

    for (i=0; i<3; i++) {
        curr_dist = FIRST_DIST;
        prevRad = 0.0;
        for (j=0; j<NUM_BUCKETS; j++) {
            currRad = (j+1) * BUCKET_DIST * 10;
            myArea =  (pow(currRad,2) - pow(prevRad, 2)) * M_PI;
            prevRad = currRad;
            fprintf(fp_out, "%f ", dmpg_SAM_xy_reg[j][i]/(myArea*num_steps*2));
        }
        fprintf(fp_out, "\n");
    }
    fclose(fp_out);
//
//  Val DEN for each PEP xy Dist PG OPP side
    snprintf(xy_dist_file, MAX_FILENAME_LEN, "%s/PG_OPP_pepxy_hist_den_%s_%s_%s-%s.dat", argv[1], argv[2], argv[3], argv[5], argv[6]);
    fp_out = fopen(xy_dist_file,"w");

    for (i=0; i<3; i++) {
        curr_dist = FIRST_DIST;
        prevRad = 0.0;
        for (j=0; j<NUM_BUCKETS; j++) {
            currRad = (j+1) * BUCKET_DIST * 10;
            myArea =  (pow(currRad,2) - pow(prevRad, 2)) * M_PI;
            prevRad = currRad;
            fprintf(fp_out, "%f ", dmpg_OPP_xy_reg[j][i]/(myArea*num_steps*2));
        }
        fprintf(fp_out, "\n");
    }
    fclose(fp_out);





    printf("DONE get_lip_den\n");
}

