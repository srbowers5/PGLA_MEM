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
#include "get_hel_ulr.h"
#include "get_tilt.h"
#include <errno.h>
#include <unistd.h>

int Num_rem;

int num_steps;
int First_aa;
int Last_aa;


#define MAX_STRUCT  300000

char helixAAUp[MAX_STRUCT][NUM_ANGLES+1];
char helixAALo[MAX_STRUCT][NUM_ANGLES+1];

float rot_ca_coord[NUM_AA][3];
float ro_angles[NUM_AA];


void get_pep_com_ulr(bool is_upper, int first_aa, int last_aa, float * xValCOM, float *yValCOM, float *zValCOM, ReadPdb * pdb, ReadDcd * dcd) {
    int atom_num;
    float xVal, yVal, zVal;
    int num_pep_atoms;
    int i;
    int CA_atom, CB_atom, N_atom, O_atom;
    int numAtoms;


    xVal = 0.0;
    yVal = 0.0;
    zVal = 0.0;
    for (i=first_aa-1; i<last_aa; i++) {
        if (is_upper == true) {
            CA_atom =  pdb->pep_up_CA[i];
            CB_atom =  pdb->pep_up_CB[i];
            N_atom =  pdb->pep_up_N[i];
            O_atom =  pdb->pep_up_O[i];
        } else {
            CA_atom =  pdb->pep_low_CA[i];
            CB_atom =  pdb->pep_low_CB[i];
            N_atom =  pdb->pep_low_N[i];
            O_atom =  pdb->pep_low_O[i];
        }
        xVal += dcd->xValArray[0][CA_atom-1];
        yVal += dcd->yValArray[0][CA_atom-1];
        zVal += dcd->zValArray[0][CA_atom-1];
        xVal += dcd->xValArray[0][CB_atom-1];
        yVal += dcd->yValArray[0][CB_atom-1];
        zVal += dcd->zValArray[0][CB_atom-1];
        xVal += dcd->xValArray[0][N_atom-1];
        yVal += dcd->yValArray[0][N_atom-1];
        zVal += dcd->zValArray[0][N_atom-1];
        xVal += dcd->xValArray[0][O_atom-1];
        yVal += dcd->yValArray[0][O_atom-1];
        zVal += dcd->zValArray[0][O_atom-1];
    }
    numAtoms = ((last_aa+1-first_aa) * 4);
    xVal = xVal / numAtoms;
    yVal = yVal / numAtoms;
    zVal = zVal / numAtoms;

    *xValCOM = xVal;
    *yValCOM = yVal;
    *zValCOM = zVal;
}


void rotate0(float xVal, float yVal, float angle, float* newX, float* newY) {

    *newX = (cos(angle) * xVal) - (sin(angle) * yVal);
    *newY = (sin(angle) * xVal) + (cos(angle) * yVal);
}

void rotate_ca(bool isUpper, float * axis_start, float *axis_end, ReadPdb * pdb, ReadDcd * dcd) {

    int i;
    int CA_atom;
    double rotAngle;
    float rotX, rotY, rotZ;
//
//  First subtract axis_start from all CA and axis end.
    axis_end[0] -= axis_start[0];
    axis_end[1] -= axis_start[1];
    axis_end[2] -= axis_start[2];

    for (i=0; i<NUM_AA; i++) {
       if (isUpper == true) {
            CA_atom =  pdb->pep_up_CA[i];
            printf("PEP UP CA %d\n", pdb->pep_up_CA[i]);
        } else {
            CA_atom =  pdb->pep_low_CA[i];
            printf("PEP LO CA %d\n", pdb->pep_low_CA[i]);
        }
        printf("X CA %f = %f - %f\n", rot_ca_coord[i][0], dcd->xValArray[0][CA_atom-1], axis_start[0]);
        rot_ca_coord[i][0] = dcd->xValArray[0][CA_atom-1] - axis_start[0];
        rot_ca_coord[i][1] = dcd->yValArray[0][CA_atom-1] - axis_start[1];
        rot_ca_coord[i][2] = dcd->zValArray[0][CA_atom-1] - axis_start[2];
        if (((i+1) >= First_aa) && ((i+1)<=Last_aa)) {
            printf("AA0 %d (%f %f %f)\n", (i+1), rot_ca_coord[i][0],  rot_ca_coord[i][1], rot_ca_coord[i][2]);
        }
    }
    axis_start[0] = axis_start[1] = axis_start[2] = 0.0;
    printf("AXIS END0 (%f %f %f)\n",  axis_end[0],  axis_end[1], axis_end[2]);

//
//  Next rotate the atoms around Z axis so axis_end[1] = 0
//
//  Get angle to rotate
    rotAngle = atan(axis_end[1]/ axis_end[0]);
    rotate0 (axis_end[0], axis_end[1], -rotAngle, &rotX, &rotY);
    axis_end[0] = rotX;
    axis_end[1] = rotY;
#ifdef DO_DEBUG
    printf("AXIS END1 (%f %f %f)\n",  axis_end[0],  axis_end[1], axis_end[2]);
#endif
    for (i=0; i<NUM_AA; i++) {
        rotate0 (rot_ca_coord[i][0], rot_ca_coord[i][1], -rotAngle, &rotX, &rotY);
        rot_ca_coord[i][0] = rotX;
        rot_ca_coord[i][1] = rotY;
#ifdef DO_DEBUG
        if (((i+1) >= First_aa) && ((i+1)<=Last_aa)) {
            printf("AA1 %d (%f %f %f)\n", (i+1), rot_ca_coord[i][0],  rot_ca_coord[i][1], rot_ca_coord[i][2]);
        }
#endif
    }
//
//  Now rotate around Y access to put helical axis on X axis.
    rotAngle = atan(axis_end[2]/ axis_end[0]);
    rotate0 (axis_end[0], axis_end[2], -rotAngle, &rotX, &rotZ);
    axis_end[0] = rotX;
    axis_end[2] = rotZ;
#ifdef DO_DEBUG
    printf("AXIS END2 (%f %f %f)\n",  axis_end[0],  axis_end[1], axis_end[2]);
#endif
    for (i=0; i<NUM_AA; i++) {
        rotate0 (rot_ca_coord[i][0], rot_ca_coord[i][2], -rotAngle, &rotX, &rotZ);
        rot_ca_coord[i][0] = rotX;
        rot_ca_coord[i][2] = rotZ;
#ifdef DO_DEBUG
        if (((i+1) >= First_aa) && ((i+1)<=Last_aa)) {
            printf("AA2 %d (%f %f %f)\n", (i+1), rot_ca_coord[i][0],  rot_ca_coord[i][1], rot_ca_coord[i][2]);
        }
#endif
    }
    if (isUpper == false) {
        rotate0 (axis_end[1], axis_end[2], M_PI, &rotY, &rotZ);
        axis_end[1] = rotY;
        axis_end[2] = rotZ;
        for (i=0; i<NUM_AA; i++) {
            rotate0 (rot_ca_coord[i][1], rot_ca_coord[i][2], M_PI, &rotY, &rotZ);
            rot_ca_coord[i][1] = rotY;
            rot_ca_coord[i][2] = rotZ;
#ifdef DO_DEBUG
            if (((i+1) >= First_aa) && ((i+1)<=Last_aa)) {
                printf("AA3 %d (%f %f %f)\n", (i+1), rot_ca_coord[i][0],  rot_ca_coord[i][1], rot_ca_coord[i][2]);
            }
#endif
        }

    }


    if ( axis_end[0] > 0) {
        rotate0 (axis_end[0], axis_end[1], M_PI, &rotX, &rotY);
        axis_end[0] = rotX;
        axis_end[1] = rotY;
#ifdef DO_DEBUG
        printf("AXIS END3 (%f %f %f)\n",  axis_end[0],  axis_end[1], axis_end[2]);
#endif
        for (i=0; i<NUM_AA; i++) {
            rotate0 (rot_ca_coord[i][0], rot_ca_coord[i][1], M_PI, &rotX, &rotY);
            rot_ca_coord[i][0] = rotX;
            rot_ca_coord[i][1] = rotY;
#ifdef DO_DEBUG
            if (((i+1) >= First_aa) && ((i+1)<=Last_aa)) {
                printf("AA3 %d (%f %f %f)\n", (i+1), rot_ca_coord[i][0],  rot_ca_coord[i][1], rot_ca_coord[i][2]);
            }
#endif
        }
    }
        
}




float get_ro_angles(bool is_upper, int first_aa, int last_aa) {
    int aa;
    float angle;
    for (aa=first_aa; aa<=last_aa; aa++) {
        angle = atan(rot_ca_coord[aa-1][2]/rot_ca_coord[aa-1][1]) * (180.0/M_PI);
        printf("UP %d (%f, %f, %f)   %f\n", aa, rot_ca_coord[aa-1][0],rot_ca_coord[aa-1][1],rot_ca_coord[aa-1][2], angle);
        if (rot_ca_coord[aa-1][1] > 0) {
            angle = 180 - angle;
        } else if (rot_ca_coord[aa-1][2] > 0){
            angle = -angle;
        } else {
            angle = 360 - angle;
        }
        printf("Angle end %f\n", angle);
        ro_angles[aa-1] = angle;

    }
    return 0;
}


float get_tilt(bool isUpper, char * helixStr, float * axis_start, float *axis_end, int first_aa, int last_aa, ReadPdb * pdb, ReadDcd * dcd) {
//    float xNValCOM, yNValCOM, zNValCOM;
//    float xCValCOM, yCValCOM, zCValCOM;
    int i;
    int num_aa;
    int nt_end;
    int ct_start;
    float xyDistSq, xyDist, zDist, tilt;

    num_aa = (last_aa - first_aa) + 1;
    if (num_aa & 0x01) {
        nt_end = first_aa + (num_aa/2);
        ct_start = first_aa + (num_aa/2);
    } else {
        nt_end = first_aa + (num_aa/2) - 1;
        ct_start = first_aa + (num_aa/2);
    }

    get_pep_com_ulr(isUpper, first_aa, nt_end, &(axis_start[0]), &(axis_start[1]), &(axis_start[2]), pdb, dcd);
    printf("COM1 = (%f %f %f)\n", axis_start[0], axis_start[1], axis_start[2]);

    get_pep_com_ulr(isUpper, ct_start, last_aa, &(axis_end[0]), &(axis_end[1]), &(axis_end[2]), pdb, dcd);
    printf("COM2 = (%f %f %f)\n", axis_end[0], axis_end[1], axis_end[2]);


    for (i=first_aa-2; i<last_aa-1; i++) {
        if (helixStr[i] != 'H') {
            return -361.0;
        }
    }


    xyDistSq = pow((axis_start[0] - axis_end[0]), 2) + pow((axis_start[1] - axis_end[1]), 2);
    xyDist = pow(xyDistSq, 0.5);
    zDist = axis_start[2] -  axis_end[2];
    printf("XY %f %f %f\n", xyDistSq, xyDist, zDist);
    if (isUpper) {
        tilt = (atan(zDist/xyDist) * (180.0/M_PI)) + 90.0;
    printf("TILT UP %f \n", tilt);
    } else {
        tilt = 90 - (atan(zDist/xyDist) * (180.0/M_PI));
    printf("TILT LO %f \n", tilt);
    }
    printf("TILT %f \n", tilt);


    return tilt;

}


#ifdef OLD


float get_tilt(bool isUpper, char * helixStr, int first_aa, int last_aa, ReadPdb * pdb, ReadDcd * dcd) {

    float xNValCOM, yNValCOM, zNValCOM;
    float xCValCOM, yCValCOM, zCValCOM;
    int i;
    int num_aa;
    int nt_end;
    int ct_start;
    float xyDistSq, xyDist, zDist, tilt;

    for (i=first_aa-2; i<last_aa-1; i++) {
        if (helixStr[i] != 'H') {
            return -361.0;
        }
    }

    num_aa = (last_aa - first_aa) + 1;
    if (num_aa & 0x01) {
        nt_end = first_aa + (num_aa/2);
        ct_start = first_aa + (num_aa/2);
    } else {
        nt_end = first_aa + (num_aa/2) - 1;
        ct_start = first_aa + (num_aa/2);
    }

    get_pep_com_ulr(isUpper, helixStr, first_aa, nt_end, &xNValCOM, &yNValCOM, &zNValCOM, pdb, dcd);
    get_pep_com_ulr(isUpper, helixStr, ct_start, last_aa, &xCValCOM, &yCValCOM, &zCValCOM, pdb, dcd);

    xyDistSq = pow((xNValCOM - xCValCOM), 2) + pow((yNValCOM - yCValCOM), 2);
    xyDist = pow(xyDistSq, 0.5);
    zDist = zNValCOM -  zCValCOM;
    if (isUpper) {
        tilt = (atan(zDist/xyDist) * (180.0/M_PI)) + 90.0;
    } else {
        tilt = 90 - (atan(zDist/xyDist) * (180.0/M_PI));
    }

    return tilt;
}

#endif

int read_helix_file(bool isUp, char * fileName, int numStr) {
    FILE * inFd;
    char line[128];
    char * lineData;
    size_t  lenField = 128;
    int i;
    int nRead;
    int errCnt;

    inFd = fopen(fileName, "r");
    printf("READ %s\n", fileName);
    if (inFd == 0) {
        printf("Cannot open file (%s)\n", fileName);
        exit(90);
    }

    lineData = &(line[0]);
    for (i=0; i<numStr; i++) {
        lineData = &(line[0]);
        nRead = getline(&lineData, &lenField, inFd);
        errCnt = 0;
        while (nRead < 0) {
            printf("GET LINE ERROR %d\n", errno);
            sleep(1);
            errCnt++;
            if (errCnt > 0) {
                return -1;
            }
            nRead = getline(&lineData, &lenField, inFd);
        }

        if (nRead < NUM_ANGLES) {
            printf("CHK %d numStr %d, Read %d, line(%s)\n", i, numStr, nRead, lineData);
            printf("Error - not enough fields\n");
            fflush(stdout);
            return -1;
        }
        if (isUp) {
            memcpy(&(helixAAUp[i][0]), line, NUM_ANGLES);
            helixAAUp[i][NUM_ANGLES] = 0;
        } else {
            memcpy(&(helixAALo[i][0]), line, NUM_ANGLES);
            helixAALo[i][NUM_ANGLES] = 0;
        }
    }
    return 0;
}

//
//ρ =   (Sum ((φ + (ref_aa − i) · 100°)%360)) /N
//
float get_ro_avg(int first_aa, int last_aa, int ref_aa) {
    int aa;
    float val;
    float totVal;

    totVal = 0.0;
    for (aa=first_aa; aa<last_aa+1; aa++) {
        val = (ro_angles[aa-1] + ((ref_aa - aa) * 100.0));
        while (val >= 360.0) {
            val -= 360.0;
        }
        while (val < 0.0) {
            val += 360.0;
        }
        totVal += val;
    }
    totVal = totVal / ((last_aa + 1) - first_aa);
    return totVal;
}



void write_ro_angles(FILE * outFpRo, bool isUpper, float * axis_start, float *axis_end, int first_aa, int last_aa, 
  ReadPdb * pdb, ReadDcd * dcd) {

    int aa;
    float ro_avg;

    rotate_ca(isUpper, axis_start, axis_end, pdb, dcd);
    get_ro_angles(isUpper, first_aa, last_aa);
    for (aa=first_aa; aa<last_aa+1; aa++) {
        fprintf(outFpRo, "%f ", ro_angles[aa-1]);
    }
    if ((first_aa <= LYS12) && (last_aa >= LYS12)) {
        ro_avg = get_ro_avg(first_aa, last_aa, LYS12);
        fprintf(outFpRo, "%f ", ro_avg);
    }
    if ((first_aa <= LYS15) && (last_aa >= LYS15)) {
        ro_avg = get_ro_avg(first_aa, last_aa, LYS15);
        fprintf(outFpRo, "%f ", ro_avg);
    }
    if ((first_aa <= LYS19) && (last_aa >= LYS19)) {
        ro_avg = get_ro_avg(first_aa, last_aa, LYS19);
        fprintf(outFpRo, "%f ", ro_avg);
    }
}


int main(int argc, char * argv[]) {

    char dcd_file[MAX_FILENAME_LEN];
    char pdb_file[MAX_FILENAME_LEN];
    char cm_file[MAX_FILENAME_LEN];
    char outFile[MAX_FILENAME_LEN];
    const char * allHelix="HHHHHHHHHHHHHHHHHHHHH";
    FILE *fpOut;
    FILE *outFpUp;
    FILE *outFpLo;
    int i, j, k;
    ReadDcd * dcd;
//    char * allHelix;
    bool helix;
    int first_aa;
    int last_aa;
    int cnt;
    float tiltU, tiltL;
    float start_com[3];
    float end_com[3];

//    allHelix = allHelixStr.c_str();


//    arg[1] = Input_dir
//    arg[2] = tr
//    arg[3] = rep
//    arg[4] = prefix
//    arg[5] = first_str
//    arg[6] = last_str
//    arg[7] = Use only helix
//    arg[8] = first_aa
//    arg[9] = last_aa
//    arg[10] = upper helix file
//    arg[11] = lower helix file
//    arg[12] = output file

    int first_step, last_step, step;
    int errVal;
    first_step = atoi(argv[5]);
    last_step = atoi(argv[6]);
    helix = atoi(argv[7]);
    first_aa = atoi(argv[8]);
    First_aa = first_aa;
    last_aa = atoi(argv[9]);
    Last_aa = last_aa;

    if ((last_step + 1 - first_step) > MAX_STRUCT) {
        printf("ERROR - Need to modify #define MAX_STRUCT, now set to %d\n", MAX_STRUCT);
        exit(88);
    }

    if (helix) {
        errVal = read_helix_file(true, argv[10], (last_step + 1 - first_step));
        if (errVal < 0) {
            printf("BAD FILE %s\n", argv[10]);
            exit(87);
        }
        errVal = read_helix_file(false, argv[11], (last_step + 1 - first_step));
        if (errVal < 0) {
            printf("BAD FILE %s\n", argv[10]);
            exit(86);
        }
    }

    snprintf(pdb_file, MAX_FILENAME_LEN, "%s/spt.pdb", argv[1]);
    ReadPdb pdb(pdb_file);
    pdb.find_lipid_tails(NUM_LIPIDS_PER_LEAFLET);
    printf("PDB %s\n", pdb_file);
    pdb.find_heavy_atoms(NUM_LIPIDS_PER_LEAFLET,NUM_PEP_PER_LEAFLET);

    snprintf(outFile, MAX_FILENAME_LEN, "%s", argv[12]);
    fpOut = fopen(outFile, "w");
    printf("OPEN %s\n", outFile);

    snprintf(outFile, MAX_FILENAME_LEN, "%s", argv[13]);
    outFpUp = fopen(outFile, "w");
    printf("OPEN %s\n", outFile);

    snprintf(outFile, MAX_FILENAME_LEN, "%s", argv[14]);
    outFpLo = fopen(outFile, "w");
    printf("OPEN %s\n", outFile);
    fflush(stdout);


    cnt = 0;
    for (step=first_step; step < (last_step+1); step++) {
        snprintf(dcd_file, MAX_FILENAME_LEN, "%s/output/%s%s_%05d_%s.dcd", argv[1], argv[4], argv[2], step, argv[3]);
        printf("READ %s\n", dcd_file);
        fflush(stdout);
        dcd = new ReadDcd(dcd_file);

        if (helix) {
            tiltU = get_tilt(true, &(helixAAUp[cnt][0]), &(start_com[0]), &(end_com[0]), first_aa, last_aa, &pdb, dcd);
            tiltL = get_tilt(false, &(helixAALo[cnt][0]), &(start_com[0]), &(end_com[0]), first_aa, last_aa, &pdb, dcd);
        } else {
            tiltU = get_tilt(true, (char *) allHelix, &(start_com[0]), &(end_com[0]), first_aa, last_aa, &pdb, dcd);
            write_ro_angles(outFpUp, true, &(start_com[0]), &(end_com[0]), first_aa, last_aa, &pdb, dcd);
            fprintf(outFpUp, "\n");

            tiltL = get_tilt(false, (char *) allHelix, &(start_com[0]), &(end_com[0]), first_aa, last_aa, &pdb, dcd);
            write_ro_angles(outFpLo, false, &(start_com[0]), &(end_com[0]), first_aa, last_aa, &pdb, dcd);
            fprintf(outFpLo, "\n");
                printf("TILT UP OUT  %f\n", tiltU);
                printf("TILT LO OUT  %f\n", tiltL);

        }
        fprintf(fpOut,"%d %f %f\n", step, tiltU, tiltL);
        cnt++;
        delete dcd;
    }
    fclose(fpOut);
    fclose(outFpUp);
    fclose(outFpLo);

    printf("DONE get_tilt\n");
}
