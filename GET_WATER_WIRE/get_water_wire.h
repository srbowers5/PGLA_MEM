#ifndef GET_WATER_WIRE_H
#define GET_WATER_WIRE_H

#define INSERT_DEPTH    18.0
#define MAX_HBOND_DIST   4.0
#define MAX_HBOND_ANGLE 30.0

#define MAX_HBONDS      2000
#define MAX_CLUSTS      MAX_HBONDS

typedef struct water_mol_ {
    int mol_index;
    float ox_mol_coord[3];
    float h1_mol_coord[3];
    float h2_mol_coord[3];
} WATER_MOL, *WATER_MOL_PTR;

typedef struct hbond_ {
    int mol1;
    int mol2;
    int hVal;
} HBOND, *HBOND_PTR;

#endif /* GET_WATER_WIRE_H */
