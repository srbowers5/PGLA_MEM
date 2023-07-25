
#ifndef GET_LIPID_DEN
#define GET_LIPID_DEN


#define MAX_FILENAME_LEN 1024


#define Z_RANGE_VAL         30

enum mol_defines{IS_WATER, IS_DMPC, IS_DMPG, IS_DMPE, IS_PEP, IS_SOD_ION, IS_CLA_ION};

#define MIN_Z_DIST         -30.5
#define MAX_Z_DIST         30.5
#define BUCKETS_PER_Z_DIST   1
#define NUM_Z_BUCKETS        (int) ((MAX_Z_DIST-MIN_Z_DIST)+0.00001)
#define BUCKET_Z_DIST        1.0


#define MAX_XY_DIST           30
#define BUCKETS_PER_DIST   1
#define NUM_BUCKETS        (MAX_XY_DIST*BUCKETS_PER_DIST)
#define BUCKET_DIST        1.0
#define FIRST_DIST         0.50

#define NEAR_DIST           0.0
#define PROXIMAL_DIST       6.0
#define DISTANT_DIST       22.0
#define FAR_DIST           22.0
#define NUM_REG            4

#define NEAR_END 6
#define PROXIMAL_END 22
#define DIST_END 30
#define FAR_END 30


#define FIXED_LIPID_BOUND 20.0

#endif
