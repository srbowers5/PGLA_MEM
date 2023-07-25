
#ifndef GET_LIPID_DEN
#define GET_LIPID_DEN


#define MAX_FILENAME_LEN 1024


#define NUM_LIPID_LEAFLET   49
#define Z_RANGE_VAL         42


#define MIN_Z_DIST         -35.5
#define MAX_Z_DIST         35.5
#define BUCKETS_PER_Z_DIST   1
#define NUM_Z_BUCKETS        (int) ((MAX_Z_DIST-MIN_Z_DIST)+0.00001)
#define BUCKET_Z_DIST        1.0


#define MAX_XY_DIST           30
#define BUCKETS_PER_DIST   1
#define NUM_BUCKETS        (int) ((MAX_XY_DIST*BUCKETS_PER_DIST)+0.00001)
#define BUCKET_DIST        1.0
#define FIRST_DIST         0.50

#define PROXIMAL_DIST       6.0
#define DISTANT_DIST       21.5
#define NUM_REG_R             3

//#define INTERIOR_Z_DIST     17.5
//#define BOUND_Z_DIST        24.0
#define NUM_REG_Z              4


//       Hydrophobic		Phos			Interface	Bulk
#define HYDROPHOBIC_Z_DIST	    12.0
#define PHOS_Z_DIST                 21.0
#define INTERFACE_Z_DIST            28.0
#define BULK_Z_DIST                 MAX_Z_DIST


#define UCELL_MAX_X    29
#define UCELL_MAX_Y    29
#define UCELL_MAX_Z    41
#define NUM_AREAS      3

#define IS_DMPC        1
#define IS_DMPG        2
#define IS_ION         3
#define IS_LIPID       4

#endif
