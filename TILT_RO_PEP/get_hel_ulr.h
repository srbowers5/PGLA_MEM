#ifndef DIH_STRUCT

#define MAX_FILENAME_LEN 1024
#define NUM_PEP_PER_LEAFLET 1
#define NUM_PEP_TOTAL (NUM_PEP_PER_LEAFLET*2)
#define NUM_LIPIDS_PER_LEAFLET 49

#define DIN_STRUCT 1
#define NUM_AA     21
#define NUM_ANGLES     (NUM_AA-2)   // Ignore first and last
#define NUM_FIELDS (NUM_ANGLES*NUM_PEP_TOTAL)+1
#define MAX_PARSE_LINE_LEN 0x100000



typedef struct dih_struct_ {
    float phi_val[NUM_ANGLES];
    float psi_val[NUM_ANGLES];
} DIH_STRUCT, * DIH_STRUCT_PTR;


#endif
