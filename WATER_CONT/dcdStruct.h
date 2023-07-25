
#ifndef DCD_STRUCT_
#define DCD_STRUCT_ 

typedef struct hdr_ {
    char    junk[4];
    char    hdr[4];
    int     nstr;
    int     dumi[8];
    float   dumr;
    int     dumi2[9];
    int     dumi3[3];
    int     nrem;
} HDR, *HDR_PTR;


#define MAX_REM  200
typedef struct rem_ {
    char    rem[80];
} REM, *REM_PTR;



typedef struct aft_rem_ {
    char junk[8];
    int numAtom;
    char junk2[8];
} AFT_REM, *AFT_REM_PTR;

#define SIZE_IBOX   6
typedef struct ibox_ {
    double ibox[SIZE_IBOX];
} IBOX, *IBOX_PTR;

#endif /* end DCD_STRUCT_ */
