class ReadPsf {

#define NUM_ATOM_STRING       "NATOM"
#define MIN_PSF_ATOM_SIZE     58
#define MAX_PSF_ATOMS         50000

  public:
    typedef struct psf_atom_ {
        int atom_num;
        int mon_num;
        char mol_type[5];
#define DMPC_MOL    1
#define DMPG_MOL    2
#define PGL1_MOL    3
#define PGL2_MOL    4
#define WAT_MOL     5
#define SOD_MOL     6
#define CLA_MOL     7
#define PGL3_MOL    8
#define PGL4_MOL    9
	int mol_type_val;
        char    monomer_type[5];

        char aa_type[5];

#define NOT_PEP_MON 0
#define SIDE_MON    1
#define BB_MON      2
        int atom_type;


#define LEAF_UN     0
#define LEAF_UP     1
#define LEAF_LO     2
        int leaf;
        double atom_charge; 
    } PSF_ATOM, *PSF_ATOM_PTR;


        ReadPsf(const char * psfFile) {
            std::string line;
            int retVal;
            int prev_atom_num;
            const char * natom_ptr;

            std::ifstream inFd(psfFile);
            if (inFd.is_open() == false) {
                printf("Failed to open %s\n", psfFile);
                exit(99);
            }
            atoms.reserve(MAX_PSF_ATOMS);
            while (std::getline(inFd, line)) {
                natom_ptr = strstr(line.c_str(), NUM_ATOM_STRING);
                if (natom_ptr != 0) {
                    break;
                }
            }
            if (natom_ptr == 0) {
                printf("Cannot find NATOM %s\n", psfFile);
                exit(98);
            }

            retVal = 1;
            prev_atom_num = 0;
            while (std::getline(inFd, line) && (retVal)) {
                retVal = parse_psf_line(line.c_str(), prev_atom_num);
                if (retVal == 1) {
                    prev_atom_num++;
                }
            }
            inFd.close();
        }
        std::vector <PSF_ATOM_PTR> atoms;
        int parse_psf_line(const char * line_ptr, int prev_atom_num);

};
