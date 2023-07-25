#ifndef READ_PDB_H
#define READ_PDB_H 

#define MAX_ATOMS  50000
#define MAX_PEP    4
#define MAX_AA     40
#define NUM_LIP_SECTIONS     5


class ReadPdb {

#define ATOM_TYPE  "ATOM  "
#define BOX_TYPE   "CRYST1"
#define END_TYPE   "END   "
#define TYPE_LEN   6


    typedef struct atom_line_ {
        char    line_type[6];
        char    atom_num[5];
        char    atom_type[5];
        char    monomer_type[5];
        char    mol_type;
        char    mol_num[5];
        char    junk1[3];
        char    xVal[8];
        char    yVal[8];
        char    zVal[8];
        char    colA[6];
        char    colB[6];
        char    chain[10];
    } ATOM_LINE, *ATOM_LINE_PTR;

    typedef struct box_line_ {
        char    line_type[6];
        char    Junk;
        char    xLen[9];
        char    yLen[9];
        char    zLen[9];
    } BOX_LINE, * BOX_LINE_PTR;
    
    typedef struct box_out_ {
        double box[3];
    } PDB_BOX, *PDB_BOX_PTR;

    public:
    typedef struct atom_out_ {
        int atom_num;
        int mol_num;
        double xVal;
        double yVal;
        double zVal;
        char atom_type[6];
        char monomer_type[6];
        char mol_type;
        char chain[11];
        char colA[7];
        char colB[7];
    } PDB_ATOM, *PDB_ATOM_PTR;

    enum mol_defines{IS_WATER, IS_DMPC, IS_DMPG, IS_DMPE, IS_PEP};
    typedef struct lipid_C_tail_ {
        int chain_num;
        int chain_pos;
        int cAtom_num;
        int hxAtom_num;
        int hyAtom_num;
        int hzAtom_num;
        mol_defines mon_type;
    } LIPID_C_TAIL, *LIPID_C_TAIL_PTR;

    typedef struct pep_struct_ {
        std::vector <int> pep_heavy_atoms;
        std::vector <int> pep_CA;
        std::vector <int> pep_CB;
        std::vector <int> pep_N;
        std::vector <int> pep_O;
        std::vector <int> pep_side_ha[MAX_AA];
        int numAA;
    } PEP_STRUCT, * PEP_STRUCT_PTR;

    typedef struct lip_struct_ {
        std::vector <int> ha_section[NUM_LIP_SECTIONS];
        mol_defines mon_type;
    } LIP_STRUCT, * LIP_STRUCT_PTR;


        ReadPdb(char * pdbFile) {
            std::string line;
            int retVal;
            int i;

            std::ifstream inFd(pdbFile);
            if (inFd.is_open() == false) {
                printf("Failed to open %s\n", pdbFile);
                exit(99);
            }
            atoms.reserve(MAX_ATOMS);
            retVal = 1;
            while (std::getline(inFd, line) && (retVal)) {
                retVal = parse_pdb_line(line.c_str());
            }
            inFd.close();
            for (i=0; i<MAX_PEP; i++) {
                pep_struct[i].numAA = 0;
            }
        }
        void print_pdb();
        int find_heavy_atoms(int lipids_per_leaf, int pep_per_leaf);
        void save_pep_atom(PEP_STRUCT_PTR pep_ptr, PDB_ATOM_PTR atom_ptr);
        int find_lipid_tails(int lipids_per_leaf);
        int find_lipid_sections();


        PDB_BOX box;
        std::vector <PDB_ATOM_PTR> atoms;
        std::vector <int> lip_up_heavy_atoms;
        std::vector <int> lip_up_Phos;
        std::vector <int> dmpc_up_heavy_atoms;
        std::vector <int> dmpg_up_heavy_atoms;
        std::vector <int> dmpe_up_heavy_atoms;
        std::vector <int> dmpc_up_P_atoms;
        std::vector <int> dmpg_up_P_atoms;
        std::vector <int> dmpe_up_P_atoms;

        std::vector <int> lip_low_heavy_atoms;
        std::vector <int> lip_low_Phos;
        std::vector <int> dmpc_low_heavy_atoms;
        std::vector <int> dmpg_low_heavy_atoms;
        std::vector <int> dmpe_low_heavy_atoms;
        std::vector <int> dmpc_low_P_atoms;
        std::vector <int> dmpg_low_P_atoms;
        std::vector <int> dmpe_low_P_atoms;

        std::vector <int> pep_up_heavy_atoms;
        std::vector <int> pep_up_CA;
        std::vector <int> pep_up_CB;
        std::vector <int> pep_up_N;
        std::vector <int> pep_up_O;

        std::vector <int> pep_low_heavy_atoms;
        std::vector <int> pep_low_CA;
        std::vector <int> pep_low_CB;
        std::vector <int> pep_low_N;
        std::vector <int> pep_low_O;
        std::vector <int> sod_ion;
        std::vector <int> cl_ion;


        PEP_STRUCT pep_struct[MAX_PEP];
        std::vector <LIP_STRUCT_PTR> lipid_sections;

        std::vector <int> water_heavy_atoms;
        std::vector <LIPID_C_TAIL_PTR> lipid_tail_c_h_up;
        std::vector <LIPID_C_TAIL_PTR> lipid_tail_c_h_low;
    private:
        int parse_pdb_line(const char * line_ptr);

};
#endif /* READ_PDB_H */
