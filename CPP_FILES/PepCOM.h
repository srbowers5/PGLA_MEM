
#include "system_conf.h"
#include "ReadPdb.h"

class PepCOM {


    public:

    double pep_com[3];
    std::vector <double> aa_comX;
    std::vector <double> aa_comY;
    std::vector <double> aa_comZ;



    PepCOM(ReadPdb * pdb, ReadDcd *dcd, int pep_num ) {
        int atom_num;
        float xVal, yVal, zVal;
        int num_pep_atoms;
        std::vector<int>::iterator first_it, end_it, it;
        std::vector<int>::iterator first_aa_it, end_aa_it, aa_it;
        ReadPdb::PEP_STRUCT * pep_struct_ptr;
        int aa;

        xVal = yVal = zVal = 0.0;
        pep_struct_ptr = &(pdb->pep_struct[pep_num]);
        first_it = std::begin(pep_struct_ptr->pep_heavy_atoms);
        end_it = std::end(pep_struct_ptr->pep_heavy_atoms);
        for (it = first_it; it != end_it; ++it) {
            atom_num=*it;
            xVal += dcd->xValArray[0][atom_num-1];
            yVal += dcd->yValArray[0][atom_num-1];
            zVal += dcd->zValArray[0][atom_num-1];
            printf("GET PEP ATOM %d\n", atom_num);
            fflush(stdout);

        }
        num_pep_atoms = pep_struct_ptr->pep_heavy_atoms.size();
        pep_com[0] = xVal / num_pep_atoms;
        pep_com[1] = yVal / num_pep_atoms;
        pep_com[2] = zVal / num_pep_atoms;


        for (aa=0; aa<pep_struct_ptr->numAA; aa++) {
            printf("GET AA ATOM %d\n", aa);
            fflush(stdout);
            first_aa_it = std::begin(pep_struct_ptr->pep_ha[aa]);
            end_aa_it = std::end(pep_struct_ptr->pep_ha[aa]);

            xVal = yVal = zVal = 0.0;
            for (aa_it = first_aa_it; aa_it != end_aa_it; ++aa_it) {
                atom_num=*it;
                xVal += dcd->xValArray[0][atom_num-1];
                yVal += dcd->yValArray[0][atom_num-1];
                zVal += dcd->zValArray[0][atom_num-1];
            }
            num_pep_atoms = pep_struct_ptr->pep_ha[aa].size();
            aa_comX.push_back(xVal / num_pep_atoms);
            aa_comY.push_back(yVal / num_pep_atoms);
            aa_comZ.push_back(zVal / num_pep_atoms);
        }
        printf("DONE PepCOM %d\n", pep_num);\
        fflush(stdout);
    }   /* End PepCom() */

};
