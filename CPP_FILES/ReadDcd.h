#ifndef READ_DCD_H
#define READ_DCD_H


 
#include "dcdStruct.h"

class ReadDcd {

    public:
        ReadDcd(char * dcdFile) {
            char * inFile;
            char * outFile;
            int i,j;
            char junk16[16];
            int arrayLen;
            int Num_rem;
            IBOX_PTR iBoxPtr;

            HDR dcdHdr;
            REM_PTR nextRem;
            AFT_REM aftRemStruct;
            IBOX iBox;
            float * xArray;
            float * yArray;
            float * zArray;

            std::ifstream inFd;
            inFd.open(dcdFile,std::ifstream::binary);
//            printf("START ");
            if (inFd) {
                inFd.read((char *)&dcdHdr,sizeof(HDR));
                Num_rem = 0;
//                printf("HDR %s\n", dcdHdr.hdr);
                printf("NSTR %d\n", dcdHdr.nstr);
                nStr = dcdHdr.nstr;
//                for (i=0; i<8; i++) {
//                    printf("DUMI %d = %d\n", i, dcdHdr.dumi[i]);
//                }
//                printf("DUMF = %.20f\n", dcdHdr.dumr);
//                for (i=0; i<9; i++) {
//                   printf("DUMI2 %d = %d\n", i, dcdHdr.dumi2[i]);
//                }
//                printf("NREM %x\n", dcdHdr.nrem);
                nRem = dcdHdr.nrem;
                for (i=0; i<dcdHdr.nrem; i++) {
                    nextRem = (REM_PTR) (malloc(sizeof(REM)));
                    if (nextRem) {
                        inFd.read((char *)(nextRem), sizeof(REM));
                        if (Num_rem < MAX_REM) {
                            allRems.push_back(nextRem);
                            Num_rem++;
                        } else {
                            printf("Too many rem in DCD, Max = %d, Num = %d %d\n", MAX_REM, i+1, dcdHdr.nrem);
                            free(nextRem);
                        }

                    }
                }
//                for (i=0; i<Num_rem; i++) {
//                    printf("REM %d - %s\n", i+1, allRems[i]->rem);
//                }
                inFd.read((char *)&aftRemStruct, sizeof(AFT_REM));
                printf("NA %d\n", aftRemStruct.numAtom);
                nAtoms = aftRemStruct.numAtom;
        
                for (j=0; j < nStr; j++) {
                    iBoxPtr = (IBOX * ) malloc(sizeof(IBOX));
                    inFd.read((char *)iBoxPtr, sizeof(IBOX));
                    boxArray.push_back(iBoxPtr);
//                    for (i=0; i<SIZE_IBOX; i++) {
//                        printf("IBOX   %d %lf\n", i, iBoxPtr->ibox[i]);
//                    }
                    inFd.read((char *)&junk16, 8);
        
                    arrayLen = sizeof(float) * (nAtoms + 2);
                    xArray = (float *) malloc(arrayLen);
                    yArray = (float *) malloc(arrayLen);
                    zArray = (float *) malloc(arrayLen);
                    inFd.read((char*)xArray, arrayLen);
                    inFd.read((char*)yArray, arrayLen);
                    inFd.read((char*)zArray, arrayLen);
        
                    xValArray.push_back(xArray);
                    yValArray.push_back(yArray);
                    zValArray.push_back(zArray);
//                    for (i=0; i<aftRemStruct.numAtom; i++) {
//                        printf("coor[%d]\t%f\t%f\t%f\n", i, xArray[i], yArray[i], zArray[i]);
//                    }
                }
        
            } else {
                printf("Cannot open file (%s)\n", inFile);
                return;
            }
//            printf("DONE\n");
        }

        ~ReadDcd() {
            while (!xValArray.empty()) {
                free(xValArray.back());
                xValArray.pop_back();
            }
            while (!yValArray.empty()) {
                free(yValArray.back());
                yValArray.pop_back();
            }
            while (!zValArray.empty()) {
                free(zValArray.back());
                zValArray.pop_back();
            }
        }

        void printDcd(int index);


        float uCell[3];

        int nRem;
        std::vector <REM_PTR> allRems;;
        std::vector <IBOX_PTR> boxArray;
        std::vector <float *> xValArray;
        std::vector <float *> yValArray;
        std::vector <float *> zValArray;
        int nAtoms;
        int nStr;

    private:
};

#endif /* READ_DCD_H */
