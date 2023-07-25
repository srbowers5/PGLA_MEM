#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include <stdarg.h>
#include <iostream>
#include <fstream>
#include <arpa/inet.h>
#include <vector>
#include "ReadCM.h"
void ReadCM::addOneLine(const char * line) {
    int lipNum, aaNum, offset;
    STEP_CONT_PTR cont_ptr;

    offset = 6;
    cont_ptr = (STEP_CONT_PTR) malloc(sizeof(STEP_CONT));
    for (aaNum=0; aaNum<NUM_AAS; aaNum++) {
        for (lipNum=0; lipNum<NUM_LIPIDS; lipNum++) {
            cont_ptr->contacts[aaNum][lipNum] = line[offset];
            if (line[offset] == '1') {
                cont_ptr->num_lip_cont[lipNum]++;
            }
//            printf("ADD %d %d %d %d = %d\n", numSteps, aaNum, lipNum, offset, cont_ptr->contacts[aaNum][lipNum]);
            offset++;
        }
    }
    contacts_ptr.push_back(cont_ptr);
    numSteps++;
}


void ReadCM::addData(const char * filename) {

    std::string line;
    std::ifstream cmFd(filename);
    if (cmFd.is_open() == false) {
        printf("Failed to open %s\n", filename);
        exit(99);
    }

//    printf("ADD data %s\n", filename);
    while (std::getline(cmFd, line)) {
//        printf("ADD one line %s\n", line.c_str());
        addOneLine(line.c_str());
    }
}
