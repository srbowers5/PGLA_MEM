#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include <stdarg.h>
#include <iostream>
#include <fstream>
#include <arpa/inet.h>
#include <vector>
#include "ReadDcd.h"


void ReadDcd::printDcd(int index) {
    int i;

    printf("BOX %f %f %f %f %f %f\n", 
      boxArray[index]->ibox[0], boxArray[index]->ibox[1], boxArray[index]->ibox[2],
      boxArray[index]->ibox[3], boxArray[index]->ibox[4], boxArray[index]->ibox[5]);
    for (i=0; i<nAtoms; i++) {
         printf("coor[%d]\t%f\t%f\t%f\n", i, xValArray[index][i], yValArray[index][i], zValArray[index][i]);
    }
}
