#include "stdio.h"
#include "stdlib.h"
#include "string.h"




int main(int argc, char * argv[]) {


    int intVal1;

    int intVal2;

    float flVal1;
    float flVal2;

    flVal1 = 0.9999;
    flVal2 = -0.9999;
    intVal1 = int(flVal1);
    intVal2 = int(flVal2);

    printf("VAL1 %d %f VAL2 %d %f\n", intVal1, flVal1, intVal2, flVal2);
}

