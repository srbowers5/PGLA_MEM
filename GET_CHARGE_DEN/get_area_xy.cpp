#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include <stdarg.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>
#include <cmath>
#include <arpa/inet.h>
#include <math.h>

#define INT_STEPS   1000

/*
*   get_areas()
*   Get the area of each of the 3 areas making up the square.
*   We will always have either 2 or 3 areas.
*
*   We assume the sides of the square are 1.
*/
void get_areas(float xVal, float yVal, float * a1Val, float * a2Val, float * a3Val) {


    float r1Val;
    float r2Val; 
    float r3Val; 
    float r4Val;

    float x2Val;
    float x3Val;
    float currX;
    float currY;
    float Ydelta;
    float intDist;
    float areaVal;

    float yDelta;
    int i;
//
//   R1val is closest corner.
//   R2Val is First ring Radius > R1Val
//   R3Val is Next Ring Radius 
//   R4Val is Radius of farthest corner.
    r1Val =  pow((pow(xVal,2) + pow(yVal,2)),0.5);
    r2Val = ((int)(r1Val)) + 1;
    r3Val = r2Val + 1;
    r4Val =  pow((pow(xVal+1,2) + pow(yVal+1,2)),0.5);


//   Get X intersect of R2
//   X = sqrt( (R**2) + (Yval**2))
    x2Val = pow((pow(r2Val,2) - pow(yVal,2)),0.5);
//    printf("X2 %f %f, R2 %f, Y %f\n", xVal, x2Val, r2Val, yVal);
    if (x2Val > (xVal + 1)) {
        x2Val = xVal +1;
    }

//
//   Integrate to get area within R2 
    intDist = (x2Val - xVal) / INT_STEPS;
    currX  = xVal + (intDist/2);
    areaVal = 0.0;
//    printf("INT DIST %f\n", intDist);
    for (i=0; i<INT_STEPS; i++) {
        currY = pow((pow(r2Val,2) - pow(currX,2)),0.5);
        if (currY > (yVal +1)) {
            currY = (yVal +1);
        }
        yDelta = currY - yVal;
        areaVal += (yDelta * intDist);
        currX += intDist;
    }
    *a1Val = areaVal;
    if (*a1Val > 1.0) {
        printf("ERROR A1 to large %f %f %f\n", xVal, yVal, areaVal);
    }

//
//   Intergrate to get Area between R3 and R4.
    if (r3Val > r4Val) {
       areaVal = 0.0;
    } else {
        x3Val = pow((pow(r3Val,2) - pow(yVal+1,2)),0.5);
        if (x3Val < xVal) {
            x3Val = xVal;
        }

        intDist = ((xVal+1)-x3Val) / INT_STEPS;
        areaVal = 0.0;
        for (i=0; i<INT_STEPS; i++) {
            currX = x3Val + (intDist/2);
            currY =  pow((pow(r3Val,2) - pow(currX,2)),0.5);
            if (currY < yVal) {
                currY = yVal;
            }
            yDelta = (yVal + 1) - currY;
            areaVal += (yDelta * intDist);
            currX += intDist;
        }
    }
    *a3Val = areaVal;
    if ((*a1Val + *a3Val) > 1) {
        printf("ERROR2 sum to big (%f,%f) %f %f\n", xVal, yVal, *a1Val, *a3Val);
    }
    *a2Val = 1 - (*a1Val + *a3Val);

//    printf(" RESULT X %f Y %f,     A1 %f,     A2 %f,     A3 %f\n", xVal, yVal, *a1Val, *a2Val, *a3Val);
    if ( (*a1Val+*a2Val+*a3Val) > 1.0001) {
        printf("ERROR in Area (%f, %f), %f %f %f\n", xVal, yVal,   *a1Val, *a2Val, *a3Val);
    }
    if ( (*a1Val+*a2Val+*a3Val) < .9999) {
        printf("ERROR in Area (%f, %f), %f %f %f\n", xVal, yVal,   *a1Val, *a2Val, *a3Val);
    }
}
