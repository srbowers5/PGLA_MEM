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


float get_min_xy_dist(float x1, float y1, float x2, float y2, float uCellXY) {

    float xDiff, yDiff;
    float uCellHalf;
    float distSq;
    float myDist;

    uCellHalf = uCellXY/2;

    xDiff = x1 - x2;
    yDiff = y2 - y2;

//    printf("GET MIN %f %f %f\n", xDiff, yDiff, uCellXY);
    if (xDiff > uCellHalf) {
        xDiff = xDiff - uCellXY;
    } else if (xDiff < -uCellHalf) {
        xDiff = xDiff + uCellXY;
    }

    if (yDiff > uCellHalf) {
        yDiff = yDiff - uCellXY;
    } else if (yDiff < -uCellHalf) {
        yDiff = yDiff + uCellXY;
    }
   // printf("GET MOD MIN %f %f %f\n", xDiff, yDiff, uCellXY);

    distSq = (pow(xDiff,2) + pow(yDiff,2));
    myDist = pow(distSq, 0.5);

   // printf("GET MOD MIN %f %f %f %f %f\n", xDiff, yDiff, uCellXY, distSq, myDist);

    return myDist;
}


float get_min_xyz_dist(float x1, float y1, float z1, float x2, float y2, float z2, float uCellXY) {

    float xDiff, yDiff, zDiff;
    float uCellHalf;
    float distSq;
    float myDist;

    uCellHalf = uCellXY/2;

    xDiff = x1 - x2;
    yDiff = y1 - y2;
    zDiff = z1 - z2;

//    printf("GET MIN %f %f %f\n", xDiff, yDiff, uCellXY);
    if (xDiff > uCellHalf) {
        xDiff = xDiff - uCellXY;
    } else if (xDiff < -uCellHalf) {
        xDiff = xDiff + uCellXY;
    }

    if (yDiff > uCellHalf) {
        yDiff = yDiff - uCellXY;
    } else if (yDiff < -uCellHalf) {
        yDiff = yDiff + uCellXY;
    }
   // printf("GET MOD MIN %f %f %f\n", xDiff, yDiff, uCellXY);

    distSq = (pow(xDiff,2) + pow(yDiff,2) + pow(zDiff,2));
    myDist = pow(distSq, 0.5);

   // printf("GET MOD MIN %f %f %f %f %f\n", xDiff, yDiff, uCellXY, distSq, myDist);

    return myDist;
}


