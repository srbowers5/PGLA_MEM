#include <stdio.h>
#include <math.h>

//
//   SCD = ((3 * cosSq) - 1) / 2;
//
//   cos = x / ( (x**2 + y**2) **0.5)
//   cos**2 = z**2 / (x**2 + y**2 + z**2)

float get_one_scd(float xC, yC, zC, xH, yH, zH) {
    float xDiffSq, yDiffSq, zDiffSq;
    float cosSq, scd;

    xDiffSq = pow( (xC - xH), 2);
    yDiffSq = pow( (yC - yH), 2);
    zDiffSq = pow( (zC - zH), 2);
    cosSq = xDiffSq / (xDiffSq + yDiffSq + zDiffSq);

    scd = ((3 * cosSq) - 1)/ 2;
    
    return scd;
