// Lagrange curve arclength integrands

#include <math.h>

double quad(int n, double args[n]){
    // args elements are [t, x1, x2, x3, y1, y2, y3, z1, z2, z3]
    //           indices [0,  1,  2,  3,  4,  5,  6,  7,  8,  9]
    return sqrt(
            pow((4.0*args[0]-3.0)*args[7]+(4.0-8.0*args[0])*args[8]+(4.0*args[0]-1.0)*args[9], 2.0)+
            pow((4.0*args[0]-3.0)*args[4]+(4.0-8.0*args[0])*args[5]+(4.0*args[0]-1.0)*args[6], 2.0)+
            pow((4.0*args[0]-3.0)*args[1]+(4.0-8.0*args[0])*args[2]+(4.0*args[0]-1.0)*args[3], 2.0)
            );
}