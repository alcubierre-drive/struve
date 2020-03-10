#include <stdio.h>
#include <stdlib.h>
#include "struve.h"

int main(int argc, char** argv) {
    int NUM = argc>1 ? atoi(argv[1]) : 100;
    double MAX = argc>2 ? atof(argv[2]) : 2.0;
    for (int i=-NUM; i<=NUM; ++i) for (int j=-NUM; j<=NUM; ++j) {
        double Real = MAX/(double)NUM * (double)i;
        double Imag = MAX/(double)NUM * (double)j;
        _Complex double z = Real + 1.0I * Imag;
        _Complex double Res1;
        _Complex double Res2;
        struve_h0( &z, &Res1 );
        struve_y0( &z, &Res2 );
        printf( "%.7e %.7e %.7e %.7e %.7e %.7e\n", z, Res1, Res2 );
    }
}
