#pragma once

/* LICENSE: LGPLv3, Author: alcubierre-drive */

#ifdef __cplusplus
extern "C" {
#endif

// z is double[2], i.e. a complex number

// ouput is double[12], i.e. six complex numbers
/*
 * kode = 1: J0,J1,Y0,Y1
 *      = 2: J0,J1,H0,H1
 *      = 3: J0,J1,Y0,Y1,H0,H1
 */
int struve_cjyhbs( int kode, const void* z, void* output );

// output is double[2], i.e. a complex number
int struve_h0( const void* z, void* output );
int struve_y0( const void* z, void* output );

#ifdef __cplusplus
}
#endif
