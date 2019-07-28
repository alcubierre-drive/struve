#pragma once

/* LICENSE: LGPLv3, Author: alcubierre-drive */

#ifdef __cplusplus
extern "C" {
#endif
typedef struct {
    double r;
    double i;
} _f2c_mycdoub;
int _f2c_struvefunctions( _f2c_mycdoub*, int*,
        _f2c_mycdoub*, _f2c_mycdoub*, _f2c_mycdoub*,
        _f2c_mycdoub*, _f2c_mycdoub*, _f2c_mycdoub* );
#ifdef __cplusplus
}

#include<array>
#include<complex>

namespace _f2c_struve {
    /* kode = 1: J0,J1,Y0,Y1
     *      = 2: J0,J1,H0,H1
     *      = 3: J0,J1,Y0,Y1,H0,H1
     */
    template<int kode>
    int cjyhbs( std::complex<double> z,
            std::array<std::complex<double>,6>& output ) {
        _f2c_mycdoub mz; mz.r = z.real(); mz.i = z.imag();
        _f2c_mycdoub moutput[6];
        int k = kode;
        int result = _f2c_struvefunctions( &mz, &k, &(moutput[0]),
                &(moutput[1]), &(moutput[2]), &(moutput[3]), &(moutput[4]),
                &(moutput[5]) );
        for (int i=0; i<6; ++i) {
            output[i] = moutput[i].r + std::complex<double>(0.0,1.0) *
                moutput[i].i;
        }

        /* result = { J_0(z), J_1(z), Y_0(z), Y_1(z), H_0(z), H_1(z) } */
        return result;
    }
    inline std::complex<double> h0( std::complex<double> z ) {
        std::array<std::complex<double>,6> results;
        cjyhbs<3>(z,results);
        return results[4];
    }
    inline std::complex<double> y0( std::complex<double> z ) {
        std::array<std::complex<double>,6> results;
        cjyhbs<3>(z,results);
        return results[2];
    }
}
#endif
