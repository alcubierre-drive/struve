#include <math.h>
#include <complex.h>
#include "struve.h"

typedef struct {
    double r;
    double i;
} _f2c_mycdoub;
static int _f2c_struvefunctions( _f2c_mycdoub*, int*,
        _f2c_mycdoub*, _f2c_mycdoub*, _f2c_mycdoub*,
        _f2c_mycdoub*, _f2c_mycdoub*, _f2c_mycdoub* );

/* translated to C with f2c, got rid of dependencies by manually writing the
   necessary functions for complex variables */
/* LICENSE: LGPLv3, Author: alcubierre-drive */

/*     WRITTEN BY D.E. AMOS AND S.L. DANIEL */

/*     REFERENCES */
/*         SLA-73-0262 */

/*         NBS HANDBOOK OF MATHEMATICAL FUNCTIONS, AMS 55, BY */
/*         M. ABRAMOWITZ AND I.A. STEGUN, DECEMBER, 1955, PP. 364, 497. */

/*     ABSTRACT */
/*         CJYHBS COMPUTES BESSEL FUNCTIONS J/SUB(NU)/(Z), Y/SUB(NU)/(Z), */
/*         AND STRUVE FUNCTIONS H/SUB(NU)/(Z), FOR DOUBLE COMPLEX Z AND NU=0 OR */
/*         1. BACKWARD RECURSION IS USED FOR THE J BESSEL FUNCTIONS OF */
/*         INTEGER ORDER TO SUM THE NEUMANN SERIES FOR THE Y AND H */
/*         FUNCTIONS FOR 0.LT.ABS(Z).LT.30. FOR CABS(Z).GT.30 THE */
/*         ASYMPTOTIC EXPANSIONS ARE USED. FOR Z, ABS(Z).GT.0. */
/*         AND -PI.LT.ARG(Z).LE.PI */

/*     DESCRIPTION OF ARGUMENTS */

/*         INPUT */
/*           KODE   - A PARAMETER TO SELECT THE PROPER FUNCTION PAIRS */
/*                    KODE=1 RETURNS J0,J1,Y0,Y1 FUNCTIONS */
/*                    KODE=2 RETURNS J0,J1,H0,H1 FUNCTIONS */
/*                    KODE=3 RETURNS J0,J1,Y0,Y1,H0,H1 FUNCTIONS */
/*           Z      - DOUBLE COMPLEX ARGUMENT, Z.NE.CMPLX(0.,0.) */
/*                    AND -PI.LT.ARG(Z).LE.PI */

/*         OUTPUT */
/*           CJ0    - BESSEL FUNCTION J/SUB(0)/(Z), A DOUBLE COMPLEX NUMBER */
/*           CJ1    - BESSEL FUNCTION J/SUB(1)/(Z), A DOUBLE COMPLEX NUMBER */
/*           CY0    - BESSEL FUNCTION Y/SUB(0)/(Z), A DOUBLE COMPLEX NUMBER */
/*           CY1    - BESSEL FUNCTION Y/SUB(1)/(Z), A DOUBLE COMPLEX NUMBER */
/*           CH0    - STRUVE FUNCTION H/SUB(0)/(Z), A DOUBLE COMPLEX NUMBER */
/*           CH1    - STRUVE FUNCTION H/SUB(1)/(Z), A DOUBLE COMPLEX NUMBER */

/*     ERROR CONDITIONS */
/*         ERROR #1, Z=0  ON INPUT, A FATAL ERROR */
/*         ERROR #2, KODE NOT 1 OR 2 OR 3 ON INPUT, A FATAL ERROR */


int struve_cjyhbs( int kode, const void* z, void* output ) {
    _f2c_mycdoub* _z = (_f2c_mycdoub*) z;
    _f2c_mycdoub* _output = (_f2c_mycdoub*) output;
    int k = kode;
    return _f2c_struvefunctions( _z, &k, &(_output[0]), &(_output[1]),
            &(_output[2]), &(_output[3]), &(_output[4]), &(_output[5]) );
}

int struve_h0( const void* z, void* output ) {
    _f2c_mycdoub* _z = (_f2c_mycdoub*) z;
    _f2c_mycdoub* _output = (_f2c_mycdoub*) output;
    _f2c_mycdoub cache[6];
    int k = 3;
    int result = _f2c_struvefunctions( _z, &k, &(cache[0]), &(cache[1]),
            &(cache[2]), &(cache[3]), &(cache[4]), &(cache[5]) );
    _output->r = cache[4].r;
    _output->i = cache[4].i;
    return result;
}

int struve_y0( const void* z, void* output ) {
    _f2c_mycdoub* _z = (_f2c_mycdoub*) z;
    _f2c_mycdoub* _output = (_f2c_mycdoub*) output;
    _f2c_mycdoub cache[6];
    int k = 3;
    int result = _f2c_struvefunctions( _z, &k, &(cache[0]), &(cache[1]),
            &(cache[2]), &(cache[3]), &(cache[4]), &(cache[5]) );
    _output->r = cache[2].r;
    _output->i = cache[2].i;
    return result;
}

static double z_abs( _f2c_mycdoub* z ) {
    return sqrt(z->r*z->r + z->i*z->i);
}
static double d_imag( _f2c_mycdoub* z ) {
    return z->i;
}
static void d_cnjg( _f2c_mycdoub* z_out, _f2c_mycdoub* z ) {
    z_out->r = z->r;
    z_out->i = -z->i;
}
static void z_prod( _f2c_mycdoub* z_out, _f2c_mycdoub* z1, _f2c_mycdoub* z2 ) {
    z_out->r = z1->r * z2->r - z1->i * z2->i;
    z_out->i = z1->r * z2->i + z1->i * z2->r;
}
static void z_div( _f2c_mycdoub* z_out, _f2c_mycdoub* z1, _f2c_mycdoub* z2 ) {
    _f2c_mycdoub z2c;
    d_cnjg( &z2c, z2 );
    double z2abs = z_abs(z2);
    _f2c_mycdoub z2inv;
    z2inv.r = z2c.r / (z2abs*z2abs);
    z2inv.i = z2c.i / (z2abs*z2abs);
    z_prod(z_out,z1,&z2inv);
}
static double dabs( double x ) { return fabs(x); }
static void z_cos( _f2c_mycdoub* z_out, _f2c_mycdoub* z ) {
    complex double zcd = z->r + I*z->i;
    complex double zcr = ccos(zcd);
    z_out->r = creal(zcr);
    z_out->i = cimag(zcr);
}
static void z_sin( _f2c_mycdoub* z_out, _f2c_mycdoub* z ) {
    complex double zcd = z->r + I*z->i;
    complex double zcr = csin(zcd);
    z_out->r = creal(zcr);
    z_out->i = cimag(zcr);
}

int _f2c_struvefunctions(_f2c_mycdoub *z__, int *kode, _f2c_mycdoub 
    *cj0, _f2c_mycdoub *cj1, _f2c_mycdoub *cy0, _f2c_mycdoub *cy1, 
    _f2c_mycdoub *ch0, _f2c_mycdoub *ch1) {
    /* Initialized data */

    static double euler = .577215664901533f;
    static double pio2 = 1.5707963267949f;
    static int iset = 1;

    /* System generated locals */
    int i__1, i__2, i__3, i__4;
    double r__1, r__2;
    _f2c_mycdoub q__1;
    _f2c_mycdoub z__1, z__2, z__3, z__4, z__5, z__6;

    /* Local variables */
    static int i__, j, k, l;
    static double x, y;
    static _f2c_mycdoub c1, c2, c3, c4, c5, c6;
    static int k0, k1;
    static _f2c_mycdoub z1, ck;
    static int kh;
    static double fn;
    static _f2c_mycdoub cr;
    static double pi, ay, ax;
    static _f2c_mycdoub sk;
    static double ey;
    static _f2c_mycdoub sn, tn;
    static int nu;
    static double sgn;
    static int nup1;
    static _f2c_mycdoub capj[71], cone, capr[71], fopi;
    static double zabs;
    static _f2c_mycdoub zlam[70], csum, topi, ctwo;
    static double trpi, theta;
    static _f2c_mycdoub czero;
    static double rtopi, theta1;

    zabs = z_abs(z__);
    if (zabs == 0.f) {
    goto L90;
    }
    if (*kode < 1 || *kode > 3) {
    goto L91;
    }
    switch (iset) {
    case 1:  goto L1;
    case 2:  goto L2;
    }
L1:
    iset = 2;

    c1.r = 0.f, c1.i = -2.f;
    c2.r = -2.f, c2.i = 0.f;
    c3.r = 0.f, c3.i = 2.f;
    c4.r = 2.f, c4.i = 0.f;
    k = 1;
    for (j = 1; j <= 17; ++j) {
    i__1 = k - 1;
    zlam[i__1].r = c1.r, zlam[i__1].i = c1.i;
    i__1 = k;
    zlam[i__1].r = c2.r, zlam[i__1].i = c2.i;
    i__1 = k + 1;
    zlam[i__1].r = c3.r, zlam[i__1].i = c3.i;
    i__1 = k + 2;
    zlam[i__1].r = c4.r, zlam[i__1].i = c4.i;
/* L3: */
    k += 4;
    }
    i__1 = k - 1;
    zlam[i__1].r = c1.r, zlam[i__1].i = c1.i;
    i__1 = k;
    zlam[i__1].r = c2.r, zlam[i__1].i = c2.i;
    pi = 3.14159265358979f;
    trpi = 2.f / pi;
    q__1.r = trpi, q__1.i = 0.f;
    topi.r = q__1.r, topi.i = q__1.i;
    r__1 = 4.f / pi;
    q__1.r = r__1, q__1.i = 0.f;
    fopi.r = q__1.r, fopi.i = q__1.i;
    rtopi = sqrt(trpi);
    ctwo.r = 2.f, ctwo.i = 0.f;
    cone.r = 1.f, cone.i = 0.f;
    czero.r = 0.f, czero.i = 0.f;
L2:
    x = z__->r;
    y = d_imag(z__);
    if (x < 0.f) {
    goto L4;
    } else if (x == 0) {
    goto L5;
    } else {
    goto L6;
    }
L4:
    if (y < 0.f) {
    goto L9;
    } else if (y == 0) {
    goto L8;
    } else {
    goto L7;
    }
L9:
    theta1 = atan(y / x);
    theta = theta1 - pi;
    sgn = -2.f;
    goto L50;
L5:
    if (y < 0.f) {
    goto L11;
    } else if (y == 0) {
    goto L90;
    } else {
    goto L13;
    }
L11:
    theta = -pio2;
    theta1 = theta;
    goto L50;
L13:
    theta = pio2;
    theta1 = theta;
    goto L50;
L7:
    theta1 = atan(y / x);
    theta = theta1 + pi;
    sgn = 2.f;
    goto L50;
L8:
    theta1 = 0.f;
    theta = pi;
    sgn = 2.f;
    goto L50;
L6:
    theta1 = atan(y / x);
    theta = theta1;
L50:

/*     BACKWARD RECURSION FOR J/SUBK/(X),K=0,1,... */

    if (zabs >= 30.f) {
    goto L200;
    }
    nu = (int) zabs + 40;
    nup1 = nu + 1;
    fn = (double) nu;
    ay = dabs(y);
    q__1.r = x, q__1.i = ay;
    z1.r = q__1.r, z1.i = q__1.i;
    sn.r = czero.r, sn.i = czero.i;
    i__1 = nup1 - 1;
    capr[i__1].r = sn.r, capr[i__1].i = sn.i;
    r__1 = fn + fn;
    q__1.r = r__1, q__1.i = 0.f;
    tn.r = q__1.r, tn.i = q__1.i;
    i__1 = nu;
    for (i__ = 1; i__ <= i__1; ++i__) {
    l = nup1 - i__;
    i__2 = l;
    z__3.r = z1.r * capr[i__2].r - z1.i * capr[i__2].i, z__3.i = z1.r * 
        capr[i__2].i + z1.i * capr[i__2].r;
    z__2.r = tn.r - z__3.r, z__2.i = tn.i - z__3.i;
    z_div(&z__1, &z1, &z__2);
    cr.r = z__1.r, cr.i = z__1.i;
    i__2 = l - 1;
    capr[i__2].r = cr.r, capr[i__2].i = cr.i;
    i__2 = l - 1;
    z__2.r = zlam[i__2].r + sn.r, z__2.i = zlam[i__2].i + sn.i;
    z__1.r = cr.r * z__2.r - cr.i * z__2.i, z__1.i = cr.r * z__2.i + cr.i 
        * z__2.r;
    sn.r = z__1.r, sn.i = z__1.i;
    z__1.r = tn.r - ctwo.r, z__1.i = tn.i - ctwo.i;
    tn.r = z__1.r, tn.i = z__1.i;
/* L10: */
    }
    z__1.r = sn.r + cone.r, z__1.i = sn.i + cone.i;
    sn.r = z__1.r, sn.i = z__1.i;
    ey = exp(ay);
    r__1 = ey * cos(x);
    r__2 = -ey * sin(x);
    q__1.r = r__1, q__1.i = r__2;
    z__2.r = q__1.r, z__2.i = q__1.i;
    z_div(&z__1, &z__2, &sn);
    capj[0].r = z__1.r, capj[0].i = z__1.i;
    i__1 = nu;
    for (i__ = 1; i__ <= i__1; ++i__) {
    i__2 = i__;
    i__3 = i__ - 1;
    i__4 = i__ - 1;
    z__1.r = capj[i__3].r * capr[i__4].r - capj[i__3].i * capr[i__4].i, 
        z__1.i = capj[i__3].r * capr[i__4].i + capj[i__3].i * capr[
        i__4].r;
    capj[i__2].r = z__1.r, capj[i__2].i = z__1.i;
/* L15: */
    }
    if (y >= 0.f) {
    goto L17;
    }
    i__1 = nup1;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L16: */
    i__2 = i__ - 1;
    d_cnjg(&z__1, &capj[i__ - 1]);
    capj[i__2].r = z__1.r, capj[i__2].i = z__1.i;
    }
L17:
    cj1->r = capj[1].r, cj1->i = capj[1].i;
    cj0->r = capj[0].r, cj0->i = capj[0].i;
    if (*kode == 2) {
    goto L100;
    }

/*     NEUMANN SERIES FOR Y0 */

    k0 = (nu - 1) / 2;
    k1 = (nu - 2) / 2;
    r__1 = log(zabs * .5f) + euler;
    q__1.r = r__1, q__1.i = theta;
    c1.r = q__1.r, c1.i = q__1.i;
    z__1.r = c1.r - cone.r, z__1.i = c1.i - cone.i;
    c2.r = z__1.r, c2.i = z__1.i;
    z__1.r = -cone.r, z__1.i = -cone.i;
    sk.r = z__1.r, sk.i = z__1.i;
    ck.r = cone.r, ck.i = cone.i;
    csum.r = czero.r, csum.i = czero.i;
    i__2 = k0;
    for (i__ = 1; i__ <= i__2; ++i__) {
    i__1 = i__ + i__;
    z__3.r = sk.r * capj[i__1].r - sk.i * capj[i__1].i, z__3.i = sk.r * 
        capj[i__1].i + sk.i * capj[i__1].r;
    z_div(&z__2, &z__3, &ck);
    z__1.r = csum.r + z__2.r, z__1.i = csum.i + z__2.i;
    csum.r = z__1.r, csum.i = z__1.i;
    z__1.r = -sk.r, z__1.i = -sk.i;
    sk.r = z__1.r, sk.i = z__1.i;
    z__1.r = ck.r + cone.r, z__1.i = ck.i + cone.i;
    ck.r = z__1.r, ck.i = z__1.i;
/* L20: */
    }
    z__3.r = c1.r * capj[0].r - c1.i * capj[0].i, z__3.i = c1.r * capj[0].i + 
        c1.i * capj[0].r;
    z__4.r = csum.r + csum.r, z__4.i = csum.i + csum.i;
    z__2.r = z__3.r - z__4.r, z__2.i = z__3.i - z__4.i;
    z__1.r = topi.r * z__2.r - topi.i * z__2.i, z__1.i = topi.r * z__2.i + 
        topi.i * z__2.r;
    cy0->r = z__1.r, cy0->i = z__1.i;

/*     NEUMANN SERIES FOR Y1 */

    c3.r = ctwo.r, c3.i = ctwo.i;
    z__1.r = ctwo.r + cone.r, z__1.i = ctwo.i + cone.i;
    c4.r = z__1.r, c4.i = z__1.i;
    c5.r = cone.r, c5.i = cone.i;
    z__1.r = -cone.r, z__1.i = -cone.i;
    sk.r = z__1.r, sk.i = z__1.i;
    csum.r = czero.r, csum.i = czero.i;
    i__2 = k1;
    for (i__ = 1; i__ <= i__2; ++i__) {
    z_div(&z__5, &c4, &c3);
    z__4.r = sk.r * z__5.r - sk.i * z__5.i, z__4.i = sk.r * z__5.i + sk.i 
        * z__5.r;
    i__1 = i__ + i__ + 1;
    z__3.r = z__4.r * capj[i__1].r - z__4.i * capj[i__1].i, z__3.i = 
        z__4.r * capj[i__1].i + z__4.i * capj[i__1].r;
    z_div(&z__2, &z__3, &c5);
    z__1.r = csum.r + z__2.r, z__1.i = csum.i + z__2.i;
    csum.r = z__1.r, csum.i = z__1.i;
    z__1.r = -sk.r, z__1.i = -sk.i;
    sk.r = z__1.r, sk.i = z__1.i;
    z__1.r = c3.r + cone.r, z__1.i = c3.i + cone.i;
    c3.r = z__1.r, c3.i = z__1.i;
    z__1.r = c4.r + ctwo.r, z__1.i = c4.i + ctwo.i;
    c4.r = z__1.r, c4.i = z__1.i;
    z__1.r = c5.r + cone.r, z__1.i = c5.i + cone.i;
    c5.r = z__1.r, c5.i = z__1.i;
/* L25: */
    }
    z__5.r = -capj[0].r, z__5.i = -capj[0].i;
    z_div(&z__4, &z__5, z__);
    z__6.r = c2.r * capj[1].r - c2.i * capj[1].i, z__6.i = c2.r * capj[1].i + 
        c2.i * capj[1].r;
    z__3.r = z__4.r + z__6.r, z__3.i = z__4.i + z__6.i;
    z__2.r = z__3.r - csum.r, z__2.i = z__3.i - csum.i;
    z__1.r = topi.r * z__2.r - topi.i * z__2.i, z__1.i = topi.r * z__2.i + 
        topi.i * z__2.r;
    cy1->r = z__1.r, cy1->i = z__1.i;
    if (*kode == 1) {
    return 0;
    }

/*     NEUMANN SERIES FOR H0 AND H1 */

L100:
    kh = nu / 2;
    z__1.r = cone.r + ctwo.r, z__1.i = cone.i + ctwo.i;
    ck.r = z__1.r, ck.i = z__1.i;
    c1.r = cone.r, c1.i = cone.i;
    c4.r = czero.r, c4.i = czero.i;
    c3.r = c4.r, c3.i = c4.i;
    l = 2;
    i__2 = kh;
    for (k = 1; k <= i__2; ++k) {
    z_div(&z__2, &capj[l - 1], &c1);
    z__1.r = c3.r + z__2.r, z__1.i = c3.i + z__2.i;
    c3.r = z__1.r, c3.i = z__1.i;
    z__3.r = c1.r * ck.r - c1.i * ck.i, z__3.i = c1.r * ck.i + c1.i * 
        ck.r;
    z_div(&z__2, &capj[l], &z__3);
    z__1.r = c4.r + z__2.r, z__1.i = c4.i + z__2.i;
    c4.r = z__1.r, c4.i = z__1.i;
    c1.r = ck.r, c1.i = ck.i;
    z__1.r = ck.r + ctwo.r, z__1.i = ck.i + ctwo.i;
    ck.r = z__1.r, ck.i = z__1.i;
    l += 2;
/* L30: */
    }
    z__1.r = fopi.r * c3.r - fopi.i * c3.i, z__1.i = fopi.r * c3.i + fopi.i * 
        c3.r;
    ch0->r = z__1.r, ch0->i = z__1.i;
    z__4.r = cone.r - capj[0].r, z__4.i = cone.i - capj[0].i;
    z__3.r = z__4.r + c4.r, z__3.i = z__4.i + c4.i;
    z__2.r = z__3.r + c4.r, z__2.i = z__3.i + c4.i;
    z__1.r = topi.r * z__2.r - topi.i * z__2.i, z__1.i = topi.r * z__2.i + 
        topi.i * z__2.r;
    ch1->r = z__1.r, ch1->i = z__1.i;
    return 0;

L200:

/*     ASYMPTOTIC EXPANSIONS FOR Y0,Y1,J0,J1,ABS(Z).GE.30 */

    theta1 *= .5f;
    r__1 = cos(theta1);
    r__2 = -sin(theta1);
    q__1.r = r__1, q__1.i = r__2;
    cr.r = q__1.r, cr.i = q__1.i;
    r__1 = rtopi / sqrt(zabs);
    z__1.r = r__1 * cr.r, z__1.i = r__1 * cr.i;
    cr.r = z__1.r, cr.i = z__1.i;
/*     FOR ABS(Z).GT.30 */
    if (x >= 0.f) {
    goto L207;
    }
    y = -y;
L207:
    ax = dabs(x);
    q__1.r = ax, q__1.i = y;
    z1.r = q__1.r, z1.i = q__1.i;
    for (l = 1; l <= 2; ++l) {
    nu = l - 1;
    i__2 = (nu << 2) * nu;
    tn.r = (double) i__2, tn.i = 0.;
    r__1 = ax - ((double) nu * .5f + .25f) * pi;
    q__1.r = r__1, q__1.i = y;
    sn.r = q__1.r, sn.i = q__1.i;
    c2.r = czero.r, c2.i = czero.i;
    ck.r = cone.r, ck.i = cone.i;
    c1.r = cone.r, c1.i = cone.i;
    sk.r = cone.r, sk.i = cone.i;
    z__1.r = z1.r * 8.f, z__1.i = z1.i * 8.f;
    c4.r = z__1.r, c4.i = z__1.i;
    c3.r = c4.r, c3.i = c4.i;
    for (k = 1; k <= 16; ++k) {
        z__4.r = sk.r * sk.r - sk.i * sk.i, z__4.i = sk.r * sk.i + sk.i * 
            sk.r;
        z__3.r = tn.r - z__4.r, z__3.i = tn.i - z__4.i;
        z_div(&z__2, &z__3, &c3);
        z__1.r = ck.r * z__2.r - ck.i * z__2.i, z__1.i = ck.r * z__2.i + 
            ck.i * z__2.r;
        ck.r = z__1.r, ck.i = z__1.i;
        if (k % 2 == 0) {
        goto L42;
        }
        z__1.r = c2.r + ck.r, z__1.i = c2.i + ck.i;
        c2.r = z__1.r, c2.i = z__1.i;
        z__1.r = -ck.r, z__1.i = -ck.i;
        ck.r = z__1.r, ck.i = z__1.i;
        goto L44;
L42:
        z__1.r = c1.r + ck.r, z__1.i = c1.i + ck.i;
        c1.r = z__1.r, c1.i = z__1.i;
        goto L44;
L44:
        if (z_abs(&ck) < 1e-13f) {
        goto L45;
        }
        z__1.r = c3.r + c4.r, z__1.i = c3.i + c4.i;
        c3.r = z__1.r, c3.i = z__1.i;
        z__1.r = sk.r + ctwo.r, z__1.i = sk.i + ctwo.i;
        sk.r = z__1.r, sk.i = z__1.i;
/* L40: */
    }
L45:
    z_cos(&z__1, &sn);
    c5.r = z__1.r, c5.i = z__1.i;
    z_sin(&z__1, &sn);
    c6.r = z__1.r, c6.i = z__1.i;
    switch (l) {
        case 1:  goto L201;
        case 2:  goto L202;
    }
L201:
    z__3.r = c1.r * c6.r - c1.i * c6.i, z__3.i = c1.r * c6.i + c1.i * 
        c6.r;
    z__4.r = c2.r * c5.r - c2.i * c5.i, z__4.i = c2.r * c5.i + c2.i * 
        c5.r;
    z__2.r = z__3.r + z__4.r, z__2.i = z__3.i + z__4.i;
    z__1.r = cr.r * z__2.r - cr.i * z__2.i, z__1.i = cr.r * z__2.i + cr.i 
        * z__2.r;
    cy0->r = z__1.r, cy0->i = z__1.i;
    z__3.r = c1.r * c5.r - c1.i * c5.i, z__3.i = c1.r * c5.i + c1.i * 
        c5.r;
    z__4.r = c2.r * c6.r - c2.i * c6.i, z__4.i = c2.r * c6.i + c2.i * 
        c6.r;
    z__2.r = z__3.r - z__4.r, z__2.i = z__3.i - z__4.i;
    z__1.r = cr.r * z__2.r - cr.i * z__2.i, z__1.i = cr.r * z__2.i + cr.i 
        * z__2.r;
    cj0->r = z__1.r, cj0->i = z__1.i;
    goto L205;
L202:
    z__3.r = c1.r * c6.r - c1.i * c6.i, z__3.i = c1.r * c6.i + c1.i * 
        c6.r;
    z__4.r = c2.r * c5.r - c2.i * c5.i, z__4.i = c2.r * c5.i + c2.i * 
        c5.r;
    z__2.r = z__3.r + z__4.r, z__2.i = z__3.i + z__4.i;
    z__1.r = cr.r * z__2.r - cr.i * z__2.i, z__1.i = cr.r * z__2.i + cr.i 
        * z__2.r;
    cy1->r = z__1.r, cy1->i = z__1.i;
    z__3.r = c1.r * c5.r - c1.i * c5.i, z__3.i = c1.r * c5.i + c1.i * 
        c5.r;
    z__4.r = c2.r * c6.r - c2.i * c6.i, z__4.i = c2.r * c6.i + c2.i * 
        c6.r;
    z__2.r = z__3.r - z__4.r, z__2.i = z__3.i - z__4.i;
    z__1.r = cr.r * z__2.r - cr.i * z__2.i, z__1.i = cr.r * z__2.i + cr.i 
        * z__2.r;
    cj1->r = z__1.r, cj1->i = z__1.i;
L205:
    ;
    }
    if (x >= 0.f) {
    goto L206;
    }
/*     FORMULA FOR Z=-Z1, -PI/2..LT.ARG(Z1).LT.PI/2. */
    q__1.r = 0.f, q__1.i = sgn;
    z__2.r = q__1.r * cj0->r - q__1.i * cj0->i, z__2.i = q__1.r * cj0->i + 
        q__1.i * cj0->r;
    z__1.r = cy0->r + z__2.r, z__1.i = cy0->i + z__2.i;
    cy0->r = z__1.r, cy0->i = z__1.i;
    z__2.r = -cy1->r, z__2.i = -cy1->i;
    q__1.r = 0.f, q__1.i = sgn;
    z__3.r = q__1.r * cj1->r - q__1.i * cj1->i, z__3.i = q__1.r * cj1->i + 
        q__1.i * cj1->r;
    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
    cy1->r = z__1.r, cy1->i = z__1.i;
    z__1.r = -cj1->r, z__1.i = -cj1->i;
    cj1->r = z__1.r, cj1->i = z__1.i;
L206:
    if (*kode == 1) {
    return 0;
    }

/*     ASYMPTOTIC EXPANSIONS FOR H0,H1,ABS(Z).GE.30 */

    z_div(&z__1, &cone, z__);
    ck.r = z__1.r, ck.i = z__1.i;
    csum.r = ck.r, csum.i = ck.i;
    sk.r = cone.r, sk.i = cone.i;
    z__1.r = ck.r * ck.r - ck.i * ck.i, z__1.i = ck.r * ck.i + ck.i * ck.r;
    cr.r = z__1.r, cr.i = z__1.i;
    for (k = 1; k <= 8; ++k) {
    z__4.r = -ck.r, z__4.i = -ck.i;
    z__3.r = z__4.r * sk.r - z__4.i * sk.i, z__3.i = z__4.r * sk.i + 
        z__4.i * sk.r;
    z__2.r = z__3.r * sk.r - z__3.i * sk.i, z__2.i = z__3.r * sk.i + 
        z__3.i * sk.r;
    z__1.r = z__2.r * cr.r - z__2.i * cr.i, z__1.i = z__2.r * cr.i + 
        z__2.i * cr.r;
    ck.r = z__1.r, ck.i = z__1.i;
    z__1.r = csum.r + ck.r, z__1.i = csum.i + ck.i;
    csum.r = z__1.r, csum.i = z__1.i;
    if (z_abs(&ck) < 1e-14f) {
        goto L62;
    }
    z__1.r = sk.r + ctwo.r, z__1.i = sk.i + ctwo.i;
    sk.r = z__1.r, sk.i = z__1.i;
/* L60: */
    }
L62:
    z__2.r = topi.r * csum.r - topi.i * csum.i, z__2.i = topi.r * csum.i + 
        topi.i * csum.r;
    z__1.r = cy0->r + z__2.r, z__1.i = cy0->i + z__2.i;
    ch0->r = z__1.r, ch0->i = z__1.i;
    ck.r = cr.r, ck.i = cr.i;
    csum.r = cr.r, csum.i = cr.i;
    c1.r = cone.r, c1.i = cone.i;
    z__1.r = c1.r + ctwo.r, z__1.i = c1.i + ctwo.i;
    c2.r = z__1.r, c2.i = z__1.i;
    for (k = 1; k <= 8; ++k) {
    z__2.r = -ck.r, z__2.i = -ck.i;
    z__4.r = c1.r * c2.r - c1.i * c2.i, z__4.i = c1.r * c2.i + c1.i * 
        c2.r;
    z__3.r = z__4.r * cr.r - z__4.i * cr.i, z__3.i = z__4.r * cr.i + 
        z__4.i * cr.r;
    z__1.r = z__2.r * z__3.r - z__2.i * z__3.i, z__1.i = z__2.r * z__3.i 
        + z__2.i * z__3.r;
    ck.r = z__1.r, ck.i = z__1.i;
    z__1.r = csum.r + ck.r, z__1.i = csum.i + ck.i;
    csum.r = z__1.r, csum.i = z__1.i;
    if (z_abs(&ck) < 1e-14f) {
        goto L66;
    }
    c1.r = c2.r, c1.i = c2.i;
    z__1.r = c2.r + ctwo.r, z__1.i = c2.i + ctwo.i;
    c2.r = z__1.r, c2.i = z__1.i;
/* L64: */
    }
L66:
    z__3.r = cone.r + csum.r, z__3.i = cone.i + csum.i;
    z__2.r = topi.r * z__3.r - topi.i * z__3.i, z__2.i = topi.r * z__3.i + 
        topi.i * z__3.r;
    z__1.r = cy1->r + z__2.r, z__1.i = cy1->i + z__2.i;
    ch1->r = z__1.r, ch1->i = z__1.i;
    return 0;
L90:
    return 0;
L91:
    return 0;
}

