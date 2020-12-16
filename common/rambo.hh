/*
 * rambo.hpp
 *
 *  Created on: Sep 15, 2016
 *      Author: akozlins
 */

#ifndef MU3ESIM_UTIL_RAMBO_H_
#define MU3ESIM_UTIL_RAMBO_H_

#include "rand.hh"

/**
 * RAMBO -  RA(NDOM) M(OMENTA) BO(OSTER)
 *
 * A DEMOCRATIC MULTI-PARTICLE PHASE SPACE GENERATOR
 * AUTHORS:  S.D. ELLIS,  R. KLEISS,  W.J. STIRLING
 *
 * N  = NUMBER OF PARTICLES (>1, IN THIS VERSION <101)
 * ET = TOTAL CENTRE-OF-MASS ENERGY
 * XM = PARTICLE MASSES ( DIM=N )
 * P  = PARTICLE MOMENTA ( DIM=(4,N) )
 * WT = WEIGHT OF THE EVENT
 * LW = FLAG FOR EVENT WEIGHTING:
 *      LW = 0 WEIGHTED EVENTS
 *      LW = 1 UNWEIGHTED EVENTS ( FLAT PHASE SPACE )
 */
inline
int rambo(const int n, const double et, double* xm, double p[][4], double* wt, const int lw) {
    static const double po2log = std::log(M_PI_2);

    /* Initialized data */

    static double acc = 1e-14;
    static int itmax = 6;
    static int ibegin = 0;
    static int iwarn[5] = { 0,0,0,0,0 };


    /* System generated locals */
    double d__1, d__2, d__3, d__4, d__5;

    /* Builtin functions */
    //double atan(), log();
    //int s_wsfe(), do_fio(), e_wsfe();
    /* Subroutine */ //int s_stop();
    //double sqrt(), cos(), sin(), exp(), pow_di(), pow_dd();

    /* Local variables */
    static double accu, rmas;
    static int iter;
    static double xmax, a, b[3], e[100], g;
    static double q[100][4], r__[4], v[100], x, w,
    z__[100], f0, g0, wtmax, p2[100], x2;
    static int nm;
    static double pm2, sm2, xm2[100], wt2, wt3, wtm, xmt;

    /* Function Body */
    /*                                                                       AAFU0024 */
    /* INITIALIZATION STEP: FACTORIALS FOR THE PHASE SPACE WEIGHT            AAFU0025 */
    if (ibegin != 0) {
        goto L103;
    }
    ibegin = 1;
    z__[1] = po2log;
    for (int k = 2; k < 100; ++k) {
        /* L101: */
        z__[k] = z__[k - 1] + po2log - 2 * std::log(double(k - 1));
    }
    for (int k = 2; k < 100; ++k) {
        /* L102: */
        z__[k] -= std::log(double(k));
    }
    /*                                                                       AAFU0035 */
    /* CHECK ON THE NUMBER OF PARTICLES                                      AAFU0036 */
    L103:
    if (n > 1 && n < 101) {
        goto L104;
    }
    //s_wsfe(&io___9);
    //do_fio(&c__1, (char *)&(*n), (ftnlen)sizeof(int));
    //e_wsfe();
    //s_stop("", (ftnlen)0);
    printf("RAMBO FAILS: # OF PARTICLES = %i, IS NOT ALLOWED\n", n);
    return -1;
    /*                                                                       AAFU0040 */
    /* CHECK WHETHER TOTAL ENERGY IS SUFFICIENT; COUNT NONZERO MASSES        AAFU0041 */
    L104:
    xmt = 0;
    nm = 0;
    for (int i = 0; i < n; ++i) {
        assert(xm[i] >= 0);
        if (xm[i] > 0) {
            ++nm;
        }
        /* L105: */
        //xmt += (d__1 = xm[i], abs(d__1));
        xmt += xm[i];
        //printf("xmt %f\n", xmt);
    }
    if (xmt <= et) {
        goto L106;
    }
    //s_wsfe(&io___13);
    //do_fio(&c__1, (char *)&xmt, (ftnlen)sizeof(double));
    //do_fio(&c__1, (char *)&(*et), (ftnlen)sizeof(double));
    //e_wsfe();
    //s_stop("", (ftnlen)0);
    printf("RAMBO FAILS: TOTAL MASS = %f IS NOT SMALLER THAN TOTAL ENERGY = %f\n", xmt, et);
    return -1;
    /*                                                                       AAFU0050 */
    /* CHECK ON THE WEIGHTING OPTION                                         AAFU0051 */
    L106:
    if (lw == 1 || lw == 0) {
        goto L201;
    }
    //  s_wsfe(&io___14);
    // do_fio(&c__1, (char *)&(*lw), (ftnlen)sizeof(int));
    // e_wsfe();
    // s_stop("", (ftnlen)0);
    printf("RAMBO FAILS: LW= %i IS NOT AN ALLOWED OPTION\n", lw);
    return -1;
    /*                                                                       AAFU0055 */
    /* THE PARAMETER VALUES ARE NOW ACCEPTED                                 AAFU0056 */
    /*                                                                       AAFU0057 */
    /* GENERATE N MASSLESS MOMENTA IN INFINITE PHASE SPACE                   AAFU0058 */
    L201:
    //printf("N: %i\n", *n);
    for (int i = 0; i < n; ++i) {
        auto u = mu3e::util::rand_u3d();
        q[i][3] = -std::log(CLHEP::RandFlat::shoot() * CLHEP::RandFlat::shoot());
        q[i][2] = q[i][3] * u.z;
        q[i][1] = q[i][3] * u.y;
        /* L202: */
        q[i][0] = q[i][3] * u.x;
    }
    /*                                                                       AAFU0067 */
    /* CALCULATE THE PARAMETERS OF THE CONFORMAL TRANSFORMATION              AAFU0068 */
    for (int i = 0; i < 4; ++i) {
        /* L203: */
        r__[i] = 0;
    }
    for (int i = 0; i < n; ++i) {
        for (int k = 0; k < 4; ++k) {
            /* L204: */
            r__[k] += q[i][k];
        }
    }
    /* Computing 2nd power */
    d__1 = r__[3];
    /* Computing 2nd power */
    d__2 = r__[2];
    /* Computing 2nd power */
    d__3 = r__[1];
    /* Computing 2nd power */
    d__4 = r__[0];
    rmas = sqrt(d__1 * d__1 - d__2 * d__2 - d__3 * d__3 - d__4 * d__4);
    for (int k = 0; k < 3; ++k) {
        /* L205: */
        b[k] = -r__[k] / rmas;
    }
    g = r__[3] / rmas;
    a = 1 / (g + 1);
    x = et / rmas;

    //printf("%f, %f, %f \n",g,a,x);
    /*                                                                       AAFU0080 */
    /* TRANSFORM THE Q'S CONFORMALLY INTO THE P'S                            AAFU0081 */
    for (int i = 0; i < n; ++i) {
        double bq = b[0] * q[i][0] + b[1] * q[i][1] + b[2] * q[i][2];
        for (int k = 0; k < 3; ++k) {
            /* L206: */
            p[i][k] = x * (q[i][k] + b[k] * (q[i][3] + a * bq));
        }
        /* L207: */
        p[i][3] = x * (g * q[i][3] + bq);
        //printf("[%i][3]: %f\n", i, p[i][3]]);
    }
    /*                                                                       AAFU0087 */
    /* RETURN FOR UNWEIGHTED MASSLESS MOMENTA                                AAFU0088 */
    *wt = 1;
    if (nm == 0 && lw == 1) {
        return 0;
    }
    /*                                                                       AAFU0091 */
    /* CALCULATE WEIGHT AND POSSIBLE WARNINGS                                AAFU0092 */
    *wt = po2log;
    if (n != 2) {
        *wt = (2 * n - 4) * std::log(et) + z__[n - 1];
    }
    if (*wt >= -180) {
        goto L208;
    }
    if (iwarn[0] <= 5) {
        //s_wsfe(&io___26);
        //  do_fio(&c__1, (char *)&(*wt), (ftnlen)sizeof(double));
        //  e_wsfe();
        printf("RAMBO WARNS: WEIGHT = EXP(%f) MAY UNDERFLOW\n", *wt);
    }
    ++iwarn[0];
    L208:
    if (*wt <= 174) {
        goto L209;
    }
    if (iwarn[1] <= 5) {
        //s_wsfe(&io___27);
        //  do_fio(&c__1, (char *)&(*wt), (ftnlen)sizeof(double));
        //  e_wsfe();
        printf("RAMBO WARNS: WEIGHT = EXP(%f) MAY OVERFLOW\n", *wt);
    }
    ++iwarn[1];
    /*                                                                       AAFU0101 */
    /* RETURN FOR WEIGHTED MASSLESS MOMENTA                                  AAFU0102 */
    L209:
    if (nm != 0) {
        goto L210;
    }
    *wt = exp(*wt);
    return 0;
    /*                                                                       AAFU0106 */
    /* MASSIVE PARTICLES: RESCALE THE MOMENTA BY A FACTOR X                  AAFU0107 */
    L210:
    /* Computing 2nd power */
    d__1 = xmt / et;
    xmax = sqrt(1 - d__1 * d__1);
    for (int i = 0; i < n; ++i) {
        /* Computing 2nd power */
        d__1 = xm[i];
        xm2[i] = d__1 * d__1;
        /* L301: */
        /* Computing 2nd power */
        d__1 = p[i][3];
        p2[i] = d__1 * d__1;
    }
    iter = 0;
    x = xmax;
    accu = et * acc;
    L302:
    f0 = -et;
    g0 = 0;
    x2 = x * x;
    for (int i = 0; i < n; ++i) {
        e[i] = sqrt(xm2[i] + x2 * p2[i]);
        f0 += e[i];
        /* L303: */
        g0 += p2[i] / e[i];
    }
    if (std::abs(f0) <= accu) {
        goto L305;
    }
    ++iter;
    if (iter <= itmax) {
        goto L304;
    }
    //s_wsfe(&io___37);
    //do_fio(&c__1, (char *)&itmax, (ftnlen)sizeof(int));
    //e_wsfe();
    printf("RAMBO WARNS: %i ITERATIONS DID NOT GIVE THE DESIRED ACCURACY = %f\n", itmax, accu);
    goto L305;
    L304:
    x -= f0 / (x * g0);
    goto L302;
    L305:
    for (int i = 0; i < n; ++i) {
        v[i] = x * p[i][3];
        for (int k = 0; k < 3; ++k) {
            /* L306: */
            p[i][k] = x * p[i][k];
        }
        /* L307: */
        p[i][3] = e[i];
    }
    /*                                                                       AAFU0134 */
    /* CALCULATE THE MASS-EFFECT WEIGHT FACTOR                               AAFU0135 */
    wt2 = 1;
    wt3 = 0;
    for (int i = 0; i < n; ++i) {
        wt2 = wt2 * v[i] / e[i];
        /* L308: */
        /* Computing 2nd power */
        d__1 = v[i];
        wt3 += d__1 * d__1 / e[i];
    }
    wtm = (2 * n - 3) * std::log(x) + std::log(wt2 / wt3 * et);
    if (lw == 1) {
        goto L401;
    }
    /*                                                                       AAFU0143 */
    /* RETURN FOR  WEIGHTED MASSIVE MOMENTA                                  AAFU0144 */
    *wt += wtm;
    //printf("%f, %f\n", wtm, *wt);
    if (*wt >= -180) {
        goto L309;
    }
    if (iwarn[2] <= 5) {
        //s_wsfe(&io___42);
        //do_fio(&c__1, (char *)&(*wt), (ftnlen)sizeof(double));
        //e_wsfe();
        printf("RAMBO WARNS: MASSIVE WEIGHT = EXP(%f) MAY UNDERFLOW\n",*wt);
    }
    ++iwarn[2];
    L309:
    if (*wt <= 174) {
        goto L310;
    }
    if (iwarn[3] <= 5) {
        //s_wsfe(&io___43);
        //do_fio(&c__1, (char *)&(*wt), (ftnlen)sizeof(double));
        //e_wsfe();
        printf("RAMBO WARNS: MASSIVE WEIGHT = EXP(%f) MAY OVERFLOW\n",*wt);
    }
    ++iwarn[3];
    L310:
    *wt = exp(*wt);
    return 0;
    /*                                                                       AAFU0154 */
    /* UNWEIGHTED MASSIVE MOMENTA REQUIRED: ESTIMATE MAXIMUM WEIGHT          AAFU0155 */
    L401:
    *wt = exp(wtm);
    if (nm > 1) {
        goto L402;
    }
    /*                                                                       AAFU0158 */
    /* ONE MASSIVE PARTICLE                                                  AAFU0159 */
    wtmax = mu3e::util::pown(xmax * xmax, 2 * n - 3);
    goto L405;
    L402:
    if (nm > 2) {
        goto L404;
    }
    /*                                                                       AAFU0163 */
    /* TWO MASSIVE PARTICLES                                                 AAFU0164 */
    sm2 = 0;
    pm2 = 0;
    for (int i = 0; i < n; ++i) {
        if (xm[i] > 0) {
            sm2 += xm2[i];
            pm2 *= xm2[i];
        }
    }
    /* Computing 2nd power */
    d__3 = et;
    /* Computing 2nd power */
    d__2 = 1 - sm2 / (d__3 * d__3);
    /* Computing 4th power */
    d__4 = et, d__4 *= d__4;
    d__1 = d__2 * d__2 - pm2 * 4 / (d__4 * d__4);
    d__5 = n - 1.5;
    wtmax = pow(d__1, d__5);
    goto L405;
    /*                                                                       AAFU0174 */
    /* MORE THAN TWO MASSIVE PARTICLES: AN ESTIMATE ONLY                     AAFU0175 */
    L404:
    wtmax = mu3e::util::pown(xmax, 2 * n - 3 + nm - 2);
    /*                                                                       AAFU0177 */
    /* DETERMINE WHETHER OR NOT TO ACCEPT THIS EVENT                         AAFU0178 */
    L405:
    w = *wt / wtmax;
    if (w <= 1) {
        goto L406;
    }
    if (iwarn[4] <= 5) {
        //s_wsfe(&io___48);
        //do_fio(&c__1, (char *)&wtmax, (ftnlen)sizeof(double));
        //do_fio(&c__1, (char *)&w, (ftnlen)sizeof(double));
        //e_wsfe();
        printf("RAMBO WARNS: ESTIMATE FOR MAXIMUM WEIGHT %f EXCEEDED BY A FACTOR %f\n",wtmax,w);
    }
    ++iwarn[4];
    L406:
    if (w < CLHEP::RandFlat::shoot()) {
        goto L201;
    }
    *wt = 1;
    return 0;
} // rambo

#endif /* MU3ESIM_UTIL_RAMBO_H_ */
