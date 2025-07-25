//****************************************************************************//
//                                                                            //
//  Header file collier.h                                                     //
//  This file is part of Package-X and CollierLink by Hiren H. Patel, and     //
//  provides a C and C++ interface to the fortran library COLLIER by          //
//  A. Denner, S. Dittmaier, L. Hofer.                                        //
//                                                                            //
//  This file is to be included with (and ONLY with) functions generated by   //
//    CollierCodeGenerate[] with option setting "Language" -> "C" or "C++",   //
//    and to be dynamically linked with libcollier.so .                       //
//                                                                            //
//  Langauge: C / C++                       Creation date: 8/01/2017          //
//  Copyright (C) 2017 by Hiren H. Patel                                      //
//                                                                            //
//  Package-X is distributed under the standard CPC license agreement         //
//  CollierLink is distributed under the CC BY-NC-ND 3.0 license              //
//  COLLIER is licensed under the GNU GPL version 3, see COPYING for details  //
//                                                                            //
//****************************************************************************//

#ifndef COLLIER_H_
#define COLLIER_H_

#ifdef __cplusplus
#include <complex>
#include <cmath>
#define dcplx std::complex<double>
#else
#include <complex.h>
#include <math.h>
#define dcplx double _Complex
#endif

//COLLIER subroutines
#ifdef __cplusplus
extern "C" {
#endif
    
    //Initializiation
    void __collier_init_MOD_init_cll(int*, int*, char*);
    void __collier_cache_MOD_initevent_cll();
    
    //Calculation of tensor integrals (section 5.3)
    void __collier_coefs_MOD_a_cll(dcplx[], dcplx[], dcplx*, int*, double[]
                                   ,int*);
    void __collier_coefs_MOD_b_main_cll(dcplx[], dcplx[], dcplx*, dcplx*, dcplx*
                                        , int*, double[], int*);
    void __collier_coefs_MOD_c_main_cll(dcplx[], dcplx[], dcplx*, dcplx*, dcplx*
                                        , dcplx*, dcplx*, dcplx*, int*, double[], int*, double[]);
    void __collier_coefs_MOD_d_main_cll(dcplx[], dcplx[], dcplx*, dcplx*, dcplx*
                                        , dcplx*, dcplx*, dcplx*, dcplx*, dcplx*, dcplx*, dcplx*, int*, double[]
                                        , int*, double[]);
    void __collier_coefs_MOD_e_main_cll(dcplx[], dcplx[], dcplx*, dcplx*, dcplx*
                                        , dcplx*, dcplx*, dcplx*, dcplx*, dcplx*
                                        , dcplx*, dcplx*, dcplx*, dcplx*, dcplx*, dcplx*, dcplx*, int*, double[]
                                        , int*, double[]);
    void __collier_coefs_MOD_f_main_cll(dcplx[], dcplx[], dcplx*, dcplx*, dcplx*
                                        , dcplx*, dcplx*, dcplx*, dcplx*, dcplx*, dcplx*, dcplx*, dcplx*, dcplx*
                                        , dcplx*, dcplx*, dcplx*, dcplx*, dcplx*, dcplx*, dcplx*, dcplx*, dcplx*
                                        , int*, double[], int*, double[]);
    
    
    void __collier_coefs_MOD_t1_checked_cll(dcplx[], dcplx[], dcplx[], const int*
                                    , const int*, double[], int*);
    void __collier_coefs_MOD_tn_main_checked_cll(dcplx[], dcplx[], dcplx[], dcplx[]
                                         , const int*, const int*, double[], int*, double[]);
    
    void __collier_coefs_MOD_db_arrays_cll(dcplx[], dcplx[], dcplx[], dcplx[]
                                           , const int*, double[]);
    
    //Regularization parameters (section 5.4.1)
    void __collier_init_MOD_setdeltauv_cll(double*);
    void __collier_init_MOD_setmuuv2_cll(double*);
    void __collier_init_MOD_getdeltauv_cll(double*);
    void __collier_init_MOD_getmuuv2_cll(double*);
    
    void __collier_init_MOD_setdeltair_cll(double*, double*);
    void __collier_init_MOD_setmuir2_cll(double*);
    void __collier_init_MOD_getdeltair_cll(double*,double*);
    void __collier_init_MOD_getmuir2_cll(double*);
    
    void __collier_init_MOD_setminf2_cll(int*, dcplx*);
    void __collier_init_MOD_addminf2_cll(dcplx*);
    void __collier_init_MOD_getnminf_cll(int*);
    void __collier_init_MOD_getminf2_cll(dcplx[]);
    void __collier_init_MOD_clearminf2_cll();
    
    //Technical parameters (section 5.4.2)
    void __collier_init_MOD_setmode_cll(int*);
    void __collier_init_MOD_getmode_cll(int*);
    
    void __collier_init_MOD_setreqacc_cll(double*);
    void __collier_init_MOD_getreqacc_cll(double*);
    
    
    void __collier_init_MOD_setcritacc_cll(double*);
    void __collier_init_MOD_getcritacc_cll(double*);
    
    void __collier_init_MOD_setcheckacc_cll(double*);
    void __collier_init_MOD_getcheckacc_cll(double*);
    
    void __collier_init_MOD_setaccuracy_cll(double*, double*, double*);
    
    void __collier_init_MOD_setritmax_cll(int*);
    void __collier_init_MOD_getritmax_cll(int*);
    
    
    //Cache system (section 5.5)
    void __cache_MOD_initcachesystem_cll(int*,int*);
    void __cache_MOD_addnewcache_cll(int*,int*);
    void __cache_MOD_setcachelevel_cll(int*,int*);
    void __cache_MOD_switchoffcachesystem_cll();
    void __cache_MOD_switchoncachesystem_cll();
    void __cache_MOD_switchoffcache_cll(int*);
    void __cache_MOD_switchoncache_cll(int*);
    
    //Error treatment and output files (section 5.6)
    void __collier_init_MOD_geterrflag_cll(int*);
    void __collier_init_MOD_seterrstop_cll(int*);
    void __collier_init_MOD_geterrstop_cll(int*);
    void __collier_init_MOD_switchofferrstop_cll();
    void __collier_init_MOD_getaccflag_cll(int*);
    void __collier_init_MOD_initaccflag_cll();
    
#ifdef __cplusplus
}

//C++ function for rapid evaluation of squares
static dcplx square(dcplx z) {return z*z;}

#else

//C (only) function for square magnitude
static double cnorm(dcplx z) {
    const double r = creal(z);
    const double i = cimag(z);return r*r + i*i;}
//C function for rapid evaluation of squares
static dcplx csquare(dcplx z) {return z*z;}

#endif


#undef dcplx
#endif
