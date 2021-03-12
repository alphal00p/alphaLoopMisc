#include <tgmath.h>
# include <quadmath.h>
# include <signal.h>
# include "numerator.h"

// polynomial in 0,c4,c5
static inline int evaluate_0_0(double complex lm[], double complex params[], double complex* out) {
	double complex Z10_,Z11_,Z12_,Z13_,Z14_,Z15_,Z16_,Z9_;
    Z9_=pow(xT,-1);
    Z10_=pow(xS,-1);
   Z11_= - lm[0] + lm[23];
   Z11_=Z9_*Z11_;
   Z12_=lm[23] + lm[36] - 2*lm[0];
   Z12_=Z10_*Z12_;
   Z11_=Z11_ + Z12_;
   Z11_=2*Z11_;
   Z12_= - Z9_ - Z10_;
   Z12_=2*Z12_;
   Z13_= - lm[30] + lm[37];
   Z14_=2*lm[24];
   Z15_= - 2*lm[33];
   Z16_=2*lm[53];
   Z13_=Z16_ + 4*lm[49] + Z15_ - lm[28] + Z14_ - 2*lm[1] - 2*lm[44] + 2*Z13_ - 
   lm[41];
   Z13_=Z10_*Z13_;
   Z14_=Z16_ + 2*lm[49] + Z15_ - lm[28] - lm[1] + Z14_;
   Z14_=Z9_*Z14_;
   Z13_=Z13_ + 1 + Z14_;


	out[0] = Z13_;
	out[2] = Z11_;
	out[7] = Z12_;
	return 8; 
}

int evaluate_0(double complex lm[], double complex params[], int conf, double complex* out) {
   switch(conf) {
		case 0: return evaluate_0_0(lm, params, out);
		default: out[0] = 0.; return 1;
    }
}

int get_rank_0(int conf) {
   switch(conf) {
		case 0: return 2;
		default: return 0;
    }
}
