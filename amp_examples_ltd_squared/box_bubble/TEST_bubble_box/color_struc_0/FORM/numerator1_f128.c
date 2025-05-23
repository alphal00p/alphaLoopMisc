#include <tgmath.h>
# include <quadmath.h>
# include <signal.h>
# include "numerator.h"


static inline int evaluate_1_0_f128(__complex128 lm[], __complex128 params[], __complex128* out) {
	__complex128 Z10_,Z11_,Z12_,Z13_,Z14_,Z15_,Z9_;

    Z9_=cpowq(xT,-1);
    Z10_=cpowq(xS,-1);
   Z11_=lm[5] - lm[0];
   Z11_=Z10_*Z11_;
   Z12_=lm[8] - lm[0];
   Z12_=Z9_*Z12_;
   Z11_=Z11_ + Z12_;
   Z11_=2*Z11_;
   Z12_= - Z10_ - Z9_;
   Z12_=2*Z12_;
   Z13_=2*lm[21];
   Z14_=2*lm[25];
   Z15_=Z14_ + Z13_ - lm[1] - lm[6] - 2*lm[23];
   Z15_=Z10_*Z15_;
   Z13_=Z14_ + Z13_ - lm[1] - 2*lm[16] + 2*lm[9] - lm[13];
   Z13_=Z9_*Z13_;
   Z13_=Z13_ + 1 + Z15_;


	out[0] = Z13_;
	out[2] = Z11_;
	out[7] = Z12_;
	return 8;
}

int evaluate_1_f128(__complex128 lm[], __complex128 params[], int conf, __complex128* out) {
    switch(conf) {
		case 0: return evaluate_1_0_f128(lm, params, out);
		default: out[0] = 0.; return 1;
    }
}

