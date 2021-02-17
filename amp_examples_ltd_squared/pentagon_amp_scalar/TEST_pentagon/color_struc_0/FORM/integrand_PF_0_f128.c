#include <tgmath.h>
# include <quadmath.h>
# include <signal.h>
# include "numerator.h"


static __complex128 diag_2_f128(__complex128 lm[], __complex128 params[]) {
	__complex128 Z10_,Z11_,Z12_,Z13_,Z14_,Z15_,Z16_,Z17_,Z18_,Z19_,Z20_,Z21_,Z22_,Z23_,Z24_,Z25_,Z26_,Z27_,Z28_,Z29_,Z6_,Z7_,Z8_,Z9_;
	__complex128 E0 = csqrtq(lm[34]+1.0q);
	__complex128 E1 = csqrtq(lm[34]+2*lm[30]+lm[2]+1.0q);
	__complex128 E2 = csqrtq(lm[34]+2*lm[30]-2*lm[18]+lm[14]-2*lm[10]+lm[2]+1.0q);
	__complex128 E3 = csqrtq(lm[34]+2*lm[30]-2*lm[27]+lm[25]-2*lm[21]-2*lm[18]+2*lm[16]+lm[14]-2*lm[10]+lm[2]+1.0q);
	__complex128 E4 = csqrtq(lm[34]-2*lm[32]+lm[7]+1.0q);
	__complex128 invd0 = 1.q/(lm[5]+E4+E0);
	__complex128 invd1 = 1.q/(lm[5]+lm[0]+E4+E1);
	__complex128 invd2 = 1.q/(-lm[8]+lm[5]+lm[0]+E4+E2);
	__complex128 invd3 = 1.q/(-lm[19]-lm[8]+lm[5]+lm[0]+E4+E3);
	__complex128 invd4 = 1.q/(lm[19]+lm[8]-lm[0]+E3+E0);
	__complex128 invd5 = 1.q/(lm[19]+lm[8]+E3+E1);
	__complex128 invd6 = 1.q/(lm[19]+E3+E2);
	__complex128 invd7 = 1.q/(lm[19]+lm[8]-lm[5]-lm[0]+E4+E3);
	__complex128 invd8 = 1.q/(lm[8]-lm[0]+E2+E0);
	__complex128 invd9 = 1.q/(lm[8]+E2+E1);
	__complex128 invd10 = 1.q/(-lm[19]+E3+E2);
	__complex128 invd11 = 1.q/(lm[8]-lm[5]-lm[0]+E4+E2);
	__complex128 invd12 = 1.q/(-lm[0]+E1+E0);
	__complex128 invd13 = 1.q/(-lm[8]+E2+E1);
	__complex128 invd14 = 1.q/(-lm[19]-lm[8]+E3+E1);
	__complex128 invd15 = 1.q/(-lm[5]-lm[0]+E4+E1);
	__complex128 invd16 = 1.q/(lm[0]+E1+E0);
	__complex128 invd17 = 1.q/(-lm[8]+lm[0]+E2+E0);
	__complex128 invd18 = 1.q/(-lm[19]-lm[8]+lm[0]+E3+E0);
	__complex128 invd19 = 1.q/(-lm[5]+E4+E0);


	Z6_= - invd17*invd3;
	Z7_= - invd3*invd12;
	Z8_=Z7_ + Z6_;
	Z9_=invd2*Z8_;
	Z10_= - invd7*invd12;
	Z11_= - invd17*invd7;
	Z12_=Z10_ + Z11_;
	Z13_=invd15*Z12_;
	Z14_= - invd12 - invd17;
	Z15_=invd15*Z14_;
	Z8_=Z8_ + Z15_;
	Z8_=invd14*Z8_;
	Z8_=Z8_ + Z9_ + Z13_;
	Z8_=invd13*Z8_;
	Z9_= - invd11*invd12;
	Z13_=Z7_ + Z9_;
	Z15_=invd8*Z13_;
	Z16_= - invd3*invd16;
	Z17_= - invd11*invd16;
	Z18_=Z16_ + Z17_;
	Z19_=invd18*Z18_;
	Z20_= - invd3 - invd11;
	Z21_=invd18*Z20_;
	Z13_=Z13_ + Z21_;
	Z13_=invd14*Z13_;
	Z13_=Z13_ + Z15_ + Z19_;
	Z13_=invd10*Z13_;
	Z15_= - invd8*invd3;
	Z19_=Z16_ + Z15_;
	Z19_=invd9*Z19_;
	Z21_= - invd2*invd16;
	Z22_= - invd16 - invd8;
	Z22_=invd9*Z22_;
	Z22_=Z21_ + Z22_;
	Z22_=invd5*Z22_;
	Z23_=invd2*Z16_;
	Z24_= - invd8 - invd2;
	Z25_=invd4*invd5*Z24_;
	Z19_=Z25_ + Z22_ + Z23_ + Z19_;
	Z19_=invd1*Z19_;
	Z22_= - invd11*invd7;
	Z23_= - invd17 - invd11;
	Z23_=invd18*Z23_;
	Z11_=Z23_ + Z11_ + Z22_;
	Z11_=invd15*Z11_;
	Z25_= - invd17*invd16;
	Z17_=Z25_ + Z17_;
	Z17_=invd18*Z17_;
	Z26_= - invd7*invd16;
	Z27_=invd17*Z26_;
	Z28_=invd11*Z26_;
	Z11_=Z11_ + Z17_ + Z27_ + Z28_;
	Z11_=invd19*Z11_;
	Z17_= - invd2*invd3;
	Z24_=invd4*Z24_;
	Z15_=Z24_ + Z15_ + Z17_;
	Z15_=invd1*Z15_;
	Z17_= - invd8*invd12;
	Z24_= - invd2*invd12;
	Z17_=Z17_ + Z24_;
	Z17_=invd4*Z17_;
	Z29_=invd8*Z7_;
	Z7_=invd2*Z7_;
	Z7_=Z15_ + Z17_ + Z29_ + Z7_;
	Z7_=invd0*Z7_;
	Z15_=invd8*Z20_;
	Z15_=Z18_ + Z15_;
	Z15_=invd10*Z15_;
	Z17_=invd8*Z22_;
	Z15_=Z15_ + Z28_ + Z17_;
	Z15_=invd9*Z15_;
	Z14_=invd2*Z14_;
	Z12_=Z12_ + Z14_;
	Z12_=invd13*Z12_;
	Z14_=invd2*Z25_;
	Z12_=Z12_ + Z27_ + Z14_;
	Z12_=invd6*Z12_;
	Z14_= - invd8*invd7;
	Z17_=Z26_ + Z14_;
	Z17_=invd9*Z17_;
	Z18_=Z26_ + Z21_;
	Z18_=invd6*Z18_;
	Z17_=Z17_ + Z18_;
	Z17_=invd5*Z17_;
	Z18_=Z10_ + Z24_;
	Z18_=invd6*Z18_;
	Z20_= - invd7 - invd2;
	Z20_=invd6*Z20_;
	Z14_=Z14_ + Z20_;
	Z14_=invd5*Z14_;
	Z20_=invd8*Z10_;
	Z14_=Z14_ + Z20_ + Z18_;
	Z14_=invd4*Z14_;
	Z9_=Z9_ + Z23_;
	Z9_=invd15*Z9_;
	Z6_=invd18*Z6_;
	Z6_=Z6_ + Z9_;
	Z6_=invd14*Z6_;
	Z9_=invd11*Z10_;
	Z10_=invd8*Z9_;
	Z16_=invd17*Z16_;
	Z18_=invd2*Z16_;
	Z16_=invd18*Z16_;
	Z9_=invd15*Z9_;


	return cpowq(2.q*M_PIq*I,1)/(2.q*E0*2.q*E1*2.q*E2*2.q*E3*2.q*E4)*(Z6_ + Z7_ + Z8_ + Z9_ + Z10_ + Z11_ + Z12_ + Z13_ + Z14_ + Z15_
       + Z16_ + Z17_ + Z18_ + Z19_);
}

static __complex128 diag_1_f128(__complex128 lm[], __complex128 params[]) {
	






	return cpowq(2.q*M_PIq*I,0)/(1)*(1);
}

static inline void evaluate_PF_0_0_f128(__complex128 lm[], __complex128 params[], __complex128* out) {
	__complex128 Z6_,Z7_;


    Z6_=diag_1_f128(lm, params);
    Z7_=diag_2_f128(lm, params);


	*out = Z7_*Z6_;
}
void evaluate_PF_0_f128(__complex128 lm[], __complex128 params[], int conf, __complex128* out) {
   switch(conf) {
		case 0: evaluate_PF_0_0_f128(lm, params, out); return;
		default: *out = 0.q;
    }
}
