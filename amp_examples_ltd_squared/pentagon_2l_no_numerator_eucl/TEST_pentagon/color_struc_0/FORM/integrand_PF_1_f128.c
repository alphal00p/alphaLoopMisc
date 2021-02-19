#include <tgmath.h>
# include <quadmath.h>
# include <signal.h>
# include "numerator.h"


static __complex128 diag_2_f128(__complex128 lm[], __complex128 params[]) {
	__complex128 Z10_,Z11_,Z12_,Z13_,Z14_,Z15_,Z16_,Z17_,Z18_,Z19_,Z20_,Z6_,Z7_,Z8_,Z9_;
	__complex128 E0 = csqrtq(lm[38]+small_mass_sq);
	__complex128 E1 = csqrtq(lm[38]+2*lm[34]+lm[2]+small_mass_sq);
	__complex128 E2 = csqrtq(lm[38]+2*lm[34]-2*lm[18]+lm[14]-2*lm[10]+lm[2]+small_mass_sq);
	__complex128 E3 = csqrtq(lm[38]+2*lm[34]-2*lm[29]+lm[27]-2*lm[23]-2*lm[18]+2*lm[16]+lm[14]-2*lm[10]+lm[2]+small_mass_sq);
	__complex128 E4 = csqrtq(lm[38]-2*lm[36]+lm[7]+small_mass_sq);
	__complex128 E5 = csqrtq(lm[47]+10000.0q);
	__complex128 invd0 = 1.q/(lm[5]+E4+E0);
	__complex128 invd1 = 1.q/(lm[5]+lm[0]+E4+E1);
	__complex128 invd2 = 1.q/(-lm[8]+lm[5]+lm[0]+E4+E2);
	__complex128 invd3 = 1.q/(-lm[21]-lm[8]+lm[5]+lm[0]+E4+E3);
	__complex128 invd4 = 1.q/(2*E5);
	__complex128 invd5 = 1.q/(lm[21]+lm[8]-lm[0]+E3+E0);
	__complex128 invd6 = 1.q/(lm[21]+lm[8]+E3+E1);
	__complex128 invd7 = 1.q/(lm[21]+E3+E2);
	__complex128 invd8 = 1.q/(lm[21]+lm[8]-lm[5]-lm[0]+E4+E3);
	__complex128 invd9 = 1.q/(lm[8]-lm[0]+E2+E0);
	__complex128 invd10 = 1.q/(lm[8]+E2+E1);
	__complex128 invd11 = 1.q/(-lm[21]+E3+E2);
	__complex128 invd12 = 1.q/(lm[8]-lm[5]-lm[0]+E4+E2);
	__complex128 invd13 = 1.q/(-lm[0]+E1+E0);
	__complex128 invd14 = 1.q/(-lm[8]+E2+E1);
	__complex128 invd15 = 1.q/(-lm[21]-lm[8]+E3+E1);
	__complex128 invd16 = 1.q/(-lm[5]-lm[0]+E4+E1);
	__complex128 invd17 = 1.q/(lm[0]+E1+E0);
	__complex128 invd18 = 1.q/(-lm[8]+lm[0]+E2+E0);
	__complex128 invd19 = 1.q/(-lm[21]-lm[8]+lm[0]+E3+E0);
	__complex128 invd20 = 1.q/(-lm[5]+E4+E0);


	Z6_= - invd7*invd5;
	Z7_= - invd14*invd7;
	Z8_= - invd16*invd14;
	Z9_= - invd9*invd5;
	Z10_= - invd16 - invd9;
	Z10_=invd12*Z10_;
	Z8_=Z10_ + Z9_ + Z8_ + Z6_ + Z7_;
	Z8_=invd8*Z8_;
	Z9_= - invd5*invd0;
	Z6_=Z7_ + Z9_ + Z6_;
	Z6_=invd2*Z6_;
	Z10_= - invd15*invd11;
	Z11_= - invd16*invd15;
	Z12_= - invd9*invd11;
	Z11_=Z12_ + Z10_ + Z11_;
	Z11_=invd12*Z11_;
	Z12_= - invd11 - invd14;
	Z12_=invd15*Z12_;
	Z13_= - invd0 - invd14;
	Z13_=invd2*Z13_;
	Z14_= - invd0 - invd11;
	Z14_=invd9*Z14_;
	Z12_=Z14_ + Z12_ + Z13_;
	Z12_=invd3*Z12_;
	Z13_= - invd15*invd14;
	Z14_=invd16*Z13_;
	Z9_=invd9*Z9_;
	Z6_=Z8_ + Z12_ + Z11_ + Z9_ + Z14_ + Z6_;
	Z6_=invd13*Z6_;
	Z8_= - invd10*invd1;
	Z9_= - invd11*invd10;
	Z11_= - invd19*invd11;
	Z12_= - invd2*invd1;
	Z15_= - invd19 - invd2;
	Z15_=invd18*Z15_;
	Z11_=Z15_ + Z12_ + Z11_ + Z8_ + Z9_;
	Z11_=invd3*Z11_;
	Z12_= - invd7*invd6;
	Z15_= - invd10*invd6;
	Z16_= - invd20 - invd10;
	Z16_=invd12*Z16_;
	Z17_= - invd20 - invd7;
	Z17_=invd18*Z17_;
	Z16_=Z17_ + Z16_ + Z12_ + Z15_;
	Z16_=invd8*Z16_;
	Z17_= - invd6*invd1;
	Z12_=Z17_ + Z12_;
	Z12_=invd2*Z12_;
	Z18_= - invd20 - invd11;
	Z18_=invd19*Z18_;
	Z18_=Z9_ + Z18_;
	Z18_=invd12*Z18_;
	Z19_= - invd19*invd20;
	Z20_= - invd2*invd7;
	Z19_=Z19_ + Z20_;
	Z19_=invd18*Z19_;
	Z17_=invd10*Z17_;
	Z11_=Z16_ + Z11_ + Z19_ + Z18_ + Z17_ + Z12_;
	Z11_=invd17*Z11_;
	Z12_= - invd1*invd0;
	Z8_=Z9_ + Z12_ + Z8_;
	Z8_=invd9*Z8_;
	Z16_= - invd19*invd15;
	Z18_= - invd2*invd14;
	Z13_=Z18_ + Z13_ + Z16_;
	Z13_=invd18*Z13_;
	Z16_=invd19*Z10_;
	Z18_=invd2*Z12_;
	Z8_=Z13_ + Z8_ + Z16_ + Z18_;
	Z8_=invd3*Z8_;
	Z13_= - invd6*invd5;
	Z15_=Z13_ + Z15_;
	Z15_=invd9*Z15_;
	Z16_= - invd16*invd20;
	Z18_= - invd9*invd10;
	Z16_=Z16_ + Z18_;
	Z16_=invd12*Z16_;
	Z18_= - invd20 - invd14;
	Z18_=invd16*Z18_;
	Z18_=Z7_ + Z18_;
	Z18_=invd18*Z18_;
	Z13_=invd7*Z13_;
	Z15_=Z18_ + Z16_ + Z13_ + Z15_;
	Z15_=invd8*Z15_;
	Z16_= - invd20 - invd15;
	Z16_=invd16*Z16_;
	Z18_=invd19*Z16_;
	Z7_=invd2*Z7_;
	Z7_=Z7_ + Z14_ + Z18_;
	Z7_=invd18*Z7_;
	Z12_=invd5*Z12_;
	Z14_= - invd6*invd5*invd1;
	Z12_=Z12_ + Z14_;
	Z13_=Z12_ + Z13_;
	Z13_=invd2*Z13_;
	Z12_=Z12_ + Z17_;
	Z12_=invd9*Z12_;
	Z10_=Z10_ + Z16_;
	Z10_=invd19*Z10_;
	Z9_=invd9*Z9_;
	Z9_=Z10_ + Z9_;
	Z9_=invd12*Z9_;
	Z6_=Z11_ + Z6_ + Z15_ + Z8_ + Z7_ + Z9_ + Z13_ + Z12_;


	return cpowq(2.q*M_PIq*I,2)/(2.q*E0*2.q*E1*2.q*E2*2.q*E3*2.q*E4*cpowq(2.q*E5,2))*(2*invd4*Z6_);
}

static __complex128 diag_1_f128(__complex128 lm[], __complex128 params[]) {
	






	return cpowq(2.q*M_PIq*I,0)/(1)*(1);
}

static inline void evaluate_PF_1_0_f128(__complex128 lm[], __complex128 params[], __complex128* out) {
	__complex128 Z6_,Z7_;


    Z6_=diag_1_f128(lm, params);
    Z7_=diag_2_f128(lm, params);


	*out =  - Z7_*Z6_;
}
void evaluate_PF_1_f128(__complex128 lm[], __complex128 params[], int conf, __complex128* out) {
   switch(conf) {
		case 0: evaluate_PF_1_0_f128(lm, params, out); return;
		default: *out = 0.q;
    }
}
