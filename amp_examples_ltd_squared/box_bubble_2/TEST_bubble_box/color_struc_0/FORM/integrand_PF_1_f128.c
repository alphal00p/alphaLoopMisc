#include <tgmath.h>
# include <quadmath.h>
# include <signal.h>
# include "numerator.h"


static __complex128 diag_4_f128(__complex128 lm[], __complex128 params[]) {
	__complex128 Z10_,Z11_,Z12_,Z13_,Z14_,Z15_;
	__complex128 E0 = csqrtq(lm[53]+1e-08);
	__complex128 E1 = csqrtq(lm[53]-2*lm[51]+2*lm[44]+lm[42]-2*lm[40]+2*lm[33]+2*lm[31]+lm[29]-2*lm[27]+lm[7]+1e-08);
	__complex128 E2 = csqrtq(lm[53]-2*lm[51]+2*lm[44]+lm[42]-2*lm[40]+lm[7]+1e-08);
	__complex128 E3 = csqrtq(lm[53]-2*lm[51]+lm[7]+1e-08);
	__complex128 E4 = csqrtq(lm[62]+1e-08);
	__complex128 E5 = csqrtq(lm[62]-2*lm[60]-2*lm[58]+2*lm[46]+lm[42]-2*lm[40]-2*lm[38]+2*lm[35]+2*lm[31]+lm[29]-2*lm[27]-2*lm[25]+2*lm[22]+2*lm[18]+2*lm[16]+lm[14]-2*lm[12]-2*lm[10]+lm[7]+2*lm[4]+lm[2]+1e-08);
	__complex128 invd0 = 1.q/(lm[5]+E3+E0);
	__complex128 invd1 = 1.q/(lm[36]+lm[23]+E3+E1);
	__complex128 invd2 = 1.q/(lm[36]+E3+E2);
	__complex128 invd3 = 1.q/(-lm[36]-lm[23]-lm[8]+lm[5]+lm[0]+E5+E4);
	__complex128 invd4 = 1.q/(lm[36]+lm[23]+lm[8]-lm[5]-lm[0]+E5+E4);
	__complex128 invd5 = 1.q/(-lm[36]+lm[5]+E2+E0);
	__complex128 invd6 = 1.q/(lm[23]+E2+E1);
	__complex128 invd7 = 1.q/(-lm[36]+E3+E2);
	__complex128 invd8 = 1.q/(-lm[36]-lm[23]+lm[5]+E1+E0);
	__complex128 invd9 = 1.q/(-lm[23]+E2+E1);
	__complex128 invd10 = 1.q/(-lm[36]-lm[23]+E3+E1);
	__complex128 invd11 = 1.q/(lm[36]+lm[23]-lm[5]+E1+E0);
	__complex128 invd12 = 1.q/(lm[36]-lm[5]+E2+E0);
	__complex128 invd13 = 1.q/(-lm[5]+E3+E0);


	Z10_=invd5*invd0;
	Z11_=invd10*invd9;
	Z12_=invd0 + invd9;
	Z12_=invd2*Z12_;
	Z13_=invd5 + invd10;
	Z13_=invd7*Z13_;
	Z10_=Z13_ + Z12_ + Z10_ + Z11_;
	Z10_=invd8*Z10_;
	Z11_=invd1*invd6;
	Z12_=invd12*invd13;
	Z13_=invd1 + invd12;
	Z13_=invd2*Z13_;
	Z14_=invd6 + invd13;
	Z14_=invd7*Z14_;
	Z11_=Z14_ + Z13_ + Z11_ + Z12_;
	Z11_=invd11*Z11_;
	Z12_=invd1*invd0;
	Z13_=invd12*invd9;
	Z12_=Z12_ + Z13_;
	Z12_=invd2*Z12_;
	Z13_=invd5*invd6;
	Z14_=invd10*invd13;
	Z13_=Z13_ + Z14_;
	Z13_=invd7*Z13_;
	Z14_=invd0 + invd6;
	Z14_=invd5*invd1*Z14_;
	Z15_=invd9 + invd13;
	Z15_=invd12*invd10*Z15_;
	Z10_=Z11_ + Z10_ + Z13_ + Z12_ + Z14_ + Z15_;
	Z11_=invd3*Z10_;
	Z10_=invd4*Z10_;


	return cpowq(2.q*M_PIq*I,2)/(2.q*E0*2.q*E1*2.q*E2*2.q*E3*2.q*E4*2.q*E5)*(Z10_ + Z11_);
}

static __complex128 diag_3_f128(__complex128 lm[], __complex128 params[]) {
	__complex128 Z10_,Z11_,Z12_,Z13_,Z14_,Z15_,Z16_,Z17_,Z18_,Z19_,Z20_,Z21_;
	__complex128 E0 = csqrtq(lm[53]+1e-08);
	__complex128 E1 = csqrtq(lm[53]-2*lm[51]+2*lm[44]+lm[42]-2*lm[40]+2*lm[33]+2*lm[31]+lm[29]-2*lm[27]+lm[7]+1e-08);
	__complex128 E2 = csqrtq(lm[53]-   2*lm[51]+2*lm[44]+lm[42]-2*lm[40]+lm[7]+1e-08);
	__complex128 E3 = csqrtq(lm[53]-2*lm[51]+lm[7]+1e-08);
	__complex128 E4 = csqrtq(lm[62]+1e-08);
	__complex128 E5 = csqrtq(lm[62]-2*lm[60]-2*   lm[58]+2*lm[46]+lm[42]-2*lm[40]-2*lm[38]+2*lm[35]+2*lm[31]+lm[29]-2*lm[27]-2*lm[25]+2*   lm[22]+2*lm[18]+2*lm[16]+lm[14]-2*lm[12]-2*lm[10]+lm[7]+2*lm[4]+lm[2]+1e-08);
	__complex128 invd0 = 1.q/(lm[5]+E3+E0);
	__complex128 invd1 = 1.q/(lm[36]+lm[23]+E3+E1);
	__complex128 invd2 = 1.q/(lm[36]+E3+E2);
	__complex128 invd3 = 1.q/(-lm[36]-lm[23]-lm[8]+lm[5]+lm[0]+E5+E4);
	__complex128 invd4 = 1.q/(lm[36]+   lm[23]+lm[8]-lm[5]-lm[0]+E5+E4);
	__complex128 invd5 = 1.q/(-lm[36]+lm[5]+E2+E0);
	__complex128 invd6 = 1.q/(lm[23]+E2+E1);
	__complex128 invd7 = 1.q/(-lm[36]+E3+E2);
	__complex128 invd8 = 1.q/(-lm[36]-   lm[23]+lm[5]+E1+E0);
	__complex128 invd9 = 1.q/(-lm[23]+E2+E1);
	__complex128 invd10 = 1.q/(-lm[36]-lm[23]+E3+E1);
	__complex128 invd11 = 1.q/(lm[36]+lm[23]-lm[5]+E1+E0);
	__complex128 invd12 = 1.q/(lm[36]-   lm[5]+E2+E0);
	__complex128 invd13 = 1.q/(-lm[5]+E3+E0);


   Z10_= - invd10*invd9;
   Z11_= - invd2*invd9;
   Z11_=Z10_ + Z11_;
   Z11_=lm[36]*Z11_;
   Z12_=E1 - lm[23];
   Z13_=invd2*Z12_;
   Z14_= - lm[36]*invd2;
   Z13_=Z13_ + Z14_;
   Z13_=invd0*Z13_;
   Z14_=invd10*Z12_;
   Z15_= - lm[36]*invd10;
   Z14_=Z14_ + Z15_;
   Z14_=invd7*Z14_;
   Z15_=invd9*Z12_;
   Z16_=invd10*Z15_;
   Z15_=invd2*Z15_;
   Z11_=Z14_ + Z13_ + Z11_ + Z16_ + Z15_;
   Z13_=invd10*invd9;
   Z14_=invd2*invd9;
   Z15_=Z13_ + Z14_;
   Z16_=invd0*invd2;
   Z17_=invd7*invd10;
   Z16_=Z17_ + Z15_ + Z16_;
   Z16_=lm[5]*Z16_;
   Z11_=2*Z11_ + Z16_;
   Z11_=lm[5]*Z11_;
   Z16_=cpowq(E1,2);
   Z17_= - 2*E1 + lm[23];
   Z17_=lm[23]*Z17_;
   Z16_=Z16_ + Z17_;
   Z17_= - E1 + lm[23];
   Z18_=2*Z17_ + lm[36];
   Z18_=lm[36]*Z18_;
   Z18_=Z16_ + Z18_;
   Z19_=invd0*Z18_;
   Z18_=invd7*Z18_;
   Z12_=Z12_ - lm[36];
   Z20_=invd0*Z12_;
   Z12_=invd7*Z12_;
   Z12_=Z20_ + Z12_;
   Z20_=invd0 + invd7;
   Z20_=lm[5]*Z20_;
   Z12_=2*Z12_ + Z20_;
   Z12_=lm[5]*Z12_;
   Z12_=Z12_ + Z19_ + Z18_;
   Z12_=invd5*Z12_;
   Z18_=invd9*Z17_;
   Z19_=invd10*Z18_;
   Z18_=invd2*Z18_;
   Z18_=Z19_ + Z18_;
   Z15_=lm[36]*Z15_;
   Z15_=2*Z18_ + Z15_;
   Z15_=lm[36]*Z15_;
   Z18_=invd2*Z16_;
   Z19_=invd2*Z17_;
   Z20_=lm[36]*invd2;
   Z19_=2*Z19_ + Z20_;
   Z19_=lm[36]*Z19_;
   Z18_=Z18_ + Z19_;
   Z18_=invd0*Z18_;
   Z19_=invd10*Z16_;
   Z17_=invd10*Z17_;
   Z20_=lm[36]*invd10;
   Z17_=2*Z17_ + Z20_;
   Z17_=lm[36]*Z17_;
   Z17_=Z19_ + Z17_;
   Z17_=invd7*Z17_;
   Z16_=invd9*Z16_;
   Z19_=invd10*Z16_;
   Z16_=invd2*Z16_;
   Z11_=Z12_ + Z11_ + Z17_ + Z18_ + Z15_ + Z19_ + Z16_;
   Z11_=invd8*Z11_;
   Z12_= - invd1 - invd9;
   Z12_=invd2*Z12_;
   Z15_=invd1*E3;
   Z15_= - 1 + Z15_;
   Z15_=invd2*Z15_;
   Z15_= - invd1 + Z15_;
   Z15_=invd0*Z15_;
   Z16_= - invd6 - invd10;
   Z16_=invd7*Z16_;
   Z17_=lm[5]*invd0*invd2*invd1;
   Z18_= - invd1*invd6;
   Z10_=Z17_ + Z16_ + 2*Z15_ + Z12_ + Z18_ + Z10_;
   Z10_=lm[5]*Z10_;
   Z12_=invd1*E2;
   Z15_= - lm[36]*invd1;
   Z12_=Z15_ - 1 + Z12_;
   Z12_=invd0*Z12_;
   Z15_=invd6*E2;
   Z16_= - lm[36]*invd6;
   Z16_=Z16_ - 1 + Z15_;
   Z16_=invd7*Z16_;
   Z15_=invd1*Z15_;
   Z17_=lm[36]*Z18_;
   Z12_=Z16_ + Z12_ + Z15_ + Z17_;
   Z15_=invd1*invd6;
   Z16_=invd0*invd1;
   Z17_=invd7*invd6;
   Z16_=Z17_ + Z15_ + Z16_;
   Z16_=lm[5]*Z16_;
   Z12_=2*Z12_ + Z16_;
   Z12_=lm[5]*Z12_;
   Z16_= - invd6*E2;
   Z17_=invd1*Z16_;
   Z18_=lm[36]*Z15_;
   Z17_=2*Z17_ + Z18_;
   Z17_=lm[36]*Z17_;
   Z18_= - invd1*E2;
   Z18_=1 + Z18_;
   Z19_=lm[36]*invd1;
   Z18_=2*Z18_ + Z19_;
   Z18_=lm[36]*Z18_;
   Z19_=cpowq(E2,2);
   Z20_=invd1*Z19_;
   Z18_=Z18_ + lm[23] - E1 - E2 + Z20_;
   Z18_=invd0*Z18_;
   Z16_=1 + Z16_;
   Z20_=lm[36]*invd6;
   Z16_=2*Z16_ + Z20_;
   Z16_=lm[36]*Z16_;
   Z19_=invd6*Z19_;
   Z16_=Z16_ + lm[23] - E1 - E2 + Z19_;
   Z16_=invd7*Z16_;
   Z19_=invd1*Z19_;
   Z12_=Z12_ + Z16_ + Z18_ + Z19_ + Z17_;
   Z12_=invd5*Z12_;
   Z16_=cpowq(E0,2);
   Z17_=Z16_*invd11*invd13;
   Z18_=E0*invd11;
   Z18_= - 1 + Z18_;
   Z18_=E0*Z18_;
   Z19_= - E2 + Z18_;
   Z19_=invd6*Z19_;
   Z20_=E0*invd13;
   Z20_= - 1 + Z20_;
   Z20_=E0*Z20_;
   Z20_=lm[23] + Z20_ - E1;
   Z20_=invd10*Z20_;
   Z21_=invd6 + invd10;
   Z21_=lm[36]*Z21_;
   Z17_=Z21_ + Z20_ + Z19_ + 1 + Z17_;
   Z17_=invd7*Z17_;
   Z20_=Z16_*invd11*invd12;
   Z18_= - E3 + Z18_;
   Z18_=invd1*Z18_;
   Z21_=E0*invd12;
   Z21_= - 1 + Z21_;
   Z21_=E0*Z21_;
   Z21_=lm[23] + Z21_ - E1;
   Z21_=invd9*Z21_;
   Z18_=Z21_ + Z18_ + 1 + Z20_;
   Z18_=invd2*Z18_;
   Z13_=Z14_ + Z15_ + Z13_;
   Z13_=lm[36]*Z13_;
   Z14_= - E3 - E2;
   Z14_=invd1*Z14_;
   Z15_=invd1*cpowq(E3,2);
   Z15_=lm[23] - E1 - E3 + Z15_;
   Z15_=invd2*Z15_;
   Z20_=invd1 + invd2;
   Z20_=lm[36]*Z20_;
   Z14_=Z20_ + Z15_ + 1 + Z14_;
   Z14_=invd0*Z14_;
   Z15_=invd12*invd13;
   Z20_=Z16_*invd11*Z15_;
   Z15_=Z16_*Z15_;
   Z15_=Z15_ + Z21_;
   Z15_=invd10*Z15_;
   Z16_=1 + Z19_;
   Z16_=invd1*Z16_;
   Z10_=Z11_ + Z12_ + Z10_ + Z17_ + Z14_ + Z13_ + Z18_ + Z15_ + Z20_ + 
   Z16_;
   Z11_=invd3*Z10_;
   Z10_=invd4*Z10_;


	return cpowq(2.q*M_PIq*I,2)/(2.q*E0*2.q*E1*2.q*E2*2.q*E3*2.q*E4*2.q*E5)*(Z10_ + Z11_);
}

static __complex128 diag_2_f128(__complex128 lm[], __complex128 params[]) {
	__complex128 Z10_,Z11_,Z12_,Z13_,Z14_,Z15_,Z16_,Z17_;
	__complex128 E0 = csqrtq(lm[53]+1e-08);
	__complex128 E1 = csqrtq(lm[53]-2*lm[51]+2*lm[44]+lm[42]-2*lm[40]+2*lm[33]+2*lm[31]+lm[29]-2*lm[27]+lm[7]+1e-08);
	__complex128 E2 = csqrtq(lm[53]-   2*lm[51]+2*lm[44]+lm[42]-2*lm[40]+lm[7]+1e-08);
	__complex128 E3 = csqrtq(lm[53]-2*lm[51]+lm[7]+1e-08);
	__complex128 E4 = csqrtq(lm[62]+1e-08);
	__complex128 E5 = csqrtq(lm[62]-2*lm[60]-2*   lm[58]+2*lm[46]+lm[42]-2*lm[40]-2*lm[38]+2*lm[35]+2*lm[31]+lm[29]-2*lm[27]-2*lm[25]+2*   lm[22]+2*lm[18]+2*lm[16]+lm[14]-2*lm[12]-2*lm[10]+lm[7]+2*lm[4]+lm[2]+1e-08);
	__complex128 invd0 = 1.q/(lm[5]+E3+E0);
	__complex128 invd1 = 1.q/(lm[36]+lm[23]+E3+E1);
	__complex128 invd2 = 1.q/(lm[36]+E3+E2);
	__complex128 invd3 = 1.q/(-lm[36]-lm[23]-lm[8]+lm[5]+lm[0]+E5+E4);
	__complex128 invd4 = 1.q/(lm[36]+   lm[23]+lm[8]-lm[5]-lm[0]+E5+E4);
	__complex128 invd5 = 1.q/(-lm[36]+lm[5]+E2+E0);
	__complex128 invd6 = 1.q/(lm[23]+E2+E1);
	__complex128 invd7 = 1.q/(-lm[36]+E3+E2);
	__complex128 invd8 = 1.q/(-lm[36]-   lm[23]+lm[5]+E1+E0);
	__complex128 invd9 = 1.q/(-lm[23]+E2+E1);
	__complex128 invd10 = 1.q/(-lm[36]-lm[23]+E3+E1);
	__complex128 invd11 = 1.q/(lm[36]+lm[23]-lm[5]+E1+E0);
	__complex128 invd12 = 1.q/(lm[36]-   lm[5]+E2+E0);
	__complex128 invd13 = 1.q/(-lm[5]+E3+E0);


   Z10_=E0*invd11*invd12;
   Z11_=E0*invd12;
   Z11_= - 1 + Z11_;
   Z11_=invd9*Z11_;
   Z12_=E0*invd11;
   Z12_= - 1 + Z12_;
   Z12_=invd1*Z12_;
   Z13_=E3 + lm[5];
   Z13_=invd1*Z13_;
   Z13_= - 1 + Z13_;
   Z13_=invd0*Z13_;
   Z10_=Z13_ + Z12_ + Z10_ + Z11_;
   Z10_=invd2*Z10_;
   Z12_=lm[5] - lm[36] + E1 - lm[23];
   Z13_=invd9*Z12_;
   Z14_=invd0*Z12_;
   Z15_=Z13_ + Z14_;
   Z15_=invd2*Z15_;
   Z16_=invd7*Z12_;
   Z14_=Z14_ + Z16_;
   Z14_=invd5*Z14_;
   Z13_=invd10*Z13_;
   Z12_=invd7*invd10*Z12_;
   Z12_=Z14_ + Z12_ + Z13_ + Z15_;
   Z12_=invd8*Z12_;
   Z13_=E2 - lm[36];
   Z14_=invd6*Z13_;
   Z15_=lm[5]*invd6;
   Z16_=Z14_ + Z15_;
   Z16_=invd1*Z16_;
   Z14_=Z15_ - 1 + Z14_;
   Z14_=invd7*Z14_;
   Z13_=Z13_ + lm[5];
   Z13_=invd1*Z13_;
   Z13_= - 1 + Z13_;
   Z13_=invd0*Z13_;
   Z13_=Z14_ + Z16_ + Z13_;
   Z13_=invd5*Z13_;
   Z14_=invd12*invd13;
   Z15_=E0*Z14_;
   Z11_=Z15_ + Z11_;
   Z11_=invd10*Z11_;
   Z15_=invd11*invd13;
   Z16_=invd6*invd11;
   Z15_=Z15_ + Z16_;
   Z15_=E0*Z15_;
   Z17_=E0*invd13;
   Z17_= - 1 + Z17_;
   Z17_=invd10*Z17_;
   Z15_=Z17_ - invd6 + Z15_;
   Z15_=invd7*Z15_;
   Z14_=E0*invd11*Z14_;
   Z16_=E0*Z16_;
   Z16_= - invd6 + Z16_;
   Z16_=invd1*Z16_;
   Z17_= - invd0*invd1;
   Z10_=Z12_ + Z13_ + Z15_ + Z10_ + Z17_ + Z11_ + Z14_ + Z16_;
   Z11_=invd3*Z10_;
   Z10_=invd4*Z10_;


	return cpowq(2.q*M_PIq*I,2)/(2.q*E0*2.q*E1*2.q*E2*2.q*E3*2.q*E4*2.q*E5)*(Z10_ + Z11_);
}

static __complex128 diag_1_f128(__complex128 lm[], __complex128 params[]) {
	






	return cpowq(2.q*M_PIq*I,0)/(1)*(1);
}

static inline void evaluate_PF_1_0_f128(__complex128 lm[], __complex128 params[], __complex128* out) {
	__complex128 Z10_,Z11_,Z12_,Z13_,Z14_,Z15_,Z16_,Z17_,Z18_,Z19_,Z20_,Z21_,Z22_;


    Z10_=diag_1_f128(lm, params);
    Z11_=diag_2_f128(lm, params);
    Z12_=cpowq(xT,-1);
    Z13_=cpowq(xS,-1);
    Z14_=diag_3_f128(lm, params);
    Z15_=diag_4_f128(lm, params);
   Z16_=lm[23]*Z11_;
   Z17_=lm[24]*Z15_;
   Z16_=Z17_ + Z16_ - Z14_;
   Z17_=lm[53] - lm[33];
   Z18_=2*lm[49] - lm[1];
   Z19_=lm[37] - lm[30] - lm[44] + Z18_ + Z17_;
   Z19_=Z15_*Z19_;
   Z20_=lm[36]*Z11_;
   Z21_= - 2*lm[0]*Z11_;
   Z19_=Z21_ + Z19_ + Z20_ + Z16_;
   Z20_= - lm[41]*Z15_;
   Z19_=2*Z19_ + Z20_;
   Z19_=Z10_*Z19_;
   Z20_=Z10_*Z15_;
   Z22_= - lm[28]*Z20_;
   Z19_=Z19_ + Z22_;
   Z19_=Z13_*Z19_;
   Z17_=Z18_ + 2*Z17_;
   Z17_=Z15_*Z17_;
   Z16_=Z21_ + Z17_ + 2*Z16_;
   Z16_=Z10_*Z16_;
   Z16_=Z16_ + Z22_;
   Z16_=Z12_*Z16_;


	*out = Z16_ + Z19_ + Z20_;
}
void evaluate_PF_1_f128(__complex128 lm[], __complex128 params[], int conf, __complex128* out) {
   switch(conf) {
		case 0: evaluate_PF_1_0_f128(lm, params, out); return;
		default: *out = 0.q;
    }
}
