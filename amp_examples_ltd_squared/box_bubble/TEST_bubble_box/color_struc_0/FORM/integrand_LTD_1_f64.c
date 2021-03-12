#include <tgmath.h>
# include <quadmath.h>
# include <signal.h>
# include "numerator.h"

static inline void evaluate_LTD_1_0(double complex lm[], double complex params[], double complex* out) {
	double complex Z10_,Z11_,Z12_,Z13_,Z14_,Z15_,Z16_,Z17_,Z18_,Z19_,Z20_,Z21_,Z22_,Z23_,Z6_,Z7_,Z8_,Z9_;
	double complex E0 = sqrt(lm[25]+0.0);
	double complex E1 = sqrt(lm[25]+2*lm[21]+lm[2]+0.0);
	double complex E2 = sqrt(lm[25]+2*lm[21]-2*lm[16]+lm[14]-2*lm[10]+lm[2]+0.0);
	double complex E3 = sqrt(lm[25]-2*lm[23]+lm[7]+0.0);
	double complex E4 = sqrt(   lm[34]+1.0);
	double complex invd0 = 1./(2*E0);
	double complex invd1 = 1./(lm[0]+E1+E0);
	double complex invd2 = 1./(lm[0]-E1+E0);
	double complex invd3 = 1./(-lm[8]+lm[0]+E2+E0);
	double complex invd4 = 1./(-lm[8]+lm[0]-E2+E0);
	double complex invd5 = 1./(-lm[5]+E3+E0);
	double complex invd6 = 1./(-lm[5]-   E3+E0);
	double complex invd7 = 1./(2*E4);
	double complex invd8 = 1./(-lm[0]+E1+E0);
	double complex invd9 = 1./(-lm[0]+E1-E0);
	double complex invd10 = 1./(2*E1);
	double complex invd11 = 1./(-lm[8]+E2+E1);
	double complex invd12 = 1./(-lm[8]-E2+E1);
	double complex invd13 = 1./(-lm[5]-lm[0]+   E3+E1);
	double complex invd14 = 1./(-lm[5]-lm[0]-E3+E1);
	double complex invd15 = 1./(2*E4);
	double complex invd16 = 1./(lm[8]-lm[0]+E2+E0);
	double complex invd17 = 1./(lm[8]-lm[0]+E2-E0);
	double complex invd18 = 1./(lm[8]+E2+E1);
	double complex invd19 = 1./(lm[8]+   E2-E1);
	double complex invd20 = 1./(2*E2);
	double complex invd21 = 1./(lm[8]-lm[5]-lm[0]+E3+E2);
	double complex invd22 = 1./(lm[8]-lm[5]-lm[0]-E3+E2);
	double complex invd23 = 1./(2*E4);
	double complex invd24 = 1./(lm[5]+E3+E0);
	double complex invd25 = 1./(lm[5]+E3-   E0);
	double complex invd26 = 1./(lm[5]+lm[0]+E3+E1);
	double complex invd27 = 1./(lm[5]+lm[0]+E3-E1);
	double complex invd28 = 1./(-lm[8]+lm[5]+lm[0]+E3+E2);
	double complex invd29 = 1./(-lm[8]+lm[5]+lm[0]+E3-E2);
	double complex invd30 = 1./(   2*E3);
	double complex invd31 = 1./(2*E4);


    Z6_=pow(xT,-1);
    Z7_=pow(xS,-1);
   Z8_=2*E2 + lm[8];
   Z8_=lm[8]*Z8_;
   Z9_= - E2 - lm[8];
   Z10_=lm[5]*Z9_;
   Z9_=Z9_ + lm[5];
   Z9_=lm[0]*Z9_;
   Z11_=lm[6] + 2*lm[23];
   Z12_=2*pow(E2,2);
   Z13_= - 2*lm[21];
   Z14_= - 2*lm[25];
   Z8_=2*Z9_ + 2*Z10_ + 2*Z8_ + Z14_ + Z13_ + lm[1] + Z11_ + Z12_;
   Z9_=pow(invd23,3);
   Z8_=Z9_*invd22*invd21*invd20*invd19*invd18*invd17*invd16*Z8_;
   Z10_= - lm[5]*E0;
   Z15_=2*pow(E0,2);
   Z16_=2*lm[0]*E0;
   Z10_=Z16_ + 2*Z10_ + Z14_ + Z13_ + lm[1] + Z11_ + Z15_;
   Z17_=pow(invd7,3);
   Z10_=Z17_*invd6*invd5*invd4*invd3*invd2*invd1*invd0*Z10_;
   Z18_= - lm[5]*E1;
   Z19_= - E1 + lm[5];
   Z19_=lm[0]*Z19_;
   Z20_=2*pow(E1,2);
   Z18_=2*Z19_ + 2*Z18_ + Z14_ + Z13_ + lm[1] + Z11_ + Z20_;
   Z19_=pow(invd15,3);
   Z18_=Z19_*invd14*invd13*invd12*invd11*invd10*invd9*invd8*Z18_;
   Z21_=lm[5]*E3;
   Z22_=2*pow(E3,2);
   Z23_=E3 + lm[5];
   Z23_=2*lm[0]*Z23_;
   Z11_=Z23_ + 2*Z21_ + Z14_ + Z13_ + lm[1] + Z11_ + Z22_;
   Z21_=pow(invd31,3);
   Z11_=Z21_*invd30*invd29*invd28*invd27*invd26*invd25*invd24*Z11_;
   Z8_=Z11_ + Z8_ + Z10_ + Z18_;
   Z8_=Z7_*Z8_;
   Z10_=2*lm[16] - 2*lm[9] + lm[13];
   Z11_= - lm[8]*E3;
   Z18_=lm[5] + 2*E3 - lm[8];
   Z18_=lm[5]*Z18_;
   Z11_=Z23_ + 2*Z18_ + 2*Z11_ + Z14_ + Z13_ + lm[1] + Z10_ + Z22_;
   Z11_=Z21_*invd30*invd29*invd28*invd27*invd26*invd25*invd24*Z11_;
   Z18_= - lm[8]*E0;
   Z15_=Z16_ + 2*Z18_ + Z14_ + Z13_ + lm[1] + Z10_ + Z15_;
   Z15_=Z17_*invd6*invd5*invd4*invd3*invd2*invd1*invd0*Z15_;
   Z16_= - lm[8]*E1;
   Z18_= - E1 + lm[8];
   Z18_=lm[0]*Z18_;
   Z16_=2*Z18_ + 2*Z16_ + Z14_ + Z13_ + lm[1] + Z10_ + Z20_;
   Z16_=Z19_*invd14*invd13*invd12*invd11*invd10*invd9*invd8*Z16_;
   Z18_=lm[8]*E2;
   Z20_= - lm[0]*E2;
   Z10_=2*Z20_ + 2*Z18_ + Z14_ + Z13_ + lm[1] + Z10_ + Z12_;
   Z10_=Z9_*invd22*invd21*invd20*invd19*invd18*invd17*invd16*Z10_;
   Z10_=Z11_ + Z10_ + Z15_ + Z16_;
   Z10_=Z6_*Z10_;
   Z11_= - Z17_*invd6*invd5*invd4*invd3*invd2*invd1*invd0;
   Z12_= - Z19_*invd14*invd13*invd12*invd11*invd10*invd9*invd8;
   Z9_= - Z9_*invd22*invd21*invd20*invd19*invd18*invd17*invd16;
   Z13_= - Z21_*invd30*invd29*invd28*invd27*invd26*invd25*invd24;
   Z8_=Z10_ + Z8_ + Z13_ + Z9_ + Z11_ + Z12_;


	*out = pow(2.*pi*I,2)/(1)*(2*Z8_);
}
void evaluate_LTD_1(double complex lm[], double complex params[], int conf, double complex* out) {
   switch(conf) {
		case 0: evaluate_LTD_1_0(lm, params, out); return;
		default: *out = 0.;
    }
}
