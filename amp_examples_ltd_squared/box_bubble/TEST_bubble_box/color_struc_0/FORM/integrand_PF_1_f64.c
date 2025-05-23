#include <tgmath.h>
# include <quadmath.h>
# include <signal.h>
# include "numerator.h"

static double complex diag_4(double complex lm[], double complex params[]) {
	double complex Z10_,Z11_,Z12_,Z13_,Z14_,Z15_;
	double complex E0 = sqrt(lm[25]+small_mass_sq);
	double complex E1 = sqrt(lm[25]+2*lm[21]+lm[2]+small_mass_sq);
	double complex E2 = sqrt(lm[25]+2*lm[21]-2*lm[16]+lm[14]-2*lm[10]+lm[2]+small_mass_sq);
	double complex E3 = sqrt(lm[25]-2*lm[23]+lm[7]+small_mass_sq);
	double complex E4 = sqrt(lm[34]+1.0);
	double complex invd0 = 1./(lm[5]+E3+E0);
	double complex invd1 = 1./(lm[5]+lm[0]+E3+E1);
	double complex invd2 = 1./(-lm[8]+lm[5]+lm[0]+E3+E2);
	double complex invd3 = 1./(2*E4);
	double complex invd4 = 1./(lm[8]-lm[0]+E2+E0);
	double complex invd5 = 1./(lm[8]+E2+E1);
	double complex invd6 = 1./(lm[8]-lm[5]-lm[0]+E3+E2);
	double complex invd7 = 1./(-lm[0]+E1+E0);
	double complex invd8 = 1./(-lm[8]+E2+E1);
	double complex invd9 = 1./(-lm[5]-lm[0]+E3+E1);
	double complex invd10 = 1./(lm[0]+E1+E0);
	double complex invd11 = 1./(-lm[8]+lm[0]+E2+E0);
	double complex invd12 = 1./(-lm[5]+E3+E0);


	Z10_=invd4*invd0;
	Z11_=invd9*invd8;
	Z12_=invd0 + invd8;
	Z12_=invd2*Z12_;
	Z13_=invd4 + invd9;
	Z13_=invd6*Z13_;
	Z10_=Z13_ + Z12_ + Z10_ + Z11_;
	Z10_=invd7*Z10_;
	Z11_=invd1*invd5;
	Z12_=invd11*invd12;
	Z13_=invd1 + invd11;
	Z13_=invd2*Z13_;
	Z14_=invd5 + invd12;
	Z14_=invd6*Z14_;
	Z11_=Z14_ + Z13_ + Z11_ + Z12_;
	Z11_=invd10*Z11_;
	Z12_=invd1*invd0;
	Z13_=invd11*invd8;
	Z12_=Z12_ + Z13_;
	Z12_=invd2*Z12_;
	Z13_=invd4*invd5;
	Z14_=invd9*invd12;
	Z13_=Z13_ + Z14_;
	Z13_=invd6*Z13_;
	Z14_=invd0 + invd5;
	Z14_=invd4*invd1*Z14_;
	Z15_=invd8 + invd12;
	Z15_=invd11*invd9*Z15_;
	Z10_=Z11_ + Z10_ + Z13_ + Z12_ + Z14_ + Z15_;


	return pow(2.*pi*I,2)/(2.*E0*2.*E1*2.*E2*2.*E3*pow(2.*E4,2))*(2*invd3*Z10_);
}
static double complex diag_3(double complex lm[], double complex params[]) {
	double complex Z10_,Z11_,Z12_,Z13_,Z14_,Z15_,Z16_,Z17_,Z18_,Z19_,Z20_,Z21_;
	double complex E0 = sqrt(lm[25]+small_mass_sq);
	double complex E1 = sqrt(lm[25]+2*lm[21]+lm[2]+small_mass_sq);
	double complex E2 = sqrt(lm[25]+2*lm[21]-2*lm[16]+lm[14]-2*lm[10]+lm[2]+small_mass_sq);
	double complex E3 = sqrt(lm[25]-2*lm[23]+lm[7]+small_mass_sq);
	double complex E4 = sqrt(   lm[34]+1.0);
	double complex invd0 = 1./(lm[5]+E3+E0);
	double complex invd1 = 1./(lm[5]+lm[0]+E3+E1);
	double complex invd2 = 1./(-lm[8]+lm[5]+lm[0]+E3+E2);
	double complex invd3 = 1./(2*E4);
	double complex invd4 = 1./(lm[8]-lm[0]+E2+E0);
	double complex invd5 = 1./(lm[8]+E2+E1);
	double complex invd6 = 1./(   lm[8]-lm[5]-lm[0]+E3+E2);
	double complex invd7 = 1./(-lm[0]+E1+E0);
	double complex invd8 = 1./(-lm[8]+E2+E1);
	double complex invd9 = 1./(-lm[5]-lm[0]+E3+E1);
	double complex invd10 = 1./(lm[0]+E1+E0);
	double complex invd11 = 1./(-   lm[8]+lm[0]+E2+E0);
	double complex invd12 = 1./(-lm[5]+E3+E0);


   Z10_=invd9*invd8;
   Z11_=invd7*Z10_;
   Z12_=invd7*invd8;
   Z13_=invd2*Z12_;
   Z14_=invd7*invd9;
   Z15_=invd6*Z14_;
   Z16_=invd0*invd2*invd7;
   Z13_=Z16_ + Z15_ + Z11_ + Z13_;
   Z13_=lm[0]*Z13_;
   Z15_= - invd9*invd8;
   Z16_=E1*invd7*Z15_;
   Z17_= - E1*invd7*invd8;
   Z17_=invd8 + 2*Z17_;
   Z17_=invd2*Z17_;
   Z18_= - E1*invd7*invd9;
   Z18_=invd5 + invd9 + 2*Z18_;
   Z18_=invd6*Z18_;
   Z19_= - E1*invd7;
   Z20_=1 + 2*Z19_;
   Z20_=invd0*invd2*Z20_;
   Z21_=invd5 + invd0;
   Z21_=invd1*Z21_;
   Z10_=Z13_ + Z21_ + Z20_ + Z18_ + Z17_ + Z10_ + 2*Z16_;
   Z10_=lm[0]*Z10_;
   Z13_=1 + Z19_;
   Z16_= - E2 - lm[8];
   Z17_=invd5*Z16_;
   Z18_=Z13_ + Z17_;
   Z18_=invd6*Z18_;
   Z19_=invd0*Z16_;
   Z17_=Z17_ + Z19_;
   Z17_=invd1*Z17_;
   Z13_=invd0*Z13_;
   Z13_=Z17_ + Z18_ + Z13_;
   Z17_=invd7 + invd5;
   Z17_=invd6*Z17_;
   Z18_=invd0*invd7;
   Z17_=Z21_ + Z17_ + Z18_;
   Z17_=lm[0]*Z17_;
   Z13_=2*Z13_ + Z17_;
   Z13_=lm[0]*Z13_;
   Z17_=E1*invd7;
   Z17_= - 1 + Z17_;
   Z17_=E1*Z17_;
   Z18_=Z16_ + Z17_;
   Z19_=pow(E2,2);
   Z20_=2*E2 + lm[8];
   Z20_=lm[8]*Z20_;
   Z19_=Z19_ + Z20_;
   Z20_=invd5*Z19_;
   Z21_=Z18_ + Z20_;
   Z21_=invd6*Z21_;
   Z19_=invd0*Z19_;
   Z19_=Z20_ + Z19_;
   Z19_=invd1*Z19_;
   Z18_=invd0*Z18_;
   Z13_=Z13_ + Z19_ + Z21_ + Z18_;
   Z13_=invd4*Z13_;
   Z18_=E0*invd10;
   Z18_= - 1 + Z18_;
   Z18_=E0*Z18_;
   Z16_=Z16_ + Z18_;
   Z16_=invd5*Z16_;
   Z19_= - E3 - lm[5];
   Z18_=Z19_ + Z18_;
   Z18_=invd2*Z18_;
   Z20_=pow(E3,2);
   Z21_=2*E3 + lm[5];
   Z21_=lm[5]*Z21_;
   Z20_=Z20_ + Z21_;
   Z20_=invd2*Z20_;
   Z20_=Z20_ - lm[8] + Z19_ - E2;
   Z20_=invd0*Z20_;
   Z18_=Z20_ + Z18_ + 1 + Z16_;
   Z18_=invd1*Z18_;
   Z20_=invd10*invd12;
   Z21_=invd9*invd12;
   Z20_=Z20_ + Z21_;
   Z20_=E0*Z20_;
   Z20_= - invd9 + Z20_;
   Z20_=E0*Z20_;
   Z14_=E1*Z14_;
   Z14_= - invd9 + Z14_;
   Z14_=E1*Z14_;
   Z14_=Z16_ + Z14_ + 1 + Z20_;
   Z14_=invd6*Z14_;
   Z16_=invd11*invd12;
   Z20_=invd8*invd11;
   Z21_=Z16_ + Z20_;
   Z21_=invd9*Z21_;
   Z16_=invd10*Z16_;
   Z16_=Z16_ + Z21_;
   Z16_=E0*Z16_;
   Z16_=Z15_ + Z16_;
   Z16_=E0*Z16_;
   Z11_=E1*Z11_;
   Z11_=Z15_ + Z11_;
   Z11_=E1*Z11_;
   Z15_=invd10*invd11;
   Z15_=Z15_ + Z20_;
   Z15_=E0*Z15_;
   Z15_= - invd8 + Z15_;
   Z15_=E0*Z15_;
   Z12_=E1*Z12_;
   Z12_= - invd8 + Z12_;
   Z12_=E1*Z12_;
   Z12_=Z12_ + 1 + Z15_;
   Z12_=invd2*Z12_;
   Z15_=Z19_ + Z17_;
   Z15_=invd2*Z15_;
   Z15_=1 + Z15_;
   Z15_=invd0*Z15_;
   Z10_=Z13_ + Z10_ + Z18_ + Z15_ + Z14_ + Z12_ + Z16_ + Z11_;


	return pow(2.*pi*I,2)/(2.*E0*2.*E1*2.*E2*2.*E3*pow(2.*E4,2))*(2*invd3*Z10_);
}
static double complex diag_2(double complex lm[], double complex params[]) {
	double complex Z10_,Z11_,Z12_,Z13_,Z14_,Z15_,Z16_,Z17_,Z18_;
	double complex E0 = sqrt(lm[25]+small_mass_sq);
	double complex E1 = sqrt(lm[25]+2*lm[21]+lm[2]+small_mass_sq);
	double complex E2 = sqrt(lm[25]+2*lm[21]-2*lm[16]+lm[14]-2*lm[10]+lm[2]+small_mass_sq);
	double complex E3 = sqrt(lm[25]-2*lm[23]+lm[7]+small_mass_sq);
	double complex E4 = sqrt(   lm[34]+1.0);
	double complex invd0 = 1./(lm[5]+E3+E0);
	double complex invd1 = 1./(lm[5]+lm[0]+E3+E1);
	double complex invd2 = 1./(-lm[8]+lm[5]+lm[0]+E3+E2);
	double complex invd3 = 1./(2*E4);
	double complex invd4 = 1./(lm[8]-lm[0]+E2+E0);
	double complex invd5 = 1./(lm[8]+E2+E1);
	double complex invd6 = 1./(   lm[8]-lm[5]-lm[0]+E3+E2);
	double complex invd7 = 1./(-lm[0]+E1+E0);
	double complex invd8 = 1./(-lm[8]+E2+E1);
	double complex invd9 = 1./(-lm[5]-lm[0]+E3+E1);
	double complex invd10 = 1./(lm[0]+E1+E0);
	double complex invd11 = 1./(-   lm[8]+lm[0]+E2+E0);
	double complex invd12 = 1./(-lm[5]+E3+E0);


   Z10_=invd8*E1;
   Z11_= - lm[0]*invd8;
   Z12_=E1 - lm[0];
   Z13_=invd0*Z12_;
   Z11_=Z13_ + Z10_ + Z11_;
   Z11_=invd2*Z11_;
   Z10_=invd9*Z10_;
   Z14_= - invd9*invd8;
   Z15_=lm[0]*Z14_;
   Z10_=Z11_ + Z10_ + Z15_;
   Z10_=invd7*Z10_;
   Z11_=invd10*invd12;
   Z15_=invd9*invd12;
   Z16_=invd5*invd10;
   Z11_=Z16_ + Z11_ + Z15_;
   Z11_=E0*Z11_;
   Z15_=invd9*E1;
   Z17_= - lm[0]*invd9;
   Z15_=Z15_ + Z17_;
   Z15_=invd7*Z15_;
   Z11_=Z15_ + Z11_ - invd9 - invd5;
   Z11_=invd6*Z11_;
   Z15_= - lm[0] + E2 + lm[8];
   Z17_=invd0*Z15_;
   Z15_=invd5*Z15_;
   Z17_=Z15_ + Z17_;
   Z17_=invd1*Z17_;
   Z12_=invd7*Z12_;
   Z12_=Z12_ - 1 + Z15_;
   Z12_=invd6*Z12_;
   Z13_=invd7*Z13_;
   Z12_=Z12_ + Z17_ - invd0 + Z13_;
   Z12_=invd4*Z12_;
   Z13_=invd11*invd12;
   Z15_=invd8*invd11;
   Z17_=Z13_ + Z15_;
   Z17_=invd9*Z17_;
   Z13_=invd10*Z13_;
   Z13_=Z13_ + Z17_;
   Z13_=E0*Z13_;
   Z17_=invd10*invd11;
   Z15_=Z17_ + Z15_;
   Z15_=E0*Z15_;
   Z15_= - invd0 - invd8 + Z15_;
   Z15_=invd2*Z15_;
   Z17_=E0*invd10;
   Z18_=E3 + lm[5];
   Z18_=invd0*Z18_;
   Z17_=Z18_ - 1 + Z17_;
   Z17_=invd2*Z17_;
   Z16_=E0*Z16_;
   Z16_=Z17_ - invd0 - invd5 + Z16_;
   Z16_=invd1*Z16_;
   Z10_=Z12_ + Z11_ + Z16_ + Z10_ + Z15_ + Z14_ + Z13_;


	return pow(2.*pi*I,2)/(2.*E0*2.*E1*2.*E2*2.*E3*pow(2.*E4,2))*(2*invd3*Z10_);
}
static double complex diag_1(double complex lm[], double complex params[]) {
	






	return pow(2.*pi*I,0)/(1)*(1);
}
static inline void evaluate_PF_1_0(double complex lm[], double complex params[], double complex* out) {
	double complex Z10_,Z11_,Z12_,Z13_,Z14_,Z15_,Z16_,Z17_,Z18_;


    Z10_=diag_1(lm, params);
    Z11_=diag_2(lm, params);
    Z12_=pow(xT,-1);
    Z13_=pow(xS,-1);
    Z14_=diag_3(lm, params);
    Z15_=diag_4(lm, params);
   Z16_=lm[16] - lm[9];
   Z17_=lm[21] + lm[25];
   Z17_= - lm[1] + 2*Z17_;
   Z16_= - lm[13] + Z17_ - 2*Z16_;
   Z16_=Z12_*Z16_;
   Z16_=1 + Z16_;
   Z16_=Z10_*Z16_;
   Z18_=Z13_*Z10_;
   Z17_= - lm[6] + Z17_ - 2*lm[23];
   Z17_=Z17_*Z18_;
   Z16_=Z16_ + Z17_;
   Z16_=Z15_*Z16_;
   Z17_=lm[5] - lm[0];
   Z17_=Z11_*Z17_;
   Z17_= - Z14_ + Z17_;
   Z17_=Z17_*Z18_;
   Z18_= - lm[0] + lm[8];
   Z18_=Z11_*Z18_;
   Z18_= - Z14_ + Z18_;
   Z18_=Z10_*Z12_*Z18_;
   Z17_=Z18_ + Z17_;


	*out = Z16_ + 2*Z17_;
}
void evaluate_PF_1(double complex lm[], double complex params[], int conf, double complex* out) {
   switch(conf) {
		case 0: evaluate_PF_1_0(lm, params, out); return;
		default: *out = 0.;
    }
}
