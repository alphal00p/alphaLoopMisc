#include <tgmath.h>
# include <quadmath.h>
# include <signal.h>
# include "numerator.h"

static inline void evaluate_LTD_1_0(double complex lm[], double complex params[], double complex* out) {
	double complex Z10_,Z11_,Z12_,Z13_,Z14_,Z15_,Z16_,Z17_,Z18_,Z19_,Z20_,Z21_,Z22_,Z23_,Z24_,Z25_,Z26_,Z27_,Z28_,Z29_,Z30_,Z31_,Z32_,Z33_,Z34_,Z35_,Z36_,Z37_,Z38_,Z39_,Z40_,Z41_,Z42_,Z43_,Z44_,Z45_,Z46_,Z47_,Z6_,Z7_,Z8_,Z9_;
	double complex E0 = sqrt(lm[53]+1e-08);
	double complex E1 = sqrt(lm[53]-2*lm[51]+2*lm[44]+lm[42]-2*lm[40]+2*lm[33]+2*lm[31]+lm[29]-2*lm[27]+lm[7]+1e-08);
	double complex E2 = sqrt(lm[53]-   2*lm[51]+2*lm[44]+lm[42]-2*lm[40]+lm[7]+1e-08);
	double complex E3 = sqrt(lm[53]-2*lm[51]+lm[7]+1e-08);
	double complex E4 = sqrt(lm[62]+1e-08);
	double complex E5 = sqrt(lm[62]-2*lm[60]-2*   lm[58]+2*lm[46]+lm[42]-2*lm[40]-2*lm[38]+2*lm[35]+2*lm[31]+lm[29]-2*lm[27]-2*lm[25]+2*   lm[22]+2*lm[18]+2*lm[16]+lm[14]-2*lm[12]-2*lm[10]+lm[7]+2*lm[4]+lm[2]+1e-08);
	double complex invd0 = 1./(2*E0);
	double complex invd1 = 1./(lm[36]+lm[23]-lm[5]+E1+E0);
	double complex invd2 = 1./(lm[36]+lm[23]-lm[5]-E1+E0);
	double complex invd3 = 1./(lm[36]-lm[5]+E2+E0);
	double complex invd4 = 1./(lm[36]-lm[5]-E2+   E0);
	double complex invd5 = 1./(-lm[5]+E3+E0);
	double complex invd6 = 1./(-lm[5]-E3+E0);
	double complex invd7 = 1./(2*E4);
	double complex invd8 = 1./(lm[36]+lm[23]+lm[8]-lm[5]-lm[0]+E5+E4);
	double complex invd9 = 1./(lm[36]+lm[23]+   lm[8]-lm[5]-lm[0]-E5+E4);
	double complex invd10 = 1./(2*E0);
	double complex invd11 = 1./(lm[36]+lm[23]-lm[5]+E1+E0);
	double complex invd12 = 1./(lm[36]+lm[23]-lm[5]-E1+E0);
	double complex invd13 = 1./(lm[36]-   lm[5]+E2+E0);
	double complex invd14 = 1./(lm[36]-lm[5]-E2+E0);
	double complex invd15 = 1./(-lm[5]+E3+E0);
	double complex invd16 = 1./(-lm[5]-E3+E0);
	double complex invd17 = 1./(-lm[36]-lm[23]-lm[8]+lm[5]+   lm[0]+E5+E4);
	double complex invd18 = 1./(-lm[36]-lm[23]-lm[8]+lm[5]+lm[0]+E5-E4);
	double complex invd19 = 1./(2*E5);
	double complex invd20 = 1./(-lm[36]-lm[23]+lm[5]+E1+E0);
	double complex invd21 = 1./(-   lm[36]-lm[23]+lm[5]+E1-E0);
	double complex invd22 = 1./(2*E1);
	double complex invd23 = 1./(-lm[23]+E2+E1);
	double complex invd24 = 1./(-lm[23]-E2+E1);
	double complex invd25 = 1./(-lm[36]-lm[23]+E3+E1);
	double complex invd26 = 1./(-   lm[36]-lm[23]-E3+E1);
	double complex invd27 = 1./(2*E4);
	double complex invd28 = 1./(lm[36]+lm[23]+lm[8]-lm[5]-lm[0]+E5+E4);
	double complex invd29 = 1./(lm[36]+lm[23]+lm[8]-lm[5]-   lm[0]-E5+E4);
	double complex invd30 = 1./(-lm[36]-lm[23]+lm[5]+E1+E0);
	double complex invd31 = 1./(-lm[36]-lm[23]+lm[5]+E1-E0);
	double complex invd32 = 1./(2*E1);
	double complex invd33 = 1./(-lm[23]+E2+E1);
	double complex invd34 = 1./(   -lm[23]-E2+E1);
	double complex invd35 = 1./(-lm[36]-lm[23]+E3+E1);
	double complex invd36 = 1./(-lm[36]-lm[23]-E3+E1);
	double complex invd37 = 1./(-lm[36]-lm[23]-lm[8]+lm[5]+lm[0]+   E5+E4);
	double complex invd38 = 1./(-lm[36]-lm[23]-lm[8]+lm[5]+lm[0]+E5-E4);
	double complex invd39 = 1./(2*E5);
	double complex invd40 = 1./(-lm[36]+lm[5]+E2+E0);
	double complex invd41 = 1./(-lm[36]+lm[5]+E2-   E0);
	double complex invd42 = 1./(lm[23]+E2+E1);
	double complex invd43 = 1./(lm[23]+E2-E1);
	double complex invd44 = 1./(2*E2);
	double complex invd45 = 1./(-lm[36]+E3+E2);
	double complex invd46 = 1./(-lm[36]-E3+E2);
	double complex invd47 = 1./(2*E4);
	double complex invd48 = 1./(lm[36]+lm[23]+   lm[8]-lm[5]-lm[0]+E5+E4);
	double complex invd49 = 1./(lm[36]+lm[23]+lm[8]-lm[5]-lm[0]-E5+E4);
	double complex invd50 = 1./(-lm[36]+lm[5]+E2+E0);
	double complex invd51 = 1./(-lm[36]+   lm[5]+E2-E0);
	double complex invd52 = 1./(lm[23]+E2+E1);
	double complex invd53 = 1./(lm[23]+E2-E1);
	double complex invd54 = 1./(2*E2);
	double complex invd55 = 1./(-lm[36]+E3+E2);
	double complex invd56 = 1./(-lm[36]-E3+E2);
	double complex invd57 = 1./(-lm[36]-   lm[23]-lm[8]+lm[5]+lm[0]+E5+E4);
	double complex invd58 = 1./(-lm[36]-lm[23]-lm[8]+lm[5]+lm[0]+E5-E4);
	double complex invd59 = 1./(2*E5);
	double complex invd60 = 1./(lm[5]+E3+E0);
	double complex invd61 = 1./(   lm[5]+E3-E0);
	double complex invd62 = 1./(lm[36]+lm[23]+E3+E1);
	double complex invd63 = 1./(lm[36]+lm[23]+E3-E1);
	double complex invd64 = 1./(lm[36]+E3+E2);
	double complex invd65 = 1./(lm[36]+E3-E2);
	double complex invd66 = 1./(2*E3);
	double complex invd67 = 1./(   2*E4);
	double complex invd68 = 1./(lm[36]+lm[23]+lm[8]-lm[5]-lm[0]+E5+E4);
	double complex invd69 = 1./(lm[36]+lm[23]+lm[8]-lm[5]-lm[0]-E5+E4);
	double complex invd70 = 1./(lm[5]+E3+   E0);
	double complex invd71 = 1./(lm[5]+E3-E0);
	double complex invd72 = 1./(lm[36]+lm[23]+E3+E1);
	double complex invd73 = 1./(lm[36]+lm[23]+E3-E1);
	double complex invd74 = 1./(lm[36]+E3+E2);
	double complex invd75 = 1./(lm[36]+E3-E2);
	double complex invd76 = 1./(2*   E3);
	double complex invd77 = 1./(-lm[36]-lm[23]-lm[8]+lm[5]+lm[0]+E5+E4);
	double complex invd78 = 1./(-lm[36]-lm[23]-lm[8]+lm[5]+lm[0]+E5-E4);
	double complex invd79 = 1./(2*E5);


    Z6_=pow(xT,-1);
    Z7_=pow(xS,-1);
   Z8_= - invd69*invd68*invd67*invd66*invd65*invd64*invd63*invd62*
   invd61*invd60;
   Z9_= - invd79*invd78*invd77*invd76*invd75*invd74*invd73*invd72*
   invd71*invd70;
   Z8_=Z8_ + Z9_;
   Z9_= - invd49*invd48*invd47*invd46*invd45*invd44*invd43*invd42*
   invd41*invd40;
   Z10_= - invd59*invd58*invd57*invd56*invd55*invd54*invd53*invd52*
   invd51*invd50;
   Z11_= - invd29*invd28*invd27*invd26*invd25*invd24*invd23*invd22*
   invd21*invd20;
   Z12_= - invd39*invd38*invd37*invd36*invd35*invd34*invd33*invd32*
   invd31*invd30;
   Z13_=Z12_ + Z11_ + Z10_ + Z8_ + Z9_;
   Z13_=lm[5]*Z13_;
   Z14_=invd49*invd48*invd47*invd46*invd45*invd44*invd43*invd42*invd41*
   invd40;
   Z15_=invd59*invd58*invd57*invd56*invd55*invd54*invd53*invd52*invd51*
   invd50;
   Z16_=invd29*invd28*invd27*invd26*invd25*invd24*invd23*invd22*invd21*
   invd20;
   Z17_=invd39*invd38*invd37*invd36*invd35*invd34*invd33*invd32*invd31*
   invd30;
   Z18_=Z17_ + Z16_ + Z14_ + Z15_;
   Z18_=lm[36]*Z18_;
   Z19_=invd69*invd68*invd67*invd66*invd65*invd64*invd63*invd62*invd61*
   invd60;
   Z20_=invd79*invd78*invd77*invd76*invd75*invd74*invd73*invd72*invd71*
   invd70;
   Z21_=Z19_ + Z20_;
   Z22_=Z15_ + Z21_ + Z14_;
   Z22_=lm[23]*Z22_;
   Z23_= - invd69*invd68*invd67*invd66*invd65*invd64*invd63*invd62*
   invd61*invd60*E3;
   Z24_= - invd79*invd78*invd77*invd76*invd75*invd74*invd73*invd72*
   invd71*invd70*E3;
   Z25_=Z23_ + Z24_;
   Z8_=lm[0]*Z8_;
   Z26_= - 2*E2 - lm[0];
   Z27_=invd49*invd48*invd47*invd46*invd45*invd44*invd43*invd42*invd41*
   invd40*Z26_;
   Z26_=invd59*invd58*invd57*invd56*invd55*invd54*invd53*invd52*invd51*
   invd50*Z26_;
   Z28_=3*lm[23];
   Z29_=Z28_ - 2*E1 - lm[0];
   Z30_=invd29*invd28*invd27*invd26*invd25*invd24*invd23*invd22*invd21*
   invd20*Z29_;
   Z29_=invd39*invd38*invd37*invd36*invd35*invd34*invd33*invd32*invd31*
   invd30*Z29_;
   Z18_=Z13_ + 2*Z18_ + Z29_ + Z30_ + Z22_ + Z26_ + Z27_ + 2*Z25_ + Z8_
   ;
   Z18_=lm[5]*Z18_;
   Z9_=Z9_ + Z10_;
   Z10_=Z12_ + Z9_ + Z11_;
   Z10_=lm[36]*Z10_;
   Z9_=lm[23]*Z9_;
   Z11_=2*E2 + lm[0];
   Z12_=invd49*invd48*invd47*invd46*invd45*invd44*invd43*invd42*invd41*
   invd40*Z11_;
   Z11_=invd59*invd58*invd57*invd56*invd55*invd54*invd53*invd52*invd51*
   invd50*Z11_;
   Z26_= - 3*lm[23] + 2*E1 + lm[0];
   Z27_=invd29*invd28*invd27*invd26*invd25*invd24*invd23*invd22*invd21*
   invd20*Z26_;
   Z26_=invd39*invd38*invd37*invd36*invd35*invd34*invd33*invd32*invd31*
   invd30*Z26_;
   Z11_=Z10_ + Z26_ + Z27_ + Z9_ + Z12_ + Z11_;
   Z11_=lm[36]*Z11_;
   Z12_=2*lm[49];
   Z26_= - 2*pow(E0,2);
   Z27_=2*lm[24];
   Z29_= - 2*lm[33];
   Z30_=2*lm[53];
   Z31_=Z30_ + Z12_ + Z29_ - lm[28] + Z27_ + Z26_ - lm[1];
   Z32_=invd9*invd8*invd7*invd6*invd5*invd4*invd3*invd2*invd1*invd0*
   Z31_;
   Z31_=invd19*invd18*invd17*invd16*invd15*invd14*invd13*invd12*invd11*
   invd10*Z31_;
   Z12_=Z30_ + Z12_ + Z29_ - lm[28] - lm[1] + Z27_;
   Z33_= - 2*pow(E3,2);
   Z34_=Z12_ + Z33_;
   Z35_=invd69*invd68*invd67*invd66*invd65*invd64*invd63*invd62*invd61*
   invd60*Z34_;
   Z34_=invd79*invd78*invd77*invd76*invd75*invd74*invd73*invd72*invd71*
   invd70*Z34_;
   Z36_= - 2*pow(E2,2);
   Z37_= - lm[0]*E2;
   Z38_=2*Z37_ + Z12_ + Z36_;
   Z39_=invd49*invd48*invd47*invd46*invd45*invd44*invd43*invd42*invd41*
   invd40*Z38_;
   Z38_=invd59*invd58*invd57*invd56*invd55*invd54*invd53*invd52*invd51*
   invd50*Z38_;
   Z40_= - 2*lm[23];
   Z41_=3*E1;
   Z42_=Z40_ + Z41_ + lm[0];
   Z42_=lm[23]*Z42_;
   Z43_= - 2*pow(E1,2);
   Z44_= - lm[0]*E1;
   Z12_=2*Z42_ + 2*Z44_ + Z12_ + Z43_;
   Z42_=invd29*invd28*invd27*invd26*invd25*invd24*invd23*invd22*invd21*
   invd20*Z12_;
   Z12_=invd39*invd38*invd37*invd36*invd35*invd34*invd33*invd32*invd31*
   invd30*Z12_;
   Z45_= - invd9*invd8*invd7*invd6*invd5*invd4*invd3*invd2*invd1*invd0*
   E0;
   Z46_= - invd19*invd18*invd17*invd16*invd15*invd14*invd13*invd12*
   invd11*invd10*E0;
   Z23_=Z24_ + Z23_ + Z45_ + Z46_;
   Z23_=lm[0]*Z23_;
   Z24_=invd9*invd8*invd7*invd6*invd5*invd4*invd3*invd2*invd1*invd0*E0;
   Z45_=invd19*invd18*invd17*invd16*invd15*invd14*invd13*invd12*invd11*
   invd10*E0;
   Z46_=invd69*invd68*invd67*invd66*invd65*invd64*invd63*invd62*invd61*
   invd60*E3;
   Z47_=invd79*invd78*invd77*invd76*invd75*invd74*invd73*invd72*invd71*
   invd70*E3;
   Z24_=Z47_ + Z46_ + Z24_ + Z45_;
   Z45_=invd49*invd48*invd47*invd46*invd45*invd44*invd43*invd42*invd41*
   invd40*E2;
   Z46_=invd59*invd58*invd57*invd56*invd55*invd54*invd53*invd52*invd51*
   invd50*E2;
   Z45_=Z46_ + Z24_ + Z45_;
   Z45_=2*lm[23]*Z45_;
   Z11_=2*Z18_ + 2*Z11_ + Z12_ + Z42_ + Z45_ + Z38_ + Z39_ + 2*Z23_ + 
   Z34_ + Z35_ + Z32_ + Z31_;
   Z11_=Z6_*Z11_;
   Z12_= - lm[30] + lm[37];
   Z12_= - 2*lm[44] + 2*Z12_ - lm[41];
   Z18_= - 2*lm[1];
   Z31_=4*lm[49];
   Z26_=Z30_ + Z31_ + Z29_ - lm[28] + Z27_ + Z18_ + Z12_ + Z26_;
   Z32_=invd9*invd8*invd7*invd6*invd5*invd4*invd3*invd2*invd1*invd0*
   Z26_;
   Z26_=invd19*invd18*invd17*invd16*invd15*invd14*invd13*invd12*invd11*
   invd10*Z26_;
   Z34_=2*lm[0];
   Z35_=Z41_ + Z34_;
   Z38_=Z35_ - 4*lm[23];
   Z39_=invd29*invd28*invd27*invd26*invd25*invd24*invd23*invd22*invd21*
   invd20*Z38_;
   Z38_=invd39*invd38*invd37*invd36*invd35*invd34*invd33*invd32*invd31*
   invd30*Z38_;
   Z34_=3*E2 + Z34_;
   Z41_=invd49*invd48*invd47*invd46*invd45*invd44*invd43*invd42*invd41*
   invd40*Z34_;
   Z34_=invd59*invd58*invd57*invd56*invd55*invd54*invd53*invd52*invd51*
   invd50*Z34_;
   Z9_=2*Z10_ + Z38_ + Z39_ + Z9_ + Z34_ + Z24_ + Z41_;
   Z9_=lm[36]*Z9_;
   Z10_=3*Z17_ + 3*Z16_ + 3*Z15_ + Z21_ + 3*Z14_;
   Z10_=lm[36]*Z10_;
   Z21_= - E2 - lm[0];
   Z24_=invd49*invd48*invd47*invd46*invd45*invd44*invd43*invd42*invd41*
   invd40*Z21_;
   Z21_=invd59*invd58*invd57*invd56*invd55*invd54*invd53*invd52*invd51*
   invd50*Z21_;
   Z8_=Z21_ + Z24_ + Z25_ + Z8_;
   Z21_= - E1 - lm[0];
   Z21_=2*Z21_ + Z28_;
   Z24_=invd29*invd28*invd27*invd26*invd25*invd24*invd23*invd22*invd21*
   invd20*Z21_;
   Z21_=invd39*invd38*invd37*invd36*invd35*invd34*invd33*invd32*invd31*
   invd30*Z21_;
   Z8_=Z13_ + Z10_ + Z21_ + Z24_ + 2*Z8_ + Z22_;
   Z8_=lm[5]*Z8_;
   Z10_=Z30_ + Z31_ + Z29_ - lm[28] + Z27_ + Z12_ + Z18_;
   Z12_=Z10_ + Z33_;
   Z13_=invd69*invd68*invd67*invd66*invd65*invd64*invd63*invd62*invd61*
   invd60*Z12_;
   Z12_=invd79*invd78*invd77*invd76*invd75*invd74*invd73*invd72*invd71*
   invd70*Z12_;
   Z18_=4*Z37_ + Z10_ + Z36_;
   Z21_=invd49*invd48*invd47*invd46*invd45*invd44*invd43*invd42*invd41*
   invd40*Z18_;
   Z18_=invd59*invd58*invd57*invd56*invd55*invd54*invd53*invd52*invd51*
   invd50*Z18_;
   Z22_=Z35_ + Z40_;
   Z22_=lm[23]*Z22_;
   Z10_=2*Z22_ + 4*Z44_ + Z10_ + Z43_;
   Z22_=invd29*invd28*invd27*invd26*invd25*invd24*invd23*invd22*invd21*
   invd20*Z10_;
   Z10_=invd39*invd38*invd37*invd36*invd35*invd34*invd33*invd32*invd31*
   invd30*Z10_;
   Z8_=2*Z8_ + 2*Z9_ + Z10_ + Z22_ + Z45_ + Z18_ + Z21_ + 4*Z23_ + Z12_
    + Z13_ + Z32_ + Z26_;
   Z8_=Z7_*Z8_;
   Z9_=invd9*invd8*invd7*invd6*invd5*invd4*invd3*invd2*invd1*invd0;
   Z10_=invd19*invd18*invd17*invd16*invd15*invd14*invd13*invd12*invd11*
   invd10;


	*out = pow(2.*pi*I,2)/(1)*(Z8_ + Z9_ + Z10_ + Z11_ + Z14_ + Z15_ + Z16_ + Z17_ + Z19_ + 
      Z20_);
}
void evaluate_LTD_1(double complex lm[], double complex params[], int conf, double complex* out) {
   switch(conf) {
		case 0: evaluate_LTD_1_0(lm, params, out); return;
		default: *out = 0.;
    }
}
