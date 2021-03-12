#include <tgmath.h>
# include <quadmath.h>
# include <signal.h>
# include "numerator.h"


static inline void evaluate_LTD_0_0_f128(__complex128 lm[], __complex128 params[], __complex128* out) {
	__complex128 Z10_,Z11_,Z12_,Z13_,Z14_,Z15_,Z16_,Z17_,Z18_,Z19_,Z20_,Z21_,Z22_,Z23_,Z24_,Z25_,Z26_,Z27_,Z28_,Z29_,Z30_,Z31_,Z32_,Z33_,Z34_,Z35_,Z36_,Z37_,Z38_,Z39_,Z40_,Z41_,Z42_,Z43_,Z44_,Z45_,Z6_,Z7_,Z8_,Z9_;
	__complex128 E0 = csqrtq(lm[25]+0.0q);
	__complex128 E1 = csqrtq(lm[25]+2*lm[21]+lm[2]+0.0q);
	__complex128 E2 = csqrtq(lm[25]+2*lm[21]-2*lm[16]+lm[14]-2*lm[10]+lm[2]+0.0q);
	__complex128 E3 = csqrtq(lm[25]-2*lm[23]+lm[7]+0.0q);
	__complex128 E4 = csqrtq(   lm[34]+0.0q);
	__complex128 E5 = csqrtq(lm[34]+2*lm[27]+lm[25]+0.0q);
	__complex128 invd0 = 1.q/(E0+E5-E4);
	__complex128 invd1 = 1.q/(-E0+E5-E4);
	__complex128 invd2 = 1.q/(lm[0]+E1+E5-E4);
	__complex128 invd3 = 1.q/(lm[0]-E1+E5-E4);
	__complex128 invd4 = 1.q/(-lm[8]+lm[0]+E2+E5-E4);
	__complex128 invd5 = 1.q/(-lm[8]+lm[0]-   E2+E5-E4);
	__complex128 invd6 = 1.q/(-lm[5]+E3+E5-E4);
	__complex128 invd7 = 1.q/(-lm[5]-E3+E5-E4);
	__complex128 invd8 = 1.q/(2*E4);
	__complex128 invd9 = 1.q/(2*E5);
	__complex128 invd10 = 1.q/(-2*E0);
	__complex128 invd11 = 1.q/(lm[0]+E1-E0);
	__complex128 invd12 = 1.q/(lm[0]-   E1-E0);
	__complex128 invd13 = 1.q/(-lm[8]+lm[0]+E2-E0);
	__complex128 invd14 = 1.q/(-lm[8]+lm[0]-E2-E0);
	__complex128 invd15 = 1.q/(-lm[5]+E3-E0);
	__complex128 invd16 = 1.q/(-lm[5]-E3-E0);
	__complex128 invd17 = 1.q/(E0+E5+E4);
	__complex128 invd18 = 1.q/(   E0+E5-E4);
	__complex128 invd19 = 1.q/(2*E5);
	__complex128 invd20 = 1.q/(-lm[0]-E1+E0);
	__complex128 invd21 = 1.q/(-lm[0]-E1-E0);
	__complex128 invd22 = 1.q/(-2*E1);
	__complex128 invd23 = 1.q/(-lm[8]+E2-E1);
	__complex128 invd24 = 1.q/(-lm[8]-E2-E1);
	__complex128 invd25 = 1.q/(-lm[5]-   lm[0]+E3-E1);
	__complex128 invd26 = 1.q/(-lm[5]-lm[0]-E3-E1);
	__complex128 invd27 = 1.q/(lm[0]+E1+E5+E4);
	__complex128 invd28 = 1.q/(lm[0]+E1+E5-E4);
	__complex128 invd29 = 1.q/(2*E5);
	__complex128 invd30 = 1.q/(lm[8]-lm[0]-E2+   E0);
	__complex128 invd31 = 1.q/(lm[8]-lm[0]-E2-E0);
	__complex128 invd32 = 1.q/(lm[8]-E2+E1);
	__complex128 invd33 = 1.q/(lm[8]-E2-E1);
	__complex128 invd34 = 1.q/(-2*E2);
	__complex128 invd35 = 1.q/(lm[8]-lm[5]-lm[0]+E3-E2);
	__complex128 invd36 = 1.q/(lm[8]-lm[5]-   lm[0]-E3-E2);
	__complex128 invd37 = 1.q/(-lm[8]+lm[0]+E2+E5+E4);
	__complex128 invd38 = 1.q/(-lm[8]+lm[0]+E2+E5-E4);
	__complex128 invd39 = 1.q/(2*E5);
	__complex128 invd40 = 1.q/(lm[5]-E3+E0);
	__complex128 invd41 = 1.q/(lm[5]-E3-   E0);
	__complex128 invd42 = 1.q/(lm[5]+lm[0]-E3+E1);
	__complex128 invd43 = 1.q/(lm[5]+lm[0]-E3-E1);
	__complex128 invd44 = 1.q/(-lm[8]+lm[5]+lm[0]-E3+E2);
	__complex128 invd45 = 1.q/(-lm[8]+lm[5]+lm[0]-E3-E2);
	__complex128 invd46 = 1.q/(   -2*E3);
	__complex128 invd47 = 1.q/(-lm[5]+E3+E5+E4);
	__complex128 invd48 = 1.q/(-lm[5]+E3+E5-E4);
	__complex128 invd49 = 1.q/(2*E5);
	__complex128 invd50 = 1.q/(2*E0);
	__complex128 invd51 = 1.q/(lm[0]+E1+E0);
	__complex128 invd52 = 1.q/(lm[0]-E1+E0);
	__complex128 invd53 = 1.q/(-lm[8]+   lm[0]+E2+E0);
	__complex128 invd54 = 1.q/(-lm[8]+lm[0]-E2+E0);
	__complex128 invd55 = 1.q/(-lm[5]+E3+E0);
	__complex128 invd56 = 1.q/(-lm[5]-E3+E0);
	__complex128 invd57 = 1.q/(2*E4);
	__complex128 invd58 = 1.q/(E0+E5+E4);
	__complex128 invd59 = 1.q/(E0-E5+   E4);
	__complex128 invd60 = 1.q/(-lm[0]+E1+E0);
	__complex128 invd61 = 1.q/(-lm[0]+E1-E0);
	__complex128 invd62 = 1.q/(2*E1);
	__complex128 invd63 = 1.q/(-lm[8]+E2+E1);
	__complex128 invd64 = 1.q/(-lm[8]-E2+E1);
	__complex128 invd65 = 1.q/(-lm[5]-lm[0]+E3+E1);
	__complex128 invd66 = 1.q/(-   lm[5]-lm[0]-E3+E1);
	__complex128 invd67 = 1.q/(2*E4);
	__complex128 invd68 = 1.q/(-lm[0]+E1+E5+E4);
	__complex128 invd69 = 1.q/(-lm[0]+E1-E5+E4);
	__complex128 invd70 = 1.q/(lm[8]-lm[0]+E2+E0);
	__complex128 invd71 = 1.q/(lm[8]-lm[0]+   E2-E0);
	__complex128 invd72 = 1.q/(lm[8]+E2+E1);
	__complex128 invd73 = 1.q/(lm[8]+E2-E1);
	__complex128 invd74 = 1.q/(2*E2);
	__complex128 invd75 = 1.q/(lm[8]-lm[5]-lm[0]+E3+E2);
	__complex128 invd76 = 1.q/(lm[8]-lm[5]-lm[0]-E3+E2);
	__complex128 invd77 = 1.q/(2*   E4);
	__complex128 invd78 = 1.q/(lm[8]-lm[0]+E2+E5+E4);
	__complex128 invd79 = 1.q/(lm[8]-lm[0]+E2-E5+E4);
	__complex128 invd80 = 1.q/(lm[5]+E3+E0);
	__complex128 invd81 = 1.q/(lm[5]+E3-E0);
	__complex128 invd82 = 1.q/(lm[5]+lm[0]+E3+   E1);
	__complex128 invd83 = 1.q/(lm[5]+lm[0]+E3-E1);
	__complex128 invd84 = 1.q/(-lm[8]+lm[5]+lm[0]+E3+E2);
	__complex128 invd85 = 1.q/(-lm[8]+lm[5]+lm[0]+E3-E2);
	__complex128 invd86 = 1.q/(2*E3);
	__complex128 invd87 = 1.q/(2*E4);
	__complex128 invd88 = 1.q/(lm[5]+   E3+E5+E4);
	__complex128 invd89 = 1.q/(lm[5]+E3-E5+E4);


    Z6_=cpowq(xT,-1);
    Z7_=cpowq(xS,-1);
   Z8_= - invd19*invd18*invd17*invd16*invd15*invd14*invd13*invd12*
   invd11*invd10*E0;
   Z9_= - invd59*invd58*invd57*invd56*invd55*invd54*invd53*invd52*
   invd51*invd50*E0;
   Z10_=invd29*invd28*invd27*invd26*invd25*invd24*invd23*invd22*invd21*
   invd20*E1;
   Z11_=invd69*invd68*invd67*invd66*invd65*invd64*invd63*invd62*invd61*
   invd60*E1;
   Z8_=Z11_ + Z10_ + Z8_ + Z9_;
   Z9_=invd29*invd28*invd27*invd26*invd25*invd24*invd23*invd22*invd21*
   invd20;
   Z12_= - invd69*invd68*invd67*invd66*invd65*invd64*invd63*invd62*
   invd61*invd60;
   Z9_=Z9_ + Z12_;
   Z12_=lm[8]*Z9_;
   Z13_=invd49*invd48*invd47*invd46*invd45*invd44*invd43*invd42*invd41*
   invd40;
   Z14_= - invd89*invd88*invd87*invd86*invd85*invd84*invd83*invd82*
   invd81*invd80;
   Z15_=Z13_ + Z14_;
   Z15_=lm[5]*Z15_;
   Z16_=invd39*invd38*invd37*invd36*invd35*invd34*invd33*invd32*invd31*
   invd30*E2;
   Z17_= - invd49*invd48*invd47*invd46*invd45*invd44*invd43*invd42*
   invd41*invd40*E3;
   Z18_=invd79*invd78*invd77*invd76*invd75*invd74*invd73*invd72*invd71*
   invd70*E2;
   Z19_= - invd89*invd88*invd87*invd86*invd85*invd84*invd83*invd82*
   invd81*invd80*E3;
   Z20_=E4 - E5;
   Z20_=invd9*invd8*invd7*invd6*invd5*invd4*invd3*invd2*invd1*invd0*
   Z20_;
   Z12_=Z15_ + Z20_ + Z19_ + Z18_ + Z17_ + Z16_ + Z8_ + Z12_;
   Z12_=lm[0]*Z12_;
   Z16_=cpowq(E2,2);
   Z18_=2*Z16_;
   Z21_=2*lm[16] - 2*lm[9] + lm[13];
   Z22_= - 2*lm[21];
   Z23_= - 2*lm[25];
   Z24_= - 2*lm[8]*E2;
   Z25_=Z24_ + Z23_ + Z22_ + lm[1] + Z21_ + Z18_;
   Z25_=invd39*invd38*invd37*invd36*invd35*invd34*invd33*invd32*invd31*
   invd30*Z25_;
   Z26_=cpowq(E3,2);
   Z27_=2*Z26_;
   Z28_=2*lm[8]*E3;
   Z29_=Z28_ + Z23_ + Z22_ + lm[1] + Z21_ + Z27_;
   Z29_=invd49*invd48*invd47*invd46*invd45*invd44*invd43*invd42*invd41*
   invd40*Z29_;
   Z16_= - 2*Z16_;
   Z30_= - 2*lm[16];
   Z31_=Z30_ + 2*lm[9] - lm[13];
   Z32_=2*lm[21];
   Z33_=2*lm[25];
   Z24_=Z24_ + Z33_ + Z32_ - lm[1] + Z31_ + Z16_;
   Z24_=invd79*invd78*invd77*invd76*invd75*invd74*invd73*invd72*invd71*
   invd70*Z24_;
   Z26_= - 2*Z26_;
   Z28_=Z28_ + Z33_ + Z32_ - lm[1] + Z31_ + Z26_;
   Z28_=invd89*invd88*invd87*invd86*invd85*invd84*invd83*invd82*invd81*
   invd80*Z28_;
   Z34_=cpowq(E4,2);
   Z35_=2*E4 - E5;
   Z35_=E5*Z35_;
   Z34_= - Z34_ + Z35_;
   Z35_=Z34_ + lm[9];
   Z36_= - E4 + E5;
   Z37_=lm[8]*Z36_;
   Z30_=2*Z37_ + Z33_ + Z32_ - lm[1] + Z30_ + 2*Z35_ - lm[13];
   Z30_=invd9*invd8*invd7*invd6*invd5*invd4*invd3*invd2*invd1*invd0*
   Z30_;
   Z35_=cpowq(E1,2);
   Z37_=2*Z35_;
   Z38_=Z23_ + Z22_ + lm[1] + Z21_ + Z37_;
   Z38_=invd29*invd28*invd27*invd26*invd25*invd24*invd23*invd22*invd21*
   invd20*Z38_;
   Z35_= - 2*Z35_;
   Z39_=Z33_ + Z32_ - lm[1] + Z31_ + Z35_;
   Z39_=invd69*invd68*invd67*invd66*invd65*invd64*invd63*invd62*invd61*
   invd60*Z39_;
   Z40_=invd19*invd18*invd17*invd16*invd15*invd14*invd13*invd12*invd11*
   invd10*E0;
   Z41_=invd59*invd58*invd57*invd56*invd55*invd54*invd53*invd52*invd51*
   invd50*E0;
   Z10_=Z11_ + Z10_ + Z40_ + Z41_;
   Z11_=lm[8]*Z10_;
   Z40_= - 2*E3;
   Z41_=Z40_ - lm[8];
   Z41_=invd49*invd48*invd47*invd46*invd45*invd44*invd43*invd42*invd41*
   invd40*Z41_;
   Z40_=Z40_ + lm[8];
   Z40_=invd89*invd88*invd87*invd86*invd85*invd84*invd83*invd82*invd81*
   invd80*Z40_;
   Z15_=Z15_ + Z41_ + Z40_;
   Z15_=lm[5]*Z15_;
   Z40_=cpowq(E0,2);
   Z41_=2*Z40_;
   Z21_=Z21_ + Z41_;
   Z21_=invd19*invd18*invd17*invd16*invd15*invd14*invd13*invd12*invd11*
   invd10*Z21_;
   Z40_= - 2*Z40_;
   Z31_=Z31_ + Z40_;
   Z31_=invd59*invd58*invd57*invd56*invd55*invd54*invd53*invd52*invd51*
   invd50*Z31_;
   Z42_=invd19*invd18*invd17*invd16*invd15*invd14*invd13*invd12*invd11*
   invd10;
   Z43_= - invd59*invd58*invd57*invd56*invd55*invd54*invd53*invd52*
   invd51*invd50;
   Z42_=Z42_ + Z43_;
   Z42_=lm[1]*Z42_;
   Z43_= - invd19*invd18*invd17*invd16*invd15*invd14*invd13*invd12*
   invd11*invd10;
   Z44_=invd59*invd58*invd57*invd56*invd55*invd54*invd53*invd52*invd51*
   invd50;
   Z43_=Z43_ + Z44_;
   Z44_=2*lm[21]*Z43_;
   Z45_=2*lm[25]*Z43_;
   Z11_=2*Z12_ + 2*Z15_ + Z30_ + Z28_ + Z24_ + Z29_ + Z25_ + 2*Z11_ + 
   Z39_ + Z38_ + Z45_ + Z44_ + Z42_ + Z21_ + Z31_;
   Z11_=Z6_*Z11_;
   Z12_=invd39*invd38*invd37*invd36*invd35*invd34*invd33*invd32*invd31*
   invd30;
   Z15_= - invd79*invd78*invd77*invd76*invd75*invd74*invd73*invd72*
   invd71*invd70;
   Z9_=Z14_ + Z15_ + Z13_ + Z9_ + Z12_;
   Z9_=lm[5]*Z9_;
   Z12_=E2 - lm[8];
   Z12_=invd39*invd38*invd37*invd36*invd35*invd34*invd33*invd32*invd31*
   invd30*Z12_;
   Z13_=E2 + lm[8];
   Z13_=invd79*invd78*invd77*invd76*invd75*invd74*invd73*invd72*invd71*
   invd70*Z13_;
   Z8_=Z9_ + Z20_ + Z19_ + Z13_ + Z17_ + Z8_ + Z12_;
   Z8_=lm[0]*Z8_;
   Z9_=invd9*invd8*invd7*invd6*invd5*invd4*invd3*invd2*invd1*invd0*Z36_
   ;
   Z9_=Z9_ + Z19_ + Z13_ + Z17_ + Z10_ + Z12_;
   Z9_=lm[5]*Z9_;
   Z10_= - 2*E2;
   Z12_=Z10_ + lm[8];
   Z12_=lm[8]*Z12_;
   Z13_=lm[6] + 2*lm[23];
   Z12_=2*Z12_ + Z23_ + Z22_ + lm[1] + Z13_ + Z18_;
   Z12_=invd39*invd38*invd37*invd36*invd35*invd34*invd33*invd32*invd31*
   invd30*Z12_;
   Z10_=Z10_ - lm[8];
   Z10_=lm[8]*Z10_;
   Z14_= - 2*lm[23];
   Z15_= - lm[6] + Z14_;
   Z10_=2*Z10_ + Z33_ + Z32_ - lm[1] + Z15_ + Z16_;
   Z10_=invd79*invd78*invd77*invd76*invd75*invd74*invd73*invd72*invd71*
   invd70*Z10_;
   Z16_=Z23_ + Z22_ + lm[1] + Z13_ + Z37_;
   Z16_=invd29*invd28*invd27*invd26*invd25*invd24*invd23*invd22*invd21*
   invd20*Z16_;
   Z17_=Z33_ + Z32_ - lm[1] + Z15_ + Z35_;
   Z17_=invd69*invd68*invd67*invd66*invd65*invd64*invd63*invd62*invd61*
   invd60*Z17_;
   Z18_=Z23_ + Z22_ + lm[1] + Z13_ + Z27_;
   Z18_=invd49*invd48*invd47*invd46*invd45*invd44*invd43*invd42*invd41*
   invd40*Z18_;
   Z19_=Z33_ + Z32_ - lm[1] + Z15_ + Z26_;
   Z19_=invd89*invd88*invd87*invd86*invd85*invd84*invd83*invd82*invd81*
   invd80*Z19_;
   Z14_=Z33_ + Z32_ - lm[1] + Z14_ + 2*Z34_ - lm[6];
   Z14_=invd9*invd8*invd7*invd6*invd5*invd4*invd3*invd2*invd1*invd0*
   Z14_;
   Z13_=Z13_ + Z41_;
   Z13_=invd19*invd18*invd17*invd16*invd15*invd14*invd13*invd12*invd11*
   invd10*Z13_;
   Z15_=Z15_ + Z40_;
   Z15_=invd59*invd58*invd57*invd56*invd55*invd54*invd53*invd52*invd51*
   invd50*Z15_;
   Z8_=2*Z8_ + 2*Z9_ + Z14_ + Z19_ + Z10_ + Z18_ + Z12_ + Z17_ + Z16_
    + Z45_ + Z44_ + Z42_ + Z13_ + Z15_;
   Z8_=Z7_*Z8_;
   Z9_= - invd29*invd28*invd27*invd26*invd25*invd24*invd23*invd22*
   invd21*invd20;
   Z10_=invd69*invd68*invd67*invd66*invd65*invd64*invd63*invd62*invd61*
   invd60;
   Z12_= - invd39*invd38*invd37*invd36*invd35*invd34*invd33*invd32*
   invd31*invd30;
   Z13_= - invd49*invd48*invd47*invd46*invd45*invd44*invd43*invd42*
   invd41*invd40;
   Z14_=invd79*invd78*invd77*invd76*invd75*invd74*invd73*invd72*invd71*
   invd70;
   Z15_=invd89*invd88*invd87*invd86*invd85*invd84*invd83*invd82*invd81*
   invd80;
   Z16_=invd9*invd8*invd7*invd6*invd5*invd4*invd3*invd2*invd1*invd0;


	*out = cpowq(2.q*M_PIq*I,2)/(1)*(Z8_ + Z9_ + Z10_ + Z11_ + Z12_ + Z13_ + Z14_ + Z15_ + Z16_ + 
      Z43_);
}
void evaluate_LTD_0_f128(__complex128 lm[], __complex128 params[], int conf, __complex128* out) {
   switch(conf) {
		case 0: evaluate_LTD_0_0_f128(lm, params, out); return;
		default: *out = 0.q;
    }
}
