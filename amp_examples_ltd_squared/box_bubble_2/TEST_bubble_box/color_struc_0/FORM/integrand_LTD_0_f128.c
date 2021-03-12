#include <tgmath.h>
# include <quadmath.h>
# include <signal.h>
# include "numerator.h"


static inline void evaluate_LTD_0_0_f128(__complex128 lm[], __complex128 params[], __complex128* out) {
	__complex128 Z10_,Z11_,Z12_,Z13_,Z14_,Z15_,Z16_,Z17_,Z18_,Z19_,Z20_,Z21_,Z22_,Z23_,Z24_,Z25_,Z26_,Z27_,Z28_,Z29_,Z30_,Z31_,Z32_,Z33_,Z34_,Z35_,Z36_,Z37_,Z38_,Z39_,Z40_,Z41_,Z42_,Z43_,Z44_,Z45_,Z46_,Z47_,Z48_,Z49_,Z50_,Z51_,Z52_,Z53_,Z54_,Z55_,Z56_,Z57_,Z58_,Z59_,Z60_,Z61_,Z62_,Z63_,Z64_,Z65_,Z66_,Z6_,Z7_,Z8_,Z9_;
	__complex128 E0 = csqrtq(lm[53]+1e-08);
	__complex128 E1 = csqrtq(lm[53]-2*lm[51]+2*lm[44]+lm[42]-2*lm[40]+2*lm[33]+2*lm[31]+lm[29]-2*lm[27]+lm[7]+1e-08);
	__complex128 E2 = csqrtq(lm[53]-   2*lm[51]+2*lm[44]+lm[42]-2*lm[40]+lm[7]+1e-08);
	__complex128 E3 = csqrtq(lm[53]-2*lm[51]+lm[7]+1e-08);
	__complex128 E4 = csqrtq(lm[62]+1e-08);
	__complex128 E5 = csqrtq(lm[62]-2*lm[60]-2*   lm[58]+2*lm[55]+lm[53]-2*lm[51]-2*lm[49]+2*lm[46]+2*lm[44]+lm[42]-2*lm[40]-2*lm[38]+2*   lm[35]+2*lm[33]+2*lm[31]+lm[29]-2*lm[27]-2*lm[25]+lm[7]+2*lm[4]+lm[2]+1e-08);
	__complex128 invd0 = 1.q/(-lm[36]-lm[23]+lm[5]+lm[0]+E0+E5-E4);
	__complex128 invd1 = 1.q/(-lm[36]-lm[23]+lm[5]+lm[0]-E0+E5-E4);
	__complex128 invd2 = 1.q/(lm[0]+E1+E5-E4);
	__complex128 invd3 = 1.q/(   lm[0]-E1+E5-E4);
	__complex128 invd4 = 1.q/(-lm[23]+lm[0]+E2+E5-E4);
	__complex128 invd5 = 1.q/(-lm[23]+lm[0]-E2+E5-E4);
	__complex128 invd6 = 1.q/(-lm[36]-lm[23]+lm[0]+E3+   E5-E4);
	__complex128 invd7 = 1.q/(-lm[36]-lm[23]+lm[0]-E3+E5-E4);
	__complex128 invd8 = 1.q/(2*E4);
	__complex128 invd9 = 1.q/(2*E5);
	__complex128 invd10 = 1.q/(-2*E0);
	__complex128 invd11 = 1.q/(lm[36]+lm[23]-lm[5]+E1-E0);
	__complex128 invd12 = 1.q/(   lm[36]+lm[23]-lm[5]-E1-E0);
	__complex128 invd13 = 1.q/(lm[36]-lm[5]+E2-E0);
	__complex128 invd14 = 1.q/(lm[36]-lm[5]-E2-E0);
	__complex128 invd15 = 1.q/(-lm[5]+E3-E0);
	__complex128 invd16 = 1.q/(-lm[5]-E3-   E0);
	__complex128 invd17 = 1.q/(-lm[36]-lm[23]+lm[5]+lm[0]+E0+E5+E4);
	__complex128 invd18 = 1.q/(-lm[36]-lm[23]+lm[5]+lm[0]+E0+E5-E4);
	__complex128 invd19 = 1.q/(2*E5);
	__complex128 invd20 = 1.q/(-   lm[36]-lm[23]+lm[5]-E1+E0);
	__complex128 invd21 = 1.q/(-lm[36]-lm[23]+lm[5]-E1-E0);
	__complex128 invd22 = 1.q/(-2*E1);
	__complex128 invd23 = 1.q/(-lm[23]+E2-E1);
	__complex128 invd24 = 1.q/(-lm[23]-E2-   E1);
	__complex128 invd25 = 1.q/(-lm[36]-lm[23]+E3-E1);
	__complex128 invd26 = 1.q/(-lm[36]-lm[23]-E3-E1);
	__complex128 invd27 = 1.q/(lm[0]+E1+E5+E4);
	__complex128 invd28 = 1.q/(lm[0]+E1+E5-E4);
	__complex128 invd29 = 1.q/(2*E5);
	__complex128 invd30 = 1.q/(-   lm[36]+lm[5]-E2+E0);
	__complex128 invd31 = 1.q/(-lm[36]+lm[5]-E2-E0);
	__complex128 invd32 = 1.q/(lm[23]-E2+E1);
	__complex128 invd33 = 1.q/(lm[23]-E2-E1);
	__complex128 invd34 = 1.q/(-2*E2);
	__complex128 invd35 = 1.q/(-lm[36]+E3-   E2);
	__complex128 invd36 = 1.q/(-lm[36]-E3-E2);
	__complex128 invd37 = 1.q/(-lm[23]+lm[0]+E2+E5+E4);
	__complex128 invd38 = 1.q/(-lm[23]+lm[0]+E2+E5-E4);
	__complex128 invd39 = 1.q/(2*E5);
	__complex128 invd40 = 1.q/(lm[5]-E3+E0);
	__complex128 invd41 = 1.q/(   lm[5]-E3-E0);
	__complex128 invd42 = 1.q/(lm[36]+lm[23]-E3+E1);
	__complex128 invd43 = 1.q/(lm[36]+lm[23]-E3-E1);
	__complex128 invd44 = 1.q/(lm[36]-E3+E2);
	__complex128 invd45 = 1.q/(lm[36]-E3-E2);
	__complex128 invd46 = 1.q/(-2*   E3);
	__complex128 invd47 = 1.q/(-lm[36]-lm[23]+lm[0]+E3+E5+E4);
	__complex128 invd48 = 1.q/(-lm[36]-lm[23]+lm[0]+E3+E5-E4);
	__complex128 invd49 = 1.q/(2*E5);
	__complex128 invd50 = 1.q/(2*E0);
	__complex128 invd51 = 1.q/(lm[36]+   lm[23]-lm[5]+E1+E0);
	__complex128 invd52 = 1.q/(lm[36]+lm[23]-lm[5]-E1+E0);
	__complex128 invd53 = 1.q/(lm[36]-lm[5]+E2+E0);
	__complex128 invd54 = 1.q/(lm[36]-lm[5]-E2+E0);
	__complex128 invd55 = 1.q/(-   lm[5]+E3+E0);
	__complex128 invd56 = 1.q/(-lm[5]-E3+E0);
	__complex128 invd57 = 1.q/(2*E4);
	__complex128 invd58 = 1.q/(lm[36]+lm[23]-lm[5]-lm[0]+E0+E5+E4);
	__complex128 invd59 = 1.q/(lm[36]+lm[23]-lm[5]-   lm[0]+E0-E5+E4);
	__complex128 invd60 = 1.q/(-lm[36]-lm[23]+lm[5]+E1+E0);
	__complex128 invd61 = 1.q/(-lm[36]-lm[23]+lm[5]+E1-E0);
	__complex128 invd62 = 1.q/(2*E1);
	__complex128 invd63 = 1.q/(-lm[23]+E2+   E1);
	__complex128 invd64 = 1.q/(-lm[23]-E2+E1);
	__complex128 invd65 = 1.q/(-lm[36]-lm[23]+E3+E1);
	__complex128 invd66 = 1.q/(-lm[36]-lm[23]-E3+E1);
	__complex128 invd67 = 1.q/(2*E4);
	__complex128 invd68 = 1.q/(-lm[0]+E1+E5+E4);
	__complex128 invd69 = 1.q/(-   lm[0]+E1-E5+E4);
	__complex128 invd70 = 1.q/(-lm[36]+lm[5]+E2+E0);
	__complex128 invd71 = 1.q/(-lm[36]+lm[5]+E2-E0);
	__complex128 invd72 = 1.q/(lm[23]+E2+E1);
	__complex128 invd73 = 1.q/(lm[23]+E2-E1);
	__complex128 invd74 = 1.q/(2*   E2);
	__complex128 invd75 = 1.q/(-lm[36]+E3+E2);
	__complex128 invd76 = 1.q/(-lm[36]-E3+E2);
	__complex128 invd77 = 1.q/(2*E4);
	__complex128 invd78 = 1.q/(lm[23]-lm[0]+E2+E5+E4);
	__complex128 invd79 = 1.q/(lm[23]-lm[0]+E2-E5+E4);
	__complex128 invd80 = 1.q/(   lm[5]+E3+E0);
	__complex128 invd81 = 1.q/(lm[5]+E3-E0);
	__complex128 invd82 = 1.q/(lm[36]+lm[23]+E3+E1);
	__complex128 invd83 = 1.q/(lm[36]+lm[23]+E3-E1);
	__complex128 invd84 = 1.q/(lm[36]+E3+E2);
	__complex128 invd85 = 1.q/(lm[36]+   E3-E2);
	__complex128 invd86 = 1.q/(2*E3);
	__complex128 invd87 = 1.q/(2*E4);
	__complex128 invd88 = 1.q/(lm[36]+lm[23]-lm[0]+E3+E5+E4);
	__complex128 invd89 = 1.q/(lm[36]+lm[23]-lm[0]+E3-E5+E4);


    Z6_=cpowq(xT,-1);
    Z7_=cpowq(xS,-1);
   Z8_= - 3*lm[0];
   Z9_=E4 - E5;
   Z10_=4*Z9_ + Z8_;
   Z10_=lm[0]*Z10_;
   Z11_= - 2*lm[23];
   Z12_= - E4 + E5;
   Z13_=3*Z12_;
   Z14_=Z13_ + 5*lm[0];
   Z15_=Z14_ + Z11_;
   Z15_=lm[23]*Z15_;
   Z16_= - lm[30] + lm[37];
   Z16_= - 2*lm[44] + 2*Z16_ - lm[41];
   Z17_=cpowq(E4,2);
   Z18_=2*E4 - E5;
   Z18_=E5*Z18_;
   Z19_= - 2*lm[1];
   Z20_=2*lm[24];
   Z21_= - 2*lm[33];
   Z22_=4*lm[49];
   Z23_=2*lm[53];
   Z14_= - 2*lm[36] + Z14_ - 4*lm[23];
   Z14_=lm[36]*Z14_;
   Z10_=2*Z14_ + 2*Z15_ + 2*Z10_ + Z23_ + Z22_ + Z21_ - lm[28] + Z20_ + 
   Z19_ + 2*Z18_ + Z16_ - 2*Z17_;
   Z10_=invd9*invd8*invd7*invd6*invd5*invd4*invd3*invd2*invd1*invd0*
   Z10_;
   Z14_= - invd19*invd18*invd17*invd16*invd15*invd14*invd13*invd12*
   invd11*invd10*E0;
   Z15_= - invd59*invd58*invd57*invd56*invd55*invd54*invd53*invd52*
   invd51*invd50*E0;
   Z24_= - invd49*invd48*invd47*invd46*invd45*invd44*invd43*invd42*
   invd41*invd40*E3;
   Z25_= - invd89*invd88*invd87*invd86*invd85*invd84*invd83*invd82*
   invd81*invd80*E3;
   Z26_= - invd39*invd38*invd37*invd36*invd35*invd34*invd33*invd32*
   invd31*invd30*E2;
   Z27_= - invd79*invd78*invd77*invd76*invd75*invd74*invd73*invd72*
   invd71*invd70*E2;
   Z28_= - invd29*invd28*invd27*invd26*invd25*invd24*invd23*invd22*
   invd21*invd20*E1;
   Z29_= - invd69*invd68*invd67*invd66*invd65*invd64*invd63*invd62*
   invd61*invd60*E1;
   Z14_=Z29_ + Z28_ + Z27_ + Z26_ + Z25_ + Z24_ + Z14_ + Z15_;
   Z14_=lm[0]*Z14_;
   Z15_=invd19*invd18*invd17*invd16*invd15*invd14*invd13*invd12*invd11*
   invd10*E0;
   Z30_=invd59*invd58*invd57*invd56*invd55*invd54*invd53*invd52*invd51*
   invd50*E0;
   Z31_=invd49*invd48*invd47*invd46*invd45*invd44*invd43*invd42*invd41*
   invd40*E3;
   Z32_=invd89*invd88*invd87*invd86*invd85*invd84*invd83*invd82*invd81*
   invd80*E3;
   Z15_=Z32_ + Z31_ + Z15_ + Z30_;
   Z30_= - invd39*invd38*invd37*invd36*invd35*invd34*invd33*invd32*
   invd31*invd30;
   Z31_=invd79*invd78*invd77*invd76*invd75*invd74*invd73*invd72*invd71*
   invd70;
   Z32_= - invd29*invd28*invd27*invd26*invd25*invd24*invd23*invd22*
   invd21*invd20;
   Z33_=invd69*invd68*invd67*invd66*invd65*invd64*invd63*invd62*invd61*
   invd60;
   Z34_=Z33_ + Z32_ + Z30_ + Z31_;
   Z35_=lm[0]*Z34_;
   Z36_=invd39*invd38*invd37*invd36*invd35*invd34*invd33*invd32*invd31*
   invd30;
   Z37_= - invd79*invd78*invd77*invd76*invd75*invd74*invd73*invd72*
   invd71*invd70;
   Z38_=Z36_ + Z37_;
   Z39_=invd29*invd28*invd27*invd26*invd25*invd24*invd23*invd22*invd21*
   invd20;
   Z40_= - invd69*invd68*invd67*invd66*invd65*invd64*invd63*invd62*
   invd61*invd60;
   Z41_=4*Z40_ + Z38_ + 4*Z39_;
   Z41_=lm[23]*Z41_;
   Z42_=Z40_ + Z38_ + Z39_;
   Z42_=lm[36]*Z42_;
   Z43_=invd39*invd38*invd37*invd36*invd35*invd34*invd33*invd32*invd31*
   invd30*E2;
   Z44_=invd79*invd78*invd77*invd76*invd75*invd74*invd73*invd72*invd71*
   invd70*E2;
   Z45_=invd29*invd28*invd27*invd26*invd25*invd24*invd23*invd22*invd21*
   invd20*E1;
   Z46_=3*Z45_;
   Z47_=invd69*invd68*invd67*invd66*invd65*invd64*invd63*invd62*invd61*
   invd60*E1;
   Z48_=3*Z47_;
   Z41_=2*Z42_ + Z41_ + 2*Z35_ + Z48_ + Z46_ + 3*Z44_ + Z15_ + 3*Z43_;
   Z41_=lm[36]*Z41_;
   Z24_=Z29_ + Z28_ + Z27_ + Z26_ + Z24_ + Z25_;
   Z25_=invd49*invd48*invd47*invd46*invd45*invd44*invd43*invd42*invd41*
   invd40;
   Z26_= - invd89*invd88*invd87*invd86*invd85*invd84*invd83*invd82*
   invd81*invd80;
   Z25_=Z40_ + Z39_ + Z37_ + Z36_ + Z25_ + Z26_;
   Z26_=lm[0]*Z25_;
   Z27_=Z24_ + Z26_;
   Z28_= - invd9*invd8*invd7*invd6*invd5*invd4*invd3*invd2*invd1*invd0;
   Z25_=Z25_ + Z28_;
   Z25_=lm[5]*Z25_;
   Z28_= - invd49*invd48*invd47*invd46*invd45*invd44*invd43*invd42*
   invd41*invd40;
   Z29_=invd89*invd88*invd87*invd86*invd85*invd84*invd83*invd82*invd81*
   invd80;
   Z36_=Z28_ + Z29_;
   Z37_=3*Z32_;
   Z49_=3*Z33_;
   Z50_=Z49_ + Z37_ + Z31_ + Z36_ + Z30_;
   Z50_=lm[23]*Z50_;
   Z36_=Z49_ + Z37_ + 3*Z31_ + Z36_ + 3*Z30_;
   Z36_=lm[36]*Z36_;
   Z37_= - 2*lm[0];
   Z49_=Z9_ + Z37_;
   Z51_=3*lm[23];
   Z49_=3*lm[36] + 2*Z49_ + Z51_;
   Z49_=invd9*invd8*invd7*invd6*invd5*invd4*invd3*invd2*invd1*invd0*
   Z49_;
   Z27_=Z25_ + Z49_ + Z36_ + 2*Z27_ + Z50_;
   Z27_=lm[5]*Z27_;
   Z36_=cpowq(E0,2);
   Z49_=2*Z36_;
   Z52_=lm[30] - lm[37];
   Z52_=2*lm[44] + 2*Z52_ + lm[41];
   Z53_=2*lm[1];
   Z54_= - 2*lm[24];
   Z55_=2*lm[33];
   Z56_= - 4*lm[49];
   Z57_= - 2*lm[53];
   Z58_=Z57_ + Z56_ + Z55_ + lm[28] + Z54_ + Z53_ + Z52_ + Z49_;
   Z58_=invd19*invd18*invd17*invd16*invd15*invd14*invd13*invd12*invd11*
   invd10*Z58_;
   Z36_= - 2*Z36_;
   Z59_=Z23_ + Z22_ + Z21_ - lm[28] + Z20_ + Z19_ + Z16_ + Z36_;
   Z59_=invd59*invd58*invd57*invd56*invd55*invd54*invd53*invd52*invd51*
   invd50*Z59_;
   Z60_=cpowq(E3,2);
   Z61_=2*Z60_;
   Z62_=Z57_ + Z56_ + Z55_ + lm[28] + Z54_ + Z53_ + Z52_ + Z61_;
   Z62_=invd49*invd48*invd47*invd46*invd45*invd44*invd43*invd42*invd41*
   invd40*Z62_;
   Z60_= - 2*Z60_;
   Z63_=Z23_ + Z22_ + Z21_ - lm[28] + Z20_ + Z19_ + Z16_ + Z60_;
   Z63_=invd89*invd88*invd87*invd86*invd85*invd84*invd83*invd82*invd81*
   invd80*Z63_;
   Z52_=Z57_ + Z56_ + Z55_ + lm[28] + Z54_ + Z52_ + Z53_;
   Z53_=cpowq(E2,2);
   Z56_=2*Z53_;
   Z64_=Z52_ + Z56_;
   Z64_=invd39*invd38*invd37*invd36*invd35*invd34*invd33*invd32*invd31*
   invd30*Z64_;
   Z16_=Z23_ + Z22_ + Z21_ - lm[28] + Z20_ + Z16_ + Z19_;
   Z19_= - 2*Z53_;
   Z22_=Z16_ + Z19_;
   Z22_=invd79*invd78*invd77*invd76*invd75*invd74*invd73*invd72*invd71*
   invd70*Z22_;
   Z53_=cpowq(E1,2);
   Z65_=2*Z53_;
   Z52_=Z52_ + Z65_;
   Z52_=invd29*invd28*invd27*invd26*invd25*invd24*invd23*invd22*invd21*
   invd20*Z52_;
   Z53_= - 2*Z53_;
   Z16_=Z16_ + Z53_;
   Z16_=invd69*invd68*invd67*invd66*invd65*invd64*invd63*invd62*invd61*
   invd60*Z16_;
   Z15_=Z48_ + Z46_ + Z44_ + Z15_ + Z43_;
   Z46_=Z32_ + Z33_;
   Z46_=lm[0]*Z46_;
   Z48_=Z39_ + Z40_;
   Z48_=2*lm[23]*Z48_;
   Z66_=Z48_ + Z15_ + 2*Z46_;
   Z66_=lm[23]*Z66_;
   Z10_=2*Z27_ + Z10_ + 2*Z41_ + 2*Z66_ + 4*Z14_ + Z16_ + Z52_ + Z22_
    + Z64_ + Z63_ + Z62_ + Z58_ + Z59_;
   Z10_=Z7_*Z10_;
   Z16_= - Z17_ + Z18_;
   Z17_=3*Z9_ + Z37_;
   Z17_=lm[0]*Z17_;
   Z11_=Z11_ + Z13_ + 4*lm[0];
   Z11_=lm[23]*Z11_;
   Z13_=2*lm[49];
   Z12_= - lm[36] - 3*lm[23] + 2*Z12_ + 3*lm[0];
   Z12_=lm[36]*Z12_;
   Z11_=2*Z12_ + 2*Z11_ + 2*Z17_ + Z23_ + Z13_ + Z21_ - lm[28] + Z20_ + 2
   *Z16_ - lm[1];
   Z11_=invd9*invd8*invd7*invd6*invd5*invd4*invd3*invd2*invd1*invd0*
   Z11_;
   Z12_=lm[36]*Z34_;
   Z8_=2*lm[36] + Z51_ + 2*Z9_ + Z8_;
   Z8_=invd9*invd8*invd7*invd6*invd5*invd4*invd3*invd2*invd1*invd0*Z8_;
   Z8_=Z25_ + Z8_ + 2*Z12_ + Z50_ + 2*Z24_ + Z26_;
   Z8_=lm[5]*Z8_;
   Z9_=Z47_ + Z45_ + Z43_ + Z44_;
   Z12_=3*Z40_ + Z38_ + 3*Z39_;
   Z12_=lm[23]*Z12_;
   Z9_=Z42_ + Z12_ + 2*Z9_ + Z35_;
   Z9_=lm[36]*Z9_;
   Z12_= - 2*lm[49];
   Z16_=Z57_ + Z12_ + Z55_ + lm[28] + Z54_ + Z49_ + lm[1];
   Z16_=invd19*invd18*invd17*invd16*invd15*invd14*invd13*invd12*invd11*
   invd10*Z16_;
   Z17_=Z23_ + Z13_ + Z21_ - lm[28] + Z20_ + Z36_ - lm[1];
   Z17_=invd59*invd58*invd57*invd56*invd55*invd54*invd53*invd52*invd51*
   invd50*Z17_;
   Z18_=Z57_ + Z12_ + Z55_ + lm[28] + Z54_ + Z61_ + lm[1];
   Z18_=invd49*invd48*invd47*invd46*invd45*invd44*invd43*invd42*invd41*
   invd40*Z18_;
   Z22_=Z23_ + Z13_ + Z21_ - lm[28] + Z20_ + Z60_ - lm[1];
   Z22_=invd89*invd88*invd87*invd86*invd85*invd84*invd83*invd82*invd81*
   invd80*Z22_;
   Z15_=Z48_ + Z15_ + Z46_;
   Z15_=lm[23]*Z15_;
   Z12_=Z57_ + Z12_ + Z55_ + lm[28] + lm[1] + Z54_;
   Z24_=Z12_ + Z56_;
   Z24_=invd39*invd38*invd37*invd36*invd35*invd34*invd33*invd32*invd31*
   invd30*Z24_;
   Z13_=Z23_ + Z13_ + Z21_ - lm[28] - lm[1] + Z20_;
   Z19_=Z13_ + Z19_;
   Z19_=invd79*invd78*invd77*invd76*invd75*invd74*invd73*invd72*invd71*
   invd70*Z19_;
   Z12_=Z12_ + Z65_;
   Z12_=invd29*invd28*invd27*invd26*invd25*invd24*invd23*invd22*invd21*
   invd20*Z12_;
   Z13_=Z13_ + Z53_;
   Z13_=invd69*invd68*invd67*invd66*invd65*invd64*invd63*invd62*invd61*
   invd60*Z13_;
   Z8_=2*Z8_ + Z11_ + 2*Z9_ + 2*Z15_ + 2*Z14_ + Z13_ + Z12_ + Z19_ + 
   Z24_ + Z22_ + Z18_ + Z16_ + Z17_;
   Z8_=Z6_*Z8_;
   Z9_= - invd19*invd18*invd17*invd16*invd15*invd14*invd13*invd12*
   invd11*invd10;
   Z11_=invd59*invd58*invd57*invd56*invd55*invd54*invd53*invd52*invd51*
   invd50;
   Z12_=invd9*invd8*invd7*invd6*invd5*invd4*invd3*invd2*invd1*invd0;


	*out = cpowq(2.q*M_PIq*I,2)/(1)*(Z8_ + Z9_ + Z10_ + Z11_ + Z12_ + Z28_ + Z29_ + Z30_ + Z31_ + 
      Z32_ + Z33_);
}
void evaluate_LTD_0_f128(__complex128 lm[], __complex128 params[], int conf, __complex128* out) {
   switch(conf) {
		case 0: evaluate_LTD_0_0_f128(lm, params, out); return;
		default: *out = 0.q;
    }
}
