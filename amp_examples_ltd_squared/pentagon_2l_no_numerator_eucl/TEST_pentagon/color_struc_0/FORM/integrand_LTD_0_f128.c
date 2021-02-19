#include <tgmath.h>
# include <quadmath.h>
# include <signal.h>
# include "numerator.h"


static inline void evaluate_LTD_0_0_f128(__complex128 lm[], __complex128 params[], __complex128* out) {
	__complex128 Z10_,Z11_,Z12_,Z13_,Z14_,Z4_,Z5_,Z6_,Z7_,Z8_,Z9_;
	__complex128 E0 = csqrtq(lm[38]+0.0q);
	__complex128 E1 = csqrtq(lm[38]+2*lm[34]+lm[2]+0.0q);
	__complex128 E2 = csqrtq(lm[38]+2*lm[34]-2*lm[29]+lm[27]-2*lm[23]+lm[2]+0.0q);
	__complex128 E3 = csqrtq(lm[38]-2*lm[36]+2*   lm[18]+lm[14]-2*lm[12]+lm[7]+0.0q);
	__complex128 E4 = csqrtq(lm[38]-2*lm[36]+lm[7]+0.0q);
	__complex128 E5 = csqrtq(lm[47]+0.0q);
	__complex128 E6 = csqrtq(lm[47]+2*lm[40]+lm[38]+0.0q);
	__complex128 invd0 = 1.q/(E0+E6-E5);
	__complex128 invd1 = 1.q/(-E0+E6-E5);
	__complex128 invd2 = 1.q/(lm[0]+E1+E6-E5);
	__complex128 invd3 = 1.q/(lm[0]-E1+E6-E5);
	__complex128 invd4 = 1.q/(-lm[21]+lm[0]+E2+E6-E5);
	__complex128 invd5 = 1.q/(-lm[21]+   lm[0]-E2+E6-E5);
	__complex128 invd6 = 1.q/(lm[8]-lm[5]+E3+E6-E5);
	__complex128 invd7 = 1.q/(lm[8]-lm[5]-E3+E6-E5);
	__complex128 invd8 = 1.q/(-lm[5]+E4+E6-E5);
	__complex128 invd9 = 1.q/(-lm[5]-E4+   E6-E5);
	__complex128 invd10 = 1.q/(2*E5);
	__complex128 invd11 = 1.q/(2*E6);
	__complex128 invd12 = 1.q/(-2*E0);
	__complex128 invd13 = 1.q/(lm[0]+E1-E0);
	__complex128 invd14 = 1.q/(lm[0]-E1-E0);
	__complex128 invd15 = 1.q/(-lm[21]+lm[0]+E2-E0);
	__complex128 invd16 = 1.q/(-lm[21]+lm[0]-   E2-E0);
	__complex128 invd17 = 1.q/(lm[8]-lm[5]+E3-E0);
	__complex128 invd18 = 1.q/(lm[8]-lm[5]-E3-E0);
	__complex128 invd19 = 1.q/(-lm[5]+E4-E0);
	__complex128 invd20 = 1.q/(-lm[5]-E4-E0);
	__complex128 invd21 = 1.q/(E0+E6+E5);
	__complex128 invd22 = 1.q/(E0+   E6-E5);
	__complex128 invd23 = 1.q/(2*E6);
	__complex128 invd24 = 1.q/(-lm[0]-E1+E0);
	__complex128 invd25 = 1.q/(-lm[0]-E1-E0);
	__complex128 invd26 = 1.q/(-2*E1);
	__complex128 invd27 = 1.q/(-lm[21]+E2-E1);
	__complex128 invd28 = 1.q/(-lm[21]-E2-E1);
	__complex128 invd29 = 1.q/(lm[8]-   lm[5]-lm[0]+E3-E1);
	__complex128 invd30 = 1.q/(lm[8]-lm[5]-lm[0]-E3-E1);
	__complex128 invd31 = 1.q/(-lm[5]-lm[0]+E4-E1);
	__complex128 invd32 = 1.q/(-lm[5]-lm[0]-E4-E1);
	__complex128 invd33 = 1.q/(lm[0]+E1+   E6+E5);
	__complex128 invd34 = 1.q/(lm[0]+E1+E6-E5);
	__complex128 invd35 = 1.q/(2*E6);
	__complex128 invd36 = 1.q/(lm[21]-lm[0]-E2+E0);
	__complex128 invd37 = 1.q/(lm[21]-lm[0]-E2-E0);
	__complex128 invd38 = 1.q/(lm[21]-E2+E1);
	__complex128 invd39 = 1.q/(   lm[21]-E2-E1);
	__complex128 invd40 = 1.q/(-2*E2);
	__complex128 invd41 = 1.q/(lm[21]+lm[8]-lm[5]-lm[0]+E3-E2);
	__complex128 invd42 = 1.q/(lm[21]+lm[8]-lm[5]-lm[0]-E3-E2);
	__complex128 invd43 = 1.q/(lm[21]-   lm[5]-lm[0]+E4-E2);
	__complex128 invd44 = 1.q/(lm[21]-lm[5]-lm[0]-E4-E2);
	__complex128 invd45 = 1.q/(-lm[21]+lm[0]+E2+E6+E5);
	__complex128 invd46 = 1.q/(-lm[21]+lm[0]+E2+E6-   E5);
	__complex128 invd47 = 1.q/(2*E6);
	__complex128 invd48 = 1.q/(-lm[8]+lm[5]-E3+E0);
	__complex128 invd49 = 1.q/(-lm[8]+lm[5]-E3-E0);
	__complex128 invd50 = 1.q/(-lm[8]+lm[5]+lm[0]-E3+E1);
	__complex128 invd51 = 1.q/(-lm[8]+lm[5]+   lm[0]-E3-E1);
	__complex128 invd52 = 1.q/(-lm[21]-lm[8]+lm[5]+lm[0]-E3+E2);
	__complex128 invd53 = 1.q/(-lm[21]-lm[8]+lm[5]+lm[0]-E3-E2);
	__complex128 invd54 = 1.q/(-2*E3);
	__complex128 invd55 = 1.q/(-lm[8]+   E4-E3);
	__complex128 invd56 = 1.q/(-lm[8]-E4-E3);
	__complex128 invd57 = 1.q/(lm[8]-lm[5]+E3+E6+E5);
	__complex128 invd58 = 1.q/(lm[8]-lm[5]+E3+E6-E5);
	__complex128 invd59 = 1.q/(2*E6);
	__complex128 invd60 = 1.q/(lm[5]-E4+E0);
	__complex128 invd61 = 1.q/(   lm[5]-E4-E0);
	__complex128 invd62 = 1.q/(lm[5]+lm[0]-E4+E1);
	__complex128 invd63 = 1.q/(lm[5]+lm[0]-E4-E1);
	__complex128 invd64 = 1.q/(-lm[21]+lm[5]+lm[0]-E4+E2);
	__complex128 invd65 = 1.q/(-lm[21]+lm[5]+   lm[0]-E4-E2);
	__complex128 invd66 = 1.q/(lm[8]-E4+E3);
	__complex128 invd67 = 1.q/(lm[8]-E4-E3);
	__complex128 invd68 = 1.q/(-2*E4);
	__complex128 invd69 = 1.q/(-lm[5]+E4+E6+E5);
	__complex128 invd70 = 1.q/(-lm[5]+E4+E6-E5);
	__complex128 invd71 = 1.q/(2*E6);
	__complex128 invd72 = 1.q/(   2*E0);
	__complex128 invd73 = 1.q/(lm[0]+E1+E0);
	__complex128 invd74 = 1.q/(lm[0]-E1+E0);
	__complex128 invd75 = 1.q/(-lm[21]+lm[0]+E2+E0);
	__complex128 invd76 = 1.q/(-lm[21]+lm[0]-E2+E0);
	__complex128 invd77 = 1.q/(lm[8]-lm[5]+E3+   E0);
	__complex128 invd78 = 1.q/(lm[8]-lm[5]-E3+E0);
	__complex128 invd79 = 1.q/(-lm[5]+E4+E0);
	__complex128 invd80 = 1.q/(-lm[5]-E4+E0);
	__complex128 invd81 = 1.q/(2*E5);
	__complex128 invd82 = 1.q/(E0+E6+E5);
	__complex128 invd83 = 1.q/(E0-E6+E5);
	__complex128 invd84 = 1.q/(-lm[0]+   E1+E0);
	__complex128 invd85 = 1.q/(-lm[0]+E1-E0);
	__complex128 invd86 = 1.q/(2*E1);
	__complex128 invd87 = 1.q/(-lm[21]+E2+E1);
	__complex128 invd88 = 1.q/(-lm[21]-E2+E1);
	__complex128 invd89 = 1.q/(lm[8]-lm[5]-lm[0]+E3+E1);
	__complex128 invd90 = 1.q/(lm[8]-   lm[5]-lm[0]-E3+E1);
	__complex128 invd91 = 1.q/(-lm[5]-lm[0]+E4+E1);
	__complex128 invd92 = 1.q/(-lm[5]-lm[0]-E4+E1);
	__complex128 invd93 = 1.q/(2*E5);
	__complex128 invd94 = 1.q/(-lm[0]+E1+E6+E5);
	__complex128 invd95 = 1.q/(-lm[0]+   E1-E6+E5);
	__complex128 invd96 = 1.q/(lm[21]-lm[0]+E2+E0);
	__complex128 invd97 = 1.q/(lm[21]-lm[0]+E2-E0);
	__complex128 invd98 = 1.q/(lm[21]+E2+E1);
	__complex128 invd99 = 1.q/(lm[21]+E2-E1);
	__complex128 invd100 = 1.q/(2*E2);
	__complex128 invd101 = 1.q/(   lm[21]+lm[8]-lm[5]-lm[0]+E3+E2);
	__complex128 invd102 = 1.q/(lm[21]+lm[8]-lm[5]-lm[0]-E3+E2);
	__complex128 invd103 = 1.q/(lm[21]-lm[5]-lm[0]+E4+E2);
	__complex128 invd104 = 1.q/(   lm[21]-lm[5]-lm[0]-E4+E2);
	__complex128 invd105 = 1.q/(2*E5);
	__complex128 invd106 = 1.q/(lm[21]-lm[0]+E2+E6+E5);
	__complex128 invd107 = 1.q/(lm[21]-lm[0]+E2-E6+E5);
	__complex128 invd108 = 1.q/(-lm[8]+lm[5]+   E3+E0);
	__complex128 invd109 = 1.q/(-lm[8]+lm[5]+E3-E0);
	__complex128 invd110 = 1.q/(-lm[8]+lm[5]+lm[0]+E3+E1);
	__complex128 invd111 = 1.q/(-lm[8]+lm[5]+lm[0]+E3-E1);
	__complex128 invd112 = 1.q/(-lm[21]-lm[8]+   lm[5]+lm[0]+E3+E2);
	__complex128 invd113 = 1.q/(-lm[21]-lm[8]+lm[5]+lm[0]+E3-E2);
	__complex128 invd114 = 1.q/(2*E3);
	__complex128 invd115 = 1.q/(-lm[8]+E4+E3);
	__complex128 invd116 = 1.q/(-lm[8]-E4+E3);
	__complex128 invd117 = 1.q/(2*   E5);
	__complex128 invd118 = 1.q/(-lm[8]+lm[5]+E3+E6+E5);
	__complex128 invd119 = 1.q/(-lm[8]+lm[5]+E3-E6+E5);
	__complex128 invd120 = 1.q/(lm[5]+E4+E0);
	__complex128 invd121 = 1.q/(lm[5]+E4-E0);
	__complex128 invd122 = 1.q/(lm[5]+lm[0]+   E4+E1);
	__complex128 invd123 = 1.q/(lm[5]+lm[0]+E4-E1);
	__complex128 invd124 = 1.q/(-lm[21]+lm[5]+lm[0]+E4+E2);
	__complex128 invd125 = 1.q/(-lm[21]+lm[5]+lm[0]+E4-E2);
	__complex128 invd126 = 1.q/(lm[8]+E4+   E3);
	__complex128 invd127 = 1.q/(lm[8]+E4-E3);
	__complex128 invd128 = 1.q/(2*E4);
	__complex128 invd129 = 1.q/(2*E5);
	__complex128 invd130 = 1.q/(lm[5]+E4+E6+E5);
	__complex128 invd131 = 1.q/(lm[5]+E4-E6+E5);


   Z4_=invd120*invd121*invd122*invd123*invd124*invd125*invd126*invd127*
   invd128*invd129*invd130*invd131;
   Z5_=invd108*invd109*invd110*invd111*invd112*invd113*invd114*invd115*
   invd116*invd117*invd118*invd119;
   Z6_=invd96*invd97*invd98*invd99*invd100*invd101*invd102*invd103*
   invd104*invd105*invd106*invd107;
   Z7_=invd84*invd85*invd86*invd87*invd88*invd89*invd90*invd91*invd92*
   invd93*invd94*invd95;
   Z8_=invd72*invd73*invd74*invd75*invd76*invd77*invd78*invd79*invd80*
   invd81*invd82*invd83;
   Z9_= - invd60*invd61*invd62*invd63*invd64*invd65*invd66*invd67*
   invd68*invd69*invd70*invd71;
   Z10_= - invd48*invd49*invd50*invd51*invd52*invd53*invd54*invd55*
   invd56*invd57*invd58*invd59;
   Z11_= - invd36*invd37*invd38*invd39*invd40*invd41*invd42*invd43*
   invd44*invd45*invd46*invd47;
   Z12_= - invd24*invd25*invd26*invd27*invd28*invd29*invd30*invd31*
   invd32*invd33*invd34*invd35;
   Z13_= - invd12*invd13*invd14*invd15*invd16*invd17*invd18*invd19*
   invd20*invd21*invd22*invd23;
   Z14_=invd0*invd1*invd2*invd3*invd4*invd5*invd6*invd7*invd8*invd9*
   invd10*invd11;


	*out = cpowq(2.q*M_PIq*I,2)/(1)*(Z4_ + Z5_ + Z6_ + Z7_ + Z8_ + Z9_ + Z10_ + Z11_ + Z12_ + Z13_
       + Z14_);
}
void evaluate_LTD_0_f128(__complex128 lm[], __complex128 params[], int conf, __complex128* out) {
   switch(conf) {
		case 0: evaluate_LTD_0_0_f128(lm, params, out); return;
		default: *out = 0.q;
    }
}
