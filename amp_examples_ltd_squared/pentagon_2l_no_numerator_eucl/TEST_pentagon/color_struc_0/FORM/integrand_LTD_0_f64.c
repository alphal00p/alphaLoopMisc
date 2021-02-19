#include <tgmath.h>
# include <quadmath.h>
# include <signal.h>
# include "numerator.h"

static inline void evaluate_LTD_0_0(double complex lm[], double complex params[], double complex* out) {
	double complex Z10_,Z11_,Z12_,Z13_,Z14_,Z4_,Z5_,Z6_,Z7_,Z8_,Z9_;
	double complex E0 = sqrt(lm[38]+0.0);
	double complex E1 = sqrt(lm[38]+2*lm[34]+lm[2]+0.0);
	double complex E2 = sqrt(lm[38]+2*lm[34]-2*lm[29]+lm[27]-2*lm[23]+lm[2]+0.0);
	double complex E3 = sqrt(lm[38]-2*lm[36]+2*   lm[18]+lm[14]-2*lm[12]+lm[7]+0.0);
	double complex E4 = sqrt(lm[38]-2*lm[36]+lm[7]+0.0);
	double complex E5 = sqrt(lm[47]+0.0);
	double complex E6 = sqrt(lm[47]+2*lm[40]+lm[38]+0.0);
	double complex invd0 = 1./(E0+E6-E5);
	double complex invd1 = 1./(-E0+E6-E5);
	double complex invd2 = 1./(lm[0]+E1+E6-E5);
	double complex invd3 = 1./(lm[0]-E1+E6-E5);
	double complex invd4 = 1./(-lm[21]+lm[0]+E2+E6-E5);
	double complex invd5 = 1./(-lm[21]+   lm[0]-E2+E6-E5);
	double complex invd6 = 1./(lm[8]-lm[5]+E3+E6-E5);
	double complex invd7 = 1./(lm[8]-lm[5]-E3+E6-E5);
	double complex invd8 = 1./(-lm[5]+E4+E6-E5);
	double complex invd9 = 1./(-lm[5]-E4+   E6-E5);
	double complex invd10 = 1./(2*E5);
	double complex invd11 = 1./(2*E6);
	double complex invd12 = 1./(-2*E0);
	double complex invd13 = 1./(lm[0]+E1-E0);
	double complex invd14 = 1./(lm[0]-E1-E0);
	double complex invd15 = 1./(-lm[21]+lm[0]+E2-E0);
	double complex invd16 = 1./(-lm[21]+lm[0]-   E2-E0);
	double complex invd17 = 1./(lm[8]-lm[5]+E3-E0);
	double complex invd18 = 1./(lm[8]-lm[5]-E3-E0);
	double complex invd19 = 1./(-lm[5]+E4-E0);
	double complex invd20 = 1./(-lm[5]-E4-E0);
	double complex invd21 = 1./(E0+E6+E5);
	double complex invd22 = 1./(E0+   E6-E5);
	double complex invd23 = 1./(2*E6);
	double complex invd24 = 1./(-lm[0]-E1+E0);
	double complex invd25 = 1./(-lm[0]-E1-E0);
	double complex invd26 = 1./(-2*E1);
	double complex invd27 = 1./(-lm[21]+E2-E1);
	double complex invd28 = 1./(-lm[21]-E2-E1);
	double complex invd29 = 1./(lm[8]-   lm[5]-lm[0]+E3-E1);
	double complex invd30 = 1./(lm[8]-lm[5]-lm[0]-E3-E1);
	double complex invd31 = 1./(-lm[5]-lm[0]+E4-E1);
	double complex invd32 = 1./(-lm[5]-lm[0]-E4-E1);
	double complex invd33 = 1./(lm[0]+E1+   E6+E5);
	double complex invd34 = 1./(lm[0]+E1+E6-E5);
	double complex invd35 = 1./(2*E6);
	double complex invd36 = 1./(lm[21]-lm[0]-E2+E0);
	double complex invd37 = 1./(lm[21]-lm[0]-E2-E0);
	double complex invd38 = 1./(lm[21]-E2+E1);
	double complex invd39 = 1./(   lm[21]-E2-E1);
	double complex invd40 = 1./(-2*E2);
	double complex invd41 = 1./(lm[21]+lm[8]-lm[5]-lm[0]+E3-E2);
	double complex invd42 = 1./(lm[21]+lm[8]-lm[5]-lm[0]-E3-E2);
	double complex invd43 = 1./(lm[21]-   lm[5]-lm[0]+E4-E2);
	double complex invd44 = 1./(lm[21]-lm[5]-lm[0]-E4-E2);
	double complex invd45 = 1./(-lm[21]+lm[0]+E2+E6+E5);
	double complex invd46 = 1./(-lm[21]+lm[0]+E2+E6-   E5);
	double complex invd47 = 1./(2*E6);
	double complex invd48 = 1./(-lm[8]+lm[5]-E3+E0);
	double complex invd49 = 1./(-lm[8]+lm[5]-E3-E0);
	double complex invd50 = 1./(-lm[8]+lm[5]+lm[0]-E3+E1);
	double complex invd51 = 1./(-lm[8]+lm[5]+   lm[0]-E3-E1);
	double complex invd52 = 1./(-lm[21]-lm[8]+lm[5]+lm[0]-E3+E2);
	double complex invd53 = 1./(-lm[21]-lm[8]+lm[5]+lm[0]-E3-E2);
	double complex invd54 = 1./(-2*E3);
	double complex invd55 = 1./(-lm[8]+   E4-E3);
	double complex invd56 = 1./(-lm[8]-E4-E3);
	double complex invd57 = 1./(lm[8]-lm[5]+E3+E6+E5);
	double complex invd58 = 1./(lm[8]-lm[5]+E3+E6-E5);
	double complex invd59 = 1./(2*E6);
	double complex invd60 = 1./(lm[5]-E4+E0);
	double complex invd61 = 1./(   lm[5]-E4-E0);
	double complex invd62 = 1./(lm[5]+lm[0]-E4+E1);
	double complex invd63 = 1./(lm[5]+lm[0]-E4-E1);
	double complex invd64 = 1./(-lm[21]+lm[5]+lm[0]-E4+E2);
	double complex invd65 = 1./(-lm[21]+lm[5]+   lm[0]-E4-E2);
	double complex invd66 = 1./(lm[8]-E4+E3);
	double complex invd67 = 1./(lm[8]-E4-E3);
	double complex invd68 = 1./(-2*E4);
	double complex invd69 = 1./(-lm[5]+E4+E6+E5);
	double complex invd70 = 1./(-lm[5]+E4+E6-E5);
	double complex invd71 = 1./(2*E6);
	double complex invd72 = 1./(   2*E0);
	double complex invd73 = 1./(lm[0]+E1+E0);
	double complex invd74 = 1./(lm[0]-E1+E0);
	double complex invd75 = 1./(-lm[21]+lm[0]+E2+E0);
	double complex invd76 = 1./(-lm[21]+lm[0]-E2+E0);
	double complex invd77 = 1./(lm[8]-lm[5]+E3+   E0);
	double complex invd78 = 1./(lm[8]-lm[5]-E3+E0);
	double complex invd79 = 1./(-lm[5]+E4+E0);
	double complex invd80 = 1./(-lm[5]-E4+E0);
	double complex invd81 = 1./(2*E5);
	double complex invd82 = 1./(E0+E6+E5);
	double complex invd83 = 1./(E0-E6+E5);
	double complex invd84 = 1./(-lm[0]+   E1+E0);
	double complex invd85 = 1./(-lm[0]+E1-E0);
	double complex invd86 = 1./(2*E1);
	double complex invd87 = 1./(-lm[21]+E2+E1);
	double complex invd88 = 1./(-lm[21]-E2+E1);
	double complex invd89 = 1./(lm[8]-lm[5]-lm[0]+E3+E1);
	double complex invd90 = 1./(lm[8]-   lm[5]-lm[0]-E3+E1);
	double complex invd91 = 1./(-lm[5]-lm[0]+E4+E1);
	double complex invd92 = 1./(-lm[5]-lm[0]-E4+E1);
	double complex invd93 = 1./(2*E5);
	double complex invd94 = 1./(-lm[0]+E1+E6+E5);
	double complex invd95 = 1./(-lm[0]+   E1-E6+E5);
	double complex invd96 = 1./(lm[21]-lm[0]+E2+E0);
	double complex invd97 = 1./(lm[21]-lm[0]+E2-E0);
	double complex invd98 = 1./(lm[21]+E2+E1);
	double complex invd99 = 1./(lm[21]+E2-E1);
	double complex invd100 = 1./(2*E2);
	double complex invd101 = 1./(   lm[21]+lm[8]-lm[5]-lm[0]+E3+E2);
	double complex invd102 = 1./(lm[21]+lm[8]-lm[5]-lm[0]-E3+E2);
	double complex invd103 = 1./(lm[21]-lm[5]-lm[0]+E4+E2);
	double complex invd104 = 1./(   lm[21]-lm[5]-lm[0]-E4+E2);
	double complex invd105 = 1./(2*E5);
	double complex invd106 = 1./(lm[21]-lm[0]+E2+E6+E5);
	double complex invd107 = 1./(lm[21]-lm[0]+E2-E6+E5);
	double complex invd108 = 1./(-lm[8]+lm[5]+   E3+E0);
	double complex invd109 = 1./(-lm[8]+lm[5]+E3-E0);
	double complex invd110 = 1./(-lm[8]+lm[5]+lm[0]+E3+E1);
	double complex invd111 = 1./(-lm[8]+lm[5]+lm[0]+E3-E1);
	double complex invd112 = 1./(-lm[21]-lm[8]+   lm[5]+lm[0]+E3+E2);
	double complex invd113 = 1./(-lm[21]-lm[8]+lm[5]+lm[0]+E3-E2);
	double complex invd114 = 1./(2*E3);
	double complex invd115 = 1./(-lm[8]+E4+E3);
	double complex invd116 = 1./(-lm[8]-E4+E3);
	double complex invd117 = 1./(2*   E5);
	double complex invd118 = 1./(-lm[8]+lm[5]+E3+E6+E5);
	double complex invd119 = 1./(-lm[8]+lm[5]+E3-E6+E5);
	double complex invd120 = 1./(lm[5]+E4+E0);
	double complex invd121 = 1./(lm[5]+E4-E0);
	double complex invd122 = 1./(lm[5]+lm[0]+   E4+E1);
	double complex invd123 = 1./(lm[5]+lm[0]+E4-E1);
	double complex invd124 = 1./(-lm[21]+lm[5]+lm[0]+E4+E2);
	double complex invd125 = 1./(-lm[21]+lm[5]+lm[0]+E4-E2);
	double complex invd126 = 1./(lm[8]+E4+   E3);
	double complex invd127 = 1./(lm[8]+E4-E3);
	double complex invd128 = 1./(2*E4);
	double complex invd129 = 1./(2*E5);
	double complex invd130 = 1./(lm[5]+E4+E6+E5);
	double complex invd131 = 1./(lm[5]+E4-E6+E5);


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


	*out = pow(2.*pi*I,2)/(1)*(Z4_ + Z5_ + Z6_ + Z7_ + Z8_ + Z9_ + Z10_ + Z11_ + Z12_ + Z13_
       + Z14_);
}
void evaluate_LTD_0(double complex lm[], double complex params[], int conf, double complex* out) {
   switch(conf) {
		case 0: evaluate_LTD_0_0(lm, params, out); return;
		default: *out = 0.;
    }
}
