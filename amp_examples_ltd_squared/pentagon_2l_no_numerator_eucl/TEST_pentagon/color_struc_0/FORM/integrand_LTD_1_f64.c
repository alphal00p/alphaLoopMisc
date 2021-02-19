#include <tgmath.h>
# include <quadmath.h>
# include <signal.h>
# include "numerator.h"

static inline void evaluate_LTD_1_0(double complex lm[], double complex params[], double complex* out) {
	double complex Z4_,Z5_,Z6_,Z7_,Z8_;
	double complex E0 = sqrt(lm[38]+0.0);
	double complex E1 = sqrt(lm[38]+2*lm[34]+lm[2]+0.0);
	double complex E2 = sqrt(lm[38]+2*lm[34]-2*lm[18]+lm[14]-2*lm[10]+lm[2]+0.0);
	double complex E3 = sqrt(lm[38]+2*lm[34]-2*   lm[29]+lm[27]-2*lm[23]-2*lm[18]+2*lm[16]+lm[14]-2*lm[10]+lm[2]+0.0);
	double complex E4 = sqrt(lm[38]-2*lm[36]+lm[7]+0.0);
	double complex E5 = sqrt(lm[47]+10000.0);
	double complex invd0 = 1./(2*E0);
	double complex invd1 = 1./(lm[0]+E1+E0);
	double complex invd2 = 1./(lm[0]-E1+E0);
	double complex invd3 = 1./(-lm[8]+lm[0]+E2+E0);
	double complex invd4 = 1./(-lm[8]+lm[0]-E2+E0);
	double complex invd5 = 1./(-lm[21]-lm[8]+lm[0]+E3+   E0);
	double complex invd6 = 1./(-lm[21]-lm[8]+lm[0]-E3+E0);
	double complex invd7 = 1./(-lm[5]+E4+E0);
	double complex invd8 = 1./(-lm[5]-E4+E0);
	double complex invd9 = 1./(2*E5);
	double complex invd10 = 1./(-lm[0]+E1+E0);
	double complex invd11 = 1./(-lm[0]+E1-   E0);
	double complex invd12 = 1./(2*E1);
	double complex invd13 = 1./(-lm[8]+E2+E1);
	double complex invd14 = 1./(-lm[8]-E2+E1);
	double complex invd15 = 1./(-lm[21]-lm[8]+E3+E1);
	double complex invd16 = 1./(-lm[21]-lm[8]-E3+E1);
	double complex invd17 = 1./(-lm[5]-   lm[0]+E4+E1);
	double complex invd18 = 1./(-lm[5]-lm[0]-E4+E1);
	double complex invd19 = 1./(2*E5);
	double complex invd20 = 1./(lm[8]-lm[0]+E2+E0);
	double complex invd21 = 1./(lm[8]-lm[0]+E2-E0);
	double complex invd22 = 1./(lm[8]+E2+E1);
	double complex invd23 = 1./(   lm[8]+E2-E1);
	double complex invd24 = 1./(2*E2);
	double complex invd25 = 1./(-lm[21]+E3+E2);
	double complex invd26 = 1./(-lm[21]-E3+E2);
	double complex invd27 = 1./(lm[8]-lm[5]-lm[0]+E4+E2);
	double complex invd28 = 1./(lm[8]-lm[5]-lm[0]-   E4+E2);
	double complex invd29 = 1./(2*E5);
	double complex invd30 = 1./(lm[21]+lm[8]-lm[0]+E3+E0);
	double complex invd31 = 1./(lm[21]+lm[8]-lm[0]+E3-E0);
	double complex invd32 = 1./(lm[21]+lm[8]+E3+E1);
	double complex invd33 = 1./(lm[21]+   lm[8]+E3-E1);
	double complex invd34 = 1./(lm[21]+E3+E2);
	double complex invd35 = 1./(lm[21]+E3-E2);
	double complex invd36 = 1./(2*E3);
	double complex invd37 = 1./(lm[21]+lm[8]-lm[5]-lm[0]+E4+E3);
	double complex invd38 = 1./(lm[21]+lm[8]-   lm[5]-lm[0]-E4+E3);
	double complex invd39 = 1./(2*E5);
	double complex invd40 = 1./(lm[5]+E4+E0);
	double complex invd41 = 1./(lm[5]+E4-E0);
	double complex invd42 = 1./(lm[5]+lm[0]+E4+E1);
	double complex invd43 = 1./(lm[5]+lm[0]+E4-E1);
	double complex invd44 = 1./(-   lm[8]+lm[5]+lm[0]+E4+E2);
	double complex invd45 = 1./(-lm[8]+lm[5]+lm[0]+E4-E2);
	double complex invd46 = 1./(-lm[21]-lm[8]+lm[5]+lm[0]+E4+E3);
	double complex invd47 = 1./(-lm[21]-   lm[8]+lm[5]+lm[0]+E4-E3);
	double complex invd48 = 1./(2*E4);
	double complex invd49 = 1./(2*E5);


   Z4_=invd40*invd41*invd42*invd43*invd44*invd45*invd46*invd47*invd48*
   pow(invd49,3);
   Z5_=invd30*invd31*invd32*invd33*invd34*invd35*invd36*invd37*invd38*
   pow(invd39,3);
   Z6_=invd20*invd21*invd22*invd23*invd24*invd25*invd26*invd27*invd28*
   pow(invd29,3);
   Z7_=invd10*invd11*invd12*invd13*invd14*invd15*invd16*invd17*invd18*
   pow(invd19,3);
   Z8_=invd0*invd1*invd2*invd3*invd4*invd5*invd6*invd7*invd8*pow(
   invd9,3);
   Z4_=Z8_ + Z7_ + Z6_ + Z4_ + Z5_;


	*out = pow(2.*pi*I,2)/(1)*(2*Z4_);
}
void evaluate_LTD_1(double complex lm[], double complex params[], int conf, double complex* out) {
   switch(conf) {
		case 0: evaluate_LTD_1_0(lm, params, out); return;
		default: *out = 0.;
    }
}
