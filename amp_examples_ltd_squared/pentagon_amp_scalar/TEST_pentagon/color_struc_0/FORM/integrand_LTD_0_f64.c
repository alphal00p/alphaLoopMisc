#include <tgmath.h>
# include <quadmath.h>
# include <signal.h>
# include "numerator.h"

static inline void evaluate_LTD_0_0(double complex lm[], double complex params[], double complex* out) {
	double complex Z4_,Z5_,Z6_,Z7_,Z8_;
	double complex E0 = sqrt(lm[34]+1.0);
	double complex E1 = sqrt(lm[34]+2*lm[30]+lm[2]+1.0);
	double complex E2 = sqrt(lm[34]+2*lm[30]-2*lm[18]+lm[14]-2*lm[10]+lm[2]+1.0);
	double complex E3 = sqrt(lm[34]+2*lm[30]-2*   lm[27]+lm[25]-2*lm[21]-2*lm[18]+2*lm[16]+lm[14]-2*lm[10]+lm[2]+1.0);
	double complex E4 = sqrt(lm[34]-2*lm[32]+lm[7]+1.0);
	double complex invd0 = 1./(2*E0);
	double complex invd1 = 1./(lm[0]+E1+E0);
	double complex invd2 = 1./(lm[0]-E1+E0);
	double complex invd3 = 1./(-lm[8]+lm[0]+E2+E0);
	double complex invd4 = 1./(-lm[8]+lm[0]-E2+E0);
	double complex invd5 = 1./(-lm[19]-lm[8]+lm[0]+E3+   E0);
	double complex invd6 = 1./(-lm[19]-lm[8]+lm[0]-E3+E0);
	double complex invd7 = 1./(-lm[5]+E4+E0);
	double complex invd8 = 1./(-lm[5]-E4+E0);
	double complex invd9 = 1./(-lm[0]+E1+E0);
	double complex invd10 = 1./(-lm[0]+E1-E0);
	double complex invd11 = 1./(2*   E1);
	double complex invd12 = 1./(-lm[8]+E2+E1);
	double complex invd13 = 1./(-lm[8]-E2+E1);
	double complex invd14 = 1./(-lm[19]-lm[8]+E3+E1);
	double complex invd15 = 1./(-lm[19]-lm[8]-E3+E1);
	double complex invd16 = 1./(-lm[5]-lm[0]+E4+   E1);
	double complex invd17 = 1./(-lm[5]-lm[0]-E4+E1);
	double complex invd18 = 1./(lm[8]-lm[0]+E2+E0);
	double complex invd19 = 1./(lm[8]-lm[0]+E2-E0);
	double complex invd20 = 1./(lm[8]+E2+E1);
	double complex invd21 = 1./(lm[8]+E2-E1);
	double complex invd22 = 1./(2*   E2);
	double complex invd23 = 1./(-lm[19]+E3+E2);
	double complex invd24 = 1./(-lm[19]-E3+E2);
	double complex invd25 = 1./(lm[8]-lm[5]-lm[0]+E4+E2);
	double complex invd26 = 1./(lm[8]-lm[5]-lm[0]-E4+E2);
	double complex invd27 = 1./(lm[19]+   lm[8]-lm[0]+E3+E0);
	double complex invd28 = 1./(lm[19]+lm[8]-lm[0]+E3-E0);
	double complex invd29 = 1./(lm[19]+lm[8]+E3+E1);
	double complex invd30 = 1./(lm[19]+lm[8]+E3-E1);
	double complex invd31 = 1./(lm[19]+   E3+E2);
	double complex invd32 = 1./(lm[19]+E3-E2);
	double complex invd33 = 1./(2*E3);
	double complex invd34 = 1./(lm[19]+lm[8]-lm[5]-lm[0]+E4+E3);
	double complex invd35 = 1./(lm[19]+lm[8]-lm[5]-lm[0]-E4+E3);
	double complex invd36 = 1./(   lm[5]+E4+E0);
	double complex invd37 = 1./(lm[5]+E4-E0);
	double complex invd38 = 1./(lm[5]+lm[0]+E4+E1);
	double complex invd39 = 1./(lm[5]+lm[0]+E4-E1);
	double complex invd40 = 1./(-lm[8]+lm[5]+lm[0]+E4+E2);
	double complex invd41 = 1./(-   lm[8]+lm[5]+lm[0]+E4-E2);
	double complex invd42 = 1./(-lm[19]-lm[8]+lm[5]+lm[0]+E4+E3);
	double complex invd43 = 1./(-lm[19]-lm[8]+lm[5]+lm[0]+E4-E3);
	double complex invd44 = 1./(2*   E4);


   Z4_= - invd36*invd37*invd38*invd39*invd40*invd41*invd42*invd43*
   invd44;
   Z5_= - invd27*invd28*invd29*invd30*invd31*invd32*invd33*invd34*
   invd35;
   Z6_= - invd18*invd19*invd20*invd21*invd22*invd23*invd24*invd25*
   invd26;
   Z7_= - invd9*invd10*invd11*invd12*invd13*invd14*invd15*invd16*invd17
   ;
   Z8_= - invd0*invd1*invd2*invd3*invd4*invd5*invd6*invd7*invd8;


	*out = pow(2.*pi*I,1)/(1)*(Z4_ + Z5_ + Z6_ + Z7_ + Z8_);
}
void evaluate_LTD_0(double complex lm[], double complex params[], int conf, double complex* out) {
   switch(conf) {
		case 0: evaluate_LTD_0_0(lm, params, out); return;
		default: *out = 0.;
    }
}
