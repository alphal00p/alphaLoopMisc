#include <tgmath.h>
# include <quadmath.h>
# include <signal.h>
# include "numerator.h"


static inline void evaluate_LTD_0_0_f128(__complex128 lm[], __complex128 params[], __complex128* out) {
	__complex128 Z4_,Z5_,Z6_,Z7_,Z8_;
	__complex128 E0 = csqrtq(lm[34]+1.0q);
	__complex128 E1 = csqrtq(lm[34]+2*lm[30]+lm[2]+1.0q);
	__complex128 E2 = csqrtq(lm[34]+2*lm[30]-2*lm[18]+lm[14]-2*lm[10]+lm[2]+1.0q);
	__complex128 E3 = csqrtq(lm[34]+2*lm[30]-2*   lm[27]+lm[25]-2*lm[21]-2*lm[18]+2*lm[16]+lm[14]-2*lm[10]+lm[2]+1.0q);
	__complex128 E4 = csqrtq(lm[34]-2*lm[32]+lm[7]+1.0q);
	__complex128 invd0 = 1.q/(2*E0);
	__complex128 invd1 = 1.q/(lm[0]+E1+E0);
	__complex128 invd2 = 1.q/(lm[0]-E1+E0);
	__complex128 invd3 = 1.q/(-lm[8]+lm[0]+E2+E0);
	__complex128 invd4 = 1.q/(-lm[8]+lm[0]-E2+E0);
	__complex128 invd5 = 1.q/(-lm[19]-lm[8]+lm[0]+E3+   E0);
	__complex128 invd6 = 1.q/(-lm[19]-lm[8]+lm[0]-E3+E0);
	__complex128 invd7 = 1.q/(-lm[5]+E4+E0);
	__complex128 invd8 = 1.q/(-lm[5]-E4+E0);
	__complex128 invd9 = 1.q/(-lm[0]+E1+E0);
	__complex128 invd10 = 1.q/(-lm[0]+E1-E0);
	__complex128 invd11 = 1.q/(2*   E1);
	__complex128 invd12 = 1.q/(-lm[8]+E2+E1);
	__complex128 invd13 = 1.q/(-lm[8]-E2+E1);
	__complex128 invd14 = 1.q/(-lm[19]-lm[8]+E3+E1);
	__complex128 invd15 = 1.q/(-lm[19]-lm[8]-E3+E1);
	__complex128 invd16 = 1.q/(-lm[5]-lm[0]+E4+   E1);
	__complex128 invd17 = 1.q/(-lm[5]-lm[0]-E4+E1);
	__complex128 invd18 = 1.q/(lm[8]-lm[0]+E2+E0);
	__complex128 invd19 = 1.q/(lm[8]-lm[0]+E2-E0);
	__complex128 invd20 = 1.q/(lm[8]+E2+E1);
	__complex128 invd21 = 1.q/(lm[8]+E2-E1);
	__complex128 invd22 = 1.q/(2*   E2);
	__complex128 invd23 = 1.q/(-lm[19]+E3+E2);
	__complex128 invd24 = 1.q/(-lm[19]-E3+E2);
	__complex128 invd25 = 1.q/(lm[8]-lm[5]-lm[0]+E4+E2);
	__complex128 invd26 = 1.q/(lm[8]-lm[5]-lm[0]-E4+E2);
	__complex128 invd27 = 1.q/(lm[19]+   lm[8]-lm[0]+E3+E0);
	__complex128 invd28 = 1.q/(lm[19]+lm[8]-lm[0]+E3-E0);
	__complex128 invd29 = 1.q/(lm[19]+lm[8]+E3+E1);
	__complex128 invd30 = 1.q/(lm[19]+lm[8]+E3-E1);
	__complex128 invd31 = 1.q/(lm[19]+   E3+E2);
	__complex128 invd32 = 1.q/(lm[19]+E3-E2);
	__complex128 invd33 = 1.q/(2*E3);
	__complex128 invd34 = 1.q/(lm[19]+lm[8]-lm[5]-lm[0]+E4+E3);
	__complex128 invd35 = 1.q/(lm[19]+lm[8]-lm[5]-lm[0]-E4+E3);
	__complex128 invd36 = 1.q/(   lm[5]+E4+E0);
	__complex128 invd37 = 1.q/(lm[5]+E4-E0);
	__complex128 invd38 = 1.q/(lm[5]+lm[0]+E4+E1);
	__complex128 invd39 = 1.q/(lm[5]+lm[0]+E4-E1);
	__complex128 invd40 = 1.q/(-lm[8]+lm[5]+lm[0]+E4+E2);
	__complex128 invd41 = 1.q/(-   lm[8]+lm[5]+lm[0]+E4-E2);
	__complex128 invd42 = 1.q/(-lm[19]-lm[8]+lm[5]+lm[0]+E4+E3);
	__complex128 invd43 = 1.q/(-lm[19]-lm[8]+lm[5]+lm[0]+E4-E3);
	__complex128 invd44 = 1.q/(2*   E4);


   Z4_= - invd36*invd37*invd38*invd39*invd40*invd41*invd42*invd43*
   invd44;
   Z5_= - invd27*invd28*invd29*invd30*invd31*invd32*invd33*invd34*
   invd35;
   Z6_= - invd18*invd19*invd20*invd21*invd22*invd23*invd24*invd25*
   invd26;
   Z7_= - invd9*invd10*invd11*invd12*invd13*invd14*invd15*invd16*invd17
   ;
   Z8_= - invd0*invd1*invd2*invd3*invd4*invd5*invd6*invd7*invd8;


	*out = cpowq(2.q*M_PIq*I,1)/(1)*(Z4_ + Z5_ + Z6_ + Z7_ + Z8_);
}
void evaluate_LTD_0_f128(__complex128 lm[], __complex128 params[], int conf, __complex128* out) {
   switch(conf) {
		case 0: evaluate_LTD_0_0_f128(lm, params, out); return;
		default: *out = 0.q;
    }
}
