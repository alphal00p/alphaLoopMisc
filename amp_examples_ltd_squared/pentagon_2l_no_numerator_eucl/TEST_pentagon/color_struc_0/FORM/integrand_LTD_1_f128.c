#include <tgmath.h>
# include <quadmath.h>
# include <signal.h>
# include "numerator.h"


static inline void evaluate_LTD_1_0_f128(__complex128 lm[], __complex128 params[], __complex128* out) {
	__complex128 Z4_,Z5_,Z6_,Z7_,Z8_;
	__complex128 E0 = csqrtq(lm[38]+0.0q);
	__complex128 E1 = csqrtq(lm[38]+2*lm[34]+lm[2]+0.0q);
	__complex128 E2 = csqrtq(lm[38]+2*lm[34]-2*lm[18]+lm[14]-2*lm[10]+lm[2]+0.0q);
	__complex128 E3 = csqrtq(lm[38]+2*lm[34]-2*   lm[29]+lm[27]-2*lm[23]-2*lm[18]+2*lm[16]+lm[14]-2*lm[10]+lm[2]+0.0q);
	__complex128 E4 = csqrtq(lm[38]-2*lm[36]+lm[7]+0.0q);
	__complex128 E5 = csqrtq(lm[47]+10000.0q);
	__complex128 invd0 = 1.q/(2*E0);
	__complex128 invd1 = 1.q/(lm[0]+E1+E0);
	__complex128 invd2 = 1.q/(lm[0]-E1+E0);
	__complex128 invd3 = 1.q/(-lm[8]+lm[0]+E2+E0);
	__complex128 invd4 = 1.q/(-lm[8]+lm[0]-E2+E0);
	__complex128 invd5 = 1.q/(-lm[21]-lm[8]+lm[0]+E3+   E0);
	__complex128 invd6 = 1.q/(-lm[21]-lm[8]+lm[0]-E3+E0);
	__complex128 invd7 = 1.q/(-lm[5]+E4+E0);
	__complex128 invd8 = 1.q/(-lm[5]-E4+E0);
	__complex128 invd9 = 1.q/(2*E5);
	__complex128 invd10 = 1.q/(-lm[0]+E1+E0);
	__complex128 invd11 = 1.q/(-lm[0]+E1-   E0);
	__complex128 invd12 = 1.q/(2*E1);
	__complex128 invd13 = 1.q/(-lm[8]+E2+E1);
	__complex128 invd14 = 1.q/(-lm[8]-E2+E1);
	__complex128 invd15 = 1.q/(-lm[21]-lm[8]+E3+E1);
	__complex128 invd16 = 1.q/(-lm[21]-lm[8]-E3+E1);
	__complex128 invd17 = 1.q/(-lm[5]-   lm[0]+E4+E1);
	__complex128 invd18 = 1.q/(-lm[5]-lm[0]-E4+E1);
	__complex128 invd19 = 1.q/(2*E5);
	__complex128 invd20 = 1.q/(lm[8]-lm[0]+E2+E0);
	__complex128 invd21 = 1.q/(lm[8]-lm[0]+E2-E0);
	__complex128 invd22 = 1.q/(lm[8]+E2+E1);
	__complex128 invd23 = 1.q/(   lm[8]+E2-E1);
	__complex128 invd24 = 1.q/(2*E2);
	__complex128 invd25 = 1.q/(-lm[21]+E3+E2);
	__complex128 invd26 = 1.q/(-lm[21]-E3+E2);
	__complex128 invd27 = 1.q/(lm[8]-lm[5]-lm[0]+E4+E2);
	__complex128 invd28 = 1.q/(lm[8]-lm[5]-lm[0]-   E4+E2);
	__complex128 invd29 = 1.q/(2*E5);
	__complex128 invd30 = 1.q/(lm[21]+lm[8]-lm[0]+E3+E0);
	__complex128 invd31 = 1.q/(lm[21]+lm[8]-lm[0]+E3-E0);
	__complex128 invd32 = 1.q/(lm[21]+lm[8]+E3+E1);
	__complex128 invd33 = 1.q/(lm[21]+   lm[8]+E3-E1);
	__complex128 invd34 = 1.q/(lm[21]+E3+E2);
	__complex128 invd35 = 1.q/(lm[21]+E3-E2);
	__complex128 invd36 = 1.q/(2*E3);
	__complex128 invd37 = 1.q/(lm[21]+lm[8]-lm[5]-lm[0]+E4+E3);
	__complex128 invd38 = 1.q/(lm[21]+lm[8]-   lm[5]-lm[0]-E4+E3);
	__complex128 invd39 = 1.q/(2*E5);
	__complex128 invd40 = 1.q/(lm[5]+E4+E0);
	__complex128 invd41 = 1.q/(lm[5]+E4-E0);
	__complex128 invd42 = 1.q/(lm[5]+lm[0]+E4+E1);
	__complex128 invd43 = 1.q/(lm[5]+lm[0]+E4-E1);
	__complex128 invd44 = 1.q/(-   lm[8]+lm[5]+lm[0]+E4+E2);
	__complex128 invd45 = 1.q/(-lm[8]+lm[5]+lm[0]+E4-E2);
	__complex128 invd46 = 1.q/(-lm[21]-lm[8]+lm[5]+lm[0]+E4+E3);
	__complex128 invd47 = 1.q/(-lm[21]-   lm[8]+lm[5]+lm[0]+E4-E3);
	__complex128 invd48 = 1.q/(2*E4);
	__complex128 invd49 = 1.q/(2*E5);


   Z4_=invd40*invd41*invd42*invd43*invd44*invd45*invd46*invd47*invd48*
   cpowq(invd49,3);
   Z5_=invd30*invd31*invd32*invd33*invd34*invd35*invd36*invd37*invd38*
   cpowq(invd39,3);
   Z6_=invd20*invd21*invd22*invd23*invd24*invd25*invd26*invd27*invd28*
   cpowq(invd29,3);
   Z7_=invd10*invd11*invd12*invd13*invd14*invd15*invd16*invd17*invd18*
   cpowq(invd19,3);
   Z8_=invd0*invd1*invd2*invd3*invd4*invd5*invd6*invd7*invd8*cpowq(
   invd9,3);
   Z4_=Z8_ + Z7_ + Z6_ + Z4_ + Z5_;


	*out = cpowq(2.q*M_PIq*I,2)/(1)*(2*Z4_);
}
void evaluate_LTD_1_f128(__complex128 lm[], __complex128 params[], int conf, __complex128* out) {
   switch(conf) {
		case 0: evaluate_LTD_1_0_f128(lm, params, out); return;
		default: *out = 0.q;
    }
}
