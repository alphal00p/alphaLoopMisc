#include <tgmath.h>
# include <quadmath.h>
# include <signal.h>
# include "numerator.h"

// polynomial in 0,c3
static inline int evaluate_0_0(double complex lm[], double complex params[], double complex* out) {
	


	out[0] = 1;
	return 1; 
}

int evaluate_0(double complex lm[], double complex params[], int conf, double complex* out) {
   switch(conf) {
		case 0: return evaluate_0_0(lm, params, out);
		default: out[0] = 0.; return 1;
    }
}

int get_rank_0(int conf) {
   switch(conf) {
		case 0: return 0;
		default: return 0;
    }
}
