#include <tgmath.h>
# include <quadmath.h>
# include <signal.h>
# include "numerator.h"


static inline int evaluate_0_0_f128(__complex128 lm[], __complex128 params[], __complex128* out) {
	



	out[0] = 1;
	return 1;
}

int evaluate_0_f128(__complex128 lm[], __complex128 params[], int conf, __complex128* out) {
    switch(conf) {
		case 0: return evaluate_0_0_f128(lm, params, out);
		default: out[0] = 0.; return 1;
    }
}

