
# include <tgmath.h>
# include <quadmath.h>
# include <signal.h>

void evaluate_PF_0(double complex[], double complex[], int conf, double complex* out);
void evaluate_LTD_0(double complex[], double complex[], int conf, double complex* out);
void evaluate_PF_0_f128(__complex128[], __complex128[], int conf, __complex128* out);
void evaluate_LTD_0_f128(__complex128[], __complex128[], int conf, __complex128* out);

void evaluate_PF(double complex lm[], double complex params[], int diag, int conf, double complex* out) {
    switch(diag) {
		case 0: evaluate_PF_0(lm, params, conf, out); return;
		default: raise(SIGABRT);
    }
}

void evaluate_PF_f128(__complex128 lm[], __complex128 params[], int diag, int conf, __complex128* out) {
    switch(diag) {
		case 0: evaluate_PF_0_f128(lm, params, conf, out); return;
		default: raise(SIGABRT);
    }
}

void evaluate_LTD(double complex lm[], double complex params[], int diag, int conf, double complex* out) {
    switch(diag) {
		case 0: evaluate_LTD_0(lm, params, conf, out); return;
		default: raise(SIGABRT);
    }
}

void evaluate_LTD_f128(__complex128 lm[], __complex128 params[], int diag, int conf, __complex128* out) {
    switch(diag) {
		case 0: evaluate_LTD_0_f128(lm, params, conf, out); return;
		default: raise(SIGABRT);
    }
}

