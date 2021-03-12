
# include <tgmath.h>
# include <quadmath.h>
# include <signal.h>

int evaluate_1(double complex[], double complex[], int conf, double complex*);
int evaluate_0(double complex[], double complex[], int conf, double complex*);
int evaluate_1_f128(__complex128[], __complex128[], int conf, __complex128*);
int evaluate_0_f128(__complex128[], __complex128[], int conf, __complex128*);
int get_rank_1(int conf);
int get_rank_0(int conf);

int evaluate(double complex lm[], double complex params[], int diag, int conf, double complex* out) {
    switch(diag) {
		case 1: return evaluate_1(lm, params, conf, out);
		case 0: return evaluate_0(lm, params, conf, out);
		default: raise(SIGABRT);
    }
}

int evaluate_f128(__complex128 lm[], __complex128 params[], int diag, int conf, __complex128* out) {
    switch(diag) {
		case 1: return evaluate_1_f128(lm, params, conf, out);
		case 0: return evaluate_0_f128(lm, params, conf, out);
		default: raise(SIGABRT);
    }
}

int get_buffer_size() {
    return 16;
}

int get_rank(int diag, int conf) {
    switch(diag) {
		case 1: return get_rank_1(conf);
		case 0: return get_rank_0(conf);
		default: raise(SIGABRT);
    }
}
