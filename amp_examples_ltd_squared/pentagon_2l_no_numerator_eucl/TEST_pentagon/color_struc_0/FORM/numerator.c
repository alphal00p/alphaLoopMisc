
# include <tgmath.h>
# include <quadmath.h>
# include <signal.h>





int evaluate(double complex lm[], double complex params[], int diag, int conf, double complex* out) {
    switch(diag) {
		default: raise(SIGABRT);
    }
}

int evaluate_f128(__complex128 lm[], __complex128 params[], int diag, int conf, __complex128* out) {
    switch(diag) {
		default: raise(SIGABRT);
    }
}

int get_buffer_size() {
    return 2;
}

int get_rank(int diag, int conf) {
    switch(diag) {
		default: raise(SIGABRT);
    }
}
