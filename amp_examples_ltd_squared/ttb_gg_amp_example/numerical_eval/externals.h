#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <complex.h>
#include <signal.h>
double complex p1[4] = {CMPLX(500.0000000000000, 0), CMPLX(0, 0), CMPLX(0, 0), CMPLX(469.1172561311297, 0)};
double complex p2[4] = {CMPLX(500.0000000000000, 0), CMPLX(0, 0), CMPLX(0, 0), CMPLX(-469.1172561311297, 0)};
double complex k1[4] = {CMPLX(499.9999999999998, 0), CMPLX(110.9242844438328, 0), CMPLX(444.8307894881214, 0), CMPLX(-199.5529299308788, 0)};
double complex ceps1[4] = {CMPLX(0, 0), CMPLX(0.707106781186548, 0), CMPLX(0, 0.7071067811865475), CMPLX(0, 0)};
double complex ceps2[4] = {CMPLX(0, 0), CMPLX(0.707106781186548, 0), CMPLX(0, -0.7071067811865475), CMPLX(0, 0)};
double complex spinu[4] = {CMPLX(5.5572244752997, 0), CMPLX(0, 0), CMPLX(31.1306481803886, 0), CMPLX(0, 0)};
double complex spinvbar[4] = {CMPLX(0, 0), CMPLX(31.1306481803886, 0), CMPLX(0, 0), CMPLX(5.5572244752997, 0)};
