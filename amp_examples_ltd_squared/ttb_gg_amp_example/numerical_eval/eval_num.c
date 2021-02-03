#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <complex.h>
#include <signal.h>
#include "externals.h"

// External functions

int evaluate(double complex lm[], double complex params[], int diag, int conf, double complex *out);

int get_rank(int diag, int conf);

// Internal functions
void print_moms(double complex moms[], int moms_size);
int form_lm(double complex ext_moms[],double complex loop_moms[],  double complex pol[],   double complex cpol[],
        double complex spinor_u[], double complex spinor_ubar[],   double complex spinor_v[], double complex spinor_vbar[],
        int n_ext ,    int n_loops ,    int n_pol,    int n_cpol,    int n_su,    int n_subar ,    int n_sv,
        int n_svbar, double complex lm[]);

bool next_sign(int sign[], int sign_size, int skip_first_n);

int main(int argc, char *argv[])
{

    
    //parameters
    const int n_diags = 2;
    const int n_config = 3; // we actually only need one I think.
    const int n_ext = 2;
    const int n_loops = 1;
    const int n_pol = 0;
    const int n_cpol = 2;
    const int n_su = 1;
    const int n_subar = 0;
    const int n_sv = 0;
    const int n_svbar = 1;


    double complex ext_moms[] =
        {
            p1[0], p1[1], p1[2], p1[3],     //
            p2[0], p2[1], p2[2], p2[3]     //
        };
    
    double complex loop_moms[] =
        {
            k1[0], k1[1], k1[2], k1[3]
        };
    double complex pol[] = {};
    double complex cpol[] =
        {
            ceps1[0], ceps1[1], ceps1[2], ceps1[3],     //
            ceps2[0], ceps2[1], ceps2[2], ceps2[3]     //
        };
    double complex spinor_u[] =
        {
            spinu[0], spinu[1], spinu[2], spinu[3]     //
        };
    double complex spinor_ubar[] =
        {
            spinu[0], spinu[1], spinu[2], spinu[3]     //
        };
    double complex spinor_v[] =  {};
    double complex spinor_vbar[] =  
        {
            spinvbar[0], spinvbar[1], spinvbar[2], spinvbar[3]     //
        };
    //Compute scalar products
    // external mom engergies + scalar prods + spatial scalar props
    int size = n_ext + 2*n_ext * (n_ext + 1) / 2 ;

 
    // loop mom lms
    size += n_loops +2*n_loops*n_ext + 2*n_loops*(n_loops+1)/2;

    
    // spatial components momenta
    size += 3*n_ext + 3*n_loops;

    
    // polarization vectors
    size += 4*n_pol + n_pol*(n_pol+1)/2 + 2*n_pol*n_loops + n_pol*n_ext + n_pol*n_cpol;

    
    // conjugate pol vectors 
    size += 4*n_cpol + n_cpol*(n_cpol+1)/2+ 2*n_cpol*n_loops + n_cpol*n_ext;

    // spinors
    size += 4*(n_su + n_subar + n_sv + n_svbar);


    double complex lm[size];
    int lm_size = form_lm(ext_moms,loop_moms,pol,cpol,spinor_u,spinor_ubar, //
        spinor_v,spinor_vbar,n_ext,n_loops,n_pol,n_cpol,n_su,n_subar,n_sv,n_svbar,lm);
    // for (int i = 0; i < lm_size; i++)
    //  {
    //      printf("lm[%d]: %f+ %f*i, ",i,creal(lm[i]), cimag(lm[i]));
    //     printf("\n");
    //  }
    
    if (size != lm_size)
    {
        printf("lm_size don't match with inputs:");
        printf("\tlm_size: %d computed_size:%d, (ext,n_loops) = (%d,%d)", lm_size,size, n_ext, n_loops);
        exit(1);
    }
    printf("lm_size: %d\n", lm_size);
    
    double complex param[] = {CMPLX(1.0, 0.0), CMPLX(1.0, 0.0), CMPLX(1.0, 0.0), CMPLX(1.0, 0.0)};
    
    

    for (int dia = 0; dia<n_diags; dia++)
    {
        double complex res =CMPLX(0.0, 0.0);
        double complex sres = CMPLX(0.0, 0.0);
        
        //for (int ci = 0; ci < bsize; ci++)
        for (int con = 0; con<n_config; con++)
        {        
            evaluate(lm,param, dia, con, &res);
            printf("diag:%d, config:%d, res: %f%+fi\n",dia,con, creal(res), cimag(res));
            sres+= res;
            
        }       
        printf("\033[1m   RESULT : diag %d: %+.15f%+.15f*i\033[0m\n", dia, creal(sres), cimag(sres));
        printf("\033[1m   RATIO WITH MATHEMATICA : diag %d: %+.15f%+.15f*i\033[0m\n", dia, creal(sres)/creal(diags[dia]), cimag(sres)/cimag(diags[dia]));
        
        printf("\n");
    }

    return 0;
}

double complex dot(double complex p[], double complex q[])
{
    return q[0] * p[0] - q[1] * p[1] - p[2] * q[2] - p[3] * q[3];
}
double complex spatial_dot(double complex p[], double complex q[])
{
    return q[1] * p[1] + p[2] * q[2] + p[3] * q[3];
}

int form_lm(double complex ext_moms[],double complex loop_moms[],  double complex pol[],   double complex cpol[],
        double complex spinor_u[], double complex spinor_ubar[],   double complex spinor_v[], double complex spinor_vbar[],
        int n_ext ,    int n_loops ,    int n_pol,    int n_cpol,    int n_su,    int n_subar ,    int n_sv,
        int n_svbar, double complex lm[])
{
    int pos = 0;
    // UNWRAP FORM INPUT
    
    // Add lm for externals only
    for (int i = 0; i < n_ext; i++)
    {   
        double complex p[4] = {ext_moms[i * 4], ext_moms[i * 4 + 1], ext_moms[i * 4 + 2], ext_moms[i * 4 + 3]};
        lm[pos++] = p[0];
        for (int j = i; j < n_ext; j++)
        {
            double complex q[4] = {ext_moms[j * 4], ext_moms[j * 4 + 1], ext_moms[j * 4 + 2], ext_moms[j * 4 + 3]};
            lm[pos++] = dot(p, q);
            lm[pos++] = spatial_dot(p, q);
        }
    }

    // Add lm for loop and externals
    for (int i = 0; i < n_loops; i++)
    {   
        double complex k[4] = {loop_moms[i * 4], loop_moms[i * 4 + 1], loop_moms[i * 4 + 2], loop_moms[i * 4 + 3]};
        lm[pos++] = k[0];
        for (int j = 0; j < n_ext; j++)
        {
            double complex p[4] = {ext_moms[j * 4], ext_moms[j * 4 + 1], ext_moms[j * 4 + 2], ext_moms[j * 4 + 3]};
            lm[pos++] = dot(k, p);
            lm[pos++] = spatial_dot(k, p);
        }
        for (int j = i; j < n_loops; j++)
        {
            double complex q[4] = {loop_moms[j * 4], loop_moms[j * 4 + 1], loop_moms[j * 4 + 2], loop_moms[j * 4 + 3]};
            lm[pos++] = dot(k, q);
            lm[pos++] = spatial_dot(k, q);
        }
    }

    // Add lm for spatial components of external
    for (int i = 0; i < n_ext; i++)
    {   
        lm[pos++] = ext_moms[i * 4 + 1];
        lm[pos++] = ext_moms[i * 4 + 2];
        lm[pos++] = ext_moms[i * 4 + 3];        
    }
    // Add lm for spatial components of loop-momenta
    for (int i = 0; i < n_loops; i++)
    {   
        lm[pos++] = loop_moms[i * 4 + 1];
        lm[pos++] = loop_moms[i * 4 + 2];
        lm[pos++] = loop_moms[i * 4 + 3];        
    }
    // printf("lm after momenta: %d",pos);
    // printf("\n");
    // Add lm for polarizations
    for (int i = 0; i < n_pol; i++)
    {   
        double complex k[4] = {pol[i * 4], pol[i * 4 + 1], pol[i * 4 + 2], pol[i * 4 + 3]};
        lm[pos++] = k[0];
        lm[pos++] = k[1];
        lm[pos++] = k[2];
        lm[pos++] = k[3];
        for (int j = i; j < n_pol; j++)
        {
            double complex q[4] = {pol[j * 4], pol[j * 4 + 1], pol[j * 4 + 2], pol[j * 4 + 3]};
            lm[pos++] = dot(k, q);
        }
        for (int j = 0; j < n_loops; j++)
        {
            double complex q[4] = {loop_moms[j * 4], loop_moms[j * 4 + 1], loop_moms[j * 4 + 2], loop_moms[j * 4 + 3]};
            lm[pos++] = dot(k, q);
            lm[pos++] = spatial_dot(k, q);
        }
        for (int j = 0; j < n_ext; j++)
        {
            double complex q[4] = {ext_moms[j * 4], ext_moms[j * 4 + 1], ext_moms[j * 4 + 2], ext_moms[j * 4 + 3]};
            lm[pos++] = dot(k, q);
        }
        for (int j = 0; j < n_cpol; j++)
        {
            double complex q[4] = {cpol[j * 4], cpol[j * 4 + 1], cpol[j * 4 + 2], cpol[j * 4 + 3]};
            lm[pos++] = dot(k, q);
        }
    }
    // printf("lm after EPS: %d",pos);
    // printf("\n");
    
    // Add lm for conjugated polarizations
    for (int i = 0; i < n_cpol; i++)
    {   
        double complex k[4] = {cpol[i * 4], cpol[i * 4 + 1], cpol[i * 4 + 2], cpol[i * 4 + 3]};
        lm[pos++] = k[0];
        lm[pos++] = k[1];
        lm[pos++] = k[2];
        lm[pos++] = k[3];
        for (int j = i; j < n_cpol; j++)
        {
            double complex q[4] = {cpol[j * 4], cpol[j * 4 + 1], cpol[j * 4 + 2], cpol[j * 4 + 3]};
            lm[pos++] = dot(k, q);
        }
        for (int j = 0; j < n_loops; j++)
        {
            double complex q[4] = {loop_moms[j * 4], loop_moms[j * 4 + 1], loop_moms[j * 4 + 2], loop_moms[j * 4 + 3]};
            lm[pos++] = dot(k, q);
            lm[pos++] = spatial_dot(k, q);
        }
        for (int j = 0; j < n_ext; j++)
        {
            double complex q[4] = {ext_moms[j * 4], ext_moms[j * 4 + 1], ext_moms[j * 4 + 2], ext_moms[j * 4 + 3]};
            lm[pos++] = dot(k, q);
        }
    }

    // printf("lm after CEPS: %d",pos);
    // printf("\n");
    // Add lm for spinor components
    for (int i = 0; i < n_sv; i++)
    {   
        lm[pos++] = spinor_v[i*4];
        lm[pos++] = spinor_v[i*4+1];
        lm[pos++] = spinor_v[i*4+2];
        lm[pos++] = spinor_v[i*4+3];
    }
    for (int i = 0; i < n_svbar; i++)
    {   
        lm[pos++] = spinor_vbar[i*4];
        lm[pos++] = spinor_vbar[i*4+1];
        lm[pos++] = spinor_vbar[i*4+2];
        lm[pos++] = spinor_vbar[i*4+3];
    }
    for (int i = 0; i < n_su; i++)
    {   
        lm[pos++] = spinor_u[i*4];
        lm[pos++] = spinor_u[i*4+1];
        lm[pos++] = spinor_u[i*4+2];
        lm[pos++] = spinor_u[i*4+3];
    }
    for (int i = 0; i < n_subar; i++)
    {   
        lm[pos++] = spinor_ubar[i*4];
        lm[pos++] = spinor_ubar[i*4+1];
        lm[pos++] = spinor_ubar[i*4+2];
        lm[pos++] = spinor_ubar[i*4+3];
    }

    // printf("lm after SPINORS: %d",pos);
    // printf("\n");
    return pos;
}

void print_moms(double complex moms[], int moms_size)
{
    for (int i = 0; i < moms_size; i++)
    {
        if ((i) % 4 == 0)
            printf("\t");
        printf("%f+%fi, ", creal(moms[i]), cimag(moms[i]));
        if ((i + 1) % 4 == 0)
            printf("\n");
    }
    printf("\n");
}

