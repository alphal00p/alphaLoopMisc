#-
ON statistics;

***********************************
* COMMENTS AMPLITUDES *
* predefined varaible
* The numbers here are PDGs. 
* See also: SM.dat for relevant qgraf
***********************************
#define GLU "21"
#define PHO "22"
#define EP "-11"
#define EM "11"
#define H "25"
#define GHO "82"
#define GHOBAR "-82"
#define FERM "-6,-5,-4,-3,-2,-1,1,2,3,4,5,6,-11,11,-12,12,-13,13"
#define Q "1,2,3,4,5,6"
#define QBAR "-1,-2,-3,-4,-5,-6"
#define L "11,12,13"
#define LBAR "-11,-12,-13"

Auto S mass;
CTable masses(-30:30);
CTable charges(-30:30);

#ifndef `OPTIMITERATIONS'
    #define OPTIMITERATIONS "100"
#endif

#ifndef `HEAVYFERMIONS'
Fill masses(1) = 0;
Fill masses(2) = 0;
Fill masses(3) = 0;
Fill masses(4) = 0;
Fill masses(5) = 0;
Fill masses(-1) = 0;
Fill masses(-2) = 0;
Fill masses(-3) = 0;
Fill masses(-4) = 0;
Fill masses(-5) = 0;
Fill masses(11) = 0;
Fill masses(12) = 0;
Fill masses(13) = 0;
Fill masses(-11) = 0;
Fill masses(-12) = 0;
Fill masses(-13) = 0;
#else
Fill masses(1) = mass_d;
Fill masses(2) = mass_u;
Fill masses(3) = mass_c;
Fill masses(4) = mass_s;
Fill masses(5) = mass_b;
Fill masses(-1) = mass_d;
Fill masses(-2) = mass_u;
Fill masses(-3) = mass_c;
Fill masses(-4) = mass_s;
Fill masses(-5) = mass_b;
Fill masses(11) = mass_e;
Fill masses(12) = mass_mu;
Fill masses(13) = mass_tau;
Fill masses(-11) = mass_e;
Fill masses(-12) = mass_mu;
Fill masses(-13) = mass_tau;
#endif

Fill masses(6) = mass_t;
Fill masses(-6) = mass_t;
Fill masses(25) = mass_h;

Fill charges(1) = -1/3; * d
Fill charges(2) = 2/3; * u
Fill charges(3) = -1/3; * s
Fill charges(4) = 2/3; * c
Fill charges(5) = -1/3; * b
Fill charges(6) = 2/3; * t
Fill charges(11) = -1; * e
Fill charges(-1) = 1/3; * d
Fill charges(-2) = -2/3; * u
Fill charges(-3) = 1/3; * s
Fill charges(-4) = -2/3; * c
Fill charges(-5) = 1/3; * b
Fill charges(-6) = -2/3; * t
Fill charges(-11) = 1; * e


***************************************************************
* COMMENTS AMPLITUDES *
* By convention, we assume loop-momenta to be k, while external
* ones are denoted by p
***************************************************************
V p1,...,p40,k1,...,k40; * force this internal ordering in FORM
Auto V p,k;
Auto S lm,ext;
Auto I mu=4,s=4;
Symbol ge, gs, gy, ghhh, type, in, out, virtual;
Auto S x, idx, t;

Set dirac: s1,...,s40;
Set lorentz: mu1,...,mu40;
Set lorentzdummy: mud1,...,mud40;

CF gamma, vector,g(s),delta(s),T, counter,color, prop;
CF f,vx, vec;
CF subs, configurations, conf, der, energy, spatial(s);
CF subgraph, uvconf, uvconf1, uvprop, uv;
Set ts: t0,...,t20;
CT penergy;
Symbol ca,cf,nf,[dabc^2/n],[d4RR/n],[d4RA/n],[d4AA/n];

S  i, m, n, ALARM;

#include- diacolor.h
Set colF: cOli1,...,cOli40;
Set colA: cOlj1,...,cOlj40;
Set colAdum: cOljj1,...,cOljj40;

* Load the diagrams
#include- input_`SGID'.h

************************************************
* COMMENTS AMPLITUDES *
* Substitute the Feynman rules for the numerator
* We expect QGRAF to deliver the dummy indices
* see the .sty's for QGRAF generation
* the used one is: orientedGraphPython.sty
************************************************

* Fix a quirk where 0 does not match to a vector
* The only 0 in a propagator or vertex is when a momentum is 0
* All indices and pdgs are non-zero
repeat id prop(?a, 0, ?b) = prop(?a, pzero, ?b);
repeat id vx(?a, 0, ?b) = vx(?a, pzero, ?b);

* ONLY RELEVANT FOR INTERFERENCES, NEEDS TO BE REPLACED BY EXTERNAL STATES
* do the spin sum external particles
repeat id prop(x?{`L'}, in, p?, idx1?)*prop(x?{`L',}, out, p?, idx2?) = gamma(dirac[idx1], p, dirac[idx2]) + masses(x)*gamma(dirac[idx1], dirac[idx2]);
repeat id prop(x?{`Q'}, in, p?, idx1?)*prop(x?{`Q'}, out, p?, idx2?) = d_(colF[idx1], colF[idx2])*(gamma(dirac[idx1], p, dirac[idx2]) + masses(x)*gamma(dirac[idx1], dirac[idx2]));
repeat id prop(x?{`LBAR'}, out, p?, idx1?)*prop(x?{`LBAR'}, in, p?, idx2?) = gamma(dirac[idx1], p, dirac[idx2]) - masses(x)*gamma(dirac[idx1], dirac[idx2]);
repeat id prop(x?{`QBAR'}, out, p?, idx1?)*prop(x?{`QBAR'}, in, p?, idx2?) = d_(colF[idx1], colF[idx2])*(gamma(dirac[idx1], p, dirac[idx2]) - masses(x)*gamma(dirac[idx1], dirac[idx2]));

* virtual edges
id prop(`GLU', virtual, p?, idx1?, idx2?) = - i_ * d_(lorentz[idx1], lorentz[idx2]) * d_(colA[idx1], colA[idx2]);
id prop(x?{`GHO',`GHOBAR'}, virtual, p?, idx1?, idx2?) = - i_ *d_(colA[idx1], colA[idx2]);
id prop(`PHO', virtual, p?, idx1?, idx2?) = - i_ * d_(lorentz[idx1], lorentz[idx2]);
id prop(x?{`L'}, virtual, p?, idx1?, idx2?) = - i_ * (gamma(dirac[idx2], p, dirac[idx1]) + masses(x) * gamma(dirac[idx2], dirac[idx1]));
id prop(x?{`LBAR'}, virtual, p?, idx1?, idx2?) = - i_ * (gamma(dirac[idx1], p, dirac[idx2]) + masses(x) * gamma(dirac[idx1], dirac[idx2]));
id prop(x?{`Q'}, virtual, p?, idx1?, idx2?) = - i_ * (gamma(dirac[idx2], p, dirac[idx1]) + masses(x) * gamma(dirac[idx2], dirac[idx1])) * d_(colF[idx2], colF[idx1]);
id prop(x?{`QBAR'}, virtual, p?, idx1?, idx2?) = - i_ * (gamma(dirac[idx1], p, dirac[idx2]) + masses(x) * gamma(dirac[idx1], dirac[idx2])) * d_(colF[idx1], colF[idx2]);
id prop(`H', virtual, p?, idx1?, idx2?) = -i_;
id prop(`H', in, p?, idx1?) = 1;
id prop(`H', out, p?, idx1?) = 1;

.sort:feynman-rules-edges;

Print +s; .sort
* vertices
id vx(x1?{`QBAR'}, `GLU', x2?{`Q'}, p1?, p2?, p3?, idx1?, idx2?, idx3?) = -gs * gamma(dirac[idx1], lorentz[idx2], dirac[idx3]) * T(colF[idx1], colA[idx2], colF[idx3]);
id vx(`GHOBAR', `GLU', `GHO', p1?, p2?, p3?, idx1?, idx2?, idx3?) = -gs * i_ * cOlf(colA[idx3], colA[idx2], colA[idx1]) * (1/2) * (p3(lorentz[idx2])-p1(lorentz[idx2]));
id vx(x1?{`QBAR'}, `PHO', x2?{`Q'}, p1?, p2?, p3?, idx1?, idx2?, idx3?) = charges(x2) * ge * i_* gamma(dirac[idx1], lorentz[idx2], dirac[idx3]) * d_(colF[idx1], colF[idx3]);
id vx(x1?{`LBAR'}, `PHO', x2?{`L'}, p1?, p2?, p3?, idx1?, idx2?, idx3?) = charges(x2) * ge * i_* gamma(dirac[idx1], lorentz[idx2], dirac[idx3]);
id vx(x1?{`QBAR'}, `H', x2?{`Q'}, p1?, p2?, p3?, idx1?, idx2?, idx3?) = -gy * i_ * d_(dirac[idx1], dirac[idx3]) * d_(colF[idx1], colF[idx3]);
id vx(x1?{`LBAR'}, `H', x2?{`L'}, p1?, p2?, p3?, idx1?, idx2?, idx3?) = -gy * i_ * d_(dirac[idx1], dirac[idx3]);
id vx(`H', `H', `H', p1?, p2?, p3?, idx1?, idx2?, idx3?) = -ghhh * i_;

* TODO: use momentum conservation to reduce the number of different terms
id vx(`GLU', `GLU', `GLU', p1?, p2?, p3?, idx1?, idx2?, idx3?) = i_ * gs * cOlf(colA[idx1], colA[idx2], colA[idx3]) *(
    - d_(lorentz[idx1], lorentz[idx3]) * p1(lorentz[idx2])
    + d_(lorentz[idx1], lorentz[idx2]) * p1(lorentz[idx3])
    + d_(lorentz[idx2], lorentz[idx3]) * p2(lorentz[idx1])
    - d_(lorentz[idx1], lorentz[idx2]) * p2(lorentz[idx3])
    - d_(lorentz[idx2], lorentz[idx3]) * p3(lorentz[idx1])
    + d_(lorentz[idx1], lorentz[idx3]) * p3(lorentz[idx2]));

* For the quartic gluon vertex we need an extra dummy index
Print +s; .sort
Multiply counter(1);
repeat id vx(`GLU', `GLU', `GLU', `GLU', p1?, p2?, p3?, p4?, idx1?, idx2?, idx3?, idx4?)*counter(i?) = - counter(i + 1) * gs^2 * i_ *(
    + cOlf(colAdum[i], colA[idx1], colA[idx2]) * cOlf(colA[idx3], colA[idx4], colAdum[i])
        * (d_(lorentz[idx1], lorentz[idx4]) * d_(lorentz[idx2], lorentz[idx3]) - d_(lorentz[idx1], lorentz[idx3]) * d_(lorentz[idx2], lorentz[idx4]))
    + cOlf(colAdum[i], colA[idx1], colA[idx3]) * cOlf(colA[idx2], colA[idx4], colAdum[i])
        * (d_(lorentz[idx1], lorentz[idx4]) * d_(lorentz[idx2], lorentz[idx3]) - d_(lorentz[idx1], lorentz[idx2]) * d_(lorentz[idx3], lorentz[idx4]))
    + cOlf(colAdum[i], colA[idx1], colA[idx4]) * cOlf(colA[idx2], colA[idx3], colAdum[i])
        * (d_(lorentz[idx1], lorentz[idx3]) * d_(lorentz[idx2], lorentz[idx4]) - d_(lorentz[idx1], lorentz[idx2]) * d_(lorentz[idx3], lorentz[idx4]))
);
id counter(x?) = 1;

if (count(vx, 1, prop, 1));
    Print "Unsubstituted propagator or vertex: %t";
    exit "Critical error";
endif;

Print +s; .sort

.sort:feynman-rules-vertices;

******************
* Color evaluation
******************
repeat id T(cOli1?,?a,cOli2?)*T(cOli2?,?b,cOli3?) = T(cOli1,?a,?b,cOli3); * collect the colour string
id  T(cOli1?, ?a, cOli1?) = cOlTr(?a);
id  cOlTr(cOli1?) = 0;
id  cOlTr = cOlNR;
Multiply color(1);
repeat id cOlTr(?a)*color(x?) = color(x * cOlTr(?a));
repeat id cOlf(cOlj1?,cOlj2?,cOlj3?)*color(x?) = color(x * cOlf(cOlj1,cOlj2,cOlj3));
B+ color;
.sort:color-prep;
Keep brackets;
* TODO: do we need to set the dimension?
* Only evaluate this part per unique color by bracketing
Argument color;
    #call color
    #call simpli
    id  cOlI2R = cOlcR*cOlNR/cOlNA;
    id  cOlNR/cOlNA*cOlcR = cOlI2R;
    id  cOld33(cOlpR1,cOlpR2) = [dabc^2/n]*cOlNR;
    id  cOlNR/cOlNA = nf/cf/2;
    id  cOlcR = cf;
    id  cOlcA = ca;
    id  cOlI2R = nf/2;
	id	cOld44(cOlpA1,cOlpA2) = [d4AA/n]*cOlNA;
	id	cOld44(cOlpR1,cOlpR2) = [d4RR/n]*cOlNA;
	id	cOld44(cOlpR1,cOlpA1) = [d4RA/n]*cOlNA;

* set the SU(3) values
    id [dabc^2/n] = 13/18;
    id cOlNR = 3;
    id cOlNA = 8;
    id cf = 4 / 3;
    id ca = 3;
    id nf = 1;
EndArgument;
.sort:color;

Print +s; .sort

* set the SU(3) values
id cOlNR = 3;

* construct gamma string
repeat id gamma(s1?,?a,s2?)*gamma(s2?,?b,s3?) = gamma(s1,?a,?b,s3);

Multiply counter(1);
repeat id gamma(mu?,?a,p?,?b,mu?)*counter(i?) = vec(p,lorentzdummy[i])*gamma(mu,?a,lorentzdummy[i],?b,mu)*counter(i+1);
id counter(n?) = 1;

*********************************************************
* COMMENTS AMPLITUDES *
* WE WILL NEED TO INSERT THE PROPER CLIFFORD ALGEBRA HERE
* AND AN EXPANSION OF GAMMA CHAINS IN TERMS OF COMPONETS
*********************************************************

#do i=1,10
    id once gamma(mu?,?a,mu?) = g_(`i',?a);
#enddo

id vec(p?,mu?) = p(mu);

#do i=1,10
    trace4 `i';
#enddo

if (count(gamma, 1));
    Print "Unsubstituted gamma string: %t";
    exit "Critical error";
endif;

.sort:gamma-traces;

Print +s; .sort

id color(x?) = x;


id pzero = 0; * Substitute the 0-momentum by 0
.sort:feynman-rules-final;

* If the expression is empty (due to color), we still write a file
#if ( termsin(F) == 0 )
    #write<out_`SGID'.proto_c> "#0 due to color\n"
#endif

*************************************************
* Process different configurations 
*************************************************

* process all the configurations
id configurations(x?) = x;

* multiply the numerator contribution of derivatives
id conf(?a,p?) = conf(?a) * penergy(p); 
id conf(?a,x?) = conf(?a) * x; * note the type difference


id subs(p1?,p2?) = replace_(p1, p2);

B+ uv;
.sort:bubble-treatment;
Keep brackets;


**************************************************
* COMMENTS AMPLITUDES *
* FROM HERE ON, WE WILL NEED TO REWRITE FOR 
* COMPONENTS OF MOMENTA, SINCE WE WONT DEAL 
* WITH (energy(k1)*energy(k2)-spatial(k1,k2))
* ONLY DUE TO THE SPIN CHAINS
*************************************************

* now extract the energy components of the LTD loop variables
id k1?.k2? = g(k1, k2);
repeat id conf(x?,?a,k1?,?b)*penergy(k1?) = conf(x,?a,k1,?b)*energy(k1); * detect loop energies from the bubble derivative
repeat id conf(x?,?a,k1?,?b)*g(k1?,k1?) = conf(x,?a,k1,?b)*(energy(k1)*energy(k1)-spatial(k1,k1));
repeat id conf(x?,?a,k1?,?b,k2?,?c)*g(k1?,k2?) = conf(x,?a,k1,?b,k2,?c)*(energy(k1)*energy(k2)-spatial(k1,k2));
repeat id conf(x?,?a,k?,?b)*g(k?,p?) = conf(x,?a,k,?b)*(energy(k)*penergy(p)-spatial(k,p));
id g(p1?,p2?) = p1.p2;

repeat id energy(?a)*energy(?b) = energy(?a,?b);
symmetrize energy;
id energy(?a) = energy(f(?a));
if (count(energy,1) == 0) Multiply energy(f(k0)); * signal with k0 that we are dealing with the constant term
.sort:energy-splitoff;

Print +s; .sort

*********************************************
* Construction of optimized numerator C code
*********************************************

* split off every energy configuration into a new expression
id conf(?a) = conf(conf(?a));
argtoextrasymbol tonumber,conf,1;
#redefine oldextrasymbols "`extrasymbols_'"
B+ conf;
.sort:conf-1;
Hide F;

#$MAXK = `NFINALMOMENTA';
#$MAXP = `NINITIALMOMENTA';

#redefine energysymbolstart "`extrasymbols_'"
#do ext={`oldextrasymbols'+1}, `extrasymbols_'
    #$tmp = extrasymbol_(`ext');
    #write<out_`SGID'.proto_c> "#CONF\n%$", $tmp;
    L FF`ext' = F[conf(`ext')];
    .sort:conf-2;
    
    argtoextrasymbol tonumber, energy, 1;
    B+ energy;
    .sort:energy-1;
    Hide FF`ext';
    #do ps={`energysymbolstart'+1}, `extrasymbols_'
        #$tmp = extrasymbol_(`ps');
        L FTMP = FF`ext'[energy(`ps')];
        .sort:energy-2;

        #if ( termsin(FTMP) > 0 )
            #write<out_`SGID'.proto_c> "#NEWMONOMIAL\n%$", $tmp;

* Convert the dot products and energies to a symbol
            #$OFFSET = 0;
            #do i=1,`$MAXP'
                id penergy(p`i') = lm`$OFFSET';
                #$OFFSET = $OFFSET + 1;
                #do j=`i',`$MAXP'
                    id p`i'.p`j' = lm`$OFFSET';
                    #$OFFSET = $OFFSET + 1;
                #enddo
            #enddo

            #do i=1,`$MAXK'
                id penergy(k`i') = lm`$OFFSET';
                #$OFFSET = $OFFSET + 1;
                #do j=1,`$MAXP'
                    id k`i'.p`j' = lm`$OFFSET';
                    #$OFFSET = $OFFSET + 1;
                    id spatial(p`j', k`i') = lm`$OFFSET';
                    #$OFFSET = $OFFSET + 1;
                #enddo

                #do j=`i',`$MAXK'
                    id k`i'.k`j' = lm`$OFFSET';
                    #$OFFSET = $OFFSET + 1;
                    id spatial(k`i', k`j') = lm`$OFFSET';
                    #$OFFSET = $OFFSET + 1;
                #enddo
            #enddo
            .sort:lm-subs;

* Optimize the output
            Format C;
            Format O4,stats=off,saIter=`OPTIMITERATIONS';
            #Optimize FTMP
            #write<out_`SGID'.proto_c> "%O\n\treturn %e",FTMP
            #clearoptimize;
            .sort:optim-`ext'-`ps';
            Format O0;
            Format normal;
        #endif
    #enddo
    Drop FF`ext';
#enddo
.end
