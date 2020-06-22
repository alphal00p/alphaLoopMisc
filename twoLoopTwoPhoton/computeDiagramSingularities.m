(* ::Package:: *)

evalFeynAmp::usage = "evalFeynAmp[amplitude, externalKinematics, internalKinematics, scalingRule] evaluates the amplitude at the user-
specified numerical values of external and internal momenta. For advanced usage for computing scaling behavior in IR / UV limits,
scalingRule can be used, but otherwise can be set to the empty list {}.";


allFeynAmps::usage = "allFeynAmps is a list of expressions whose sum is the subtracted finite integrand for qqbar -> \!\(\*SuperscriptBox[\(\[Gamma]\), \(*\)]\) + \!\(\*SuperscriptBox[\(\[Gamma]\), \(*\)]\).";


(* ::Section:: *)
(*Settings*)


On[Assert];


(*subUVfudgeFactor = 1;*)
(*commend out this line to track the effect of subUV counterterms on the double-UV counterterms for
modified bubble / vertex diagrams*)


(*scheme = Input["Please choose the scheme, GS or gammaHat\n"];*)
scheme = GS;


Assert[scheme == GS || scheme == gammaHat];


$runNumericalChecks = False;


timingPrint[{a_?NumericQ, b_}] := (Print["Computation time: ", a]; Print["result: ", b]);
SetAttributes[myTiming, HoldAll];
myTiming[expr_] := timingPrint[AbsoluteTiming[expr]];


externalKinematics = {
p1 -> {1,0,0,1},
p2 -> {1,0,0,-1},
q1 -> {1, 0, 2, 1},
q2 -> {1, 0, -2, -1},
(*This is u for incoming fermion*)
spinor[-1] -> {0, -7, 3, 0}(*{{0}, {-7}, {3}, {0}}*),
(*This is uBar for incoming anti-fermion*)
spinor[-3] -> {0, 13, -9, 0}(*{{0, 13, -9, 0}}*),
(*photon q1*)
polVec[-2] -> {4, -4, 1, 2},
(*Photon q2*)
polVec[-4] -> {2, -1, -5, 8}
};

randomNumerics = Thread[
    {k[1,minus],k[1,plus],k[1,x],k[1,y],k[2,minus],k[2,plus],k[2,x],k[2,y], M} ->
    {33/17, -48/89, 21/23, 21/41, 47/23, -7/61, -37/73, -39/67, 5/3}
];


Assert[
    MatchQ[p1/.externalKinematics, {a_, 0, 0, a_}] && MatchQ[p2/.externalKinematics, {a_, 0, 0, b_}/;a+b==0]
];


(* Massive reference vector for projecting out longitudinal polarizations*)
refEta = 2(p1 + p2);


(*If set to true, perform kT \[Rule] -kT symmetrization for both loops, to further suppress various IR divergences. *)
$transvereSym = True;


(* ::Section::Closed:: *)
(*Feynman rules*)


insertFeynmanRules[diag_] := diag //. {uvSubtracted[a_]:>a, uvSubtracted[a__]:>Sequence[a]} /. gluon1 -> gluon /.
    Join[vertexBubbleRule, vertexBubbleSubtractedRule] /.
    Join[
        externalStateRules,
        propRules,
        vertexRules,
        projectorRules,
        projectOutScalarRule,
        fermionBubbleFactorRule,
        {hardBubbleRule}
    ];


(*Without modifying vertices and bubbles*)
insertFeynmanRulesNoMod[diag_] := diag //. {uvSubtracted[a_]:>a, uvSubtracted[a__]:>Sequence[a]} /.
    vertexBubbleRuleNoMod /.
    Join[
        externalStateRules,
        propRules,
        vertexRules,
        projectorRules,
        projectOutScalarFakeRule,
        fermionBubbleFactorRule,
        {hardBubbleRule}
    ];


vertexBubbleRuleNoMod = {
    modifiedVertex1[_, gluon[label_, _], _,_] :> (vertex1original /. \[Mu]->lorentzIndex[label]),
    modifiedBubble1[_,_, gluon[label_, _], _] :> (bubble1original /. \[Mu]->lorentzIndex[label]),
    modifiedVertex2[_, gluon[label_, _], _,_] :> (vertex2original /. \[Mu]->lorentzIndex[label]),
    modifiedBubble2[_,_, gluon[label_, _], _] :> (bubble2original /. \[Mu]->lorentzIndex[label])
};


vertexBubbleRule = {
    modifiedVertex1[_, gluon[label_, _], _,_] :> (vertex1modified  /. \[Mu]->lorentzIndex[label]),
    modifiedBubble1[_,_, gluon[label_, _], _] :> (bubble1modified /. \[Mu]->lorentzIndex[label]),
    modifiedVertex2[_, gluon[label_, _], _,_] :> (vertex2modified /. \[Mu]->lorentzIndex[label]),
    modifiedBubble2[_,_, gluon[label_, _], _] :> (bubble2modified /. \[Mu]->lorentzIndex[label])
};


vertexBubbleSubtractedRule = {
    modifiedVertex1Subtracted[_, gluon[label_, _], _,_] :> (vertex1modified + subUVfudgeFactor*vertex1UV /. \[Mu]->lorentzIndex[label]),
    modifiedBubble1Subtracted[_,_, gluon[label_, _], _] :> (bubble1modified + subUVfudgeFactor*bubble1UV /. \[Mu]->lorentzIndex[label]),
    modifiedVertex2Subtracted[_, gluon[label_, _], _,_] :> (vertex2modified + subUVfudgeFactor*vertex2UV /. \[Mu]->lorentzIndex[label]),
    modifiedBubble2Subtracted[_,_, gluon[label_, _], _] :> (bubble2modified + subUVfudgeFactor*bubble2UV /. \[Mu]->lorentzIndex[label])
};


externalStateRules = {
   cpol[(quark|qbar)[number_, p_]] :> spinor[number],
   pol[photon[i_, q_]] :> 1 (*the nontrivial part will be absorbed into the definition of the photon-quark vertex*)
};


propRules = {
	prop[quark[i_, q_], qbar[j_, q_]] :> I * slash[q] * propagator[q],
	prop[gluon[i_, q_], gluon[j_, q_]] :> -I * metricTensor[lorentzIndex[i], lorentzIndex[j]] * propagator[q],
	propGluonUV[k_, {a1_, a2_}] :> -I * metricTensor[lorentzIndex[a1], lorentzIndex[a2]] * propagator[k, M],
	propFermionUV[k_] :> I * slash[k] * propagator[k, M]
};


(*3-gluon and 4-gluon vertices not yet implemented*)
vertexRules = {
	v3[qbar[i_, q1_], photon[j_, q2_], quark[k_, q3_]] :>
		-I * lambda * slash[polVec[j]],
	v3[qbar[i_, q1_], gluon[j_, q2_], quark[k_, q3_]] :> 
    	I * g * gamma[lorentzIndex[j]]
};


projectorRules = {
    projectorP2P1 -> 1/s * fermionChain[slash[p2], slash[p1]],
    projectorP1P2 -> 1/s * fermionChain[slash[p1], slash[p2]]
};


projectOutScalarRule = {projectOutScalar[i1_, i2_, momentum_] :>
    KroneckerDelta[lorentzIndex[i1], lorentzIndex[i2]] -
        component[refEta, lorentzIndex[i1]] * component[lower[momentum], lorentzIndex[i2]] / lorentzDot[refEta, momentum]
};


projectOutScalarFakeRule = {projectOutScalar[i1_, i2_, momentum_] :>
    KroneckerDelta[lorentzIndex[i1], lorentzIndex[i2]]
};


fermionBubbleFactorRule = {fermionBubbleFactor -> I * (2 *2(1-eps))/(3-2eps) * lambda^2 *
    (propagator[k1] * propagator[k1+k2] - propagator[k1, M]^2)
};


(* ::Section::Closed:: *)
(*Modified sub-diagrams*)


symmetrize[expr_, rule_] := (1/2 * expr) + (1/2 * expr /. rule);


(*Eq. 17 of "short singlet-IR-fact-20.01.14.pdf", with additional factor (-I)*I^6 = +I, adding a missing overall factor of -2,
and with a sign correction for the last term*)
vertex1original = -2*I*g^3 propagator[-k+p1] *
(
(1-eps) * (
    2 component[l,\[Mu]] * propagator[l] * propagator[-l+p1] *
    (
        1 * id -
        (lorentzDot[k] * id - fermionChain[slash[l],slash[k]]) * propagator[k-l]
    ) -
    2 component[l,\[Mu]] * id * propagator[k-l] * propagator[l] -
    (2 * (component[p1,\[Mu]] - component[k,\[Mu]]) * id + fermionChain[gamma[\[Mu]], slash[k]]) * propagator[k-l] * propagator[-l+p1]
) -
fermionChain[
    slash[p1] - slash[k],
    fermionChain[slash[l], gamma[\[Mu]], slash[k]] - eps * fermionChain[slash[k], gamma[\[Mu]], slash[l]]
] * propagator[k-l] * propagator[l] * propagator[-l+p1]
) /. {k->-k2, l->p1+k1};


(* Must use id to denote 4*4 Dirac Matrix!*)
vertex1Term1Extra := id*(I*g^3) * propagator[k2+p1] * 2 * (2-2*eps) *
       (1 + 2 * lorentzDot[p2, k2] / s) *
       (
          (
            2 * component[k1, \[Mu]] + 
            component[k2, \[Mu]]
          ) *
          lorentzDot[k1, p1] +
          component[p1, \[Mu]] *
          (
            -(1/2) lorentzDot[k1] +
            lorentzDot[k1, p1]
          )
       ) * propagator[k1] * propagator[k1+p1] * propagator[k1+k2+p1];


vertex1modifiedGS := -2*I* propagator[-k+p1]*
(
(1-eps) * (
    propagator[l] * propagator[-l+p1] *
    (
        component[p1, \[Mu]] * id -
        component[l, \[Mu]] * (lorentzDot[k] * id - fermionChain[slash[l], slash[k]]) * propagator[k-l] -
        (component[p1, \[Mu]] - component[l, \[Mu]]) * (lorentzDot[k] * id - fermionChain[slash[p1] - slash[l], slash[k]]) * propagator[p1-k-l]
    ) -
    component[k,\[Mu]] * id * propagator[k-l] * propagator[l] -
    (2 (component[p1,\[Mu]] - component[k,\[Mu]]) * id + fermionChain[gamma[\[Mu]], slash[k]]) * propagator[k-l] * propagator[-l+p1]
) -
fermionChain[
    slash[p1] - slash[k],
    fermionChain[slash[l], gamma[\[Mu]], slash[k]] - eps * fermionChain[slash[k], gamma[\[Mu]], slash[l]]
] * propagator[k-l] * propagator[l] * propagator[-l+p1]
) /. {k->-k2, l->p1+k1};


vertex1modified = Switch[scheme,
    gammaHat,
    vertex1original + vertex1Term1Extra,
    GS,
    vertex1modifiedGS,
    _,
    vertex1original
];


(*Eq 20 of "short singlet-IR-fact-20.01.14.pdf, with additional factor (-I)*I^6 = +I, and with correction
slash[p1-l] --> slash[l-k] "*)
bubble1original = I*g^3 * (-2)(1-eps)* fermionChain[slash[p1-k], slash[l-k], slash[p1-k], gamma[\[Mu]]] *
    propagator[p1-k] * propagator[l-k] * propagator[p1-l] * propagator[p1-k] /. {k->-k1, l->p1+k2};


plusPart[{p0_, p1_, p2_, p3_}] := {(p0+p3)/2, 0,0, (p0+p3)/2};
minusPart[{p0_, p1_, p2_, p3_}] := {(p0-p3)/2, 0,0, -(p0-p3)/2};
(*delayVec[{p0_, p1_, p2_, p3_}] := {p0, 0,0, p3};*)
(*Partial tensor reduction in Eq. 22, explained in George's email on 2020.01.21*)
minusReducedL1 := l - minusPart[l]+ (-1/2) * minusPart[k];
minusReducedL2 := l - minusPart[l]+ (+1/2) * minusPart[k];


bubble1modifiedGS := I * (1-eps) * propagator[p1-k] *
fermionChain[
    slash[p1-k-minusReducedL1] * propagator[l] * propagator[p1-l-k] -
        slash[2p1-k-minusReducedL2] * propagator[p1-l] * propagator[l-k],
    gamma[\[Mu]]
] /. {k->-k1, l->p1+k2};


bubble1modified = Switch[scheme,
    gammaHat,
    (*symmetrize[bubble1original, k2 -> -p1 - k1 - k2],*)
    (* The above expression is correct, but has linear divergence term-by-term, making it inconveninent for
       automatic generation of UV CTs. So we'll type it by hand.*)
    I*g^3 * (-2)(1-eps)* (1/2)(*from symmetrization*) fermionChain[slash[p1-k], gamma[\[Mu]]] *
    propagator[p1-k] * propagator[l-k] * propagator[p1-l] /. {k->-k1, l->p1+k2},
    GS,
    bubble1modifiedGS,
    _,
    bubble1original
];


fermionChainReverse[a__] := Reverse[fermionChain[a]];


bubble2original = bubble1original /. {fermionChain -> fermionChainReverse, p1->-p2};


bubble2modified = bubble1modified /. {fermionChain -> fermionChainReverse, p1->-p2, minusPart -> plusPart};


vertex2original = vertex1original /. {fermionChain -> fermionChainReverse, p1->-p2, p2->-p1, k1->k2, k2->k1};


vertex2modified = vertex1modified /. {fermionChain -> fermionChainReverse, p1->-p2, p2->-p1, k1->k2, k2->k1};


(* ::Section::Closed:: *)
(*Sub-UV counterterm for hard bubbles / vertices and form-factor vertices*)


hardBubbleRule = hardBubbleCT[p_, {a1_, a2_, k_}] :>(-I) * fermionChain[
    I * g * gamma[lorentzIndex[a1]],
    (
       I * fermionChain[slash[k]] * propagator[k, M]^2 -
       I * fermionChain[slash[k], slash[p], slash[k]] * propagator[k, M]^3
    ),
    I*g* gamma[lorentzIndex[a2]]
] * metricTensor[lorentzIndex[a1], lorentzIndex[a2]];


(* ::Section::Closed:: *)
(*Sub-UV counterterm for modified diagrams*)


(*Need to be corrected*)
(*bubble1UV = I* (1-eps) fermionChain[propagator[k2, M] propagator[k2, M] slash[-k2 + minusPart[k2]] -
    propagator[k2, M] propagator[k2, M] slash[-k2 + minusPart[k2]], gamma[\[Mu]]] propagator[k1+p1];*)


(*Let's copy the expression of vertex1modified again.*)
(*vertex1modified == -2 \[ImaginaryI] (
    -fermionChain[
        -slash[-k2]+slash[p1],
        -eps fermionChain[slash[-k2],gamma[\[Mu]],slash[k1+p1]] +
            fermionChain[slash[k1+p1],gamma[\[Mu]],slash[-k2]]
    ] * propagator[-k1] propagator[-k1-k2-p1] propagator[k1+p1] +
    (1-eps) * (
        -(2 id (-component[-k2,\[Mu]]+component[p1,\[Mu]]) + fermionChain[gamma[\[Mu]],slash[-k2]]) *
            propagator[-k1] propagator[-k1-k2-p1] -
        id component[-k2,\[Mu]] propagator[-k1-k2-p1] propagator[k1+p1] +
        propagator[-k1] *(
            id component[p1,\[Mu]] -
            (component[p1,\[Mu]]-component[k1+p1,\[Mu]]) *
                (-fermionChain[slash[p1]-slash[k1+p1],slash[-k2]] + id lorentzDot[-k2]) * propagator[-k1+k2] -
            component[k1+p1,\[Mu]] (-fermionChain[slash[k1+p1],slash[-k2]]+id lorentzDot[-k2]) propagator[-k1-k2-p1]
        ) * propagator[k1+p1]
    )
) * propagator[k2+p1];*)


If[scheme == GS,
vertex1UV = (-1)* (-2 I) (
    (1-eps) * (
        -(2 id (-component[-k2,\[Mu]]+component[p1,\[Mu]]) + fermionChain[gamma[\[Mu]],slash[-k2]]) *
            propagator[k1, M] propagator[k1, M] -
        id component[-k2,\[Mu]] * propagator[k1, M] propagator[k1, M] +
        propagator[k1, M] * (
            id component[p1,\[Mu]]  -
                2 * component[k1,\[Mu]] (-fermionChain[slash[k1],slash[-k2]]) propagator[k1, M]
        ) * propagator[k1, M]
    )
) * propagator[k2+p1]
];


(*commented out old version before looking at George's correction.*)
(*If[scheme == GS,
bubble1UV = (-1)*I (1 - eps) * propagator[k1 + p1] * fermionChain[
    propagator[k2, M]^2 * slash[k1 - k2 + minusPart[-k1]/2 + minusPart[k2 + p1]] + propagator[k2, M]^3 * 2 lorentzDot[k2, k1-p1] * slash[- k2 + minusPart[k2]]  - 
    propagator[k2, M]^2 * slash[k1 - k2 + p1 + minusPart[-k1]/2 + minusPart[k2 + p1]] - propagator[k2, M]^3 * 2 lorentzDot[k2, -k1-p1] * slash[- k2 + minusPart[k2]], 
    gamma[\[Mu]]]
];*)
If[scheme == GS,
bubble1UV = (-1)*I (1 - eps) * propagator[k1 + p1] * fermionChain[
    propagator[k2, M]^2 * slash[k1 - k2 + minusPart[-k1]/2 + minusPart[k2 + p1]] + propagator[k2, M]^3 * 2 lorentzDot[k2, k1-p1] * slash[- k2 + minusPart[k2]]  - 
    propagator[k2, M]^2 * slash[k1 - k2 + p1 - minusPart[-k1]/2 + minusPart[k2 + p1]] - propagator[k2, M]^3 * 2 lorentzDot[k2, -k1-p1] * slash[- k2 + minusPart[k2]], 
    gamma[\[Mu]]]
];


If[scheme == gammaHat,
vertex1UV = 2 * I * (1 - eps) * g^3 * fermionChain[
    slash[k2] + slash[p1],
    slash[k1], 
    gamma[\[Mu]],
    slash[k1]
] * propagator[k2 + p1] * propagator[k1, M]^3 - 
(*extra term; must use id to denote 4*4 Dirac Matrix!*)
id * 2 * I * (2 - 2 eps) * g^3 (
    -(1/2) component[p1, \[Mu]] * lorentzDot[k1, k1] + 
    2 component[k1, \[Mu]] lorentzDot[k1, p1]
) * (1 + (2 lorentzDot[p2, k2]) / s) * propagator[k2 + p1] * propagator[k1, M]^3
];


If[scheme == gammaHat,
bubble1UV = I * (1-eps) * g^2 * fermionChain[slash[k1+p1], gamma[\[Mu]]] * propagator[k2,M]^2 * propagator[k1+p1]
];


bubble2UV = bubble1UV /. {fermionChain -> fermionChainReverse, p1->-p2, minusPart -> plusPart};


vertex2UV = vertex1UV /. {fermionChain -> fermionChainReverse, p1->-p2, p2->-p1, k1->k2, k2->k1};


(* ::Section:: *)
(*Load diagrams*)


diags = Get["qgrafProcessed/twoLoopTwoPhotonSubtracted.m"];


Print["Loaded diagrams including categories:"];
Print[Keys[diags]];


(* ::Section:: *)
(*Numerical evaluation*)


(* ::Subsection::Closed:: *)
(*Routines for evaluating lorentz and Dirac expression*)


dot[___, 0, ___]=0;
dot[a___] := Dot[a];
TrDot[a___]:=Tr[dot[a]];


signature={1,-1,-1,-1};


component[p:{_,_,_,_}, i_Integer]:=p[[i+1]];


(* Lowering Lorentz index of a four-vector*)
lower[p:{_,_,_,_}]:=p*signature;
upper=lower;


(*Lorentz dot product*)
lorentzDot[p1:{_,_,_,_}, p2:{_,_,_,_}]:=lower[p1].p2;
lorentzDot[p1:{_,_,_,_}]:=lorentzDot[p1,p1];


propagator[p1:{_,_,_,_}]:=1 / lorentzDot[p1];
propagator[p1:{_,_,_,_}, M_]:=1 / (lorentzDot[p1] - M^2);


metricTensor[\[Mu]_/;0<=\[Mu]<=3, \[Nu]_/;0<=\[Nu]<=3]:=KroneckerDelta[\[Mu], \[Nu]]*signature[[1+\[Mu]]];


pauli[1]=({
 {0, 1},
 {1, 0}
});
pauli[2]=({
 {0, -I},
 {I, 0}
});
pauli[3]=({
 {1, 0},
 {0, -1}
});
pauli[0]=({
 {1, 0},
 {0, 1}
});


(* Alternatively, pauliVec[[i]] = pauli[i] for 1\[LessEqual]i\[LessEqual]3 *)
pauliVec = Array[pauli, 3];


(*gamma matrices in chiral representation *)
gamma[n_/;0<=n<=3]:=Join[
Join[({
 {0, 0},
 {0, 0}
}), pauli[n], 2],
Join[signature[[1+n]]*pauli[n], ({
 {0, 0},
 {0, 0}
}), 2]
];


gammaFourVec=Table[gamma[n], {n,0,3}];


slash[mom:{_,_,_,_}]:=lorentzDot[mom, gammaFourVec];


(*gamma matrice with lower Lorentz index*)
gammaLower[n_/;0<=n<=3]:=gamma[n] * signature[[n+1]];


gammaLowerFourVec=Table[gammaLower[n], {n,0,3}];


s := lorentzDot[p1+p2];


Assert[(lorentzDot[q1, polVec[-2]]/.externalKinematics) == (lorentzDot[q2, polVec[-4]]/.externalKinematics) == 0];


Assert[
    Flatten[((slash[p1]. spinor[-1])/.externalKinematics)] ==
    Flatten[((spinor[-3] . slash[p2])/.externalKinematics)] ==
    {0,0,0,0}
];


(*SetAttributes[reverseDot, HoldAll];
reverseDot[a__]:=dot @@ Reverse[Hold[a]];*)


loopParam = {
    k1 -> k[1, plus] * {1,0,0,1} + k[1, minus] * {1,0,0,-1} + {0, k[1,x], k[1,y], 0},
    k2 -> k[2, plus] * {1,0,0,1} + k[2, minus] * {1,0,0,-1} + {0, k[2,x], k[2,y], 0}
};


symmetrizeK1K2[expr_] := 1/2 * expr + (1/2 * expr /. {k[1, a_] :> k[2, a], k[2, a_] :> k[1, a]});


symmetrizeTransverse[expr_] := 1/2 * expr + (1/2 * expr /. {k[1, a:(x|y)] :> -k[1, a], k[2, a:(x|y)] :> -k[2, a]});


(*eval[diag_, externalKinematics_] := Module[{feynamp, repeatedIndices, uniqueIndices, sumLimits, summand, levels},
    feynamp = insertFeynmanRules[diag];
    repeatedIndices = Cases[feynamp, lorentzIndex[_], Infinity] // DeleteDuplicates;
    If[repeatedIndices === {},
    (
    feynamp /.
            {lambda->1, g->1, eps->0} /.
            Join[externalKinematics, loopParam]  /.
            {fermionChain->dot, fermionTrace[a__]:>Tr[dot[a]], id->IdentityMatrix[4]}
    ),
    (
    sumLimits = {#, 0, 3} & /@ repeatedIndices;
    levels = Length[repeatedIndices];
    Table[
        feynamp /.
            {lambda->1, g->1, eps->0} /.
            Join[externalKinematics, loopParam] /.
            {fermionChain->dot, fermionTrace[a__]:>Tr[dot[a]], id->IdentityMatrix[4]},
        Evaluate[Sequence @@ sumLimits]
    ] // Flatten[#, levels-1]& // Total
    )
    ]
];*)


eval[diag_, externalKinematics_List, internalNumerics_List, internalScaling_List] :=
evalFeynAmp[insertFeynmanRules[diag], externalKinematics, internalNumerics, internalScaling];


(*e.g. internalScaling can be {k[1, i_] \[RuleDelayed] \[Delta]*k[1,i]} for the single-soft limit, while internalNumerics should give all components of both k1 and k2.*)
evalFeynAmp[feynamp_, externalKinematics_List, internalNumerics_List, internalScaling_List] :=
Module[{repeatedIndices, uniqueIndices, sumLimits, summand, levels, loopParamNum, extraNumerics},
    extraNumerics = Select[internalNumerics, !MatchQ[#[[1]], k[1|2, plus|minus|x|y]]&];
    loopParamNum = Join[
        Table[loopParam[[i, 1]] -> (loopParam[[i, 2]] /. internalScaling /. internalNumerics), {i, Length[loopParam]}],
        extraNumerics
    ];
    repeatedIndices = Cases[feynamp, lorentzIndex[_], Infinity] // DeleteDuplicates;
    If[repeatedIndices === {},
    (
    feynamp /.
            {lambda->1, g->1, eps->0} /.
            Join[externalKinematics, loopParamNum,
                {fermionChain->dot, fermionTrace->TrDot, id->IdentityMatrix[4]}
            ]
    ),
    (
    sumLimits = {#, 0, 3} & /@ repeatedIndices;
    levels = Length[repeatedIndices];
    Table[
        feynamp /.
            {lambda->1, g->1, eps->0} /.
            Join[externalKinematics, loopParamNum,
                {fermionChain->dot, fermionTrace->TrDot, id->IdentityMatrix[4]}
            ],
        Evaluate[Sequence @@ sumLimits]
    ] // Flatten[#, levels-1]& // Total
    )
    ]
];


eval[diag_, externalKinematics_] := eval[diag, externalKinematics, {k[1, plus], k[1, minus], k[1,x], k[1,y], k[2, plus], k[2, minus], k[2,x], k[2,y], M}];


eval[diag_, externalKinematics_, numericalValues_List /; FreeQ[numericalValues, k[1|2, plus|minus|x|y] | M] && Length[numericalValues] == 9] :=
    eval[diag, externalKinematics, Thread[{k[1, plus], k[1, minus], k[1,x], k[1,y], k[2, plus], k[2, minus], k[2,x], k[2,y], M} -> numericalValues], {}];


(* ::Subsection:: *)
(*Prepare all diagrams for evaluation*)


Print["inserting numerical external kinematics"];
diagsParam = Association[];
Do[
    (
        Print[key];
        diagsParam[key] = eval[#, externalKinematics] & /@ diags[key]
    ),
    {key, Keys[diags]}
];


Print["Symmetrizing photonic diagrams."];
diagsParamSym = Association[];
Do[
    (
        Print[key];
        diagsParamSym[key] = If[scheme === GS,
            symmetrizeTransverse /@ symmetrizeK1K2 /@ diagsParam[key],
            symmetrizeK1K2 /@ diagsParam[key]
        ];
    ),
    {key, Select[Keys[diags], StringMatchQ[#[[1]], "photonic"~~___]&]}
];


diagsSymTotal = Total[diagsParamSym[{"photonicHardSubtracted", "mod"}]] +
    Total[diagsParamSym[{"photonicHardSubtracted", "double"}]] +
    Total[diagsParamSym[{"photonicHardSubtracted", "single"}]];


diagsFbubbleTotal = Total[diagsParam[{"fbubbleHardSubtracted"}]];


diagsFboxTotal = Total[diagsParam[{"fboxHardSubtracted"}]];


(* ::Section::Closed:: *)
(*Consistency checks*)


subUVfudgeFactor=1;


(* ::Subsection:: *)
(*Consistency check - if we replace modified vertices and bubbles by original expressions, we should agree*)
(*with automated evaluation of original diagrams*)


If[KeyExistsQ[diags, {"photonic", "orig"}],
Assert[
    numDiffOrigMod[n_Integer] := Block[{insertFeynmanRules = insertFeynmanRulesNoMod},
    eval[diags[{"photonic", "orig"}][[n]], externalKinematics] -
        eval[diags[{"photonic", "mod"}][[n]], externalKinematics] /.
        randomNumerics
    ];
    modifiedPositions = Position[diags[{"photonic", "mod"}],
        modifiedVertex1 | modifiedVertex2 | modifiedBubble1 | modifiedBubble2
    ][[All,1]];
    Union[(numDiffOrigMod /@ modifiedPositions)] == {0}
]
];


(* ::Subsection:: *)
(*Consistency check for gammaHat scheme - modified bubbles & vertices are suppressed in region-2*)


If[scheme == gammaHat,
Assert[
Module[{position},
position = Intersection[
    Position[diags[{"photonic", "mod"}], modifiedVertex1[__]][[All,1]],
    Position[diags[{"photonic", "mod"}], projectOutScalar[__]][[All,1]]
][[1]];
    0 == (
    diagsParam[{"photonic", "mod"}][[position]] /.
    {k[2, plus] :> k[2, plus] * \[Delta]^2, k[2, i:(x|y)]:> k[2, i] * \[Delta]} /. randomNumerics // Series[#, {\[Delta],0,-4}]& // Normal
    )
]
]
];


If[scheme == gammaHat,
Assert[
Module[{position},
position = Intersection[
    Position[diags[{"photonic", "mod"}], modifiedBubble1[__]][[All,1]],
    Position[diags[{"photonic", "mod"}], projectOutScalar[__]][[All,1]]
][[1]];
    0 == (
    diagsParam[{"photonic", "mod"}][[position]] /.
    {k[1, plus] :> k[1, plus] * \[Delta]^2, k[1, i:(x|y)]:> k[1, i] * \[Delta]} /. randomNumerics // Series[#, {\[Delta],0,-4}]& // Normal
    )
]
]
];


If[scheme == gammaHat,
Assert[
Module[{position},
position = Intersection[
    Position[diags[{"photonic", "mod"}], modifiedVertex2[__]][[All,1]],
    Position[diags[{"photonic", "mod"}], projectOutScalar[__]][[All,1]]
][[1]];
    0 == (
    diagsParam[{"photonic", "mod"}][[position]] /.
    {k[1, minus] :> k[1, minus] * \[Delta]^2, k[1, i:(x|y)]:> k[1, i] * \[Delta]} /. randomNumerics // Series[#, {\[Delta],0,-4}]& // Normal
    )
]
]
];


If[scheme == gammaHat,
Assert[
Module[{position},
position = Intersection[
    Position[diags[{"photonic", "mod"}], modifiedBubble2[__]][[All,1]],
    Position[diags[{"photonic", "mod"}], projectOutScalar[__]][[All,1]]
][[1]];
    0 == (
    diagsParam[{"photonic", "mod"}][[position]] /.
    {k[1, minus] :> k[1, minus] * \[Delta]^2, k[1, i:(x|y)]:> k[1, i] * \[Delta]} /. randomNumerics // Series[#, {\[Delta],0,-4}]& // Normal
    )
]
]
];


(* ::Subsection:: *)
(*Consistency check for gammaHat scheme - modified bubbles & vertices have single-soft divergences canceled by the double-IR form factor CT*)


If[scheme == gammaHat,
Assert[
Module[{position},
position = Position[diags[{"photonic", "mod"}], projectOutScalar[__]][[All,1]];
{0} == Union[
    diagsParam[{"photonic", "mod"}][[position]] + diagsParam[{"photonic", "double"}][[position]] /.
    {k[2, i_]:> k[2, i] * \[Delta]} /. randomNumerics // Series[#, {\[Delta],0,-4}]& // Normal
    ]
]
]
];


(* ::Subsection:: *)
(*Consistency check for gammaHat scheme - Lack of "loop polarizations" in modified vertices - O(1/\[Delta]^2) terms should cancel*)


(*Print["Checking O(1/\[Delta]^2) cancellation from lack of loop polarizations in modified vertex 1: "];*)
Assert[{0,0,0,0} == (fermionChain[vertex1modified * component[lower[p1], \[Mu]], spinor[-1]] /. \[Mu] -> lorentzIndex[1] //
    eval[#, externalKinematics] /. {k[2, minus] -> k[2, minus] * \[Delta]^2, k[2, i:(x|y)] :> k[2, i] * \[Delta]} /. randomNumerics& //
    Series[#, {\[Delta], 0, -2}]& // Normal)
];


(*Print["Checking O(1/\[Delta]^2) cancellation from lack of loop polarizations in modified vertex 2: "];*)
Assert[{0,0,0,0} == (fermionChain[spinor[-3], vertex2modified * component[lower[p2], \[Mu]]] /. \[Mu] -> lorentzIndex[1] //
    eval[#, externalKinematics] /. {k[1, plus] -> k[1, plus] * \[Delta]^2, k[1, i:(x|y)] :> k[1, i] * \[Delta]} /. randomNumerics& //
    Series[#, {\[Delta], 0, -2}]& // Normal)
];


(* ::Subsection:: *)
(*Consistency check for GS scheme - Ward identities relating modified vertices to modified bubbles in region-2*)


(*Print["Check region-2 O(\[Delta]^0) cancellation between modified vertex 1 and modifed bubble 1: "];*)
If[scheme == GS,
Assert[{0,0,0,0} == Normal[
(fermionChain[vertex1modified * component[lower[p2], \[Mu]], spinor[-1]] /. \[Mu] -> lorentzIndex[1] //
    eval[#, externalKinematics] /. {k[2, plus] -> k[2, plus]* \[Delta]^2, k[2, i:(x|y)] :> k[2, i] * \[Delta]} /. randomNumerics& //
    Series[#, {\[Delta], 0, 0}]&) +
(fermionChain[bubble1modified * component[lower[p2], \[Mu]], spinor[-1]] /. \[Mu] -> lorentzIndex[1] //
    eval[#, externalKinematics] /. {k[2, a_]:>k[1,a], k[1, a_]:>k[2,a]} /. {k[2, plus] -> k[2, plus]* \[Delta]^2, k[2, i:(x|y)] :> k[2, i] * \[Delta]} /. randomNumerics& //
    Series[#, {\[Delta], 0, 0}]&)
]
]
];


(*Print["Check region-2 O(\[Delta]^0) cancellation between modified vertex 2 and modifed bubble 2: "];*)
If[scheme == GS,
Assert[{0,0,0,0} == Normal[
(fermionChain[vertex2modified * component[lower[p1], \[Mu]], spinor[-1]] /. \[Mu] -> lorentzIndex[1] //
    eval[#, externalKinematics] /. {k[1, minus] -> k[1, minus]* \[Delta]^2, k[1, i:(x|y)] :> k[1, i] * \[Delta]} /. randomNumerics& //
    Series[#, {\[Delta], 0, 0}]&) +
(fermionChain[bubble2modified * component[lower[p1], \[Mu]], spinor[-1]] /. \[Mu] -> lorentzIndex[1] //
    eval[#, externalKinematics] /. {k[1, minus] -> k[1, minus]* \[Delta]^2, k[1, i:(x|y)] :> k[1, i] * \[Delta]} /. randomNumerics& //
    Series[#, {\[Delta], 0, 0}]&)
]
]
];


(* ::Subsection:: *)
(*Consistency check for GS scheme - Ward identities relating modified vertices to modified bubbles in soft region*)


(*Print["Check soft region O(\[Delta]^-1) cancellation between modified vertex 1 and modifed bubble 1: "];*)
If[scheme == GS,
Assert[{0,0,0,0} == Normal[
(fermionChain[vertex1modified * component[lower[p2], \[Mu]], spinor[-1]] /. \[Mu] -> lorentzIndex[1] //
    eval[#, externalKinematics] /. {k[2, i_] :> k[2, i] * \[Delta]} /. randomNumerics& //
    Series[#, {\[Delta], 0, -1}]&) +
(fermionChain[bubble1modified * component[lower[p2], \[Mu]], spinor[-1]] /. \[Mu] -> lorentzIndex[1] //
    eval[#, externalKinematics] /. {k[2, a_]:>k[1,a], k[1, a_]:>k[2,a]} /. {k[2, i_] :> k[2, i] * \[Delta]} /. randomNumerics& //
    Series[#, {\[Delta], 0, -1}]&)
]
]
];


(*Print["Check soft region O(\[Delta]^-1) cancellation between modified vertex 2 and modifed bubble 2: "];*)
If[scheme == GS,
Assert[{0,0,0,0} == Normal[
(fermionChain[vertex2modified * component[lower[p1], \[Mu]], spinor[-1]] /. \[Mu] -> lorentzIndex[1] //
    eval[#, externalKinematics] /. {k[1, i_] :> k[1, i] * \[Delta]} /. randomNumerics& //
    Series[#, {\[Delta], 0, -1}]&) +
(fermionChain[bubble2modified * component[lower[p1], \[Mu]], spinor[-1]] /. \[Mu] -> lorentzIndex[1] //
    eval[#, externalKinematics] /. {k[1, minus] -> k[1, minus]* \[Delta]^2, k[1, i_] :> k[1, i] * \[Delta]} /. randomNumerics& //
    Series[#, {\[Delta], 0, -1}]&)
]
]
];


(* ::Subsection:: *)
(*Consistency check - cancellation of sub-UV divergence in diagrams w/ hard bubble*)


Assert[ 0 ==
(diags[{"photonicHardSubtracted","mod"}][[13]] // eval[#, externalKinematics] /. {k[2, i_] :> k[2, i] / \[Delta]} /.
    randomNumerics & // Series[#, {\[Delta],0,4}]& // Normal)
];


(* ::Subsection:: *)
(*Consistency check - cancellation of sub-UV divergence in diagrams w/ hard photon vertex*)


Assert[ 0 ==
(diags[{"photonicHardSubtracted","mod"}][[3]] // eval[#, externalKinematics] /. {k[2, i_] :> k[2, i] / \[Delta]} /.
    randomNumerics & // Series[#, {\[Delta],0,4}]& // Normal)
];


(* ::Subsection:: *)
(*Consistency check - cancellation of sub-UV divergence in diagrams w/ modified vertex / bubble*)


Position[diags[{"photonicHardSubtracted","mod"}], modifiedBubble1Subtracted] // Print;


Assert[
    diags[{"photonicHardSubtracted","mod"}][[11]] // eval[#, externalKinematics] /. {k[2, i_] :> k[2, i] / \[Delta]} /.
    randomNumerics & // Series[#, {\[Delta],0,4}]& // Normal // #==0&
];


Assert[
    diags[{"photonicHardSubtracted","mod"}][[30]] // eval[#, externalKinematics] /. {k[2, i_] :> k[2, i] / \[Delta]} /.
    randomNumerics & // Series[#, {\[Delta],0,4}]& // Normal // #==0&
];


Position[diags[{"photonicHardSubtracted","mod"}], modifiedVertex1Subtracted] // Print;


Assert[
    diags[{"photonicHardSubtracted","mod"}][[9]] // eval[#, externalKinematics] /. {k[1, i_] :> k[1, i] / \[Delta]} /.
    randomNumerics & // Series[#, {\[Delta],0,4}]& // Normal // #==0&
];


Position[diags[{"photonicHardSubtracted","mod"}], modifiedBubble2Subtracted] // Print;


Assert[
    diags[{"photonicHardSubtracted","mod"}][[10]] // eval[#, externalKinematics] /. {k[2, i_] :> k[2, i] / \[Delta]} /.
    randomNumerics & // Series[#, {\[Delta],0,4}]& // Normal // #==0&
];


Assert[
    diags[{"photonicHardSubtracted","mod"}][[29]] // eval[#, externalKinematics] /. {k[2, i_] :> k[2, i] / \[Delta]} /.
    randomNumerics & // Series[#, {\[Delta],0,4}]& // Normal // #==0&
];


Position[diags[{"photonicHardSubtracted","mod"}], modifiedVertex2Subtracted] // Print;


Assert[
    diags[{"photonicHardSubtracted","mod"}][[7]] // eval[#, externalKinematics] /. {k[2, i_] :> k[2, i] / \[Delta]} /.
    randomNumerics & // Series[#, {\[Delta],0,4}]& // Normal // #==0&
];


Assert[
    diags[{"photonicHardSubtracted","mod"}][[26]] // eval[#, externalKinematics] /. {k[2, i_] :> k[2, i] / \[Delta]} /.
    randomNumerics & // Series[#, {\[Delta],0,4}]& // Normal // #==0&
];


(* ::Subsection:: *)
(*Consistency check - modified fermion bubble diagrams have no singularity associated with UV / soft /collinear limits of the gluon momentum*)


Assert[
    N[(1/\[Delta]^4* Total[diagsParam[{"fbubbleHardSubtracted"}]] /. {k[2, i_] :> k[2, i] / \[Delta]} /.
        Join[randomNumerics, {\[Delta]->1/10^(120/4)}]), 2*$MachinePrecision]//Chop//#==0&
];


Assert[
    N[(\[Delta]^4* Total[diagsParam[{"fbubbleHardSubtracted"}]] /. {k[2, x_] :> k[2,x] * \[Delta]} /. Join[randomNumerics, {\[Delta]->1/10^(120/4)}]), 2*$MachinePrecision]//Chop//#==0&
];


Assert[
    N[(\[Delta]^2*Total[diagsParam[{"fbubbleHardSubtracted"}]]  /. {k[2, minus] -> k[2, minus] * \[Delta], k[2, i:(x|y)]:> k[2, i] * Sqrt[\[Delta]]} /.
        Join[randomNumerics, {\[Delta]->1/10^(120/2)}]), 2*$MachinePrecision]//Chop//#==0&
];


(* ::Section:: *)
(*Add double-UV CTs*)


Get["buildUVfeynamp.m"];


doubleUVlimit[diagram_]:=diagram//insertFeynmanRules//removeNestedFermionChain// insertMassRegulatorAndUVexpand//breakUpFermionChain //
    uvScalingDiagram// Normal[Series[#, {delta, 0, 8}]] /. delta->1 & // Expand[#, diracGamma[__]]& // rebuildFermionChain//#/.id[_]->id&;


doubleUVpos[diags_] := Select[diags //Length //Range,
0 != (N[
  (1/\[Delta]^8* eval[diags[[#]], externalKinematics]/. {k[a_, i_]:> k[a, i] / \[Delta] /. randomNumerics } /.
        Join[randomNumerics, {\[Delta]->1/10^(120/8)}]), 2*$MachinePrecision
]//Chop) &
];


Print["Finding double-UV divergent diagrams. ", DateString[]];
doubleUVpositions = <|
    {"photonicHardSubtracted", "mod"} -> doubleUVpos[diags[{"photonicHardSubtracted", "mod"}]],
    {"photonicHardSubtracted", "single"} -> doubleUVpos[diags[{"photonicHardSubtracted", "single"}]],
    {"photonicHardSubtracted", "double"} -> doubleUVpos[diags[{"photonicHardSubtracted", "double"}]]
|>;
Print["Finished finding double-UV divergent diagrams. ", DateString[]];


Clear[subUVfudgeFactor];


doubleUVCTs = Table[
    key -> Table[(*Print[i];*) If[MemberQ[doubleUVpositions[key], i], (-1)*doubleUVlimit[diags[key][[i]]], 0],
        {i, Length[diags[key]]}
    ],
    {key, Keys[doubleUVpositions]}
] // Association;


Print["inserting numerical external kinematics into double-UV CTs."];
diagsParamDoubleUV = Association[];
Do[
    (
        Print[key];
        diagsParamDoubleUV[key] = eval[#, externalKinematics] & /@ doubleUVCTs[key]
    ),
    {key, Keys[doubleUVCTs]}
];


Print["Symmetrizing double-UV CTs."];
diagsParamSymDoubleUV = Association[];
Do[
    (
        Print[key];
        diagsParamSymDoubleUV[key] = If[scheme === GS,
            symmetrizeTransverse /@ symmetrizeK1K2 /@ diagsParamDoubleUV[key],
            symmetrizeK1K2 /@ diagsParamDoubleUV[key]
        ];
    ),
    {key, Keys[diagsParamDoubleUV]}
];


diagsSymDoubleUVTotal = Total[diagsParamSymDoubleUV[{"photonicHardSubtracted", "mod"}]] +
    Total[diagsParamSymDoubleUV[{"photonicHardSubtracted", "double"}]] +
    Total[diagsParamSymDoubleUV[{"photonicHardSubtracted", "single"}]];


diagsAllTotal = diagsSymTotal + diagsSymDoubleUVTotal;


Print["Finished evaluating diagrams and constructing double-UV CTs. Time used by MathKernel: ", TimeUsed[]];


(* ::Section:: *)
(*Export (unsymmetrized) diagrams for numerical integration*)


subUVfudgeFactor=1;


cleanFermionChain[expr_] := expr//removeNestedFermionChain // breakUpFermionChain // Expand[#, diracGamma[__]]& //
    rebuildFermionChain// # /. id[_]->id  /. fermionChain[a___, id, b___] :> fermionChain[a, b] &;


allFeynAmps = Join[
    diagsParam[{"photonicHardSubtracted", "mod"}],
    diagsParam[{"photonicHardSubtracted", "double"}],
    diagsParam[{"photonicHardSubtracted", "single"}],
    diagsParamDoubleUV[{"photonicHardSubtracted", "mod"}],
    diagsParamDoubleUV[{"photonicHardSubtracted", "double"}],
    diagsParamDoubleUV[{"photonicHardSubtracted", "single"}]] /.
    eval[a_, b___] :> insertFeynmanRules[a] // DeleteCases[#, 0]&;


allFeynAmps = cleanFermionChain /@ allFeynAmps;


Print[Length[allFeynAmps], " integrand expressions to be exported."];


(*Put[allFeynAmps, "twoPhotonIntegrand.m"];*)


If[$runNumericalChecks =!= True,
    Print["$runNumericalChecks is not set to True. Stop further evaluations by Abort[]."];
    Abort[]
];


(* ::Section:: *)
(*Check cancellations in IR / UV limits - photonic diagrams*)


(* The first argument is a diagram expression with numerical external kinematics but analytic internal kinematics.
   The second argument is the scaling rule of +/-/x/y components of x and y, using \[Delta] as a small parameter*)
limitScaling[diag_, limit:{(HoldPattern[Rule[__] | RuleDelayed[__]]) ..}, randomNumerics_] := Module[
    {measure, measureExponent, diagScaled, \[Delta]1, \[Delta]2, result1, result2, result},
    Print["Computing scaling behavior in limit: ", limit];
    measure = {k[1, plus], k[1, minus], k[1, x], k[1, y], k[2, plus], k[2, minus], k[2, x], k[2, y]};
    measureExponent = Exponent[
        Det[Table[D[measure[[i]]/.limit, measure[[j]]], {i, Length[measure]}, {j, Length[measure]}]],
        \[Delta]
    ];
    (*measureExponent = Exponent[measure/.limit, \[Delta]];*)
    Print["measureExponent = ", measureExponent];
    \[Delta]1 = 1 / 10^(120 / Abs[measureExponent] - 1);
    \[Delta]2 = 1 / 10^(120 / Abs[measureExponent] + 1);
    diagScaled = \[Delta]^measureExponent * diag /. limit;
    result1 = N[diagScaled /. Join[randomNumerics, {\[Delta] -> \[Delta]1}], 2*$MachinePrecision];
    result2 = N[diagScaled /. Join[randomNumerics, {\[Delta] -> \[Delta]2}], 2*$MachinePrecision];
    result = Log[result1/result2] / Log[\[Delta]1/\[Delta]2];
    Print["unrounded result: ", result];
    Round[result]
];    


limitScaling[diag_, transformation:{(HoldPattern[Rule[__] | RuleDelayed[__]]) ..},
    limit:{(HoldPattern[Rule[__] | RuleDelayed[__]]) ..}, randomNumerics_] :=
Module[{limit1},
    limit1 = #[[1]] -> (#[[2]] /. limit) & /@ transformation;
    limitScaling[diag, limit1, randomNumerics]
];


Print["Computing limits of photonic diagrams."];


Print["single soft (k1 soft) limit"];
myTiming[
    limitScaling[diagsAllTotal,
    {k[1, x_] :> k[1,x] * \[Delta]},
    randomNumerics]
];


Print["single collinear (k1 // p1) limit"];
myTiming[
    limitScaling[diagsAllTotal,
    {k[1, minus] :> k[1, minus] * \[Delta]^2, k[1, i:(x|y)]:> k[1, i] * \[Delta]},
    randomNumerics]
];


Print["double-soft limit"];
myTiming[
    limitScaling[diagsAllTotal, {k[i_, x_] :> k[i,x] * \[Delta]}, randomNumerics]
];


Print["soft k1, k2 // p1 limit"];
myTiming[
    limitScaling[diagsAllTotal, {k[1, x_] :> k[1,x] * \[Delta]^2, k[2, minus] :> k[2, minus] * \[Delta]^2, k[2, i:(x|y)]:> k[2, i] * \[Delta]}, randomNumerics]
];


Print["k1 // p1, k2 // p2 limit"];
myTiming[
    limitScaling[diagsAllTotal,
    {k[1, minus] :> k[1, minus] * \[Delta]^2, k[1, i:(x|y)]:> k[1, i] * \[Delta],
    k[2, plus] :> k[2, plus] * \[Delta]^2, k[2, i:(x|y)]:> k[2, i] * \[Delta]},
    randomNumerics]
];


Print["k1 & k2 // p1 limit"];
myTiming[
    limitScaling[diagsAllTotal,
    {k[1, minus] :> k[1, minus] * \[Delta]^2, k[1, i:(x|y)]:> k[1, i] * \[Delta],
    k[2, minus] :> k[2, minus] * \[Delta]^2, k[2, i:(x|y)]:> k[2, i] * \[Delta]},
    randomNumerics]
];


Print["double UV"];
myTiming[
    limitScaling[diagsAllTotal, {k[a_, i_] :> k[a, i] / \[Delta]}, randomNumerics]
]


Print["single UV"];
myTiming[
    limitScaling[diagsAllTotal, {k[1, i_] :> k[1, i] / \[Delta]}, randomNumerics]
];


(* ::Section::Closed:: *)
(*Check the lack of divergence in "potential region" due to the massive reference vector \[Eta] - photonic diagrams (only for gammaHat scheme)*)


Print["Computing potential spurious limits of photonic diagrams."];


If[scheme === gammaHat,
myTiming[
    limitScaling[diagsAllTotal,
    {k[1, minus]->k[1, plus]-k[1, minus], k[1,plus]->k[1,plus]+k[1,minus]},
    {k[1, plus]->\[Delta]^2 k[1,plus], k[1, i:(minus|x|y)] :>\[Delta] k[1, i]},
    randomNumerics]
];
];


If[scheme === gammaHat,
myTiming[
    limitScaling[diagsAllTotal,
    {k[1, minus]->k[1, plus]-k[1, minus], k[1,plus]->k[1,plus]+k[1,minus]},
    {k[a_, plus] :> \[Delta]^2 k[a,plus], k[a_, i:(minus|x|y)] :>\[Delta] k[a, i]},
    randomNumerics]
];
]


(* ::Section::Closed:: *)
(*Check cancellations in IR / UV limits - fermion bubble diagrams*)


Print["Computing limits of fermion box diagrams."];


Print["l // p2 limit"];
myTiming[
    limitScaling[diagsFbubbleTotal,
    {k[2, plus] :> k[2, plus] * \[Delta]^2, k[2, i:(x|y)]:> k[2, i] * \[Delta]},
    randomNumerics]
];


Print["l soft limit limit"];
myTiming[
    limitScaling[diagsFbubbleTotal,
    {k[2, x_] :> k[2,x] * \[Delta]},
    randomNumerics]
];


Print["double-soft limit"];
myTiming[
    limitScaling[diagsFbubbleTotal,
    {k[i_, x_] :> k[i,x] * \[Delta]},
    randomNumerics]
];


Print["k, l large limit"];
myTiming[
    limitScaling[diagsFbubbleTotal,
    {k[a_, i_] :> k[a, i] / \[Delta]},
    randomNumerics]
];


Print["k large limit"];
myTiming[
    limitScaling[diagsFbubbleTotal,
    {k[1, i_] :> k[1, i] / \[Delta]},
    randomNumerics]
];


(* ::Section::Closed:: *)
(*Check cancellations in IR / UV limits - fermion box diagrams*)


(*In the paper, k1 ~ l, k2 ~ k*)


Print["Computing limits of fermion box diagrams."];


Print["k1 soft limit"];
myTiming[
    limitScaling[diagsFboxTotal,
    {k[1, i_] :> k[1, i] * \[Delta]},
    randomNumerics]
];


Print["k1 // p1 limit"];
myTiming[
    limitScaling[diagsFboxTotal,
    {k[1, minus] :> k[1, minus]*\[Delta]^2, k[1, i:x|y] :> k[1, i] * \[Delta]},
    randomNumerics]
];


Print["k1 // p2 limit"];
myTiming[
    limitScaling[diagsFboxTotal,
    {k[1, plus] :> k[1, plus]*\[Delta]^2, k[1, i:x|y] :> k[1, i] * \[Delta]},
    randomNumerics]
];


Print["k2 large limit"];
myTiming[
    limitScaling[diagsFboxTotal,
    {k[2, i_] :> k[2, i] / \[Delta]},
    randomNumerics]
];


Print["Checked all limits. Time used by MathKernel: ", TimeUsed[]];


(* ::Section:: *)
(*Exporting double-UV counterterms for analytic integration*)


doubleUVlimit1[diagram_]:=diagram/.projectorP1P2->projectorHere//insertFeynmanRules//
    #/.projectorHere->projectorP1P2&//
    removeNestedFermionChain// insertMassRegulatorAndUVexpand//breakUpFermionChain //
    uvScalingDiagram// SeriesCoefficient[#, {delta, 0, 8}]& // Expand[#, diracGamma[__]]& // rebuildFermionChain//#/.id[_]->id&;


doubleUVCTs1 = Table[
    key -> Table[(*Print[i];*) If[MemberQ[doubleUVpositions[key], i], (-1)*doubleUVlimit1[diags[key][[i]]], 0],
        {i, Length[diags[key]]}
    ],
    {key, Keys[doubleUVpositions]}
] // Association;


Export["doubleUVCTs.txt", 
doubleUVCTs1[{"photonicHardSubtracted","mod"}]/.lambda->e]


Put[doubleUVCTs1[{"photonicHardSubtracted","mod"}]/.lambda->e, "doubleUVCTs.m"];


Export["doubleUVCTs-FF.txt", 
doubleUVCTs1[{"photonicHardSubtracted","double"}]/.lambda->e]


Put[doubleUVCTs1[{"photonicHardSubtracted","double"}]/.lambda->e, "doubleUVCTs-FF.m"];
